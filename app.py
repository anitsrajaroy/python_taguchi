from flask import Flask, render_template, request
import io, csv, math
from decimal import Decimal, InvalidOperation
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import base64

app = Flask(__name__)

PAPER_FACTORS = [
    {"name": "Load (N)", "levels": ["10", "30", "40", "50"]},
    {"name": "Speed (m/s)", "levels": ["0.6", "0.9", "1.2", "1.5"]},
    {"name": "Distance (km)", "levels": ["0.25", "0.5", "0.75", "1"]},
]
PAPER_ROWS = [
[10,0.6,0.25,0.003041367],[10,0.9,0.5,0.006429936],[10,1.2,0.75,0.019119821],[10,1.5,1,0.054649643],
[30,0.6,0.5,0.089785191],[30,0.9,0.25,0.117041974],[30,1.2,1,0.166156451],[30,1.5,0.75,0.180918178],
[40,0.6,0.75,0.215269354],[40,0.9,1,0.227123522],[40,1.2,0.25,0.26218299],[40,1.5,0.5,0.292906082],
[50,0.6,1,0.306220215],[50,0.9,0.75,0.304216223],[50,1.2,0.5,0.297652719],[50,1.5,0.25,0.30806682],
]

def parse_factor_config(raw):
    factors = []
    for line in raw.strip().splitlines():
        if not line.strip():
            continue
        if "|" not in line:
            raise ValueError("Each factor line must use: Factor name | level1, level2, level3")
        name, levels = line.split("|", 1)
        vals = [v.strip() for v in levels.split(",") if v.strip()]
        if not name.strip() or len(vals) < 2:
            raise ValueError("Each factor needs a name and at least two levels.")
        if len(set(vals)) != len(vals):
            raise ValueError(f"Duplicate levels found for factor {name.strip()}.")
        factors.append({"name": name.strip(), "levels": vals})
    if not factors:
        raise ValueError("Enter at least one factor.")
    return factors

def _normalise_column_name(value):
    """Normalise headers so BOMs, extra spaces and optional units do not break matching."""
    import re
    value = str(value).replace("\ufeff", "").strip().lower()
    value = re.sub(r"\([^)]*\)", "", value)  # e.g. 'Load (N)' -> 'Load'
    value = re.sub(r"\s+", " ", value)
    return value.strip()

def _canonical_level(value):
    """Canonicalise numeric levels so values like 1, 1.0 and 1.000 match."""
    text = str(value).replace("\ufeff", "").strip()
    try:
        number = Decimal(text)
        if number.is_finite():
            return format(number.normalize(), "f")
    except (InvalidOperation, ValueError):
        pass
    return text.casefold()

def parse_doe(raw, factors, response_name):
    raw = (raw or "").strip()
    if not raw:
        raise ValueError("DOE table is empty. Paste a CSV table with a header row and experiment data.")
    try:
        # Excel commonly copies data as tab-delimited text; some regional settings use semicolons.
        # Let pandas detect the delimiter rather than assuming commas only.
        try:
            df = pd.read_csv(io.StringIO(raw), sep=None, engine="python")
        except Exception:
            df = None
        if df is None or len(df.columns) == 1:
            candidates = []
            for delimiter in [",", "\t", ";", "|"]:
                try:
                    candidate = pd.read_csv(io.StringIO(raw), sep=delimiter)
                    candidates.append(candidate)
                except Exception:
                    continue
            if candidates:
                df = max(candidates, key=lambda item: len(item.columns))
        if df is None:
            raise ValueError("Could not detect the table delimiter.")
    except Exception as exc:
        raise ValueError("Could not read the DOE table. Paste comma-, tab-, or semicolon-delimited data. Details: " + str(exc))

    # Clean headers and map them to the configured names. Matching ignores case,
    # surrounding spaces, BOM characters, and optional units such as '(N)' or '(g)'.
    df.columns = [str(c).replace("\ufeff", "").strip() for c in df.columns]
    expected = [f["name"] for f in factors] + [response_name]
    header_lookup = {}
    for col in df.columns:
        header_lookup.setdefault(_normalise_column_name(col), col)
    rename = {}
    missing = []
    for wanted in expected:
        actual = header_lookup.get(_normalise_column_name(wanted))
        if actual is None:
            missing.append(wanted)
        else:
            rename[actual] = wanted
    if missing:
        detected = ", ".join(repr(c) for c in df.columns) or "(no columns detected)"
        raise ValueError(
            "DOE table is missing column(s): " + ", ".join(missing) +
            ". Detected column(s): " + detected +
            ". Ensure the first row contains the factor names and response name, or click 'Load paper data'."
        )
    df = df.rename(columns=rename)[expected].copy()
    if df.empty:
        raise ValueError("DOE table must contain at least one experiment.")
    for f in factors:
        column = f["name"]
        configured = {_canonical_level(level): level for level in f["levels"]}
        observed = df[column].astype(str).str.replace("\ufeff", "", regex=False).str.strip()
        bad = sorted({value for value in observed if _canonical_level(value) not in configured})
        if bad:
            raise ValueError(
                f"DOE values for '{column}' are not listed in its levels: {', '.join(bad)}. "
                f"Configured levels: {', '.join(f['levels'])}. Numeric equivalents such as 1 and 1.0 are accepted."
            )
        # Store the exact configured level label after validating numeric equivalence.
        # This keeps grouping, optimum-level lookup and plot labels consistent.
        df[column] = observed.map(lambda value: configured[_canonical_level(value)])
    df[response_name] = pd.to_numeric(df[response_name], errors="coerce")
    if df[response_name].isna().any():
        raise ValueError("Response values must all be numeric.")
    if (df[response_name] == 0).any() or (df[response_name] < 0).any():
        raise ValueError("S/N calculations require positive, non-zero response values.")
    return df

def calculate_sn(y, condition):
    y = np.asarray(y, dtype=float)
    if condition == "maximum":
        return -10.0 * np.log10(np.mean(1.0 / (y ** 2)))
    return -10.0 * np.log10(np.mean(y ** 2))

def make_plot(df, factors, response_name, condition):
    fig, axes = plt.subplots(1, len(factors), figsize=(5 * len(factors), 4), squeeze=False)
    for ax, factor in zip(axes[0], factors):
        means = df.groupby(factor["name"], sort=False)["S/N ratio (dB)"].mean()
        levels = factor["levels"]
        vals = [means.get(level, np.nan) for level in levels]
        ax.plot(range(1, len(levels)+1), vals, marker="o", linewidth=2.5, markersize=7)
        ax.set_xticks(range(1, len(levels)+1), levels)
        ax.set_xlabel(factor["name"])
        ax.set_ylabel("Mean S/N ratio (dB)")
        ax.set_title(factor["name"])
        ax.grid(True, alpha=.25)
    fig.suptitle(f"Main effects plot — {condition.replace('_',' ').title()} response", fontsize=13, fontweight="bold")
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("ascii")

def analyze(df, factors, response_name, condition):
    df = df.copy()
    df["S/N ratio (dB)"] = df[response_name].apply(lambda y: float(calculate_sn([y], condition)))
    grand_mean = df["S/N ratio (dB)"].mean()
    optimum = []
    level_tables = {}
    for f in factors:
        means = df.groupby(f["name"])["S/N ratio (dB)"].mean()
        ordered = [(level, float(means[level])) for level in f["levels"] if level in means.index]
        if not ordered:
            raise ValueError(f"No experimental rows found for factor '{f['name']}'.")
        best = (max if condition == "maximum" else min)(ordered, key=lambda item: item[1])
        optimum.append({"factor": f["name"], "level": best[0], "mean_sn": best[1], "levels": f["levels"],
                        "means": [{"level": lv, "mean_sn": val} for lv, val in ordered]})
        level_tables[f["name"]] = means
    predicted_sn = grand_mean + sum(item["mean_sn"] - grand_mean for item in optimum)
    if condition == "maximum":
        predicted_response = 10 ** (predicted_sn / 20.0)
    else:
        predicted_response = 10 ** (-predicted_sn / 20.0)

    # Classical one-factor-at-a-time ANOVA on S/N ratios, with pooled residual error.
    n = len(df)
    total_ss = float(((df["S/N ratio (dB)"] - grand_mean) ** 2).sum())
    rows = []
    model_ss = 0.0
    total_df = n - 1
    factor_dfs = 0
    for f in factors:
        name = f["name"]
        grouped = df.groupby(name)["S/N ratio (dB)"]
        ss = float(sum(len(g) * (g.mean() - grand_mean) ** 2 for _, g in grouped))
        df_factor = max(len(grouped) - 1, 0)
        ms = ss / df_factor if df_factor else float("nan")
        model_ss += ss
        factor_dfs += df_factor
        rows.append({"source": name, "df": df_factor, "ss": ss, "ms": ms, "f": None, "p": None, "contribution": 0.0})
    residual_ss = max(total_ss - model_ss, 0.0)
    residual_df = total_df - factor_dfs
    residual_ms = residual_ss / residual_df if residual_df > 0 else float("nan")
    for r in rows:
        if residual_df > 0 and residual_ms > 0 and np.isfinite(r["ms"]):
            r["f"] = r["ms"] / residual_ms
            r["p"] = float(stats.f.sf(r["f"], r["df"], residual_df))
        r["contribution"] = (r["ss"] / total_ss * 100.0) if total_ss else 0.0
    rows.append({"source": "Residual error", "df": residual_df, "ss": residual_ss,
                 "ms": residual_ms, "f": None, "p": None,
                 "contribution": residual_ss / total_ss * 100 if total_ss else 0.0})
    rows.append({"source": "Total", "df": total_df, "ss": total_ss, "ms": None, "f": None, "p": None, "contribution": 100.0 if total_ss else 0.0})

    # Predicted S/N for all observed trials and selected optimum setting.
    df["Predicted S/N (additive)"] = grand_mean
    for item in optimum:
        df["Predicted S/N (additive)"] += df[item["factor"]].map(level_tables[item["factor"]]) - grand_mean
    df["Prediction residual"] = df["S/N ratio (dB)"] - df["Predicted S/N (additive)"]
    plot = make_plot(df, factors, response_name, condition)
    return {
        "df": df, "grand_mean": grand_mean, "optimum": optimum, "predicted_sn": predicted_sn,
        "predicted_response": predicted_response, "anova": rows, "plot": plot,
        "n": n, "response_name": response_name, "condition": condition
    }

def default_csv():
    header = ["Load (N)", "Speed (m/s)", "Distance (km)", "Mass loss (g)"]
    rows = [header] + PAPER_ROWS
    out = io.StringIO()
    writer = csv.writer(out)
    writer.writerows(rows)
    return out.getvalue().strip()

@app.route("/", methods=["GET", "POST"])
def index():
    default_factor_text = "\n".join(f'{f["name"]} | {", ".join(f["levels"])}' for f in PAPER_FACTORS)
    default_response = "Mass loss (g)"
    default_condition = "minimum"
    doe_text = default_csv()
    result = None
    error = None
    if request.method == "POST":
        default_factor_text = request.form.get("factor_config", default_factor_text)
        default_response = request.form.get("response_name", default_response).strip() or "Response"
        default_condition = request.form.get("condition", "minimum")
        doe_text = request.form.get("doe_csv", doe_text)
        try:
            factors = parse_factor_config(default_factor_text)
            df = parse_doe(doe_text, factors, default_response)
            result = analyze(df, factors, default_response, default_condition)
            result["table_html"] = result["df"].to_html(index=False, classes="data-table", border=0, float_format=lambda x: f"{x:.6g}")
        except Exception as exc:
            error = str(exc)
    return render_template("index.html", factor_config=default_factor_text, response_name=default_response,
                           condition=default_condition, doe_csv=doe_text, result=result, error=error)

if __name__ == "__main__":
    app.run(debug=True)
