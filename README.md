# TaguchiLab — Flask DOE Optimizer

A mobile-responsive Flask application for Taguchi-style factor-level optimization. It includes editable factor configuration and DOE CSV inputs, minimum/maximum objectives, S/N calculations, a main-effects plot, ANOVA on S/N ratios, optimum factor levels, and additive-model predictions.

## Run locally

1. Install Python 3.10 or later.
2. Extract the ZIP and open a terminal in the project folder.
3. Create and activate a virtual environment (recommended).
4. Install packages:

   ```bash
   pip install -r requirements.txt
   ```

5. Start the app:

   ```bash
   python app.py
   ```

6. Open `http://127.0.0.1:5000`.

## Entering your own experiment

### Factor configuration
One factor per line, formatted as `Factor name | level1, level2, level3`. Factor names must match the DOE CSV column headers exactly. The levels in the DOE table must match the configured labels.

### DOE table
Paste comma-separated values with a header row. Include each factor column and one response column. Set the response column name in the response input. Response values must be positive and non-zero for the implemented S/N formulas.

## Implemented calculations

- **Smaller-the-better:** `S/N = -10 log10(mean(y²))`
- **Larger-the-better:** `S/N = -10 log10(mean(1/y²))`
- **Optimum level:** highest mean S/N at each factor.
- **Predicted optimum S/N:** `grand mean + sum(optimum factor-level mean − grand mean)`.
- **Predicted response:** single-observation equivalent obtained by inverting the predicted S/N.
- **ANOVA:** factor sums of squares on S/N, pooled residual error, F and p values when estimable, and percentage contributions.
- **Main effects plot:** mean S/N at each configured level.

## Important statistical notes

- The supplied example uses the 16 experimental runs and mass-loss values transcribed from the uploaded paper, “Sliding wear behavior of Ti6Al4V implant alloy in the Simulated body Environment” (2017).
- The paper reports a predicted S/N of 45.61 dB and predicted mass loss of 0.005205 g at Load = 10 N, Speed = 0.6 m/s, Distance = 0.25 km. The app recomputes results from the table; small differences may occur because the paper's printed S/N table and ANOVA contain rounding or method-specific differences.
- The app computes ANOVA on S/N values. This is not guaranteed to reproduce every MINITAB output unless the same response transformation, pooling method, and degrees-of-freedom conventions are used.
- The additive prediction is a Taguchi main-effects estimate. Run a confirmation experiment at the recommended settings before adopting them.
- The app accepts any DOE table whose rows and factor-level combinations are appropriate for the design. It does not automatically construct or validate orthogonal arrays, detect aliasing, or fit interaction effects.
- This is an educational analysis utility, not a replacement for statistical review of the experiment design.

## Project structure

```text
taguchi_flask_project/
├── app.py
├── requirements.txt
├── README.md
├── templates/
│   └── index.html
└── static/
    └── style.css
```
