# Competition models for ENSCI 3120

R/Shiny app using the **same Shinylive + GitHub Actions + GitHub Pages delivery architecture as the SIR app**. Each browser runs its own R session. No shinyapps.io hosting or shared R server is required.

## Publish on GitHub Pages

1. Extract the ZIP and open `competition-shiny-app`. Keep the hidden `.github` folder.
2. Create a **new, empty repository**, for example `competitionShinyApp`. Do not initialize it with a README when using the commands below. Keep your SIR repository unchanged.
3. **Enable Pages first:** repository **Settings → Pages → Build and deployment → Source → GitHub Actions**. This avoids the “Get Pages site failed / Not Found” error encountered with the SIR deployment. No additional workflow template is needed.
4. Open Git Bash in the extracted folder and run:

```bash
git init -b main
git add .
git commit -m "Add competition teaching app"
git remote add origin https://github.com/loeckete/competitionShinyApp.git
git push -u origin main
```

Change the URL if your repository has another name. These commands are for an empty repository; do not force-push over existing work. For a repository already initialized on GitHub, clone it, copy these files into its root, then commit and push.

5. Watch **Actions → Test and deploy competition app**. R model and Shiny interaction tests run before export and deployment.
6. Use the link from the successful deployment. With the example name, it will be `https://loeckete.github.io/competitionShinyApp/`. This is a proposed URL, **not an app already published by this package**.
7. Check the browser behavior below before class. Have students open the URL early to download the R runtime.

If Pages configuration fails, verify step 3 and rerun all jobs. If no workflow appears, ensure `.github/workflows/deploy.yml` is at repository root, and the branch is `main`.

## Model and controls

This original implementation follows the carrying-capacity formulation of the EcoEvoApps two-species Lotka–Volterra competition model:

```
dN1/dt = r1*N1 * [1 - (N1 + alpha*N2)/K1]
dN2/dt = r2*N2 * [1 - (N2 + beta*N1)/K2]
```

alpha measures species 2's effect on species 1; beta measures species 1's effect on species 2. **beta is not the SIR transmission coefficient.** Competition coefficients are relative to within-species effects. Time units are arbitrary and rates are per time unit. Resources are implicit; this is not a mechanistic model of resource concentrations.

Use positive growth rates, nonnegative competition coefficients, and a constant environment. N1 and N2 may be zero, but an absent species cannot reappear without immigration (not modeled). If either competitive coefficient is zero, its corresponding direct effect is neutral rather than negative. Counts are continuous, and tiny populations are not artificially rounded to extinction.

The default model is the formulation in the Competition SoftChalk lesson, not EcoEvoApps' alternative absolute-coefficient parameterization. No EcoEvoApps source code or images were copied; attribution is in the app.

## Classroom presets and answer guidance

All presets use r1=r2=0.2, K1=K2=100, and duration 150.

| Preset | N1 | N2 | alpha | beta | Expected long-term outcome |
|---|---:|---:|---:|---:|---|
| A | 20 | 20 | 0.5 | 0.5 | Stable coexistence, approximately 66.67 each |
| B | 20 | 20 | 0.5 | 1.5 | Species 1 excludes species 2 |
| C | 20 | 20 | 1.5 | 0.5 | Species 2 excludes species 1 |
| D | 30 | 20 | 1.5 | 1.5 | Species 1 wins |
| E | 20 | 30 | 1.5 | 1.5 | Species 2 wins |

Student-facing preset names do not announce the outcome. **Load settings** changes controls but does not run the model. Students predict before pressing **Run simulation**. Predictions remain in the browser session; the app does not collect or send responses to Top Hat.

**20-minute sequence for November 3:**

- 0–5: Load A, predict, and run. Show the independent no-competitor curves. Ask why coexistence can still represent a negative/negative interaction.
- 5–10: Keep A and load B. Explain the meaning of the coefficients before running. Can a winner still be negatively affected? C reverses the asymmetry.
- 10–16: Compare D and E, keeping D before running E. Only initial abundances change. Ask students to explain what changed and what stayed fixed.
- 16–20: Individual Top Hat explanation and feedback. Return to the interaction table: signs describe effects relative to absence, not simply who wins.

Extension: In D, set N1=N2=20. The perfectly symmetric trajectory approaches (40,40), an unstable interior equilibrium, despite the priority-effects regime. Then change N1 to 21 and extend time to 500. Numerical symmetry is not evidence of stable coexistence. Ask whether both species can recover when rare.

## Graphs and interpretation

- Population curves: solid lines for the completed run, dashed for one kept run, dotted for independent logistic growth without competitors (optional).
- Phase plane: current trajectory, kept trajectory, and the **current run's** zero-growth isoclines. The kept run does not contribute isoclines. Arrows are normalized directions, not growth-speed magnitudes. The phase-plane range focuses on carrying capacities and observed abundances, so distant intercepts can fall outside the plot.
- Reveal panel: outcome classification from invasion growth, not merely the populations at the last displayed time.
- CSV: current trajectory, no-competitor trajectories, and the complete run settings. Controls edited without Run do not change the exported run.

Species 1 invasion growth has sign `1-alpha*K2/K1`; species 2 has sign `1-beta*K1/K2`. Both positive means stable coexistence, both negative means priority effects. One positive and the other negative gives exclusion. Zero invasion growth is identified as a boundary case. Coincident isoclines have a continuum of equilibria. These classifications assume both species are initially present. For priority effects, exact starts on the separating trajectory are exceptional.

A winning species can still be harmed by competition: compare its growth with its own no-competitor curve. Do not infer effects from a comparison of the two species' abundances alone.

## Run locally

From the repository root:

```r
install.packages(c("shiny", "shinylive"))
shiny::runApp("app")
```

Test and export:

```bash
Rscript tests/model-tests.R
Rscript tests/shiny-tests.R
```

```r
shinylive::export("app", "site")
httpuv::runStaticServer("site")
```

Only `app/` is exported; instructor guidance and tests stay outside the student app. The solver uses base-R RK4 with rate-bounded substeps. Extreme parameter combinations exceeding the calculation budget produce an informative message rather than freezing the browser.

## Verification and browser checklist

Completed locally: five high-accuracy SciPy reference scenarios; a Python reproduction of the RK4 scheme agreed within 4e-9 individuals. R is unavailable in the build environment, so **native R tests, Shiny execution, Shinylive export, and browser behavior have not been executed here**. GitHub Actions runs the included R tests before publishing. The tests cover independent references, logistic limits, absent species, equilibrium classification, invalid inputs, step refinement, and Run/Keep/error behavior.

After the workflow succeeds:

1. Open the deployed site in a fresh browser and run all five presets.
2. Confirm A approaches 66.67 each and D/E reverse the winner.
3. Show no-competitor curves and inspect phase-plane arrows/isoclines.
4. Keep a run, change controls, and verify the old plots remain until Run.
5. Set alpha=beta=0; confirm independent logistic growth and vertical/horizontal isoclines. Set one starting abundance to zero and check it stays absent.
6. Download CSV and verify its settings. Enter K1=0 and confirm validation preserves the prior successful plots.
7. Check on a phone and a student laptop on classroom Wi-Fi.

There is no shared R-server session bottleneck, but simultaneous first-load downloads still depend on network capacity, hosting, and browser/device performance. This app has not been tested with 400 simultaneous clients. Refreshing clears session predictions and comparisons.

## Sources

- EcoEvoApps, Kenji Hayashi, Lotka–Volterra competition app: https://ecoevoapps.gitlab.io/apps/lotka-volterra-competition/
- Model equations and documentation: https://ecoevoapps.gitlab.io/docs/articles/lotka-volterra-competition.html
- Shinylive export and static hosting: https://posit-dev.github.io/r-shinylive/
