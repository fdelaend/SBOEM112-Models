# Shiny apps for SBOEM112 (part: linear statistical modeling)

| App | Directory | Status |
|---|---|---|
| GLM simulator | `glm/` | v1 |
| LMM simulator | `lmm/` | planned |

## Running locally

```r
install.packages(c("shiny", "ggplot2", "DT"))
shiny::runApp("shiny/glm")
```

## Publishing to Posit Connect Cloud

Connect Cloud pulls from a GitHub repo. It needs two things in the app
directory: the primary file (`app.R`) and a `manifest.json`.

**1. Generate the manifest** — must be done with `shiny/glm` as the working
directory, not the repo root:

```r
install.packages("rsconnect")
setwd("shiny/glm")
rsconnect::writeManifest()
setwd("../..")
```

Commit the resulting `manifest.json`. Regenerate it whenever the package
dependencies change (not needed for ordinary edits to `app.R`).

**2. Publish** — on connect.posit.cloud, choose Shiny (R), pick the repo and
branch, and set the primary file to `shiny/glm/app.R`. Automatic republishing
is on by default, so every push to the connected branch redeploys.

### Which repo?

The free Connect Cloud tier can only pull from **public** repos; private repos
need a paid plan. `SBOEM112` currently tracks copyrighted PDFs in `resources/`,
so it should not be made public as it stands. Options:

- **Separate public repo** for the apps (recommended on the free tier). Keeps
  the course repo private, and stops every slide commit from triggering a
  redeploy. Both apps live there as sibling directories.
- **Keep it here**, if you either have a paid plan, or first remove
  `resources/` from the repo *and its history*.

Either way, add the student evaluation PDF to `.gitignore` — it is currently
untracked but not ignored, so a `git add .` would commit it.

## Extending to LMMs

`app.R` is deliberately split into a simulation engine (section 1) and the
Shiny layer (sections 2–3). The engine knows nothing about Shiny.

- `DISTRIBUTIONS` holds everything distribution-specific — link, inverse link,
  variance function, RNG — in the same order as the tables in the slides. A new
  `D` is one new entry.
- For LMMs the only structural change is in `simulate_glm()`: build a `Z` matrix
  of group indicators, draw `u ~ N(0, Sigma_u)`, and compute
  `eta <- X %*% beta + Z %*% u`. The plot and table code needs no changes beyond
  showing the group and `u_i` columns.
