# Within-sample negativity-asymmetry check

WIP exploratory re-analysis of Paolini et al. (2024), asking whether negative-contact associations are larger than positive-contact associations **within the same sample and coded outcome**.

Read the rendered report: `docs/index.html` (published via GitHub Pages).

To reproduce locally:

```bash
curl -L -o data/F2.xlsx 'https://osf.io/download/2vctu/'
npm install
node scripts/analyse.js
quarto render report.qmd --output-dir docs
cp docs/report.html docs/index.html
```

The source data remain on the authors' [OSF project](https://osf.io/38rpj/); they are not redistributed here.
