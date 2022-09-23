# Parallel and fixed mutations in SARS-CoV-2 evolution

## Download code and denpendency

Download the repository from GitHub. Raw data (not included in the repo due to file size and GISAID terms of use) can be acquired on request.

This project uses [Bioconda](https://bioconda.github.io/) to manage package dependency.

```bash
conda env update
conda activate ncov_paraFix
```

## Run the analysis

To prepare the data for `phangorn` (for homoplasyFinder) and `hyphy`, run the following commands.

```bash
python scripts/prep_homoplasyFinder.py
python scripts/prep_hyphy.py
```

To run the analysis, use the following commands.

```bash
Rscript scripts/run_homoplasyFinder.R # homoplasy
bash scripts/run_hyphy_meme.sh # episodic positive selection
bash scripts/run_hyphy_fubar.sh # pervasive  positive selection
Rscript scripts/run_sitePath.R # paraFix sites
```

To prepare for plotting the spatial-temporal distribution of mutated sites, run the following.

```bash
python scripts/prep_mutation_num.py
```

## Visualize the result

Launch the jupyter notebook from the root directory of the project

```bash
jupyter-notebook
```

Open `parse_hyphy.ipynb` and run the notebook to parse the output of `hyphy`.

Open `plot_result_comparison.ipynb` and run the notebook to visualize the comparison between paraFix and episodic positive selection sites.

Open `plot_mutation_num.ipynb` and run the notebook to visualize the spatial-temporal distribution of mutated sites
