## Install GitHub-only R packages (ALRA, scType) into the active R library.
## Run once after: conda env create -f envs/env_r.yaml && conda activate env_r
install-r-github:
	Rscript -e "remotes::install_github('nalab-stanford/ALRA',              upgrade = 'never')"
	Rscript -e "remotes::install_github('IanevskiAleksandr/sc-type',        upgrade = 'never')"
