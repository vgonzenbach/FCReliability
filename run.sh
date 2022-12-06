cd $(dirname $0)

# install packages
Rscript -e "renv::restore()"
python -m venv .venv \
    && source .venv/bin/activate \
    && pip install -r requirements.txt 

# render notebook -> save plots in figures/
Rscript -e "rmarkdown::render('notebooks/viz.Rmd')"

# join plots
 .venv/bin/python3 scripts/join_plots.py