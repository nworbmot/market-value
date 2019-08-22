
from shutil import copy

files = ["config.yaml",
         "Snakefile",
         "solve_network.py"]

for f in files:
    copy(f,snakemake.config['results_dir'] + '/' + snakemake.config['run'] + '/configs/')
