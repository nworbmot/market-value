configfile: "config.yaml"

wildcard_constraints:
    policy="[\-a-zA-Z0-9\.]+",
    parameter="[0-9]*"

rule plot_summary:
    input:
        summary=config['results_dir'] + "/" + config['run'] + "/csvs/summary.csv"
    output: 'paper_graphics/' + config['run'] + "/pen-compare-" + config["run"] + ".pdf"
    threads: 2
    resources: mem_mb=2000
    script:
        'plot_summary.py'

rule solve_network:
    input:
        config=config['results_dir'] + '/' + config['run'] + '/configs/config.yaml'
    output:
        network=config['results_dir'] + "/" + config['run'] + "/networks/{policy}_{parameter}.nc"
    log:
        solver="logs/{policy}_{parameter}_solver.log",
        python="logs/{policy}_{parameter}_python.log",
        memory="logs/{policy}_{parameter}_memory.log"
    threads: 4
    resources: mem=6000
    script: "solve_network.py"


rule make_summary:
    input:
        expand(config['results_dir'] + "/" + config['run'] + "/networks/{policy}_{parameter}.nc",
               **config['scenario'])
    output:
        summary=config['results_dir'] + "/" + config['run'] + "/csvs/summary.csv"
    threads: 2
    resources: mem_mb=2000
    script: 'make_summary.py'

rule copy_config:
    input:
    output:
        config=config['results_dir'] + '/' + config['run'] + '/configs/config.yaml'
    threads: 1
    resources: mem_mb=1000
    script:
        'copy_config.py'
