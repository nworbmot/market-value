configfile: "config.yaml"

wildcard_constraints:
    policy="[\-a-zA-Z0-9\.]+",
    parameter="[0-9]*",
    assumptions="[\-a-zA-Z0-9\.]+"

rule plot_summary:
    input:
        summary=config['results_dir'] + "/" + config['run'] + "/csvs/summary.csv"
    output: config['results_dir'] + "/" + config['run'] + "/graphs/rmvs-res_" + config["scenario"]["assumptions"][0] + ".pdf"
    threads: 2
    resources: mem_mb=2000
    script:
        'plot_summary.py'

rule solve_network:
    input:
    output:
        network=config['results_dir'] + "/" + config['run'] + "/networks/{policy}_{parameter}_{assumptions}.nc"
#    shadow: "shallow"
    log:
        solver="logs/{policy}_{parameter}_{assumptions}_solver.log",
        python="logs/{policy}_{parameter}_{assumptions}_python.log",
        memory="logs/{policy}_{parameter}_{assumptions}_memory.log"
#    benchmark: "benchmarks/solve_network/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}"
    threads: 4
    resources: mem=6000
    script: "solve_network.py"


rule make_summary:
    input:
        expand(config['results_dir'] + "/" + config['run'] + "/networks/{policy}_{parameter}_{assumptions}.nc",
               **config['scenario'])
    output:
        summary=config['results_dir'] + "/" + config['run'] + "/csvs/summary.csv"
    threads: 2
    resources: mem_mb=2000
    script: 'make_summary.py'
