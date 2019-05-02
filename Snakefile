configfile: "config.yaml"

wildcard_constraints:
    policy="[\-a-z0-9\.]+",
    parameter="[0-9]*",


rule solve_network:
    input:
    output:
        network=config['results_dir'] + "/" + config['run'] + "/networks/{policy}_{parameter}.nc"
#    shadow: "shallow"
    log:
        solver="logs/{policy}_{parameter}_solver.log",
        python="logs/{policy}_{parameter}_python.log",
        memory="logs/{policy}_{parameter}_memory.log"
#    benchmark: "benchmarks/solve_network/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}"
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


rule plot_summary:
    input:
        costs=config['results_dir'] + '/' + config['run'] + '/csvs/costs.csv',
        energy=config['results_dir'] + '/' + config['run'] + '/csvs/energy.csv'
    output:
        costs=config['results_dir'] + '/' + config['run'] + '/graphs/costs.pdf',
        energy=config['results_dir'] + '/' + config['run'] + '/graphs/energy.pdf'
    threads: 2
    resources: mem_mb=10000
    script:
        'scripts/plot_results.py'
