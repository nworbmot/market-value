
import pypsa

import numpy as np

import pandas as pd

config = snakemake.config

index = pd.MultiIndex.from_product([config["scenario"]["policy"],config["scenario"]["parameter"],config["scenario"]["assumptions"]],names=["policy","parameter","asusmptions"])

df = pd.DataFrame(index=index,dtype=float)

for i in df.index:
    network_name = config['results_dir'] + "/" + config['run'] + "/networks/{policy}_{parameter}_{assumptions}.nc".format(policy=i[0],parameter=i[1],assumptions=i[2])
    network = pypsa.Network(network_name)

    elec_buses = network.buses.index[network.buses.carrier == "AC"]

    load = network.loads_t.p.sum().sum()

    #load-weighted
    df.at[i,"mp"] = (network.buses_t.marginal_price[elec_buses].mean()*network.loads_t.p_set.sum()).sum()/load
    df.at[i,"dual"] = network.penetration_dual
    df.at[i,"objective"] = network.objective

    if network.global_constraints.empty:
        df.at[i,"co2_shadow"] = 0.
    else:
        df.at[i,"co2_shadow"] = network.global_constraints.at["CO2Limit","mu"]

    for techs in [["solar"],["wind"],["solar","wind"],["nucl"]]:
        gens = network.generators.index[network.generators.carrier.isin(techs)]

        gen_by_bus = network.generators_t.p[gens].groupby(network.generators.bus,axis=1).sum()

        mv_by_bus = (gen_by_bus*network.buses_t.marginal_price[elec_buses]).sum()/gen_by_bus.sum()

        tech_name = "-".join(techs)

        # load-weighted
        df.at[i,tech_name+"-mv"] = (mv_by_bus*network.loads_t.p_set.sum()).sum()/load
        df.at[i,tech_name+"-rmv"] = df.at[i,tech_name+"-mv"]/df.at[i,"mp"]

        df.at[i,tech_name+"-penetration"] = network.generators_t.p[gens].sum().sum()/load

        df.at[i,tech_name+"-mc"] = (network.generators.loc[gens,"capital_cost"]*network.generators.loc[gens,"p_nom_opt"] +  network.generators_t["p"][gens].multiply(network.generators.loc[gens,"marginal_cost"]).multiply(network.snapshot_weightings,axis=0).sum()).sum()/network.generators_t.p[gens].multiply(network.snapshot_weightings,axis=0).sum().sum()

df.to_csv(snakemake.output["summary"])
