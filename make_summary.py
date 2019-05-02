
import pypsa

import numpy as np

import pandas as pd

config = snakemake.config

index = pd.MultiIndex.from_product([config["scenario"]["policy"],config["scenario"]["parameter"]],names=["policy","parameter"])

df = pd.DataFrame(index=index,dtype=float)

for policy in config["scenario"]["policy"]:
    for parameter in config["scenario"]["parameter"]:
        network_name = config['results_dir'] + "/" + config['run'] + "/networks/{policy}_{parameter}.nc".format(policy=policy,parameter=parameter)
        network = pypsa.Network(network_name)

        elec_buses = network.buses.index[network.buses.carrier == "AC"]

        df.at[(policy,parameter),"mp"] = network.buses_t.marginal_price[elec_buses].mean().mean()
        df.at[(policy,parameter),"dual"] = network.penetration_dual
        df.at[(policy,parameter),"objective"] = network.objective

        for tech in ["solar","wind"]:
            gens = network.generators.index[network.generators.carrier == tech]

            gen_by_bus = network.generators_t.p[gens].groupby(network.generators.bus,axis=1).sum()

            mv_by_bus = (gen_by_bus*network.buses_t.marginal_price[elec_buses]).sum()/gen_by_bus.sum()

            df.at[(policy,parameter),tech+"-mv"] = mv_by_bus.mean()
            df.at[(policy,parameter),tech+"-rmv"] = df.at[(policy,parameter),tech+"-mv"]/df.at[(policy,parameter),"mp"]
        #df[tech][penetration]["rmv"] = df[tech][penetration]["mv"]/df[tech][penetration]["mp"]

        #ignore variable costs here
        #df[tech][penetration]["mc"] = (network.generators.loc[gens,"capital_cost"]*network.generators.loc[gens,"p_nom_opt"]/network.generators_t.p[gens].multiply(network.snapshot_weightings,axis=0).sum()).rename(network.generators.bus)
        #mc = (network.generators.capital_cost*network.generators.p_nom_opt + network.generators_t.p.multiply(network.generators.marginal_cost,axis=1).multiply(network.snapshot_weightings,axis=0).sum())/(network.generators_t.p.multiply(network.snapshot_weightings,axis=0)).sum()

df.to_csv(snakemake.output["summary"])
