
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

        load = network.loads_t.p.sum().sum()

        df.at[(policy,parameter),"mp"] = network.buses_t.marginal_price[elec_buses].mean().mean()
        df.at[(policy,parameter),"dual"] = network.penetration_dual
        df.at[(policy,parameter),"objective"] = network.objective

        if network.global_constraints.empty:
            df.at[(policy,parameter),"co2_shadow"] = 0.
        else:
            df.at[(policy,parameter),"co2_shadow"] = network.global_constraints.at["CO2Limit","mu"]

        for techs in [["solar"],["wind"],["solar","wind"]]:
            gens = network.generators.index[network.generators.carrier.isin(techs)]

            gen_by_bus = network.generators_t.p[gens].groupby(network.generators.bus,axis=1).sum()

            mv_by_bus = (gen_by_bus*network.buses_t.marginal_price[elec_buses]).sum()/gen_by_bus.sum()

            tech_name = "-".join(techs)
            df.at[(policy,parameter),tech_name+"-mv"] = mv_by_bus.mean()
            df.at[(policy,parameter),tech_name+"-rmv"] = df.at[(policy,parameter),tech_name+"-mv"]/df.at[(policy,parameter),"mp"]

            df.at[(policy,parameter),tech_name+"-penetration"] = network.generators_t.p[gens].sum().sum()/load

            df.at[(policy,parameter),tech_name+"-mc"] = (network.generators.loc[gens,"capital_cost"]*network.generators.loc[gens,"p_nom_opt"]).sum()/network.generators_t.p[gens].multiply(network.snapshot_weightings,axis=0).sum().sum()

df.to_csv(snakemake.output["summary"])
