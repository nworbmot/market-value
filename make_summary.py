
import pypsa

import numpy as np

import pandas as pd

config = snakemake.config

index = pd.MultiIndex.from_product([config["scenario"]["policy"],config["scenario"]["parameter"]],names=["policy","parameter"])

df = pd.DataFrame(index=index,dtype=float)

for i in df.index:
    network_name = config['results_dir'] + "/" + config['run'] + "/networks/{policy}_{parameter}.nc".format(policy=i[0],parameter=i[1])
    network = pypsa.Network(network_name)

    elec_buses = network.buses.index[network.buses.carrier == "AC"]

    load_by_bus = network.loads_t.p_set.sum()

    df.at[i,"load"] = load_by_bus.sum()

    prices = {"regular" : network.buses_t.marginal_price[elec_buses],
              "zeroed" : network.buses_t.marginal_price[elec_buses].copy()}

    prices["zeroed"][prices["zeroed"] < 0.] = 0.


    #load-weighted
    df.at[i,"mp"] = (prices["regular"]*network.loads_t.p_set).sum().sum()/df.at[i,"load"]
    df.at[i,"mp-zeroed"] = (prices["zeroed"]*network.loads_t.p_set).sum().sum()/df.at[i,"load"]
    df.at[i,"dual"] = network.penetration_dual
    df.at[i,"objective"] = network.objective

    if network.global_constraints.empty:
        df.at[i,"co2_shadow"] = 0.
    else:
        df.at[i,"co2_shadow"] = network.global_constraints.at["CO2Limit","mu"]

    for techs in [["solar"],["wind"],["wind","solar"],["nucl"]]:
        gens = network.generators.index[network.generators.carrier.isin(techs)]

        #total generation
        gen_by_bus = network.generators_t.p[gens].groupby(network.generators.bus,axis=1).sum()

        mv_by_bus = (gen_by_bus*prices["regular"]).sum()/gen_by_bus.sum()
        mv_by_bus_zeroed = (gen_by_bus*prices["zeroed"]).sum()/gen_by_bus.sum()

        #LCOE
        mc_by_bus = (network.generators.loc[gens,"capital_cost"]*network.generators.loc[gens,"p_nom_opt"] +  network.generators_t["p"][gens].multiply(network.generators.loc[gens,"marginal_cost"]).multiply(network.snapshot_weightings,axis=0).sum()).groupby(network.generators.bus).sum()/gen_by_bus.sum()

        tech_name = "-".join(techs)

        # generation-weighted
        df.at[i,tech_name+"-mv"] = (mv_by_bus*gen_by_bus.sum()).sum()/gen_by_bus.sum().sum()
        df.at[i,tech_name+"-rmv"] = df.at[i,tech_name+"-mv"]/df.at[i,"mp"]
        df.at[i,tech_name+"-mv-zeroed"] = (mv_by_bus_zeroed*gen_by_bus.sum()).sum()/gen_by_bus.sum().sum()
        df.at[i,tech_name+"-rmv-zeroed"] = df.at[i,tech_name+"-mv-zeroed"]/df.at[i,"mp"]
        df.at[i,tech_name+"-mc"] = (mc_by_bus*gen_by_bus.sum()).sum()/gen_by_bus.sum().sum()

        # load-weighted
        df.at[i,tech_name+"-mv-load"] = (mv_by_bus*load_by_bus).sum()/load_by_bus.sum()
        df.at[i,tech_name+"-rmv-load"] = df.at[i,tech_name+"-mv-load"]/df.at[i,"mp"]
        df.at[i,tech_name+"-mv-load-zeroed"] = (mv_by_bus_zeroed*load_by_bus).sum()/load_by_bus.sum()
        df.at[i,tech_name+"-rmv-load-zeroed"] = df.at[i,tech_name+"-mv-load-zeroed"]/df.at[i,"mp"]
        df.at[i,tech_name+"-mc-load"] = (mc_by_bus*load_by_bus).sum()/load_by_bus.sum()

        df.at[i,tech_name+"-penetration"] = network.generators_t.p[gens].sum().sum()/df.at[i,"load"]

    for ln in ["generators","stores","links"]:
        ldf = getattr(network,ln)
        pnl = getattr(network,ln+"_t")
        ldf["carrier"] = ldf.index.str[4:]
        trans = ldf.index[ldf.index.str.contains("->")]
        ldf.loc[trans,"carrier"] = "AC"

        attr = "e" if ln == "stores" else "p"
        pattr = "p0" if ln == "links" else "p"

        cost = (ldf["capital_cost"]*ldf[attr + "_nom_opt"] + pnl[pattr].multiply(ldf["marginal_cost"]).multiply(network.snapshot_weightings,axis=0).sum()).groupby(ldf.carrier).sum()

        for carrier in cost.index:
            df.at[i,carrier + "-cost"] = cost.at[carrier]


df.to_csv(snakemake.output["summary"])
