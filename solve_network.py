## Reproduce Lion Hirth market value factors
#
#This script reproduces the EMMA model presented in the paper
#
#Lion Hirth, "The market value of variable renewables: The effect of solar wind power variability on their relative price", Energy Economics, 2013, https://doi.org/10.1016/j.eneco.2013.02.004
#
#using the PyPSA modelling framework, then makes some changes to the model setup.
#
### Library requirements
#
#To run this script you need the following Python libraries:
#
#ipython, pandas, pypsa, numpy, matplotlib, pyomo
#
#which are available from [PyPI](https://pypi.org/) with pip.
#
### Data requirements
#
#The datatables (saved as Excel .xls files) of the EMMA model is available as supplmentary material, downloadable from https://doi.org/10.1016/j.eneco.2013.02.004 as a .zip file. Download the file 1-s2.0-S0140988313000285-mmc2.zip to the directory emma_folder (adjust this path as you see fit) and unzip it.
#
#In addition you need the file assumptions-mv.csv for the storage assumptions.

import pypsa

import pandas as pd
import logging
logger = logging.getLogger(__name__)


# Suppress logging of the slack bus choices
pypsa.pf.logger.setLevel(logging.WARNING)

from vresutils.benchmark import memory_logger


idx = pd.IndexSlice

from pyomo.environ import Constraint

import numpy as np

# locations of data
emma_folder = 'emma/'
assumptions_file = 'assumptions-mv.csv'

### Read in data

wind_pu = pd.read_excel(emma_folder + "data_ts.xls",
              sheet_name="wind",
              index_col=0,
              skiprows=5,
              #usecols=list(range(43)), for some reason this is broken in pandas 0.24, had to edit file
              header=[0,1],
              convert_float=False)

solar_pu = pd.read_excel(emma_folder + "data_ts.xls",
              sheet_name="sola",
              index_col=0,
              skiprows=5,
              #usecols=42,
              header=[0,1],
              convert_float=False)

load = pd.read_excel(emma_folder + "data_ts.xls",
              sheet_name="loa",
              index_col=0,
              skiprows=5,
              #usecols=30,
              header=[0,1],
             convert_float=False)

ntcs = pd.read_excel(emma_folder + "data.xls",
                     sheet_name="NTC",
                     index_col=0,
                     skiprows=1,
                     usecols=list(range(1,32)),
                     convert_float=False)

costs = pd.read_excel(emma_folder + "data.xls",
                     sheet_name="cost",
                     index_col=0,
                     skiprows=1,
                     usecols=list(range(0,7)),
                     convert_float=False,
                      skipfooter=8)
costs = costs.rename(columns={"Unnamed: 4" : "fuel"}).fillna(0.)

costs

### Major settings (see config.yaml)

year_start=2010

year_end = 2010

Nyears = year_end - year_start + 1

#years
lifetime = 25

#per unit
discount_rate=0.07

#EUR/tCO2
co2_price = 20


convs = ["nucl","coal","lign","OCGT","CCGT"]

### Required functions

def annuity(lifetime,rate):
    if rate == 0.:
        return 1/lifetime
    else:
        return rate/(1. - 1. / (1. + rate)**lifetime)

assumptions = costs.copy()

#1e3 is kW to MW
assumptions["fixed"] = 1e3*Nyears*(annuity(lifetime,discount_rate)*assumptions["invest"] + assumptions["qfixcost"])


assumptions['variable'] = assumptions['varcost'] + (assumptions['fuel'] + co2_price*assumptions["co2int"])/ assumptions['eff']

usd_to_eur=1/1.2

assumptions_year=2030

assumptions_prev = pd.read_csv(assumptions_file,index_col=list(range(3))).sort_index()
assumptions_prev.loc[assumptions_prev.unit.str.contains("/kW"),"value"]*=1e3
assumptions_prev.loc[assumptions_prev.unit.str.contains("USD"),"value"]*=usd_to_eur

assumptions_prev = assumptions_prev.loc[idx[:,assumptions_year,:],"value"].unstack(level=2).groupby(level="technology").sum(min_count=1)

#fill defaults
assumptions_prev = assumptions_prev.fillna({"FOM" : assumptions_prev.at["default","FOM"],
                                      "discount rate" : discount_rate,
                                      "lifetime" : lifetime,
                                      "CO2 intensity" : 0,
                                      "VOM" : 0,
                                      "efficiency" : 1,
                                      "fuel" : 0,
                                      "investment" : 0})

#annualise investment costs, add FOM
assumptions_prev["fixed"] = [(annuity(v["lifetime"],v["discount rate"])+v["FOM"]/100.)*v["investment"]*Nyears for i,v in assumptions_prev.iterrows()]
st_techs = ["H2 CCGT","H2 electrolysis","H2 steel tank storage","H2 underground storage","battery inverter","battery storage"]

assumptions = assumptions.reindex(assumptions.index.append(pd.Index(st_techs)))

for attr in ["investment","lifetime","discount rate","FOM","fixed","efficiency"]:
    assumptions.loc[st_techs,attr] = assumptions_prev.loc[st_techs,attr]

print(assumptions)

cts = ["GER","FRA","BEL","NLD","POL"]

def prepare_network(allow_transmission_expansion=False):
    #technologies to remove

    network = pypsa.Network()

    full_snapshots = pd.date_range("{}-01-01".format(year_start),"{}-12-31 23:00".format(year_end),
                              freq="1H")

    snapshots = pd.date_range("{}-01-01".format(year_start),"{}-12-31 23:00".format(year_end),
                              freq=str(frequency)+"H")

    network.set_snapshots(snapshots)

    network.snapshot_weightings = pd.Series(float(frequency),index=network.snapshots)

    network.madd("Carrier",
                 convs,
                 co2_emissions=assumptions.loc[convs,"co2int"])

    for ct in cts:

        network.add("Bus",ct)
        network.add("Load",ct,
                bus=ct,
                p_set=pd.Series(index=full_snapshots,data=load[ct,year_start].values))

        network.add("Generator",ct+" solar",
                bus=ct,
                p_max_pu = pd.Series(index=full_snapshots,data=solar_pu[ct,year_start].values),
                p_nom_extendable = True,
                carrier="solar",
                marginal_cost = assumptions.at["sola","variable"],
                capital_cost = assumptions.at['sola','fixed'])

        network.add("Generator",ct+" wind",
                bus=ct,
                p_max_pu = pd.Series(index=full_snapshots,data=wind_pu[ct,year_start].values),
                carrier="wind",
                p_nom_extendable = True,
                marginal_cost = assumptions.at["sola","variable"],
                capital_cost = assumptions.at['wind','fixed'])

        for conv in convs:
            network.add("Generator",ct+" " + conv,
                        bus=ct,
                        p_nom_extendable = True,
                        carrier=conv,
                        efficiency=assumptions.at[conv,'eff'],
                        marginal_cost = assumptions.at[conv,'variable'],
                        capital_cost = assumptions.at[conv,'fixed'])
    for ct1 in cts:
        for ct2 in cts:
            if not pd.isnull(ntcs.at[ct1,ct2]):
                print("adding link",ct1,ct2,ntcs.at[ct1,ct2])
                network.add("Link","{}->{}".format(ct1,ct2),
                            bus0=ct1,
                            bus1=ct2,
                            p_nom_extendable=allow_transmission_expansion,
                            p_nom=ntcs.at[ct1,ct2])
    for ct in cts:
        if add_battery:
            network.add("Bus",ct + " battery",
                        carrier="battery")

            network.add("Store",ct + " battery storage",
                bus = ct + " battery",
                e_nom_extendable = True,
                e_cyclic=True,
                capital_cost=assumptions.at['battery storage','fixed'])

            network.add("Link",ct + " battery charge",
                bus0 = ct,
                bus1 = ct + " battery",
                efficiency = assumptions.at['battery inverter','efficiency'],
                p_nom_extendable = True,
                capital_cost=assumptions.at['battery inverter','fixed'])

            network.add("Link",ct + " battery discharge",
                bus0 = ct + " battery",
                bus1 = ct,
                p_nom_extendable = True,
                efficiency = assumptions.at['battery inverter','efficiency'])


        if add_hydrogen:

            network.add("Bus",
                     ct + " H2",
                     carrier="H2")

            network.add("Link",
                    ct + " H2 electrolysis",
                    bus1=ct + " H2",
                    bus0=ct,
                    p_nom_extendable=True,
                    efficiency=assumptions.at["H2 electrolysis","efficiency"],
                    capital_cost=assumptions.at["H2 electrolysis","fixed"])

            network.add("Link",
                     ct + " H2 to power",
                     bus0=ct + " H2",
                     bus1=ct,
                     p_nom_extendable=True,
                     efficiency=assumptions.at["H2 CCGT","efficiency"],
                     capital_cost=assumptions.at["H2 CCGT","fixed"]*assumptions.at["H2 CCGT","efficiency"])  #NB: fixed cost is per MWel

            network.add("Store",
                     ct + " H2 storage",
                     bus=ct + " H2",
                     e_nom_extendable=True,
                     e_cyclic=True,
                     capital_cost=assumptions.at["H2 underground storage","fixed"])

    return network

def solve_network(network,penetration,load,techs,emissions):

    network.add("GlobalConstraint", "CO2Limit",
                carrier_attribute="co2_emissions", sense="<=",
                constant=emissions*load)


    def extra_functionality(network,snapshots):

        if add_battery:
            def battery(model,ct):
                return model.link_p_nom[ct + " battery charge"] == model.link_p_nom[ct + " battery discharge"]*network.links.at[ct + " battery charge","efficiency"]

            network.model.battery = Constraint(cts,rule=battery)

        if penetration is not None:
            network.model.penetration = Constraint(expr=sum([network.model.generator_p[gen, sn]*network.snapshot_weightings.at[sn] for gen in network.generators.index if network.generators.at[gen,"carrier"] in techs for sn in snapshots]) == penetration*load)

    if solver_name == "gurobi":
        solver_options = {"threads" : 4,
                          "method" : 2,
                          "crossover" : 0,
                          "BarConvTol": 1.e-5,
                          "FeasibilityTol": 1.e-6 }
    else:
        solver_options = {}


    def extra_postprocessing(n, snapshots, duals):
        if penetration is None:
            n.penetration_dual =  0.
        else:
            index = list(n.model.penetration.keys())
            cdata = pd.Series(list(n.model.penetration.values()),
                              index=index)
            n.penetration_dual =  -cdata.map(duals).sum()
        print("penetration dual:",n.penetration_dual)



    network.consistency_check()

    network.lopf(solver_name=solver_name,
                 solver_logfile=snakemake.log.solver,
                 solver_options=solver_options,
                 extra_functionality=extra_functionality,
                 extra_postprocessing=extra_postprocessing)



if __name__ == "__main__":
    # Detect running outside of snakemake and mock snakemake for testing
    if 'snakemake' not in globals():
        from vresutils.snakemake import MockSnakemake, Dict
        snakemake = MockSnakemake(
            wildcards=dict(assumptions='wind1100-sola750-nucl4000', policy="OCGT",parameter="101"),
            #input=dict(network="networks/{network}_s{simpl}_{clusters}_lv{lv}_{opts}.nc"),
            #output=["results/networks/s{simpl}_{clusters}_lv{lv}_{opts}-test.nc"],
            #log=dict(gurobi="logs/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_gurobi-test.log",
            #         python="logs/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_python-test.log")
        )
        import yaml
        with open('config.yaml') as f:
            snakemake.config = yaml.load(f)

    #change to cbc or glpk for open-source solvers
    solver_name = snakemake.config["solver"]["name"]

    #1 is hourly, 3 is 3-hourly
    frequency = snakemake.config["frequency"]

    if "storage" in snakemake.wildcards.policy:
        add_hydrogen = True
        add_battery = True
    else:
        add_hydrogen = False
        add_battery = False


    techs_to_remove=[]
    allow_transmission_expansion=False

    changed_assumptions = snakemake.wildcards.assumptions

    for opt in changed_assumptions.split("-"):
        for tech in assumptions.index:
            if tech in opt:
                if opt[len(tech):] == "None":
                    print(tech,"None")
                    techs_to_remove.append(tech)
                else:
                    print(tech,float(opt[len(tech):]))
                    assumptions.at[tech,"invest"] = float(opt[len(tech):])
                    assumptions.at[tech,"fixed"] = 1e3*Nyears*(annuity(lifetime,discount_rate)*assumptions.at[tech,"invest"] + assumptions.at[tech,"qfixcost"])
        if opt == "trans":
            allow_transmission_expansion=True

            ###
    #assumptions["fixed"] = 1e3*Nyears*(annuity(lifetime,discount_rate)*assumptions["invest"] + assumptions["qfixcost"])
    #assumptions_prev["fixed"] = [(annuity(v["lifetime"],v["discount rate"])+v["FOM"]/100.)*v["investment"]*Nyears for i,v in assumptions_prev.iterrows()]


    techs=[]
    if "co2" in snakemake.wildcards.policy:
        emissions = float(snakemake.wildcards.parameter)/1e3 #tCO2/Mwh_el on average
        penetration = None
        #remove nuclear to force in wind and solar
        for tech in techs_to_remove:
            convs.remove(tech)
    else:
        for tech in ["solar","wind","OCGT","CCGT"]:
            if tech in snakemake.wildcards.policy:
                techs.append(tech)
        penetration = float(snakemake.wildcards.parameter)/1e3
        emissions = 2.


    print("solving network for policy {} and penetration {} and emissions {} and assumptions {}".format(snakemake.wildcards.policy,penetration,emissions,changed_assumptions))


    print(assumptions)

    with memory_logger(filename=getattr(snakemake.log, 'memory', None), interval=30.) as mem:
        network = prepare_network(allow_transmission_expansion=allow_transmission_expansion)

        total_load = (network.loads_t.p_set.multiply(network.snapshot_weightings,axis=0)).sum().sum()

        solve_network(network,penetration,total_load,techs,emissions)

        network.export_to_netcdf(snakemake.output[0])

    logger.info("Maximum memory usage: {}".format(mem.mem_usage))
