import pandas as pd

#allow plotting without Xwindows
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

import numpy as np

colors={"wind" : "b",
        "solar" : "y",
        "nucl" : "r",
        "wind-solar" : "g",
        "wind+solar" : "g",
        "CCGT" : "orange",
        "OCGT" : "wheat",
        "battery storage" : "gray",
        "coal" : "k",
        "lign" : "brown",
        "hydrogen storage" : "m",
        "shed" : "pink"}

def choose_color(s):
    if "system cost" in s:
        return "k"
    #force wind-solar first
    for tech in ["wind-solar","wind+solar"] + list(colors.keys()):
        if tech in s:
            return colors[tech]
    return "c"


def choose_style(s):
    if "MV" in s:
        return "-"
    if "market price" in s:
        return "-."
    if "LCOE" in s:
        return "--"
    if "penetration" in s:
        return "--"
    elif "FiP" in s or "CO$_2$ price" in s:
        return ":"
    if "system cost" in s:
        return ":"
    else:
        return "-"

ret_color = "crimson"
co2_color = "darkblue"


def plot_re_penetration(policy):

    if "nuclNone" in policy:
        suffix = ""
    else:
        suffix = "-load"
    #suffix = ""  #"" for gen-weighted, -load for load-weighted

    fig, ax = plt.subplots()
    fig.set_size_inches((4,3))

    for tech in ["solar","wind","wind-solar"]:

        p = policy.format(tech.replace("-",""))
        s = df.loc[p][tech+"-rmv"+suffix].copy()
        s.index = s.index/20*float(pen)
        s[ (s>1.3)] = np.nan
        s.plot(ax=ax,label=tech,color=colors[tech],linewidth=2)

        print(tech)
        for x in [0,15,30]:
            print(x,np.interp(x,s.index,s.values))

        ax.legend()

        ax.set_ylim([0,1.3])
        ax.set_xlim([0,30])
        ax.set_xlabel("penetration [%]")
        ax.set_ylabel("relative market value [per unit]")

        ax.grid()

        fig.tight_layout()

    fig.savefig("paper_graphics/{}/rmv-{}-{}{}.pdf".format(scenario,scenario,policy.format(""),suffix),transparent=True)


def ret_tech_fig(tech,policy):
    fig, ax = plt.subplots()
    fig.set_size_inches((4,3))

    plot_df = df.loc[policy].copy()

    plot_df.index = plot_df.index/20*float(policy[3:6])

    plot_df["dual"] *= -1

    plot_df.rename(columns={tech + "-mc" : tech.replace("-","+") + " LCOE",
                            "mp" : "market price",
                            tech + "-mv" : tech.replace("-","+") + " MV",
                            "dual" : tech.replace("-","+") + " FiP"},inplace=True)


    plot_df["system cost"] = plot_df[[c for c in plot_df.columns if "-cost" in c]].sum(axis=1)/plot_df["load"]
    sel = [tech.replace("-","+")+" MV",tech.replace("-","+") +" LCOE",tech.replace("-","+") + " FiP","market price"]#,"system cost"]


    plot_df[sel].plot(ax=ax,linewidth=2,color=choose_color(tech) if tech != "wind-solar" else ret_color,
                      style=[choose_style(s) for s in sel])
        #s[ (s>1.3)] = np.nan
        #s.plot(ax=ax,label=tech,color=colors[tech])

    ax.legend(prop={'size': 9},ncol=2)
    xlim = {"wind" : 65,
            "solar" : 34,
            "wind-solar" : 70,
            "nucl" : 100}
    ax.set_xlim([0,xlim[tech]])
    ax.set_ylim([0,ylim_comparison])
    ax.set_xlabel("{} penetration [%]".format(tech.replace("-","+")))
    ax.set_ylabel("energy price [€/MWh]")

    #ax.grid()

    fig.tight_layout()

    fig.savefig("paper_graphics/{}/mwh-{}-{}-{}.pdf".format(scenario,tech,scenario,policy),transparent=True)



def ret_tech_fig_clean(tech,policy):
    fig, ax = plt.subplots()
    fig.set_size_inches((4,3))

    plot_df = df.loc[policy].copy()

    plot_df.index = plot_df.index/20*float(policy[3:6])

    plot_df["dual"] *= -1

    plot_df.rename(columns={tech + "-mc" : tech.replace("-","+") + " LCOE",
                            "mp" : "market price",
                            tech + "-mv" : tech.replace("-","+") + " MV",
                            "dual" : tech.replace("-","+") + " FiP"},inplace=True)


    plot_df["system cost"] = plot_df[[c for c in plot_df.columns if "-cost" in c]].sum(axis=1)/plot_df["load"]
    sel = [tech.replace("-","+")+" MV",tech.replace("-","+") +" LCOE",tech.replace("-","+") + " FiP"]#,"system cost"]


    plot_df[sel].plot(ax=ax,linewidth=2,color=choose_color(tech) if tech != "wind-solar" else ret_color,
                      style=[choose_style(s) for s in sel])
        #s[ (s>1.3)] = np.nan
        #s.plot(ax=ax,label=tech,color=colors[tech])

    ax.legend(prop={'size': 9},ncol=2)
    xlim = {"wind" : 65,
            "solar" : 34,
            "wind-solar" : 70,
            "nucl" : 100}
    ax.set_xlim([0,xlim[tech]])
    ax.set_ylim([0,ylim_comparison])
    ax.set_xlabel("{} penetration [%]".format(tech.replace("-","+")))
    ax.set_ylabel("energy price [€/MWh]")

    #ax.grid()

    fig.tight_layout()

    fig.savefig("paper_graphics/{}/mwh-{}-{}-{}-clean.pdf".format(scenario,tech,scenario,policy),transparent=True)


def pen_compare():

    tech = "wind+solar"
    #tech = "wind"
    #tech = "solar"

    fig, ax = plt.subplots()
    fig.set_size_inches((4,3))
    for policy in ["pen{}{}-{}-nuclNone-lCCSNone".format(pen,tech.replace("+",""),assumptions),"pen100{}-{}-nuclNone-lCCSNone-trans-storage".format(tech.replace("+",""),assumptions)]:

        plot_df = df.loc[policy].copy()

        plot_df.index = plot_df.index/20*float(policy[3:6])

        policy_text = "with flexibility" if "trans-storage" in policy else "without flexibility"

        plot_df.rename(columns={"mp" : "market price "+policy_text,
                                tech.replace("+","-") + "-mv" : tech + " MV "+policy_text},inplace=True)


        sel = [tech+" MV "+policy_text,"market price "+policy_text]
        print(sel)
        plot_df[sel].plot(ax=ax,linewidth=2,color=ret_color,
                          style=["-","-."],
                          alpha=1. if "trans-storage" in policy else 0.4)

            #s[ (s>1.3)] = np.nan
            #s.plot(ax=ax,label=tech,color=colors[tech])


        ax.set_xlim([0,100])
        ax.set_xlabel(tech + " penetration [%]")
        ax.set_ylim([0,70])
        ax.set_ylabel("energy price [€/MWh]")
        ax.legend(loc="upper right",prop={'size': 9})

        fig.tight_layout()

        fig.savefig("paper_graphics/{}/pen-compare-{}.pdf".format(scenario,scenario),transparent=True)




def co2_as_em(policy):

    fig, ax = plt.subplots()
    fig.set_size_inches((4,3))


    plot_df = df.loc[policy].copy()

    plot_df.index = (plot_df.index/20)*float(policy[3:6])/100.


    plot_df["dual"] *= -1

    plot_df["wind-penetration"] *= 100
    plot_df["solar-penetration"] *= 100

    plot_df.rename(columns={"wind-mc" : "wind LCOE  [€/MWh]",
                            "mp" : "average market price [€/MWh]",
                            "wind-mv" : "wind MV = LCOE [€/MWh]",
                            "solar-mv" : "solar MV = LCOE [€/MWh]",
                            "wind-penetration" : "wind penetration [%]",
                            "solar-penetration" : "solar penetration [%]",
                            "co2_shadow" : "CO$_2$ price [€/tCO$_2$]"},inplace=True)

    sel = ["average market price [€/MWh]","wind MV = LCOE [€/MWh]","solar MV = LCOE [€/MWh]","CO$_2$ price [€/tCO$_2$]",
             "wind penetration [%]", "solar penetration [%]"
            ]
    plot_df[sel].plot(ax=ax,linewidth=2,color=[choose_color(s) for s in sel],
                      style=[choose_style(s) for s in sel])
        #s[ (s>1.3)] = np.nan
        #s.plot(ax=ax,label=tech,color=colors[tech])

    if "trans-storage" in policy:
        xlim = 0.
    else:
        xlim = 0.2
    ax.set_xlim([1.2,xlim])
    ax.set_xlabel("average emissions [tCO2/MWhel]")
    ax.set_ylim([0,200])
    #ax.set_ylabel("energy cost [€/MWh]")
    ax.legend(loc="upper left",prop={'size': 9})

    fig.tight_layout()

    fig.savefig("paper_graphics/{}/mwh-co2-{}-{}.pdf".format(scenario,scenario,policy),transparent=True)





def co2_as_pen(policy):
    fig, ax = plt.subplots()
    fig.set_size_inches((4,3))

    plot_df = df.loc[policy].copy()

    plot_df.index = plot_df["wind-solar-penetration"].values*100.
    plot_df.drop(plot_df.index[plot_df.index < 0.01],inplace=True)

    plot_df["dual"] *= -1

    plot_df["wind-penetration"] *= 100
    plot_df["solar-penetration"] *= 100

    plot_df.rename(columns={"wind-mc" : "wind LCOE  [€/MWh]",
                            "mp" : "average market price [€/MWh]",
                            "wind-mv" : "wind MV = LCOE [€/MWh]",
                            "solar-mv" : "solar MV = LCOE [€/MWh]",
                            "wind-penetration" : "wind penetration [%]",
                            "solar-penetration" : "solar penetration [%]",
                            "co2_shadow" : "CO$_2$ price [€/tCO$_2$]"},inplace=True)

    sel = ["average market price [€/MWh]","wind MV = LCOE [€/MWh]","solar MV = LCOE [€/MWh]","CO$_2$ price [€/tCO$_2$]",
             "wind penetration [%]", "solar penetration [%]"
            ]
    plot_df[sel].plot(ax=ax,linewidth=2,color=[choose_color(s) for s in sel],
                      style=[choose_style(s) for s in sel])
        #s[ (s>1.3)] = np.nan
        #s.plot(ax=ax,label=tech,color=colors[tech])

    if "trans-storage" in policy:
        xlim = 100.
    else:
        xlim = 70.

    ax.set_xlim([0,xlim])
    ax.set_xlabel("wind+solar penetration [%]")
    ax.set_ylim([0,200])
    #ax.set_ylabel("energy cost [€/MWh]")
    ax.legend(loc="upper left",prop={'size': 9})

    fig.tight_layout()

    fig.savefig("paper_graphics/{}/mwh-pen-{}-{}.pdf".format(scenario,scenario,policy),transparent=True)




def cot_flex_compare():
    fig, ax = plt.subplots()
    fig.set_size_inches((4,3))
    for policy in ["co2120-{}-nuclNone-lCCSNone".format(assumptions),"co2120-trans-storage-{}-nuclNone-lCCSNone".format(assumptions)]:

        plot_df = df.loc[policy].copy()

        plot_df.index = plot_df["wind-solar-penetration"].values*100.
        plot_df.drop(plot_df.index[plot_df.index < 0.01],inplace=True)

        plot_df["dual"] *= -1

        plot_df["wind-solar-penetration"] *= 100

        policy_text = "with flexibility" if "trans-storage" in policy else "without flexibility"

        plot_df.rename(columns={"mp" : "market price "+policy_text,
                                "wind-solar-mv" : "wind+solar MV "+policy_text},inplace=True)


        sel = ["wind+solar MV "+policy_text,"market price "+policy_text]

        plot_df[sel].plot(ax=ax,linewidth=2,color=co2_color,
                          style=["-","-."],
                          alpha=1. if "trans-storage" in policy else 0.4)

        ax.set_xlim([0,100])
        ax.set_xlabel("wind+solar penetration [%]")
        ax.set_ylim([0,140])
        ax.set_ylabel("energy price [€/MWh]")
        ax.legend(loc="lower right",prop={'size': 9})

        #ax.grid()

        fig.tight_layout()

        fig.savefig("paper_graphics/{}/co2-compare-{}.pdf".format(scenario,scenario),transparent=True)


def cot_flex_compare_clean():
    fig, ax = plt.subplots()
    fig.set_size_inches((4,3))
    for policy in ["co2120-{}-nuclNone-lCCSNone".format(assumptions),"co2120-trans-storage-{}-nuclNone-lCCSNone".format(assumptions)]:

        plot_df = df.loc[policy].copy()

        plot_df.index = plot_df["wind-solar-penetration"].values*100.
        plot_df.drop(plot_df.index[plot_df.index < 0.01],inplace=True)

        plot_df["dual"] *= -1

        plot_df["wind-solar-penetration"] *= 100

        policy_text = "with flexibility" if "trans-storage" in policy else "without flexibility"

        plot_df.rename(columns={"mp" : "market price "+policy_text,
                                "wind-solar-mv" : "wind+solar MV "+policy_text},inplace=True)


        sel = ["wind+solar MV "+policy_text]

        plot_df[sel].plot(ax=ax,linewidth=2,color=co2_color,
                          style=["-","-."],
                          alpha=1. if "trans-storage" in policy else 0.4)

        ax.set_xlim([0,100])
        ax.set_xlabel("wind+solar penetration [%]")
        ax.set_ylim([0,140])
        ax.set_ylabel("energy price [€/MWh]")
        ax.legend(loc="lower right",prop={'size': 9})

        #ax.grid()

        fig.tight_layout()

        fig.savefig("paper_graphics/{}/co2-compare-{}-clean.pdf".format(scenario,scenario),transparent=True)



def cot_pen_mu(policy):

    fig, ax = plt.subplots()
    fig.set_size_inches((4.5,3))

    plot_df = df.loc[policy].copy()

    plot_df.index = plot_df["wind-solar-penetration"].values*100.

    plot_df["dual"] *= -1

    plot_df["wind-solar-penetration"] *= 100

    plot_df.drop(plot_df.index[plot_df.index < 0.01],inplace=True)

    policy_text = "with flexibility" if "trans-storage" in policy else "without flexibility"

    plot_df.rename(columns={"mp" : "market price",
                            "wind-solar-mv" : "wind+solar MV = LCOE"},inplace=True)

    sel = ["wind+solar MV = LCOE","market price"]
    plot_df[sel].plot(ax=ax,linewidth=2,color=co2_color,
                      style=[choose_style(s) for s in sel],legend=False)
        #s[ (s>1.3)] = np.nan
        #s.plot(ax=ax,label=tech,color=colors[tech])


    ax2 = ax.twinx()

    plot_df["co2_shadow"].plot(ax=ax2,label="CO$_2$ price (right axis)",color=co2_color,style=":")

    ax.set_xlim([0,70])
    ax.set_xlabel("wind+solar penetration [%]")
    ax.set_ylim([0,ylim_comparison])
    ax.set_ylabel("energy price [€/MWh]")
    fig.legend(prop={'size': 9},
              loc='upper left', bbox_to_anchor=(0.15, 0.93))#loc="upper left",

    ax2.set_ylim([0,250])
    ax2.set_ylabel("CO$_2$ price [€/tCO$_2$]")
    fig.tight_layout()

    fig.savefig("paper_graphics/{}/mwh-pen-co2-{}-{}.pdf".format(scenario,scenario,policy),transparent=True)





def cot_pen_mu_clean(policy):

    fig, ax = plt.subplots()
    fig.set_size_inches((4.5,3))

    plot_df = df.loc[policy].copy()

    plot_df.index = plot_df["wind-solar-penetration"].values*100.

    plot_df["dual"] *= -1

    plot_df["wind-solar-penetration"] *= 100

    plot_df.drop(plot_df.index[plot_df.index < 0.01],inplace=True)

    policy_text = "with flexibility" if "trans-storage" in policy else "without flexibility"

    plot_df.rename(columns={"mp" : "market price",
                            "wind-solar-mv" : "wind+solar MV = LCOE"},inplace=True)

    sel = ["wind+solar MV = LCOE"]
    plot_df[sel].plot(ax=ax,linewidth=2,color=co2_color,
                      style=[choose_style(s) for s in sel],legend=False)
        #s[ (s>1.3)] = np.nan
        #s.plot(ax=ax,label=tech,color=colors[tech])


    ax2 = ax.twinx()

    plot_df["co2_shadow"].plot(ax=ax2,label="CO$_2$ price (right axis)",color=co2_color,style=":")

    ax.set_xlim([0,70])
    ax.set_xlabel("wind+solar penetration [%]")
    ax.set_ylim([0,ylim_comparison])
    ax.set_ylabel("energy price [€/MWh]")
    fig.legend(prop={'size': 9},
              loc='upper left', bbox_to_anchor=(0.15, 0.93))#loc="upper left",

    ax2.set_ylim([0,250])
    ax2.set_ylabel("CO$_2$ price [€/tCO$_2$]")
    fig.tight_layout()

    fig.savefig("paper_graphics/{}/mwh-pen-co2-{}-{}-clean.pdf".format(scenario,scenario,policy),transparent=True)






def sys_cost(policy):

    fig, ax = plt.subplots()
    fig.set_size_inches((4,3))

    plot_df = df.loc[policy].copy()

    sdf = plot_df[plot_df.columns[plot_df.columns.str.contains("-cost")]].rename(columns=(lambda col:col.replace("-cost","")))


    sdf = sdf.divide(plot_df["load"],axis=0)

    rename = pd.Series(sdf.columns,sdf.columns)
    rename.loc[rename.index.str.contains("H2")] = "hydrogen storage"
    rename.loc[rename.index.str.contains("battery")] = "battery storage"
    rename.loc["lign"] = "lignite"
    rename.loc["shed"] = "load-shedding"
    plot_df = sdf.groupby(rename,axis=1).sum()


    plot_df.index = (plot_df.index/20)*float(policy[3:6])/100.

    plot_df["CO$_2$ price"] = plot_df.index*df.loc[policy,"co2_shadow"].values


    plot_df.drop(columns=plot_df.columns[(plot_df.abs() < 2.).all()],inplace=True)

    preferred_order = pd.Index(["nucl","lignite","coal","CCGT","OCGT","wind","solar","battery storage","hydrogen storage","load-shedding"])

    new_index = (preferred_order&plot_df.columns).append(plot_df.columns.difference(preferred_order))

    plot_df[new_index].plot(kind="area",stacked=True,ax=ax,linewidth=0,
                            color=[choose_color(s) for s in new_index])

    #s = df.loc[policy,"mp"]

    #s.index = (s.index/20)*limit

    #s.plot(ax=ax,color="k")

    ax.set_xlim([1.2,0])

    ax.set_ylim([0,150])

    handles,labels = ax.get_legend_handles_labels()

    handles.reverse()
    labels.reverse()

    ax.legend(handles,labels,loc="upper left",ncol=2,prop={'size': 8})

    ax.set_xlabel("average emissions [tCO2/MWhel]")

    ax.set_ylabel("average system cost [€/MWhel]")

    fig.tight_layout()

    fig.savefig("paper_graphics/{}/sys_cost-{}-{}.pdf".format(scenario,scenario,policy),transparent=True)






def comparison():
    tech = "wind-solar"
    fig, ax = plt.subplots()
    fig.set_size_inches((4,3))
    if True:
        policy = "pen{}windsolar-{}-nuclNone-lCCSNone".format(pen,assumptions)
        s = df.loc[policy][tech+"-mv"].copy()
        s.index = float(policy[3:6])*s.index/20
        s[ (s> 200)] = np.nan
        s.plot(ax=ax,label="VRE policy",color=ret_color,linewidth=2)

        policy = "co2120-{}-nuclNone-lCCSNone".format(assumptions)

        s = df.loc[policy][tech+"-mv"].copy()
        s.index = df.loc[policy][tech+"-penetration"].values*100.
        s.drop(s.index[s.index < 0.01],inplace=True)
        s.plot(ax=ax,label="CO$_2$ policy",color=co2_color,linewidth=2)

        ax.legend(loc="upper left")

        ax.set_ylim([0,110])
        ax.set_xlim([0,70])
        ax.set_xlabel("wind+solar penetration [%]")
        ax.set_ylabel("market value [€/MWh]")

        fig.tight_layout()

        fig.savefig("paper_graphics/{}/comparison-{}.pdf".format(scenario,scenario),transparent=True)



def syscost_v_mv():
    tech = "wind-solar"

    policy = "pen{}windsolar-{}-nuclNone-lCCSNone".format(pen,assumptions)

    fig, ax = plt.subplots()
    fig.set_size_inches((4,3))

    plot_df = df.loc[policy].copy()

    costs = plot_df.columns[plot_df.columns.str.contains("-cost")]

    s = plot_df[costs].sum(axis=1)


    s = s.divide(plot_df["load"],axis=0)

    s.index = (s.index/20)*float(policy[3:6])

    s.plot(ax=ax,linewidth=2,color=ret_color,label="VRE system cost")

    s = plot_df["mp"]

    s.index = (s.index/20)*float(policy[3:6])

    s.plot(ax=ax,linewidth=2,color=ret_color,label="VRE market price",style="-.")





    policy = "co2120-{}-nuclNone-lCCSNone".format(assumptions)


    plot_df = df.loc[policy].copy()

    s = plot_df[costs].sum(axis=1)


    s = s.divide(plot_df["load"],axis=0)

    s.index = plot_df[tech+"-penetration"].values*100.
    s.plot(ax=ax,label="CO$_2$ system cost",color=co2_color,linewidth=2)


    s = plot_df["mp"]

    s.index = plot_df[tech+"-penetration"].values*100.

    s.plot(ax=ax,linewidth=2,color=co2_color,label="CO$_2$ market price",style="-.")


    plot_df

    #s = df.loc[policy,"mp"]

    #s.index = (s.index/20)*limit

    #s.plot(ax=ax,color="k")

    ax.set_xlim([0,70])

    ax.set_ylim([0,170])

    ax.legend(loc="upper right",ncol=2,prop={'size': 9})

    ax.set_xlabel("wind+solar penetration [%]")

    ax.set_ylabel("energy price [€/MWh]")

    fig.tight_layout()

    fig.savefig("paper_graphics/{}/compare-sys_cost-{}.pdf".format(scenario,scenario),transparent=True)


def compare_with_co2():

    tech = "wind-solar"

    policy = "pen{}windsolar-{}-nuclNone-lCCSNone".format(pen,assumptions)

    fig, ax = plt.subplots()
    fig.set_size_inches((4,3))

    plot_df = df.loc[policy].copy()

    costs = plot_df.columns[plot_df.columns.str.contains("-cost")]

    s = plot_df[costs].sum(axis=1)


    s = s.divide(plot_df["load"],axis=0)

    s.index = plot_df["emissions"]

    s.plot(ax=ax,linewidth=2,color=ret_color,label="VRE policy")

    xmax = s.index.max()
    xmin = s.index.min()


    #s = plot_df["mp"]

    #s.index = (s.index/20)*float(policy[3:6])

    #s.plot(ax=ax,linewidth=2,color=ret_color,label="VRE market price",style="-.")

    policy = "co2120-{}-nuclNone-lCCSNone".format(assumptions)

    plot_df = df.loc[policy].copy()

    s = plot_df[costs].sum(axis=1)

    s = s.divide(plot_df["load"],axis=0)

    s.index = plot_df["emissions"]
    s.plot(ax=ax,label="CO$_2$ policy",color=co2_color,linewidth=2)

    ax.set_xlim([xmax,xmin])

    ax.legend(loc="upper left",ncol=1,prop={'size': 9})

    ax.set_ylim([0,100])

    ax.set_xlabel("average emissions [tCO2/MWhel]")

    ax.set_ylabel("average system cost [€/MWh]")

    fig.tight_layout()

    fig.savefig("paper_graphics/{}/compare-sys_cost-co2-{}.pdf".format(scenario,scenario),transparent=True)


if __name__ == "__main__":
    # Detect running outside of snakemake and mock snakemake for testing
    if 'snakemake' in globals():
        fn = snakemake.input["summary"]
        config = snakemake.config
    else:
        from pypsa.descriptors import Dict
        config = Dict()
        config["run"] = "190919-final"
        config["results_dir"] = "results"
        fn = "{}/{}/csvs/summary.csv".format(config["results_dir"],config["run"])

    scenario =  config["run"]

    df = pd.read_csv(fn,index_col=[0,1])

    assumptions = "wind1040-sola510"
    pen = "075"
    ylim_comparison = 140.


    for policy in ["pen{}{}-{}-nuclNone-lCCSNone".format(pen,"{}",assumptions),"pen{}{}-co2price20".format(pen,"{}"),"availpen{}{}-co2price20".format(pen,"{}")]:
        plot_re_penetration(policy)

    for tech in ["wind","solar","wind-solar"]:
        for policy in ["pen{}{}-{}-nuclNone-lCCSNone".format(pen,tech.replace("-",""),assumptions),"pen100{}-{}-nuclNone-lCCSNone-trans-storage".format(tech.replace("-",""),assumptions),"pen100{}-{}-nuclNone-lCCSNone-battery".format(tech.replace("-",""),assumptions)]:
            if tech == "wind-solar" and "battery" in policy:
                continue
            ret_tech_fig(tech,policy)
            ret_tech_fig_clean(tech,policy)

    tech = "nucl"
    for policy in ["pen100{}-{}-nucl6000-lCCSNone".format(tech,assumptions),"pen100{}-{}-nucl10000-lCCSNone-trans-storage".format(tech,assumptions),f"pen100{tech}-{assumptions}-nucl10000-lCCSNone-coalNone-lignNone-OCGTNone-CCGTNone-trans-storage"]:
        ret_tech_fig(tech,policy)


    pen_compare()

    for policy in ["co2120-trans-storage-{}-nuclNone-lCCSNone".format(assumptions),"co2120-{}-nuclNone-lCCSNone".format(assumptions)]:
        co2_as_em(policy)
        co2_as_pen(policy)

    cot_flex_compare()
    cot_flex_compare_clean()
    cot_pen_mu("co2120-{}-nuclNone-lCCSNone".format(assumptions))
    cot_pen_mu_clean("co2120-{}-nuclNone-lCCSNone".format(assumptions))
    sys_cost("co2120-trans-storage-{}-nuclNone-lCCSNone".format(assumptions))
    comparison()
    syscost_v_mv()
    compare_with_co2()
