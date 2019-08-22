import pandas as pd

#allow plotting without Xwindows
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

colors={"wind" : "b",
        "solar" : "y",
        "solar-wind" : "g"}


def plot_res(df, suffix,assumptions):
    fig, ax = plt.subplots()
    fig.set_size_inches((6,4))

    for tech in ["solar","wind","solar-wind"]:

        s = df.loc[pd.IndexSlice[tech+suffix,:,assumptions]][tech+"-rmv"]
        s.index = s.index/1e3
        s.plot(ax=ax,label=tech,color=colors[tech])

        ax.legend()

        ax.set_ylim([0,1.3])
        ax.set_xlabel("penetration [pu]")
        ax.set_ylabel("relative market value [pu]")

        fig.tight_layout()

        fig.savefig(config['results_dir'] + "/" + config['run'] + "/graphs/rmvs-res{}_{}.pdf".format(suffix,assumptions),transparent=True)

def plot_co2(df, suffix, assumptions):

    fig, ax = plt.subplots()
    fig.set_size_inches((6,4))

    for tech in ["solar","wind","solar-wind"]:

        s = df.loc[pd.IndexSlice["co2"+suffix,:,assumptions]][tech+"-rmv"]
        s.index = s.index/1e3


        s.plot(ax=ax,label=tech,color=colors[tech])

        ax.legend()

        ax.set_ylim([0,1.5])
        ax.set_xlim([0.7,0])
        ax.set_xlabel("average_emissions [tCO2/MWhel]")
        ax.set_ylabel("relative market value [pu]")

        fig.tight_layout()

        fig.savefig(config['results_dir'] + "/" + config['run'] + "/graphs/rmvs-co2{}_{}.pdf".format(suffix,assumptions),transparent=True)


def plot_co2_penetration(df, suffix,assumptions):
    """Plot CO2 scenario results, but against penetration, not CO2."""

    fig, ax = plt.subplots()
    fig.set_size_inches((6,4))

    for tech in ["solar","wind","solar-wind"]:

        s = df.loc[pd.IndexSlice["co2"+suffix,:,assumptions]][tech+"-rmv"]
        s.index = df.loc[pd.IndexSlice["co2"+suffix,:,assumptions]][tech+"-penetration"].values



        s.plot(ax=ax,label=tech,color=colors[tech])

        ax.legend()

        ax.set_ylim([0,1.5])
        ax.set_xlabel("penetration [pu]")
        ax.set_ylabel("relative market value [pu]")

        fig.tight_layout()

        fig.savefig(config['results_dir'] + "/" + config['run'] + "/graphs/rmvs-co2-penetration{}_{}.pdf".format(suffix,assumptions),transparent=True)

if __name__ == "__main__":
    # Detect running outside of snakemake and mock snakemake for testing
    if 'snakemake' not in globals():
        from pypsa.descriptors import Dict

    config = snakemake.config

    fn = snakemake.input["summary"]
    df = pd.read_csv(fn,index_col=[0,1,2])

    for assumptions in df.index.levels[2]:
        for suffix in ["","-storage"]:
            plot_res(df, suffix,assumptions)
            plot_co2(df, suffix,assumptions)
            plot_co2_penetration(df, suffix,assumptions)
