import pandas as pd

#allow plotting without Xwindows
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

colors={"wind" : "b",
        "solar" : "y",
        "solar-wind" : "g"}


def plot_res(df, suffix):
    fig, ax = plt.subplots()
    fig.set_size_inches((6,4))

    for tech in ["solar","wind","solar-wind"]:

        s = df.loc[tech+suffix][tech+"-rmv"]
        s.index = s.index/1e3
        s.plot(ax=ax,label=tech,color=colors[tech])

        ax.legend()

        ax.set_ylim([0,1.3])
        ax.set_xlabel("penetration [pu]")
        ax.set_ylabel("relative market value [pu]")

        fig.tight_layout()

        fig.savefig(config['results_dir'] + "/" + config['run'] + "/graphs/rmvs-res{}.pdf".format(suffix),transparent=True)

def plot_co2(df, suffix):

    fig, ax = plt.subplots()
    fig.set_size_inches((6,4))

    for tech in ["solar","wind","solar-wind"]:

        s = df.loc["co2"+suffix][tech+"-rmv"]
        s.index = s.index/1e3


        s.plot(ax=ax,label=tech,color=colors[tech])

        ax.legend()

        ax.set_ylim([0,1.5])
        ax.set_xlim([0.7,0])
        ax.set_xlabel("average_emissions [tCO2/MWhel]")
        ax.set_ylabel("relative market value [pu]")

        fig.tight_layout()

        fig.savefig(config['results_dir'] + "/" + config['run'] + "/graphs/rmvs-co2{}.pdf".format(suffix),transparent=True)


def plot_co2_penetration(df, suffix):
    """Plot CO2 scenario results, but against penetration, not CO2."""

    fig, ax = plt.subplots()
    fig.set_size_inches((6,4))

    for tech in ["solar","wind","solar-wind"]:

        s = df.loc["co2"+suffix][tech+"-rmv"]
        s.index = df.loc["co2"+suffix][tech+"-penetration"].values



        s.plot(ax=ax,label=tech,color=colors[tech])

        ax.legend()

        ax.set_ylim([0,1.5])
        ax.set_xlabel("penetration [pu]")
        ax.set_ylabel("relative market value [pu]")

        fig.tight_layout()

        fig.savefig(config['results_dir'] + "/" + config['run'] + "/graphs/rmvs-co2-penetration{}.pdf".format(suffix),transparent=True)

if __name__ == "__main__":
    # Detect running outside of snakemake and mock snakemake for testing
    if 'snakemake' not in globals():
        from pypsa.descriptors import Dict

    config = snakemake.config

    fn = snakemake.input["summary"]
    df = pd.read_csv(fn,index_col=[0,1])

    for suffix in ["","-storage"]:
        plot_res(df, suffix)
        plot_co2(df, suffix)
        plot_co2_penetration(df, suffix)
