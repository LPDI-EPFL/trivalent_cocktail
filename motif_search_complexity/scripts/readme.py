# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
import gzip
import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import rstoolbox


def load_master(dfs):
    """Load the master search processed data.
    """
    for k in dfs:
        if not isinstance(dfs[k], int):
            continue
        dfs[k] = pd.read_csv(os.path.join("data", k, "master_search.csv.gz"))
        if os.path.isfile(os.path.join("data", k, "ddg_match.csv.gz")):
            tmp = pd.read_csv(os.path.join("data", k, "ddg_match.csv.gz"))
            dfs[k] = dfs[k].merge(tmp, how="left", on=["cluster", "str"])
    return dfs


def parse_pfam(f):
    """Read and process the PFAM file
    """
    def pfam_length(pfam_range):
        r = [int(x) for x in pfam_range.split("-")]
        return r[1] - r[0] + 1
    data = {"pdb": [], "chain": [], "pfamrange": [], "pfam": []}
    with gzip.open(f, 'rt') as fd:
        for line in fd:
            ln = [x.strip() for x in line.strip().split(";")]
            data["pdb"].append(ln[0].lower())
            data["chain"].append(ln[1])
            data["pfamrange"].append(ln[-2])
            data["pfam"].append(ln[3])

    df = pd.DataFrame(data)
    df["length"] = df.apply(lambda row: pfam_length(row["pfamrange"]), axis=1)
    return df


def pfam2master(dfs, pfam):
    """Assign pfam data to master searches data.
    """
    def fit_in_range(motif_range, pfam_range):
        """As a single protein can have multiple PFAM domains,
        this function evaluates if the motif match falls
        inside the assigned PFAM domain.
        """
        if pfam_range == np.nan or isinstance(pfam_range, float):
            return False
        pfam_range = [int(x) for x in pfam_range.split("-")]
        motif_range = (motif_range.replace(")(", "),(")
                                  .replace("(", "")
                                  .replace(")", "")
                                  .strip("[]"))
        cmatch = 0
        for r in motif_range.split(","):
            r = [int(x) for x in r.split("-")]
            if r[0] >= pfam_range[0] and r[1] <= pfam_range[1]:
                cmatch += 1
            else:
                return False
        return cmatch == len(motif_range.split(","))

    for k in dfs:
        dfs[k] = dfs[k].merge(pfam, how="left", on=["pdb", "chain"])
        dfs[k]["inrange"] = dfs[k].apply(
            lambda row: fit_in_range(row["range"],
                                     row["pfamrange"]),
            axis=1)
    return dfs


def plot_all(fig, dfs, total_master_list, min_domain_size,
             max_domain_size, top_limit, mode):
    """Plot data for all cases.
    """
    grid = (2, 3)
    ax11 = plt.subplot2grid(grid, (0, 0), fig=fig)
    ax12 = plt.subplot2grid(grid, (0, 1), fig=fig)
    ax13 = plt.subplot2grid(grid, (0, 2), fig=fig)
    ax21 = plt.subplot2grid(grid, (1, 0), fig=fig)
    ax22 = plt.subplot2grid(grid, (1, 1), fig=fig)
    ax23 = plt.subplot2grid(grid, (1, 2), fig=fig)
    if isinstance(dfs["2fx7"], pd.DataFrame):
        plot(dfs["2fx7"], "2fx7", total_master_list, "2FX7: HIV 4e10 (H)",
             ax11, [0.89], ["1Z6N"], kcolor=["black"], annotate=True,
             top=top_limit, min_length=min_domain_size,
             max_length=max_domain_size)

    if isinstance(dfs["3ixt"], pd.DataFrame):
        if mode == 0:
            plot(dfs["3ixt"], "3ixt", total_master_list,
                 "3IXT: RSV mota (HLH)", ax12, [1.10], ["3LHP"],
                 kcolor=["black"], annotate=True, top=top_limit,
                 min_length=min_domain_size, max_length=max_domain_size)
        else:
            plot(dfs["3ixt"], "3ixt", total_master_list,
                 "3IXT: RSV mota (HLH)", ax12, [1.10, 2.37], ["3LHP", "1KX8"],
                 kcolor=["black", "blue"], min_length=min_domain_size,
                 max_length=max_domain_size)

    if isinstance(dfs["5tpn"], pd.DataFrame):
        plot(dfs["5tpn"], "5tpn", total_master_list, "5TPN: RSV RSV90 (HLE)",
             ax13, annotate=True, top=top_limit, min_length=min_domain_size,
             max_length=max_domain_size)

    if isinstance(dfs["3o41"], pd.DataFrame):
        if mode == 0:
            plot(dfs["3o41"], "3o41", total_master_list, "3O41: RSV 101F (E)",
                 ax21, annotate=True, top=top_limit,
                 min_length=min_domain_size, max_length=max_domain_size)
        else:
            plot(dfs["3o41"], "3o41", total_master_list, "3O41: RSV 101F (E)",
                 ax21, [2.30], ["TOP7"], kcolor=["blue"],
                 min_length=min_domain_size, max_length=max_domain_size)

    if isinstance(dfs["3vtt"], pd.DataFrame):
        plot(dfs["3vtt"], "3vtt", total_master_list, "3VTT: DEN3 ED3 (ELE)",
             ax22, [0.68], ["3WEI"], kcolor=["black"], annotate=True,
             top=top_limit, min_length=min_domain_size,
             max_length=max_domain_size)

    if isinstance(dfs["4jhw"], pd.DataFrame):
        plot(dfs["4jhw"], "4jhw", total_master_list, "4JHW: RSV D25 (HxL)",
             ax23, annotate=True, top=top_limit, min_length=min_domain_size,
             max_length=max_domain_size)


def plot(df, selfname, maxim, title, ax, known=None, knames=None, kcolor=None,
         annotate=False, top=10, rmsd_lim=5, min_length=50, max_length=100):
    """
    Full plot
    """
    def data_plot(df, ax, maxim, rmsd_lim, linestyle, color):
        """Calculate true cumulative curves
        (seaborn does an aproximation that does not work well with very
        small recovery)
        """
        allvalues = (df.sort_values("rmsd")
                       .groupby(["pdb", "chain"])
                       .head(1)[["rmsd"]].values)
        raw, y, x = rstoolbox.analysis.cumulative(allvalues, max_count=maxim,
                                                  upper_limit=rmsd_lim)
        ax.plot(x, y, color=sns.color_palette()[color], lw=4,
                linestyle=linestyle)
        return raw, x, y

    def first_marker_plot(raw, x, ax, top, shape, color, linestyle):
        """Top x marker
        """
        idx = (np.abs(np.array(raw)-top)).argmin()
        ax.axvline(x=x[idx], ymin=0, ymax=0.2, c=sns.color_palette()[color],
                   linewidth=4, linestyle=linestyle, zorder=10)
        ax.plot([x[idx]], [0.2], shape, c=sns.color_palette()[color],
                markersize=12)

    # 1. Filter self-hits and error assignments
    df = df[(df["pdb"] != "eeee") & (df["pdb"] != selfname)]

    # 2. Make cumulative plot for best-hit/pdb-chain for all dataset
    raw, x, y = data_plot(df, ax, maxim, rmsd_lim, "solid", 0)
    if annotate:
        first_marker_plot(raw, x, ax, top, 'o', 0, "solid")

    # 3. Make cumulative plot for best-hit/pdb-chain for proteins
    # of size "max_length" or smaller
    sdf = df[(df["length"] >= min_length) &
             (df["length"] <= max_length) &
             (df["inrange"])]
    raw, x, y = data_plot(sdf, ax, maxim, rmsd_lim, "dashed", 0)
    if annotate:
        first_marker_plot(raw, x, ax, top, 's', 0, "dashed")

    # 4. Add clashes data
    if "ddg" in df:
        # 4.1. Make cumulative plot for best-hit/pdb-chain for
        # all non-clashing dataset
        raw, x, y = data_plot(df[(df["ddg"] <= 0)], ax, maxim,
                              rmsd_lim, "solid", 1)
        ax.fill_between(x, 0, y, color=sns.color_palette()[1], alpha=0.3)
        if annotate:
            first_marker_plot(raw, x, ax, top, 'o', 1, "solid")

        # 4.2. Make cumulative plot for best-hit/pdb-chain for all non-clashing
        # dataset for proteins of size "max_length" or smaller
        sdf = df[(df["length"] >= min_length) &
                 (df["ddg"] <= 0) & (df["length"] <= max_length) &
                 (df["inrange"])]
        raw, x, y = data_plot(sdf, ax, maxim, rmsd_lim, "dashed", 1)
        if annotate:
            first_marker_plot(raw, x, ax, top, 's', 1, "dashed")

    # 5. Label where known designs match
    if known is not None:
        for i, k in enumerate(known):
            ax.axvline(x=k, ymin=0, ymax=1, c=kcolor[i], linewidth=2,
                       linestyle="dashed", zorder=10)

    # 6. Label where own designs mathc
    if knames is not None:
        for i, k in enumerate(knames):
            ax.text(known[i] - 0.3, 0.8 - (0.06 * i), k, color=kcolor[i],
                    horizontalalignment='center', size=25, weight='semibold')

    # 7. Show the max expected number to retrieve defined as "maxim";
    # value = 1 in y scale
    ax.text(4.02, 0.15, "n = {}".format(maxim), color="black",
            horizontalalignment='left', size=25, weight='semibold')

    # 8. Format axis
    ax.set_ylim(0, 1)
    ax.set_xlim(0, 5)
    ax.set_xlabel("RMSD")
    ax.set_ylabel("")
    rstoolbox.utils.add_top_title(ax, title, size=30)


def plot_all_usable(fig, dfs, total_master_list, min_domain_size,
                    max_domain_size, top_limit):
    """Plot data for all cases.
    """
    grid = (2, 3)
    ax11 = plt.subplot2grid(grid, (0, 0), fig=fig)
    ax12 = plt.subplot2grid(grid, (0, 1), fig=fig)
    ax13 = plt.subplot2grid(grid, (0, 2), fig=fig)

    ax22 = plt.subplot2grid(grid, (1, 1), fig=fig)
    ax23 = plt.subplot2grid(grid, (1, 2), fig=fig)

    linethickness = 2

    if isinstance(dfs["2fx7"], pd.DataFrame):
        plot_usable(dfs["2fx7"], "2fx7", total_master_list,
                    "2FX7: HIV 4e10 (H)", ax11,
                    top=top_limit, min_length=min_domain_size,
                    max_length=max_domain_size)
        ax11.set_xticklabels([0, 1, 2, 3, 4, ''])
        for axis in ['top', 'bottom', 'left', 'right']:
            ax11.spines[axis].set_linewidth(linethickness)
            ax11.spines[axis].set_color('black')

    if isinstance(dfs["3ixt"], pd.DataFrame):
        plot_usable(dfs["3ixt"], "3ixt", total_master_list,
                    "3IXT: RSV mota (HLH)", ax12, top=top_limit,
                    min_length=min_domain_size, max_length=max_domain_size)
        ax12.set_xticklabels([])
        ax12.set_xlabel('')
        ax12.set_yticklabels([])
        ax12.set_ylabel('')
        for axis in ['top', 'bottom', 'left', 'right']:
            ax12.spines[axis].set_linewidth(linethickness)
            ax12.spines[axis].set_color('black')

    if isinstance(dfs["5tpn"], pd.DataFrame):
        plot_usable(dfs["5tpn"], "5tpn", total_master_list,
                    "5TPN: RSV RSV90 (HLE)", ax13, top=top_limit,
                    min_length=min_domain_size, max_length=max_domain_size)
        ax13.set_xticklabels([])
        ax13.set_xlabel('')
        ax13.set_yticklabels([])
        ax13.set_ylabel('')
        for axis in ['top', 'bottom', 'left', 'right']:
            ax13.spines[axis].set_linewidth(linethickness)
            ax13.spines[axis].set_color('black')

    if isinstance(dfs["3o41"], pd.DataFrame):
        plot_usable(dfs["3o41"], "3o41", total_master_list,
                    "3O41: RSV 101F (E)", ax22, top=top_limit,
                    min_length=min_domain_size, max_length=max_domain_size)
        ax22.set_yticklabels([0.0, 0.2, 0.4, 0.6, 0.8, ''])
        for axis in ['top', 'bottom', 'left', 'right']:
            ax22.spines[axis].set_linewidth(linethickness)
            ax22.spines[axis].set_color('black')

    if isinstance(dfs["4jhw"], pd.DataFrame):
        plot_usable(dfs["4jhw"], "4jhw", total_master_list,
                    "4JHW: RSV D25 (HxL)", ax23, top=top_limit,
                    min_length=min_domain_size, max_length=max_domain_size)
        ax23.set_yticklabels([])
        ax23.set_ylabel('')
        for axis in ['top', 'bottom', 'left', 'right']:
            ax23.spines[axis].set_linewidth(linethickness)
            ax23.spines[axis].set_color('black')

def plot_usable(df, selfname, maxim, title, ax,
                top=10, rmsd_lim=5, min_length=50, max_length=100):
    """
    Full plot
    """
    def data_plot(df, ax, maxim, rmsd_lim, linestyle, color):
        """Calculate true cumulative curves
        (seaborn does an aproximation that does not work well with very
        small recovery)
        """
        allvalues = (df.sort_values("rmsd")
                       .groupby(["pdb", "chain"])
                       .head(1)[["rmsd"]].values)
        raw, y, x = rstoolbox.analysis.cumulative(allvalues, max_count=maxim,
                                                  upper_limit=rmsd_lim)
        ax.plot(x, y, color=sns.color_palette()[color], lw=4,
                linestyle=linestyle)
        return raw, x, y

    def first_marker_plot(raw, x, ax, top, shape, color, linestyle):
        """Top x marker
        """
        idx = (np.abs(np.array(raw)-top)).argmin()
        ax.axvline(x=x[idx], ymin=0, ymax=0.2, c=sns.color_palette()[color],
                   linewidth=4, linestyle=linestyle, zorder=10)
        ax.plot([x[idx]], [0.2], shape, c=sns.color_palette()[color],
                markersize=12)

    # 1. Filter self-hits and error assignments
    df = df[(df["pdb"] != "eeee") & (df["pdb"] != selfname)]

    # 4.2. Make cumulative plot for best-hit/pdb-chain for all non-clashing
    # dataset for proteins of size "max_length" or smaller
    sdf = df[(df["length"] >= min_length) &
             (df["ddg"] <= 0) & (df["length"] <= max_length) &
             (df["inrange"])]
    raw, x, y = data_plot(sdf, ax, maxim, rmsd_lim, "solid", 1)
    ax.fill_between(x, 0, y, color=sns.color_palette()[1], alpha=0.3)
    first_marker_plot(raw, x, ax, top, 's', 1, "solid")

    # 7. Show the max expected number to retrieve defined as "maxim";
    # value = 1 in y scale
    ax.text(4.9, 0.3, "n = {}".format(maxim), color="black",
            horizontalalignment='right', size=25, weight='semibold')

    # 8. Format axis
    ax.set_ylim(0, 1)
    ax.set_xlim(0, 5)
    ax.set_xlabel("RMSD")
    ax.set_ylabel("")
    ax.text(4.9, 0.90, title, color="black",
            horizontalalignment='right', size=30, weight='semibold')
