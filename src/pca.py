#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PCA plots

Created on  Nov 15  2022
@author: johanna
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import seaborn as sns
from scipy import stats
import argparse
from matplotlib.patches import Ellipse
import functions_general as fg


def pca_args():
    parser = argparse.ArgumentParser(
        prog="python -m DIMet.src.pca",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('config', type=str,
                        help="configuration file in absolute path")

    parser.add_argument('--save_pca_tables',
                        action=argparse.BooleanOptionalAction, default=False,
                        help="save pca results to .csv files")

    parser.add_argument('--draw_ellipses',
                        action=argparse.BooleanOptionalAction, default=False,
                        help="draw ellipse for each group in the pca plot")

    parser.add_argument('--run_iris_demo',
                        action=argparse.BooleanOptionalAction, default=False,
                        help="draws full pca, with ellipses, in iris data")

    return parser


def make_ellipse(mean, cov, level=0.95, color=None):
    """Support function for scatter_ellipse.
       Use make_ellipse OR eigsorted, but not both"""

    v, w = np.linalg.eigh(cov)
    u = w[0] / np.linalg.norm(w[0])
    angle = np.arctan(u[1] / u[0])
    angle = 180 * angle / np.pi  # convert to degrees
    v = 2 * np.sqrt(
        v * stats.chi2.ppf(level, 2))  # get size corresponding to level
    ell = Ellipse(mean[:2], v[0], v[1], 180 + angle, facecolor='none',
                  edgecolor="gray", alpha=0.4,
                  linestyle='-',
                  lw=1.5)
    return ell


def eigsorted(cov):
    """
    used for calculating ellipses
    many thanks to :
    https://rayblick.gitbooks.io/my-python-scrapbook/content/
    analysis/plotting/scatterplot_ellipse.html

    use make_ellipse OR eigsorted, but not both
    """
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    return vals[order], vecs[:, order]


def clean_reduce_datadf_4pca(df, metadatasub):
    df = df[metadatasub['name_to_plot']]
    df = df.loc[~(df == 0).all(axis=1)]  # drop 'zero all' rows
    df = df.fillna(df[df > 0].min().min())
    df = df.replace(0, df[df > 0].min().min())
    df_red = fg.give_reduced_df(df, ddof=0)  # reduce rows
    if np.isinf(np.array(df_red)).any():  # avoid Inf error
        return df
    else:
        return df_red


def compute_pca(mymat, metadata):
    dims = min(metadata.shape[0],
               mymat.shape[0])  # min(nb_samples, nb_features)

    X = np.transpose(np.array(mymat))

    pca = PCA(n_components=dims)

    pc = pca.fit_transform(X)

    pc_df = pd.DataFrame(data=pc,
                         columns=['PC' + str(i) for i in range(1, dims + 1)])
    pc_df = pc_df.assign(name_to_plot=mymat.columns)

    pc_df = pd.merge(pc_df, metadata, on='name_to_plot')

    var_explained_df = pd.DataFrame({
        'var': pca.explained_variance_ratio_ * 100,
        'PC': ['PC' + str(i) for i in range(1, dims + 1)]})

    return pc_df, var_explained_df


def save_pca_plots(title, pc_df, var_explained_df, col1, col2, pointlabels,
                   out_plot_dir, *args) -> None:
    """
       *args :  options for advanced PCA plotting:
             args[0] is column name for ellipses
    """

    col_ellipses = None
    try:
        col_ellipses = args[0]
    except IndexError:
        pass  # ellipses args not set, not drawn by default
    except Exception as e:
        print("unknown error! ", e)

    # barplot
    plt.figure()
    sns.barplot(x='PC', y="var", data=var_explained_df, color="cadetblue")
    plt.title("Principal components: explained variance")
    plt.ylabel("variance (%)")
    plt.savefig(f'{out_plot_dir}variance_pca_{title.replace(" ", "-")}.pdf')
    plt.close()

    # scatterplot

    # sns.set_style("whitegrid")
    fig, ax = plt.subplots()
    sns.scatterplot(x="PC1", y="PC2",
                    ax=ax,
                    data=pc_df,
                    hue=col1,
                    style=col2,
                    legend=True,
                    s=80, zorder=3)
    ax.axhline(0, ls="--", color="gray", zorder=1)
    ax.axvline(0, ls="--", color="gray", zorder=1)

    yesnolabel = "no"

    if pointlabels != "":
        yesnolabel = "yes"
        for i, row in pc_df.iterrows():
            ax.text(pc_df.at[i, 'PC1'] + 0.2, pc_df.at[i, 'PC2'],
                    pc_df.at[i, pointlabels],
                    size='x-small')
    # end if

    row_xlab = var_explained_df.iloc[0, :]
    row_ylab = var_explained_df.iloc[1, :]

    plt.xlabel(
        f"{row_xlab['PC']} {round(row_xlab['var'], 2)} %")
    plt.ylabel(
        f'{row_ylab["PC"]} {round(row_ylab["var"], 2)} %')
    plt.title(title)

    # ellipses, if true
    if col_ellipses is not None:
        myellipsesnames = pc_df[col_ellipses].unique()
        for lab in myellipsesnames:
            xdata = pc_df.loc[pc_df[col_ellipses] == lab, 'PC1']
            ydata = pc_df.loc[pc_df[col_ellipses] == lab, 'PC2']

            # get values to build the ellipse
            cov = np.cov(xdata, ydata)
            # ell = make_ellipse((np.mean(xdata), np.mean(ydata)),cov,
            #                  level= 0.95, color=None)
            vals, vecs = eigsorted(cov)
            theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
            # print(vals, vecs, theta)
            w, h = 2 * 2 * np.sqrt(vals)

            # create the ellipse
            ell = Ellipse(xy=(np.mean(xdata), np.mean(ydata)),
                          width=w, height=h,
                          angle=theta,
                          edgecolor="lightgray",
                          linestyle='-', facecolor='none')

            # ell.set_facecolor() # reference the colour for each factor
            ax.add_artist(ell)
    # end  if ellipses

    plt.savefig(
        f'{out_plot_dir}pca_{title.replace(" ", "-")}_label{yesnolabel}.pdf',
        format="pdf")


def run_steps_pca(type_measure: str, table_prefix: str,
                  metadatadf: pd.DataFrame,
                  out_plot_dir: str, confidic: dict, args) -> None:
    out_path = os.path.expanduser(confidic['out_path'])
    suffix = confidic['suffix']
    compartments = metadatadf['short_comp'].unique().tolist()
    # dynamically open the file based on prefix, compartment and suffix:
    for co in compartments:
        meta_sub = metadatadf.loc[metadatadf['short_comp'] == co, :]
        fn0 = f'{out_path}results/prepared_tables/'
        fn1 = f'{table_prefix}--{co}--{suffix}.tsv'
        fn = fn0 + fn1
        measurements = pd.read_csv(fn, sep='\t', header=0, index_col=0)

        timepoints = meta_sub['timepoint'].unique().tolist()

        if len(timepoints) >= 2:  # worthed
            mat = clean_reduce_datadf_4pca(measurements, meta_sub)
            pc_df, dfvare = compute_pca(mat, meta_sub)
            pc_df = pc_df.assign(
                col_label=pc_df['name_to_plot'].str.replace(co, ""))
            title = f'{type_measure} {co} {suffix}'
            if args.draw_ellipses:
                # with label
                save_pca_plots(title, pc_df, dfvare, "timepoint", "condition",
                               "col_label", out_plot_dir, "timepoint")

                # no label
                save_pca_plots(title, pc_df, dfvare, "timepoint", "condition",
                               "", out_plot_dir, "timepoint")

            else:
                # with label
                save_pca_plots(title, pc_df, dfvare, "timepoint", "condition",
                               "col_label", out_plot_dir)

                # no label
                save_pca_plots(title, pc_df, dfvare, "timepoint", "condition",
                               "", out_plot_dir)

            if args.save_pca_tables:
                pc_df.to_csv(f'{out_plot_dir}{title.replace(" ","-")}_pc.csv',
                             sep='\t')
                dfvare.to_csv(f'{out_plot_dir}{title.replace(" ","-")}_var.csv',
                             sep='\t')
        # end if

        # valid for any timepoints length:
        for ti in timepoints:
            title_ti = f'{type_measure} {co} {ti} {suffix}'
            meta_ti = meta_sub.loc[(meta_sub['short_comp'] == co) & (
                        meta_sub['timepoint'] == ti), :]
            mat_ti = clean_reduce_datadf_4pca(measurements, meta_ti)
            pc_df, dfvare = compute_pca(mat_ti, meta_ti)


            # always with labels
            pc_df = pc_df.assign(
                col_label=pc_df['name_to_plot'].str.replace(co, ""))
            if args.draw_ellipses:
                save_pca_plots(title_ti, pc_df, dfvare, "condition",
                               "condition", "col_label",
                               out_plot_dir, "condition")
            else:
                save_pca_plots(title_ti, pc_df, dfvare, "condition",
                               "condition", "col_label",
                               out_plot_dir)

            if args.save_pca_tables:
                pc_df.to_csv(f'{out_plot_dir}{title_ti.replace(" ","-")}_pc.csv',
                             sep='\t')
                dfvare.to_csv(f'{out_plot_dir}{title_ti.replace(" ","-")}_var.csv',
                              sep='\t')
    # end for


def run_pca_in_iris(outdir) -> None:
    iris = sns.load_dataset("iris")
    sns.relplot(data=iris, x="sepal_width", y="petal_width",
                hue="species")
    iris = iris.assign(name_to_plot=[str(i) for i in iris.index])
    fakemeta = iris[['name_to_plot', "species"]]
    iris = iris.drop(columns=['name_to_plot', "species"])
    fakedf = iris.T  # variables rows, samples columns
    fakedf = fakedf.div(fakedf.std(axis=1, ddof=0), axis=0)

    fakedf.columns = [str(i) for i in iris.index]
    pc_df, dfvare = compute_pca(fakedf, fakemeta)
    save_pca_plots("Iris", pc_df, dfvare,
                   "species", "species", "", outdir, "species")


if __name__ == "__main__":

    parser = pca_args()
    args = parser.parse_args()
    configfile = os.path.expanduser(args.config)
    confidic = fg.open_config_file(configfile)
    fg.auto_check_validity_configuration_file(confidic)
    out_path = os.path.expanduser(confidic['out_path'])
    meta_path = os.path.expanduser(confidic['metadata_path'])
    clean_tables_path = out_path + "results/prepared_tables/"

    if args.run_iris_demo:
        run_pca_in_iris(out_path + "results/plots/")

    # tpd : tables prefixes dictionary
    tpd = fg.clean_tables_names2dict(
        f'{out_path}results/prepared_tables/TABLESNAMES.csv')
    metadatadf = fg.open_metadata(meta_path)

    # separately using fracContrib and abundance for the pca plots

    # pca's for abund

    abund_tab_prefix = tpd['name_abundance']
    out_plot_dir_pca_abun = out_path + "results/plots/pca_Abundance/"
    fg.detect_and_create_dir(out_plot_dir_pca_abun)
    run_steps_pca("Abundance", abund_tab_prefix, metadatadf,
                  out_plot_dir_pca_abun, confidic, args)

    # pca's for fraccon

    fraccon_tab_prefix = tpd['name_meanE_or_fracContrib']
    out_plot_dirpca_fc = out_path + "results/plots/pca_fracCorME/"
    fg.detect_and_create_dir(out_plot_dirpca_fc)
    run_steps_pca("fracContrib", fraccon_tab_prefix, metadatadf,
                  out_plot_dirpca_fc, confidic, args)

# end

# calculate Ellipses , many thanks to :
# https://rayblick.gitbooks.io/my-python-scrapbook/content/analysis/plotting/scatterplot_ellipse.html

# ex4:
# https://www.programcreek.com/python/example/61396/matplotlib.patches.Ellipse
