import matplotlib;

matplotlib.use('Agg')  # forces matplotlib to not launch X11 window
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import seaborn as sns
import pandas as pd
import cv2

from scipy import stats

__all__ = ["focus_plot"]


def make_scatter(twas_df, scale="logp"):
    """
    Make a scatterplot of zscore values with gene names as xtick labels.

    :param twas_df: pandas.DataFrame containing at least zscores and gene-names
    :param scale: string the scale to plot on. values are 'logp' or 'zscore'. Default logp

    :return: numpy.ndarray (RGB) formatted scatterplot of zscores
    """
    mpl.rcParams["figure.figsize"] = [6.4, 4.8]
    fig = plt.figure()

    if scale == "zscore":
        ax = sns.stripplot(x="mol_name", y="twas_z", data=twas_df, color='black', size=5)
        # make ymin, ymax symmetric
        zmin = np.min(twas_df['twas_z'])
        zmax = np.max(twas_df['twas_z'])
        abs_bound = round(max(np.abs(zmin), np.abs(zmax)))

        if abs_bound % 2 != 0:
            abs_bound += 1

        ax.set_ylim([abs_bound * -1, abs_bound])
        ylabel = 'Z-score'
    elif scale == "logp":
        size_arr = []
        color_arr = []
        custom_palette = ["#e4f1fe", "#89c4f4", "#2574a9", "#013243"]
        size_palette = [4, 8, 10, 12]

        for i, row in twas_df.iterrows():
            pip = row['pip']
            if pip >= 0 and pip < .20:
                color_arr.append(custom_palette[0])
                size_arr.append(size_palette[0])
            elif pip >= .20 and pip < .40:
                color_arr.append(custom_palette[1])
                size_arr.append(size_palette[1])
            elif pip >=.40 and pip <.60:
                color_arr.append(custom_palette[2])
                size_arr.append(size_palette[2])
            elif pip >= .60 and pip <.80:
                color_arr.append(custom_palette[3])
                size_arr.append(size_palette[3])
            else:
                color_arr.append(custom_palette[3])
                size_arr.append(size_palette[3])

        twas_df = twas_df.assign(logp=-stats.chi2.logsf(twas_df["twas_z"].values ** 2, 1))
        n_rows = len(twas_df.index)
        x_values = np.arange(1,n_rows + 1)
        ax = sns.scatterplot(x=x_values, y="logp", data=twas_df, legend=False, size=size_arr, color=color_arr, edgecolor='black')
        ylabel = "-log10 p-value"
        plt.legend(loc="best", fontsize=10, scatterpoints=1)
        from matplotlib.patches import Patch
        from matplotlib.lines import Line2D
        legend_elements = [
                   Line2D([0], [0], marker='o', color='w', label='[0.0, 0.2)',
                          markerfacecolor=custom_palette[0], markersize=size_palette[0], markeredgecolor='k'),
                  Line2D([1], [1], marker='o', color='w', label='[0.2, 0.4)',
                         markerfacecolor=custom_palette[1], markersize=size_palette[1], markeredgecolor='k'),
               Line2D([2], [2], marker='o', color='w', label='[0.4, 0.6)',
                      markerfacecolor=custom_palette[2], markersize=size_palette[2], markeredgecolor='k'),
               Line2D([3], [3], marker='o', color='w', label='[0.8, 1.0]',
                      markerfacecolor=custom_palette[3], markersize=size_palette[3], markeredgecolor='k')]
        ax.legend(handles=legend_elements, loc='best', title="PIP")

    else:
        raise ValueError("Invalid scale for scatter-plot")

    n_rows = len(twas_df.index)
    plt.xticks(np.arange(1,n_rows + 1, 1.0))
    gene_names = twas_df['mol_name'].values
    ax.set_xticklabels(gene_names, rotation=90, ha="center")

    # thresh = 1.0
    # ax.axhline(thresh, linestyle='--', color='red', linewidth=.75, dashes=(5, 7))
    # ax.axhline(-1*thresh, linestyle='--', color='red', linewidth=.75, dashes=(5, 7))

    plt.gcf().subplots_adjust(bottom=0.25)  # make room for xlabel
    plt.ylabel(ylabel, fontsize=18)
    plt.xlabel("")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.tight_layout()
    fig.canvas.draw()

    # save as numpy.ndarray format
    data = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
    scatter_plot = data.reshape(fig.canvas.get_width_height()[::-1] + (3,))

    return scatter_plot


def heatmap(wcor):
    """
    Make a scatterplot of zscore values with gene names as xtick labels.

    :param wcor: numpy.ndarray matrix of sample correlation structure for predicted expression

    :return: numpy.ndarray (RGB) formatted heatmap of correlation structure
    """
    mpl.rcParams["figure.figsize"] = [6.4, 6.4]
    fig = plt.figure()
    fig.subplots_adjust(bottom=0.20, left=0.28)
    mask = np.zeros_like(wcor, dtype=np.bool)
    mask[np.triu_indices_from(mask)] = True
    ax = sns.heatmap(wcor, mask=mask, cmap="RdBu_r", square=True,
                     linewidths=0, cbar=False, xticklabels=False, yticklabels=False, ax=None,
                     vmin=-1, vmax=1)
    ax.margins(2)
    ax.set_aspect('equal', 'box')
    fig.canvas.draw()

    # save image as numpy array
    data = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
    img = data.reshape(fig.canvas.get_width_height()[::-1] + (3,))

    # rotate heatmap to make upside-down triangle shape
    rows, cols, ch = img.shape
    M = cv2.getRotationMatrix2D((cols / 2, rows / 2), 45, 1)
    dst = cv2.warpAffine(img, M, (cols, rows), borderMode=cv2.BORDER_CONSTANT,
                         borderValue=(255, 255, 255))

    # trim extra whitespace
    crop_img = dst[int(dst.shape[0] / 2.5):int(dst.shape[0] / 1.1)]

    return crop_img


def heatmap_colorbar():
    """
    Make a colorbar legend for correlation heatmap for range [-1,1].

    :param: None

    :return: numpy.ndarray (RGB) formatted colorbar for range [-1,1]
    """
    fig = plt.figure(figsize=(3.0, 1.0))
    ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
    norm = mpl.colors.Normalize(vmin=-1, vmax=1)
    mpl.colorbar.ColorbarBase(ax1, cmap="RdBu_r", norm=norm, orientation='horizontal', ticks=[-1, -.50, 0, .50, 1])

    # convert to numpy array format
    fig.canvas.draw()
    data = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
    colorbar = data.reshape(fig.canvas.get_width_height()[::-1] + (3,))

    # reshape and fill to match width of heatmap
    new_size = (colorbar.shape[0], colorbar.shape[1])
    desired_size_w = 640
    desired_size_h = 100
    delta_w = desired_size_w - new_size[1]
    delta_h = desired_size_h - new_size[0]
    top, bottom = delta_h // 2, delta_h - (delta_h // 2)
    left, right = delta_w // 2, delta_w - (delta_w // 2)
    left += 18  # make colorbar line up with heatmap
    right -= 18

    # fill in extra space with whitespace
    color = [255, 255, 255]
    colorbar = cv2.copyMakeBorder(colorbar, top, bottom, left, right, cv2.BORDER_CONSTANT, value=color)

    return colorbar


def focus_plot(wcor, twas_df, scale="logp"):
    """
    Plot zscores and local correlation structure for TWAS fine-mapping.

    :param twas_df: pandas.DataFrame containing at least zscores and gene-names
    :param wcor: numpy.ndarray matrix of sample correlation structure for predicted expression
    :param scale: string the scale to plot on. values are 'logp' or 'zscore'. Default logp

    :return: matplotlib.figure.Figure object containing plot of zscores and local correlation heatmap
    """

    # filter out the null model
    twas_df = twas_df[twas_df["ens_gene_id"] != "NULL.MODEL"]

    scatter_plot = make_scatter(twas_df, scale="logp")
    crop_img = heatmap(wcor)
    colorbar = heatmap_colorbar()

    # combine plots
    numpy_vertical_concat = np.concatenate((scatter_plot, crop_img), axis=0)
    numpy_vertical_concat = np.concatenate((numpy_vertical_concat, colorbar), axis=0)
    numpy_vertical_concat = cv2.resize(numpy_vertical_concat, (0, 0), fx=2.5, fy=2.5)

    fig = plt.figure()
    plt.imshow(numpy_vertical_concat)
    plt.title("")
    plt.axis('off')
    plot_arr = [fig]
    plt.close('all')

    return plot_arr
