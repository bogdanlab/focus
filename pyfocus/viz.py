import matplotlib

matplotlib.use("Agg")  # forces matplotlib to not launch X11 window
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import seaborn as sns
import pandas as pd
import cv2

from scipy import stats
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

__all__ = ["focus_plot"]


def make_scatter(twas_df):
    """
    Make a scatterplot of zscore values with gene names as xtick labels.

    :param twas_df: pandas.DataFrame containing at least zscores and gene-names

    :return: numpy.ndarray (RGB) formatted scatterplot of zscores
    """
    mpl.rcParams["figure.figsize"] = [6.4, 4.8]
    fig, ax = plt.subplots()

    size_arr = []
    color_arr = []
    custom_palette = ["#e4f1fe", "#89c4f4", "#2574a9", "#013243"]
    size_palette = [4, 8, 10, 12]

    for i, row in twas_df.iterrows():
        pip = row["pip"]
        if pip < 0.2:
            color_arr.append(custom_palette[0])
            size_arr.append(size_palette[0])
        elif 0.2 <= pip < 0.40:
            color_arr.append(custom_palette[1])
            size_arr.append(size_palette[1])
        elif 0.40 <= pip < 0.60:
            color_arr.append(custom_palette[2])
            size_arr.append(size_palette[2])
        elif 0.60 <= pip < 0.80:
            color_arr.append(custom_palette[3])
            size_arr.append(size_palette[3])
        else:
            color_arr.append(custom_palette[3])
            size_arr.append(size_palette[3])

    n_rows = len(twas_df.index)
    x_values = np.arange(1, n_rows + 1)
    ax.scatter(x=x_values, y=twas_df["logp"].values, s=size_arr, c=color_arr, edgecolor="black")

    # create legend
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", label="[0.0, 0.2)",
               markerfacecolor=custom_palette[0], markersize=size_palette[0], markeredgecolor="k"),
        Line2D([1], [1], marker="o", color="w", label="[0.2, 0.4)",
               markerfacecolor=custom_palette[1], markersize=size_palette[1], markeredgecolor="k"),
        Line2D([2], [2], marker="o", color="w", label="[0.4, 0.6)",
               markerfacecolor=custom_palette[2], markersize=size_palette[2], markeredgecolor="k"),
        Line2D([3], [3], marker="o", color="w", label="[0.8, 1.0]",
               markerfacecolor=custom_palette[3], markersize=size_palette[3], markeredgecolor="k")]
    plt.legend(handles=legend_elements, loc="best", title="PIP")

    n_rows = len(twas_df.index)
    gene_names = twas_df["mol_name"].values

    plt.xticks(np.arange(1, n_rows + 1, 1.0))
    plt.xticks(x_values, labels=gene_names, rotation="vertical")
    plt.subplots_adjust(bottom=0.25)  # make room for xlabel

    plt.ylabel("$-\log_{10}(p)$", fontsize=18)
    plt.xlabel("")

    # drop right/top axis bars
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    fig.tight_layout()
    fig.canvas.draw()

    # save as numpy.ndarray format
    data = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep="")
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
    ax.set_aspect("equal", "box")
    fig.canvas.draw()

    # save image as numpy array
    data = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep="")
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
    mpl.colorbar.ColorbarBase(ax1, cmap="RdBu_r", norm=norm, orientation="horizontal", ticks=[-1, -.50, 0, .50, 1])

    # convert to numpy array format
    fig.canvas.draw()
    data = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep="")
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


def focus_plot(wcor, twas_df):
    """
    Plot zscores and local correlation structure for TWAS fine-mapping.

    :param twas_df: pandas.DataFrame containing at least zscores and gene-names
    :param wcor: numpy.ndarray matrix of sample correlation structure for predicted expression

    :return: matplotlib.figure.Figure object containing plot of zscores and local correlation heatmap
    """

    # filter out the null model
    twas_df = twas_df[twas_df["ens_gene_id"] != "NULL.MODEL"]

    # add p-value to compute -log10 p
    twas_df = twas_df.assign(logp=-stats.chi2.logsf(twas_df["twas_z"].values ** 2, 1))

    scatter_plot = make_scatter(twas_df)
    crop_img = heatmap(wcor)
    colorbar = heatmap_colorbar()

    # combine plots
    numpy_vertical_concat = np.concatenate((scatter_plot, crop_img), axis=0)
    numpy_vertical_concat = np.concatenate((numpy_vertical_concat, colorbar), axis=0)
    numpy_vertical_concat = cv2.resize(numpy_vertical_concat, (0, 0), fx=2.5, fy=2.5)

    fig = plt.figure()
    plt.imshow(numpy_vertical_concat)
    plt.title("")
    plt.axis("off")
    plot_arr = [fig]
    plt.close("all")

    return plot_arr
