import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import seaborn as sns
import pandas as pd
import logging
import pdb
import cv2

__all__ = ["focus_plot"]


def zscore_scatter(twas_df):
    """
    Make a scatterplot of zscore values with gene names as xtick labels.

    :param twas_df: pandas.DataFrame containing at least zscores and gene-names

    :return: numpy.ndarray (RGB) formatted scatterplot of zscores
    """
    mpl.rcParams["figure.figsize"] = [6.4, 4.8]
    fig = plt.figure()

    ax = sns.stripplot(x="ens_gene_id", y="twas_z", data=twas_df, color='black', size=5)

    gene_names = twas_df.sort_index(0)['ens_gene_id'].values
    N = len(gene_names)
    ax.set_xticklabels(gene_names,rotation=90, ha="right")

    # make ymin,ymax symmetric
    zmin = np.min(twas_df['twas_z'])
    zmax = np.max(twas_df['twas_z'])
    abs_bound = round(max(np.abs(zmin), np.abs(zmax)))

    if abs_bound % 2 != 0:
        abs_bound += 1

    ax.set_ylim([abs_bound*-1, abs_bound])

    #thresh = 1.0
    #ax.axhline(thresh, linestyle='--', color='red', linewidth=.75, dashes=(5, 7))
    #ax.axhline(-1*thresh, linestyle='--', color='red', linewidth=.75, dashes=(5, 7))
    plt.gcf().subplots_adjust(bottom=0.25) # make room for xlabel
    plt.ylabel('Z-score', fontsize=18)
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
    ax=sns.heatmap(wcor, mask=mask, cmap="RdBu_r", square=True,
                linewidths=0, cbar=False, xticklabels=False, yticklabels=False, ax=None,
                vmin=-1, vmax=1)
    ax.margins(2)
    ax.set_aspect('equal', 'box')
    fig.canvas.draw()

    # save image as numpy array
    data = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
    img = data.reshape(fig.canvas.get_width_height()[::-1] + (3,))

    # rotate heatmap to make upside-down triangle shape
    rows,cols, ch = img.shape
    M = cv2.getRotationMatrix2D((cols/2,rows/2),45,1)
    dst = cv2.warpAffine(img,M,(cols,rows), borderMode=cv2.BORDER_CONSTANT,
                               borderValue=(255,255,255))

    # trim extra whitespace
    crop_img = dst[int(dst.shape[0]/2.5):int(dst.shape[0]/1.1)]

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
    mpl.colorbar.ColorbarBase(ax1, cmap="RdBu_r", norm=norm, orientation='horizontal', ticks=[-1,-.50,0,.50, 1])

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
    top, bottom = delta_h//2, delta_h-(delta_h//2)
    left, right = delta_w//2, delta_w-(delta_w//2)
    left += 18 # make colorbar line up with heatmap
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
    scatter_plot = zscore_scatter(twas_df)
    crop_img = heatmap(wcor)
    colorbar = heatmap_colorbar()

    # combine plots
    numpy_vertical_concat = np.concatenate((scatter_plot, crop_img), axis=0)
    numpy_vertical_concat = np.concatenate((numpy_vertical_concat, colorbar), axis=0)
    numpy_vertical_concat = cv2.resize(numpy_vertical_concat, (0,0), fx=2.5, fy=2.5)

    fig = plt.figure()
    plt.imshow(numpy_vertical_concat)
    plt.title("")
    plt.axis('off')
    plot_arr = [fig]
    plt.close('all')
    return plot_arr


def credible_set(twas_df, thresh=.90):
    """
    Plot zscores and local correlation structure for TWAS fine-mapping.

    :param twas_df: pandas.DataFrame containing at least pip and gene-names
    :param thresh: credible set threshold (0,1]; default=.90

    :return: list of gene names in specified thresh credible set
    """
    sorted_twas_df = twas_df.sort_values(by="pip", ascending=False)
    credible_gene_set = []
    running_sum = 0

    for index, row in sorted_twas_df.iterrows():
        if running_sum <= thresh:
            running_sum += row['pip']
            credible_gene_set.append(row['ens_gene_id'])
        else:
            break

    return credible_gene_set
