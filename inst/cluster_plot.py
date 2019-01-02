#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import plotting
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Ellipse
from scipy.stats import linregress
import math
import warnings

def main():
    parser = argparse.ArgumentParser(description="This script is for clustering visulization of RNA-seq paired DGE files")
    parser.add_argument('-d', '--out_dir', required=True, help='output directory for resulting plots')
    parser.add_argument('-f', '--file_paths', nargs='+', help='a list of filepaths. e.g. ./multidimension.py -l path1 path2 path3')
    parser.add_argument('-m', '--MCA_result', default=0, type=int, help='whether or not the input is a MCA result. 1 as True, 0 as False')
    parser.add_argument('-s', '--subplots', default=0, type=int, help='whether to generate subplots. 1 as True, 0 as False')
    parser.add_argument('-n1', '--x_file_number', default=0, type=int, help='file number in list for x axis. index starts at 0')
    parser.add_argument('-n2', '--y_file_number', default=1, type=int, help='file number in list for y axis. index starts at 0')
    parser.add_argument('-g', '--gene_col', required=True, type=str, help='gene ID column name') 
    parser.add_argument('-r', '--clustering_result', required=True, help='file path of the clustering result')
    parser.add_argument('-a', '--adj_pvalue', default=1, type=int, help='whether to use adjusted pvalue or pvalue. 1 as True, 0 as False')
    parser.add_argument('-c', '--color', default='brg', help='cmap color for visualization')
    parser.add_argument('-z', '--dot_size', default=20, type=int, help='size of data points')
    args = parser.parse_args()
    
    warnings.filterwarnings('ignore') # ignore runtime warnings
    
    clustered_dat = pd.read_table(args.clustering_result)   
    
    if args.MCA_result == 0:
        dataset_1 = pd.read_table(args.file_paths[args.x_file_number])
        dataset_2 = pd.read_table(args.file_paths[args.y_file_number])
        merged_set = dataset_1.merge(dataset_2, left_on=args.gene_col, right_on=args.gene_col)
        merged_set = merged_set.merge(clustered_dat, left_on=args.gene_col, right_on=clustered_dat.columns[0], how='inner')
        x_axis = 'log2FoldChange_x'
        y_axis = 'log2FoldChange_y'
    else:
        MCA = pd.read_table(args.file_paths[0])
        merged_set = MCA.merge(clustered_dat, left_on=MCA.columns[0], right_on=clustered_dat.columns[0], how='inner')
        x_axis = 'V1'
        y_axis = 'V2'

    # split dataset
    i = 1
    group_list = list()
    while ('Group.' + str(i)) in merged_set['ind'].tolist():
        group_list.append(merged_set[merged_set['ind'] == 'Group.' + str(i)].reset_index(drop=True))
        i += 1   
    cmap = cm.get_cmap(args.color)(np.linspace(0, 1.0, len(group_list)))

    # -- main plot --
    # prepare plotting axes
    plt.close()
    fig = plt.figure(figsize=(18, 18))
    ax = fig.add_subplot(111)
    if args.MCA_result == 0:
        ax = plotting.scatter_plot(file_paths=args.file_paths, gene_col=args.gene_col, x_file_number=args.x_file_number, y_file_number=args.y_file_number, adj_pvalue=args.adj_pvalue, for_cluster_plot=True, dot_size=args.dot_size)
    else:
        title = 'MCA plot'
        xtitle = 'Dim1'
        ytitle = 'Dim2'
        # ax.axvline(x=0, linestyle='dotted', color='grey')
        # ax.axhline(y=0, linestyle='dotted', color='grey')
        ax.set_title(title, fontweight='bold', fontsize=36, y=1.02)
        ax.set_xlabel(xtitle, fontsize=35)
        ax.set_ylabel(ytitle, fontsize=35)
        ax.tick_params(labelsize=22)
    
    # plotting
    groups = list()
    group_names = list()  
    for index, group in enumerate(group_list):
        color = cmap[index]
        x = group[x_axis]
        y = group[y_axis]

        # plot the centroid for each cluster
        cx = np.mean(x)
        cy = np.mean(y)
        ax.scatter(cx, cy, s=12, facecolors='none', edgecolors='black', alpha=0.5)
        ax.annotate(str(index + 1), xy=(cx, cy), xytext=(cx + 0.01, cy + 0.01))

        # plot data points
        groups.append(ax.scatter(x, y, s=args.dot_size, c=([color,] * group.shape[0]), alpha=0.6))
        group_names.append('cluster ' + str(index + 1) + ' (' + str(group.shape[0]) + ')')

        # plot lines from centroid to points
        for row in range(0, group.shape[0]):
            ax.plot([cx, x[row]], [cy, y[row]], 
                    linestyle='dashed', linewidth=0.5, c=color, alpha=0.3)
        
        # encircle the group points
        slope, intercept, r_value, p_value, std_err = linregress(x,y)
        angle = math.degrees(math.atan(slope))
        if np.std(x) < np.std(y):
                    angle = angle + 90
        ell = Ellipse(xy=[cx, cy], width=np.std(x) * 5, height=np.std(y) * 5, angle=angle, ec='k', fc=color, alpha=0.2)
        ax.add_patch(ell)
    
    ax.legend(groups, group_names, markerscale=3, prop={'size': 26})
    plt.savefig(args.out_dir + '/cluster_all.png')

    if args.subplots == 0:
        return
    # -- subplots --
    for j in range(0, len(group_list)):
        # prepare plotting axes
        plt.close()
        fig = plt.figure(figsize=(18, 18))
        ax = fig.add_subplot(111)
        if args.MCA_result == 0:
            ax = plotting.scatter_plot(file_paths=args.file_paths, gene_col=args.gene_col, x_file_number=args.x_file_number, y_file_number=args.y_file_number, adj_pvalue=args.adj_pvalue, for_cluster_plot=True, dot_size=args.dot_size)
        else:
            title = 'MCA plot'
            xtitle = 'Dim1'
            ytitle = 'Dim2'
            ax.set_title(title, fontweight='bold', fontsize=36, y=1.02)
            ax.set_xlabel(xtitle, fontsize=35)
            ax.set_ylabel(ytitle, fontsize=35)
            ax.tick_params(labelsize=22)
        
        # plotting
        for index, group in enumerate(group_list):
            color = cmap[index]
            x = group[x_axis]
            y = group[y_axis]

            if index != j:
                # plot the centroid for each cluster
                cx = np.mean(x)
                cy = np.mean(y)
                ax.scatter(cx, cy, s=8, facecolors='none', edgecolors='grey', alpha=0.3)
                ax.annotate(str(index + 1), xy=(cx, cy), xytext=(cx + 0.01, cy + 0.01), color='grey')

                # plot lines from centroid to points
                for row in range(0, group.shape[0]):
                    ax.plot([cx, x[row]], [cy, y[row]],
                            linestyle='dashed', linewidth=0.5, c='grey', alpha=0.3)
                 
                # plot data points
                ax.scatter(x, y, s=args.dot_size, c=('grey'), alpha=0.3)

                # encircle the group points
                slope, intercept, r_value, p_value, std_err = linregress(x,y)
                angle = math.degrees(math.atan(slope))
                if np.std(x) < np.std(y):
                    angle = angle + 90
                ell = Ellipse(xy=[cx, cy], width=np.std(x) * 5, height=np.std(y) * 5,
                        angle=angle, ec='#7f7f7f', fc='grey', alpha=0.1)
                ax.add_patch(ell)
        
        # plot the centroid for cluster j
        color = cmap[j]
        x = group_list[j][x_axis]
        y = group_list[j][y_axis]
        
        cx = np.mean(x)
        cy = np.mean(y)
        ax.scatter(cx, cy, s=12, facecolors='none', edgecolors='grey')
        ax.annotate(str(j + 1), xy=(cx, cy), xytext=(cx + 0.01, cy + 0.01))
                
        # plot lines from centroid to cluster j points
        for row in range(0, group_list[j].shape[0]):
            ax.plot([cx, x[row]], [cy, y[row]], linestyle='dashed', linewidth=0.5, c=color)
        
        # plot data points
        g1 = ax.scatter(x, y, s=args.dot_size, c=([color,] * group.shape[0]), edgecolors='black', linewidth=0.8)
        g1_name = 'cluster ' + str(j + 1) + ' (' + str(group_list[j].shape[0]) + ')'
                 
        # encircle the cluster j group points
        slope, intercept, r_value, p_value, std_err = linregress(x,y)
        angle = math.degrees(math.atan(slope))
        if np.std(x) < np.std(y):
             angle = angle + 90 
        ell = Ellipse(xy=[cx, cy], width=np.std(x) * 5, height=np.std(y) * 5, angle=angle, ec='k', fc=color, alpha=0.2)
        ax.add_patch(ell)

        ax.legend((g1,), (g1_name,), markerscale=3, prop={'size': 26})
        plt.savefig(args.out_dir + '/cluster_' + str(j + 1) + '.png')

if __name__ == '__main__':
    main()

