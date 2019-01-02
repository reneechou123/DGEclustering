#!/usr/bin/env python3

import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from scipy.stats.mstats import mquantiles
from scipy.stats import beta

def scatter_plot(file_paths, gene_col, x_file_number=0, y_file_number=1, out_dir='./', x_threshold=0.05, y_threshold=0.05, adj_pvalue=True, for_cluster_plot=False, return_sig_plot=False, out_file_name=None, dot_size=20):
    # check input file names for regex search
    for index, file_path in enumerate(file_paths):
        if re.search(r".+\/(.+).tsv", file_path) == None:
            file_paths[index] = './' + file_path

    datasets = list()
    for file_path in file_paths:
        datasets.append(pd.read_table(file_path))
    keys = np.arange(len(datasets)).astype(str)
    merged_set = pd.concat([x.set_index(gene_col) for x in datasets], axis=1, keys=keys, join='inner', ignore_index=False)
    merged_set.columns = merged_set.columns.map('_'.join)
    
    # create output filepath
    filename_1 = re.search(r".+\/(.+).tsv", file_paths[x_file_number]).group(1)
    filename_2 = re.search(r".+\/(.+).tsv", file_paths[y_file_number]).group(1)
    if out_file_name == None:
        out = out_dir + '/' + filename_1 + '_vs_' + filename_2
    else:
        out = out_dir + '/' + out_file_name

    # create subsets
    if adj_pvalue == True:
        padj_x = str(x_file_number) + '_padj'
        padj_y = str(y_file_number) + '_padj'
        sig_vs_sig = merged_set[(merged_set[padj_x] < x_threshold) & (merged_set[padj_y] < y_threshold)]
        sig_vs_NS = merged_set[(merged_set[padj_x] < x_threshold) & (merged_set[padj_y] >= y_threshold)]
        NS_vs_sig = merged_set[(merged_set[padj_x] >= x_threshold) & (merged_set[padj_y] < y_threshold)]
        NS_vs_NS = merged_set[(merged_set[padj_x] >= x_threshold) & (merged_set[padj_y] >= y_threshold)]
        non_NA_set = merged_set[(merged_set.isna()[padj_x] == False) & (merged_set.isna()[padj_y] == False)]
    
    else:
        pvalue_x = str(x_file_number) + '_pvalue'
        pvalue_y = str(y_file_number) + '_pvalue'
        sig_vs_sig = merged_set[(merged_set[pvalue_x] < x_threshold) & (merged_set[pvalue_y] < y_threshold)]
        sig_vs_NS = merged_set[(merged_set[pvalue_x] < x_threshold) & (merged_set[pvalue_y] >= y_threshold)]
        NS_vs_sig = merged_set[(merged_set[pvalue_x] >= x_threshold) & (merged_set[pvalue_y] < y_threshold)]
        NS_vs_NS = merged_set[(merged_set[pvalue_x] >= x_threshold) & (merged_set[pvalue_y] >= y_threshold)]
        non_NA_set = merged_set[(merged_set.isna()[pvalue_x] == False) & (merged_set.isna()[pvalue_y] == False)]  

    log2FoldChange_x = str(x_file_number) + '_log2FoldChange'
    log2FoldChange_y = str(y_file_number) + '_log2FoldChange'
    
    # -- plotting section --
    xtitle = filename_1.replace('_', ' ').replace('.', ' ')
    ytitle = filename_2.replace('_', ' ').replace('.', ' ')

    # general scatter plot (for the two specified files only)
    fig = plt.figure(figsize=(18, 18))
    ax = fig.add_subplot(111)
    if for_cluster_plot == False:
        plt.close()
        fig = plt.figure(figsize=(18, 18))
        ax = fig.add_subplot(111)
        g4 = ax.scatter(NS_vs_NS[log2FoldChange_x], NS_vs_NS[log2FoldChange_y], s=dot_size, c='grey', alpha=0.3)
        g3 = ax.scatter(NS_vs_sig[log2FoldChange_x], NS_vs_sig[log2FoldChange_y], s=dot_size, c=(31 / 255., 119 / 255., 180 / 255.), alpha=0.6)
        g2 = ax.scatter(sig_vs_NS[log2FoldChange_x], sig_vs_NS[log2FoldChange_y], s=dot_size, c=(255 / 255., 127 / 255., 14 / 255.), alpha=0.6)
        g1 = ax.scatter(sig_vs_sig[log2FoldChange_x], sig_vs_sig[log2FoldChange_y], s=dot_size+8, c=(214 / 255., 39 / 255., 40 / 255.), alpha=1.0)

        ax.legend((g1, g2, g3, g4),(
            'sig vs sig (' + str(sig_vs_sig.shape[0]) + ')',
            'sig vs NS (' + str(sig_vs_NS.shape[0]) + ')',
            'NS vs sig (' + str(NS_vs_sig.shape[0]) + ')',
            'NS vs NS (' + str(NS_vs_NS.shape[0]) + ')'),
            markerscale=2, prop={'size': 26})
        ax.set_xlim(min(non_NA_set[log2FoldChange_x].min(), non_NA_set[log2FoldChange_y].min()) - 0.5, 
                max(non_NA_set[log2FoldChange_x].max(), non_NA_set[log2FoldChange_y].max()) + 0.5)
        ax.set_ylim(min(non_NA_set[log2FoldChange_x].min(), non_NA_set[log2FoldChange_y].min()) - 0.5, 
                max(non_NA_set[log2FoldChange_x].max(), non_NA_set[log2FoldChange_y].max()) + 0.5)
        ax.axvline(x=0, linestyle='dotted', color='grey')
        ax.axhline(y=0, linestyle='dotted', color='grey')

        title = '(' + xtitle + ') vs (' + ytitle + ')\n(gene number=' + str(merged_set.shape[0]) + ')'
        sig_discordant = sig_vs_sig[((sig_vs_sig[log2FoldChange_x] < 0) & (sig_vs_sig[log2FoldChange_y] > 0)) |
                   ((sig_vs_sig[log2FoldChange_x] > 0) & (sig_vs_sig[log2FoldChange_y] < 0))]
        anchored_text = AnchoredText('# of sig vs sig in II and IV: ' + str(sig_discordant.shape[0]), loc=3, prop={'size': 26})
        anchored_text.patch.set(color='red', alpha=0.3)

        ax.set_title(title, fontweight='bold', fontsize=36, y=1.02)
        ax.set_xlabel(xtitle + u' (log\u2082 fold change)', fontsize=35)
        ax.set_ylabel(ytitle + u' (log\u2082 fold change)', fontsize=35)
        ax.add_artist(anchored_text)
        ax.tick_params(labelsize=22)

        fig.savefig(out + '_scatter_plot.png')
    
    # significant scatter plot
    if return_sig_plot == True:
        # prepare output significant dataset(s) and set plotting parameters
        if len(datasets) == 2:
            sig_discordant = sig_vs_sig[((sig_vs_sig[log2FoldChange_x] < 0) & (sig_vs_sig[log2FoldChange_y] > 0)) |
                      ((sig_vs_sig[log2FoldChange_x] > 0) & (sig_vs_sig[log2FoldChange_y] < 0))]
            sig_concordant = sig_vs_sig[((sig_vs_sig[log2FoldChange_x] >= 0) & (sig_vs_sig[log2FoldChange_y] >=0)) | 
                      ((sig_vs_sig[log2FoldChange_x] <= 0) & (sig_vs_sig[log2FoldChange_y] <= 0))]
            if sig_discordant.shape[0] > 0:
                sig_discordant.to_csv(out + '_disagreeing_genes.tsv', sep='\t', index=True)
            if sig_concordant.shape[0] > 0:
                sig_concordant.to_csv(out + '_agreeing_genes.tsv', sep='\t', index=True)
                
            title = '(' + xtitle + ') vs (' + ytitle + ')\n(gene number=' + str(merged_set.shape[0]) + ')'
            anchored_text = AnchoredText('# of sig vs sig in II and IV: ' + str(sig_discordant.shape[0]), loc=3, prop={'size': 26})
            anchored_text.patch.set(color='red', alpha=0.3)

        else: # the significant threshold for multiple files will be the smaller one between x_ and y_threshold
            if adj_pvalue == True:
                temp = pd.concat((sig_vs_sig['%i_padj'%i] < min(x_threshold, y_threshold) for i in np.arange(len(datasets))), axis=1).all(axis=1)
            else:
                temp = pd.concat((sig_vs_sig['%i_pvalue'%i] < min(x_threshold, y_threshold) for i in np.arange(len(datasets))), axis=1).all(axis=1)
            all_sig = sig_vs_sig[temp]
            if all_sig.shape[0] > 0:
                all_sig.to_csv(out + '_all_sig_genes.tsv', sep='\t', index=True)
            
            title = '(' + xtitle + ') vs (' + ytitle + ') (multi-dimensional)'
            anchored_text = None

        plt.close()
        fig = plt.figure(figsize=(18, 18))
        ax = fig.add_subplot(111)
       
        g2 = ax.scatter(non_NA_set[log2FoldChange_x], non_NA_set[log2FoldChange_y], s=dot_size, c='grey', alpha=0.3)
        if len(datasets) == 2:
            g1 = ax.scatter(sig_vs_sig[log2FoldChange_x], sig_vs_sig[log2FoldChange_y], s=dot_size+8, c=(214 / 255., 39 / 255., 40 / 255.), alpha=1.0)
            ax.legend((g1,), ('sig vs sig (' + str(sig_vs_sig.shape[0]) + ')',), markerscale=2, prop={'size': 26})
        else:
            g1 = ax.scatter(all_sig[log2FoldChange_x], all_sig[log2FoldChange_y], s=dot_size, c=(214 / 255., 39 / 255., 40 / 255.), alpha=1.0)
            ax.legend((g1,), ('all sig (' + str(all_sig.shape[0]) + ')',), markerscale=2, prop={'size': 26})    
        
        ax.set_xlim(min(non_NA_set[log2FoldChange_x].min(), non_NA_set[log2FoldChange_y].min()) - 0.5,
                max(non_NA_set[log2FoldChange_x].max(), non_NA_set[log2FoldChange_y].max()) + 0.5)
        ax.set_ylim(min(non_NA_set[log2FoldChange_x].min(), non_NA_set[log2FoldChange_y].min()) - 0.5,
                max(non_NA_set[log2FoldChange_x].max(), non_NA_set[log2FoldChange_y].max()) + 0.5)
        ax.axvline(x=0, linestyle='dotted', color='grey')
        ax.axhline(y=0, linestyle='dotted', color='grey')
        
        ax.set_title(title, fontweight='bold', fontsize=36, y=1.02)
        ax.set_xlabel(xtitle + u' (log\u2082 fold change)', fontsize=35)
        ax.set_ylabel(ytitle + u' (log\u2082 fold change)', fontsize=35)
        if len(datasets) == 2:
            ax.add_artist(anchored_text)
        ax.tick_params(labelsize=22)
        
        fig.savefig(out + '_sig_plot.png')

    # plotting frame for cluster plot
    if for_cluster_plot == True:
        if len(datasets) == 2:
            title = '(' + xtitle + ') vs (' + ytitle + ')'
        else:
            title = '(' + xtitle + ') vs (' + ytitle + ') (multi-dimensional)'
        
        plt.close()
        fig = plt.figure(figsize=(18, 18))
        ax = fig.add_subplot(111)
        ax.set_xlim(min(non_NA_set[log2FoldChange_x].min(), non_NA_set[log2FoldChange_y].min()) - 0.5,
                max(non_NA_set[log2FoldChange_x].max(), non_NA_set[log2FoldChange_y].max()) + 0.5)
        ax.set_ylim(min(non_NA_set[log2FoldChange_x].min(), non_NA_set[log2FoldChange_y].min()) - 0.5,
                max(non_NA_set[log2FoldChange_x].max(), non_NA_set[log2FoldChange_y].max()) + 0.5)
        ax.axvline(x=0, linestyle='dotted', color='grey')
        ax.axhline(y=0, linestyle='dotted', color='grey')
        ax.set_title(title, fontweight='bold', fontsize=36, y=1.02)
        ax.set_xlabel(xtitle + u' (log\u2082 fold change)', fontsize=35)
        ax.set_ylabel(ytitle + u' (log\u2082 fold change)', fontsize=35)
        ax.tick_params(labelsize=22)

        return ax


def fish_plot(file_path_1, file_path_2, gene_col, output_dir):
    # check file names
    if re.search(r".+\/(.+).tsv", file_path_1) == None:
        file_path_1 = './' + file_path_1
    if re.search(r".+\/(.+).tsv", file_path_2) == None:
        file_path_2 = './' + file_path_2 
    
    dataset = pd.read_table(file_path_1)
    dataset_2 = pd.read_table(file_path_2)
    merged_set = dataset.merge(dataset_2, left_on=gene_col, right_on=gene_col)
    merged_set['-log10_pvalue_x'] = - np.sign(merged_set['log2FoldChange_x']) * np.sign(merged_set['log2FoldChange_y']) * np.log10(merged_set['pvalue_x'])
    merged_set['-log10_pvalue_y'] = - np.sign(merged_set['log2FoldChange_x']) * np.sign(merged_set['log2FoldChange_y']) * np.log10(merged_set['pvalue_y'])

    plt.close()
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    ax.scatter(merged_set['-log10_pvalue_x'], merged_set['-log10_pvalue_y'], s=3, c=(31 / 255., 119 / 255., 180 / 255.), alpha=0.5)

    filename_1 = re.search(r".+\/(.+).tsv", file_path_1).group(1)
    filename_2 = re.search(r".+\/(.+).tsv", file_path_2).group(1)
    xtitle = filename_1.replace('_', ' ')
    xtitle = xtitle.replace('.', ' ')
    ytitle = filename_2.replace('_', ' ')
    ytitle = ytitle.replace('.', ' ')
    ax.set_xlim(-40, 40)
    ax.set_ylim(-40, 40)
    ax.axvline(x=0, linestyle='dotted', color='grey')
    ax.axhline(y=0, linestyle='dotted', color='grey')
    ax.set_title('(' + xtitle + ') vs (' + ytitle + ')\n(gene number=' + str(merged_set.shape[0]) + ')', fontweight='bold', fontsize=24, y=1.02)
    ax.set_xlabel(xtitle + u' -log\u2081\u2080 pvalue', fontsize=22)
    ax.set_ylabel(ytitle + u' -log\u2081\u2080 pvalue', fontsize=22)
    ax.tick_params(labelsize=12)

    fig.savefig(output_dir + '/' + filename_1 + '_vs_' + filename_2 + '_fish_plot.png')


def qq_plot(output_dir, file_path=None, dataset=None):    
    
    if dataset is None:
        dataset = pd.read_table(file_path)
    dataset = dataset[dataset.isna()['pvalue'] == False]
    dataset['-log10_pvalue'] = - np.log10(dataset['pvalue'])
    gene_size = dataset.shape[0]

    exp = np.concatenate([np.arange(100) / gene_size, np.logspace(-np.log10(gene_size) + 2, 0, 200)])
    obs = mquantiles(dataset['pvalue'], prob=exp, alphap=0, betap=1)

    lower = list()
    upper = list()
    for i in range(0, len(exp)):
        CI_values = beta.interval(0.95, gene_size * exp[i], gene_size - gene_size * exp[i])
        lower.append(CI_values[0])
        upper.append(CI_values[1])

    exp = -np.log10(exp)
    obs = -np.log10(obs)
    up = -np.log10(lower)
    low = -np.log10(upper)

    plt.close()
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    ax.fill_between(exp, up, low, color='grey', alpha=0.5)
    ax.set_xlim(np.nanmin(exp[exp != -np.inf]), np.nanmax(exp[exp != np.inf]) + 0.1)
    ax.set_ylim(np.nanmin(obs[obs != -np.inf]), max(np.nanmax(obs[obs != np.inf]), np.nanmax(up[up != np.inf])) + 0.5)
    ax.plot(ax.get_xlim(), ax.get_xlim(), linestyle='--', color='black')
    ax.scatter(exp, obs, s=3, c=(31 / 255., 119 / 255., 180 / 255.))

    if file_path is not None:
        # check file names
        if re.search(r".+\/(.+).tsv", file_path) == None:
            file_path = './' + file_path
        filename = re.search(r".+\/(.+).tsv", file_path).group(1)
        title = filename.replace('_', ' ')
        title = title.replace('.', ' ')
    else:
        filename = 'null'
        title = 'Null'
    ax.set_title(title + ' QQ-Plot', fontweight='bold', fontsize=24, y=1.02)
    ax.set_xlabel('expected -log\u2081\u2080 pvalue', fontsize=22)
    ax.set_ylabel('observed -log\u2081\u2080 pvalue', fontsize=22)
    ax.tick_params(labelsize=12)

    fig.savefig(output_dir + '/' + filename + '_qq_plot.png')

