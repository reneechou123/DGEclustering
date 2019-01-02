#!/usr/bin/env python3

import argparse
import os
import re
import pandas as pd
import numpy as np
import sqlite3 as sql
import plotting
from matplotlib.colors import ListedColormap
import warnings
import random

def main():
    parser = argparse.ArgumentParser(description="This script is for RNA-seq DGE downstream workflow automation")
    parser.add_argument('-d', '--root_dir', required=True, help='starting directory')
    parser.add_argument('-g', '--gene_col', required=True, help='gene ID column name')
    parser.add_argument('-x', '--x_threshold', default=0.05, type=float, help='(adjusted) pvalue for scatter plot x axis')
    parser.add_argument('-y', '--y_threshold', default=0.05, type=float, help='(adjusted) pvalue for scatter plot y axis')
    parser.add_argument('-a', '--adj_pvalue', default=1, type=int, help='whether to use adjusted pvalue or pvalue. 1 as True, 0 as False')
    parser.add_argument('-q', '--qq_plot', default=1, type=int, help='generate Q-Q plots. 1 as True, 0 as False')
    parser.add_argument('-f', '--fish_plot', default=1, type=int, help='generate fish plots. 1 as True, 0 as False')
    parser.add_argument('-s', '--scatter_plot', default=1, type=int, help='generate scatter plots. 1 as True, 0 as False')
    parser.add_argument('-z', '--dot_size', default=20, type=int, help='size of data points')
    args = parser.parse_args()

    warnings.filterwarnings('ignore') # ignore runtime warnings
    if args.qq_plot == 1:
        os.system('mkdir -p ' + args.root_dir + '/qq_plots')
        os.system('rm -f ' + args.root_dir + '/qq_plots/*')
    if args.fish_plot == 1:
        os.system('mkdir -p ' + args.root_dir + '/fish_plots')
        os.system('rm -f ' + args.root_dir + '/fish_plots/*')
    if args.scatter_plot == 1:
        os.system('mkdir -p ' + args.root_dir + '/scatter_plots')
        os.system('rm -f ' + args.root_dir + '/scatter_plots/*')

    # recursively search for candidate tsv and Rmd files
    os.system('find ' + args.root_dir + ' -name "rnaseq*.Rmd" > rmd_filepaths')
    os.system('find ' + args.root_dir + ' -name "*.tsv" > tsv_filepaths')
    print("searching for candidate Rmd and tsv files, done!")

    # check if Rmd files are rnaseq*.rmd files and cache the rmd and tsv file paths in db
    os.system('echo > ' + args.root_dir + '/rnaseq.db')
    conn = sql.connect(args.root_dir + '/rnaseq.db')
    c = conn.cursor()
    c.execute('''CREATE TABLE rnaseq_rmd_files (
            file_id INTEGER PRIMARY KEY AUTOINCREMENT,
            file_path VARCHAR(255) NOT NULL
            );''')
    c.execute('''CREATE TABLE deseq2_output_files (
            file_id INTEGER PRIMARY KEY AUTOINCREMENT,
            file_path VARCHAR(255) NOT NULL,
            rmd_file_id INT NOT NULL,
            FOREIGN KEY (rmd_file_id) REFERENCES rnaseq_rmd_files(file_id)
            );''')
    conn.commit()

    for line in open('rmd_filepaths'):
        rmd_filepath = line.rstrip()
        dir_path = re.search("(\/).+.Rmd", rmd_filepath).group(1) # assume Rmd and deseq2 output files in same directory
        lib_deseq2 = False
        lib_ggplot2 = False
        deseq2_function = False
        res_list = False

        for line in open(rmd_filepath):
            if res_list == True:
                m = re.search(r"\s+([^=]+)=", line.rstrip())
                if not m and re.search(r"\s+\)", line.rstrip()):
                    break                
                if not m:
                    continue
                c.execute('''INSERT INTO deseq2_output_files (file_path, rmd_file_id) 
                        VALUES (?, (SELECT file_id FROM rnaseq_rmd_files WHERE file_path=?))''', 
                        (dir_path + m.group(1) + ".tsv", rmd_filepath))
            elif (lib_deseq2 == True and lib_ggplot2 == True and deseq2_function == True) and re.search(r"res.list <- list\(", line.rstrip()):
                res_list = True
                c.execute("INSERT INTO rnaseq_rmd_files (file_path) VALUES (?)", (rmd_filepath,))
            elif re.search(r"library\(DESeq2\)", line.rstrip()):
                lib_deseq2 = True
            elif re.search(r"library\(ggplot2\)", line.rstrip()):
                lib_ggplot2 = True
            elif re.search(r"DESeqDataSetFrom", line.rstrip()):
                deseq2_function = True
    
    conn.commit()
    conn.close()
    os.system('rm rmd_filepaths')
    print("checking Rmd files and chaching tsv and Rmd files in db, done!")

    # check if tsv files are deseq2 output files and cache the correct file paths
    deseq2_outputs = list()
    for tsv_filepath in open('tsv_filepaths'):
        for line in open(tsv_filepath.rstrip()):
            if re.search(r"{0}.+{1}.+{2}".format('baseMean', 'log2FoldChange', 'lfcSE'), line.rstrip()):
                deseq2_outputs.append(tsv_filepath.rstrip())
                break
    os.system('rm tsv_filepaths')
    print("checking deseq2 output files, done!")

    # align a deseq2 output file with another one based on column #1
    paired_files = list()
    column_name = args.gene_col
    for index, filepath in enumerate(deseq2_outputs):
        for filepath_2 in deseq2_outputs[index+1:len(deseq2_outputs)]:
            dataset = pd.read_table(filepath, sep='\t')
            dataset_2 = pd.read_table(filepath_2, sep='\t')
            if dataset.shape[0] == dataset_2.shape[0] and column_name in dataset.columns and column_name in dataset_2.columns:
                comparison = pd.DataFrame(columns=['result'])
                comparison['result'] = np.where(dataset.loc[:,column_name] == dataset_2.loc[:,column_name], True, False)
                c_series = pd.Series(comparison['result'].tolist())
                if not 'False' in c_series:    
                    paired_files.append({'file_1': filepath, 'file_2': filepath_2})
    print("pairing files based on gene IDs, done!")

    # generate pairwaise comparison scatter plots
    """
    conn = sql.connect(args.root_dir + '/rnaseq.db')
    c = conn.cursor()
    """

    for paired_file in paired_files:
        # search the source Rmd files in db
        """
        c.execute('''SELECT r.file_path FROM rnaseq_rmd_files AS r JOIN deseq2_output_files AS d ON r.file_id = d.rmd_file_id 
                WHERE d.file_path=? OR d.file_path=?''', (paired_file['file_1'], paired_file['file_2']))
        result = c.fetchall()
        if not len(result) == 2:
            continue
        rmd_path = result[0][0]
        rmd_path_2 = result[1][0]
        
        conn.commit()
        conn.close()
        """

        # emit plots and diagnostics
        if args.qq_plot == 1:
            plotting.qq_plot(output_dir=args.root_dir+'/qq_plots', file_path=paired_file['file_1'])
            plotting.qq_plot(output_dir=args.root_dir+'/qq_plots', file_path=paired_file['file_2'])
        if args.fish_plot == 1:
            plotting.fish_plot(paired_file['file_1'], paired_file['file_2'], gene_col=args.gene_col, output_dir=args.root_dir+'/fish_plots')
        if args.scatter_plot == 1:
            file_paths = [paired_file['file_1'], paired_file['file_2']]
            temp = plotting.scatter_plot(file_paths, gene_col=args.gene_col, out_dir=args.root_dir+'/scatter_plots', x_threshold=args.x_threshold, y_threshold=args.y_threshold, adj_pvalue=args.adj_pvalue, dot_size=args.dot_size)

    # generate null Q-Q plot
    if args.qq_plot == 1:
        random.seed(123)
        dataset = pd.DataFrame(data=np.random.uniform(low=0, high=1, size=17000), columns=['pvalue']) 
        plotting.qq_plot(output_dir=args.root_dir+'/qq_plots', dataset=dataset)

if __name__ == '__main__':
    main()

