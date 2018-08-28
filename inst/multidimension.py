#!/usr/bin/env python3

import argparse
import plotting

def main():
    parser = argparse.ArgumentParser(description="This script is for plotting and generating datasets for multi-dimensional clustering")
    parser.add_argument('-l', '--list', nargs='+', required=True, help='a list of filepaths e.g. ./multidimension.py -l path1 path2 path3')
    parser.add_argument('-p', '--plot_dir', required=True, help='plotting directory')
    parser.add_argument('-d', '--data_dir', required=True, help='significant data directory')
    parser.add_argument('-x', '--x_threshold', default=0.05, type=float, help='(adjusted) pvalue for scatter plot x axis')
    parser.add_argument('-y', '--y_threshold', default=0.05, type=float, help='(adjusted) pvalue for scatter plot y axis')
    parser.add_argument('-a', '--adj_pvalue', default=1, type=int, help='whether to use adjusted pvalue or pvalue, 1 as True, 0 as False'



if __name__ == '__main__':
    main()

