import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser(description='Process PLINK --r2 output and plot LD decay.')
    parser.add_argument('--input', required=True, help='Path to the PLINK LD output file (gzipped).')
    parser.add_argument('--bin_size_bp', type=int, default=100, help='Bin size for distance (in base pairs).')
    parser.add_argument('--output_prefix', required=True, help='Prefix for the output plot files.')

    args = parser.parse_args()

    # Read the input file
    ld_data = pd.read_csv(args.input, delim_whitespace=True, compression='gzip')

    # Calculate pairwise physical distance
    ld_data['dist'] = ld_data['BP_B'] - ld_data['BP_A']
    ld_data['dist_kbp'] = ld_data['dist'] / 1000.0  # Convert distance to Kbp

    # Bin distances and calculate mean R2 per bin
    bins = np.arange(0, ld_data['dist'].max() + args.bin_size_bp, args.bin_size_bp)
    ld_data['bin'] = pd.cut(ld_data['dist'], bins=bins, labels=bins[:-1])
    mean_r2_per_bin = ld_data.groupby('bin')['R2'].mean().reset_index()
    mean_r2_per_bin['bin'] = mean_r2_per_bin['bin'].astype(float) / 1000.0  # Convert bin labels to Kbp

    # Plot mean R2 by distance bin in Kbp
    plt.figure(figsize=(10, 6))
    plt.plot(mean_r2_per_bin['bin'], mean_r2_per_bin['R2'], marker='o')
    plt.xlabel('Distance between SNPs (Kbp)')
    plt.ylabel('Mean R2')
    plt.title('LD Decay')
    plt.grid(True)
    plt.savefig(f'{args.output_prefix}_ld_decay.svg')

    # Plot with log10 x-axis in Kbp
    plt.figure(figsize=(10, 6))
    plt.plot(mean_r2_per_bin['bin'], mean_r2_per_bin['R2'], marker='o')
    plt.xscale('log')
    plt.xlabel('Distance between SNPs (Kbp)')
    plt.ylabel('Mean R2')
    plt.title('LD Decay (Log Scale)')
    plt.grid(True)
    plt.savefig(f'{args.output_prefix}_ld_decay_log.svg')

    # Histogram of distances
    plt.figure(figsize=(10, 6))
    plt.hist(ld_data['dist_kbp'], bins=50, edgecolor='black')
    plt.xlabel('Distance between SNPs (Kbp)')
    plt.ylabel('Frequency')
    plt.title('Distribution of Pairwise SNP Distances')
    plt.grid(True)
    plt.savefig(f'{args.output_prefix}_distance_histogram.svg')

if __name__ == '__main__':
    main()

