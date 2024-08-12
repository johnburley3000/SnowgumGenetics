import pandas as pd
import matplotlib.pyplot as plt
import gzip
import os
import argparse

def plot_ld_decay(data_dir, bin_size, output_dir):
    files = [f for f in os.listdir(data_dir) if f.endswith('_LD.stat.gz')]
    files = sorted(files, key=lambda x: int(x.split('_')[2][3:]))
    ld_data = {}

    for file in files:
        chromosome = file.split('_')[2]
        print(chromosome)
        file_path = os.path.join(data_dir, file)
        with gzip.open(file_path, 'rt') as f:
            # Reading the file and skipping comment lines
            df = pd.read_csv(f, sep='\t')
            df.columns = ['Dist', 'Mean_r^2', 'Mean_D\'', 'Sum_r^2', 'Sum_D\'', 'NumberPairs']
            df['Dist_kbp'] = df['Dist'] * bin_size / 1000  # Convert distance to Kbp
            ld_data[chromosome] = df

    # Plotting
    plt.figure(figsize=(14, 7))

    # Plot with linear x-axis
    plt.subplot(1, 2, 1)
    for chromosome, df in ld_data.items():
        plt.plot(df['Dist_kbp'], df['Mean_r^2'], label=chromosome)
    plt.xlabel('Distance (Kbp)')
    plt.ylabel('Mean r^2')
    plt.title('LD Decay (Linear Scale)')
    plt.legend()
    plt.grid(True)

    # Plot with log x-axis
    plt.subplot(1, 2, 2)
    for chromosome, df in ld_data.items():
        plt.plot(df['Dist_kbp'], df['Mean_r^2'], label=chromosome)
    plt.xlabel('Distance (Kbp)')
    plt.ylabel('Mean r^2')
    plt.xscale('log')
    plt.title('LD Decay (Log Scale)')
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    output_path = os.path.join(output_dir, 'LD_decay_plots.png')
    plt.savefig(output_path)
    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot LD decay for multiple chromosomes.')
    parser.add_argument('--bin', type=int, required=True, help='Size of the distance bin in base pairs.')
    parser.add_argument('--data_dir', type=str, required=True, help='Directory containing the LD.stat.gz files.')
    parser.add_argument('--output_dir', type=str, required=True, help='Directory to save the output plot.')
    args = parser.parse_args()

    plot_ld_decay(args.data_dir, args.bin, args.output_dir)

