import pandas as pd
import matplotlib.pyplot as plt
import os
import logging

def main():

    # Main function to generate a depth coverage histogram plot.
    # Reads depth data from a file, processes it, and saves the plot as an SVG file.

    # Input and output paths from Snakemake
    depth_file = snakemake.input.depth
    output_plot = snakemake.output.plot
    
    FIXED_DPI = 350  # Fixed DPI for the output plot

    # Create the output directory if it does not exist
    output_dir = os.path.dirname(output_plot)
    
    print(f"Output path: {output_plot}")  # Debugging: Check if the output path is correct
    print(f"Output directory: {output_dir}")  # Debugging: Check if the directory path is correct

    os.makedirs(output_dir, exist_ok=True)

    # Load depth data
    df = pd.read_csv(depth_file, sep='\t', 
                     names=['Chromosome', 'Position', 'Depth'],  # Column names for the depth file
                     dtype={'Depth': 'uint16'})  # Optimize memory usage for depth values
    
    # Plot configuration
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Main histogram configuration
    ax.hist(df['Depth'], bins=50, alpha=0.75, 
            color='#2c7bb6', edgecolor='black',  # Set histogram color and edge color
            density=True)  # Normalize the histogram to show frequency
    
    # Annotate statistical metrics
    median_depth = df['Depth'].median()  # Calculate the median depth
    ax.axvline(median_depth, color='#d7191c', 
               linestyle='--', label=f'Median: {median_depth}X')  # Add a vertical line for the median
    
    # Set axis labels
    ax.set_xlabel('Sequencing Depth (X)', fontweight='bold')  # X-axis label
    ax.set_ylabel('Normalized Frequency', fontweight='bold')  # Y-axis label
    plt.legend()  # Add a legend to the plot

    # Save the plot to the specified output file
    try:
        plt.savefig(output_plot, dpi=FIXED_DPI, 
                    bbox_inches='tight', format='svg')  # Save as an SVG file
        plt.close()  # Close the plot to free memory
    except Exception as e:
        logging.error(f"Failed to save plot: {str(e)}")  # Log any errors during saving
        raise

if __name__ == "__main__":
    main()