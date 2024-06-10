import os
import pickle
import pandas as pd
import numpy as np
from lightgbm import LGBMClassifier

# Get the directory path of the current script
script_dir = '/home/zmk214/zmk214/zmk214workingspace/genomewide/'
os.chdir(script_dir)

threshold = 0.5
lib_name = "poolSub_tcs_K562_N"
input_profiles_subtnorm_path = "data/results/poolSub_k562_tc_profiles/profiles_subtnorm/tcs/profiles_subtnorm_tcs_K562_N.csv"
input_profiles_count_path = "data/results/poolSub_k562_tc_profiles/profiles/tcs/profiles_count_tcs_K562_N.csv"
input_model_path = 'data/resources/model_Meena_v3_gm12878_fullModel.sav'
model = pickle.load(open(input_model_path, 'rb'))

# Specify the directory path
output_directory = 'data/results/poolSub_k562_tc_profiles/prediction/'
# Check if the directory exists
if not os.path.exists(output_directory):
    # If it doesn't exist, create the directory
    os.makedirs(output_directory)
    print(f"Directory '{output_directory}' created successfully.")
else:
    print(f"Directory '{output_directory}' already exists.")

output_all_results = output_directory + 'pred_all_' + lib_name + '.bed'
output_slt_results = output_directory + 'pred_slt_' + lib_name + '.bed'

# Read CSV file
df_subtnorm = pd.read_csv(input_profiles_subtnorm_path, header=0, index_col=0)
df_count = pd.read_csv(input_profiles_count_path, header=0, index_col=0)

# Functions
def extract_ranges(df):

    """
    Extracts chromosome, start, end, and strand from DataFrame index.

    Args:
    df (pandas.DataFrame): DataFrame with index containing strings in the format "chrom:start-end;strand".

    Returns:
    pandas.DataFrame: DataFrame with 'chrom', 'chromStart', 'chromEnd', and 'strand' columns.
    """
    # Split index and expand into separate columns
    split_index = df.index.str.split(':|-|;', expand=True)

    # Untuple the split_index
    chrom, chromStart, chromEnd, strand = zip(*split_index.values)

    # Create a DataFrame with column names
    result_df = pd.DataFrame({
        'chrom': chrom,
        'chromStart': chromStart,
        'chromEnd': chromEnd,
        'strand': strand
    })

    return result_df

def extract_profiles(df):
    """
    Extract profile vectors from a DataFrame.

    Parameters:
    df (pandas.DataFrame): DataFrame containing profile data.

    Returns:
    numpy.ndarray: NumPy array representing profile vectors.
    """
    # Convert DataFrame to NumPy array
    numpy_array = df.values
    
    # Convert NaN values to 0
    numpy_array[np.isnan(numpy_array)] = 0

    return numpy_array


# Extract ranges and profiles
ranges = extract_ranges(df_subtnorm)

profiles_subtnorm = extract_profiles(df_subtnorm)
y_proba = model.predict_proba(profiles_subtnorm)
ranges['score'] = y_proba[:, 1]

profiles_count = extract_profiles(df_count)
ranges['countsum'] = profiles_count.sum(axis=1)

ranges['name'] = lib_name

ranges = ranges.reindex(columns=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'countsum'])

selected_ranges = ranges[ranges['score'] >= threshold]



# Save results to .bed files
ranges.to_csv(output_all_results, sep='\t', header=True, index=False)
selected_ranges.to_csv(output_slt_results, sep='\t', header=True, index=False)
