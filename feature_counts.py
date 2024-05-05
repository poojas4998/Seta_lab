import pandas as pd

def filter_feature_counts(file_path):
    # Read the GTF file
    df = pd.read_csv(file_path, sep='\t', comment='#')
    
    # Remove the suffix from column names
    df.columns = df.columns.str.replace('.Aligned.sortedByCoord.out.bam', '')

    # Selecting the first column and columns starting from the 7th position onwards
    selected_columns = df.columns[[0] + list(range(6, len(df.columns)))]

    # Filtering the DataFrame to keep only the selected columns
    filtered_df = df[selected_columns]

    # Printing the filtered DataFrame
    print(filtered_df)
    
    # Saving the filtered DataFrame as a text file
    filtered_df.to_csv('featurecounts_matrix.txt', sep='\t', index=False)

# Provide the file path to the function
file_path = '/projectnb/setagrp/pooja/data/featureCounts.txt'
filter_feature_counts(file_path)

