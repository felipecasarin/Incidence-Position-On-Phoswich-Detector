import uproot
import pandas as pd
import numpy as np
import awkward as ak

x_span=np.array([-6,-4,-2,0,2,4,6])
y_span=np.array([-6,-4,-2,0,2,4])

# Function to convert a ROOT file to a CSV file with specific columns and value range
def root_to_csv(root_file_path, tree_name, csv_file_path,x_value,y_value):
    # Open the ROOT file
    with uproot.open(root_file_path) as file:
        # Access the tree
        tree = file[tree_name]
        # Convert the tree to a pandas DataFrame
        #df = tree.arrays(library='pd')
        ak_array = tree.arrays(library='ak')
        df = ak.to_dataframe(ak_array)
        
        # Select only the desired columns
        columns = ['C1', 'C2', 'C3', 'C4', 'L1', 'L2', 'L3', 'L4']
        df = df[columns]
        df = df.apply(pd.to_numeric, errors='coerce')
        
        # Filter rows where values are between 0 and 200 for each column
        #mask = (df['C1'] > 0) & (df['C1'] < 200) & (df['C2'] > 0) & (df['C2'] < 200) & (df['C3'] > 0) & (df['C3'] < 200) & (df['C4'] > 0) & (df['C4'] < 200) & (df['L1'] > 0) & (df['L1'] < 200) & (df['L2'] > 0) & (df['L2'] < 200) & (df['L3'] > 0) & (df['L3'] < 200) & (df['L4'] > 0) & (df['L4'] < 200)
        mask = (df > 0) & (df < 200)
        df = df[mask.all(axis=1)]
        
        # Save the DataFrame to a CSV file
        df.to_csv(csv_file_path, index=False)

# Example usage
for x_value in x_span:
    for y_value in y_span:
        root_file_path = 'collection_files/' + str(x_value) +'_' + str(y_value) +'_NS.root'
        tree_name = 'tree'
        csv_file_path = 'collection_files_csv_filtered/' + str(x_value) + '_' + str(y_value) + '.csv'
        root_to_csv(root_file_path, tree_name, csv_file_path, x_value, y_value)
        print(str(x_value) + '_' + str(y_value) + '.csv is done!')