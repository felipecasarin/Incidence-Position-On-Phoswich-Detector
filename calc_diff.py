
import pandas as pd
import numpy as np
import os

x_span = np.array([-6, -4, -2, 0, 2, 4, 6])
y_span = np.array([-6, -4, -2, 0, 2, 4])

# Function to calculate the differences and add columns
def calc_diff(xlsx_file_path, x_value, y_value):
    df = pd.read_excel(xlsx_file_path)
    
    # Ensure 'X' and 'Y' columns exist in the DataFrame
    if 'X' not in df.columns or 'Y' not in df.columns:
        raise KeyError("'X' or 'Y' columns not found in the Excel file.")

    # Calculate differences and add new columns
    df['X_exp'] = x_value
    df['Y_exp'] = y_value
    df['X_diff'] = x_value - df['X']
    df['Y_diff'] = y_value - df['Y']
    
    return df[['X', 'Y', 'X_exp', 'Y_exp', 'X_diff', 'Y_diff']]  # Return only desired columns

# List to store DataFrame results
results = []

# Process each generated xlsx file
for x_value in x_span:
    for y_value in y_span:
        xlsx_file_path = f'outfiles_pos_calc/{x_value}_{y_value}_pos.xlsx'
        
        # Check if the file exists
        if os.path.exists(xlsx_file_path):
            diffs = calc_diff(xlsx_file_path, x_value, y_value)
            results.append(diffs)
        else:
            print(f"File {xlsx_file_path} does not exist.")

# Concatenate all DataFrames into a single DataFrame
results_df = pd.concat(results, ignore_index=True)

# Save the results to a new CSV file
results_df.to_csv('diff_values.csv', index=False)
print("Values differences have been calculated and saved to 'diff_values.csv'")