import pandas as pd
import numpy as np

x_span = np.array([-6, -4, -2, 0, 2, 4, 6])
y_span = np.array([-6, -4, -2, 0, 2, 4])

# List to store the results
results = []

# Function to calculate the average values of specific columns in a CSV file
def calculate_averages(csv_file_path, x_value, y_value):
    df = pd.read_csv(csv_file_path)
    columns = ['C1', 'C2', 'C3', 'C4', 'L1', 'L2', 'L3', 'L4']
    averages = df[columns].mean()
    averages['x_value'] = x_value
    averages['y_value'] = y_value
    return averages

# Process each generated CSV file
for x_value in x_span:
    for y_value in y_span:
        csv_file_path = 'collection_files_csv_filtered/' + str(x_value) + '_' + str(y_value) + '.csv'
        averages = calculate_averages(csv_file_path, x_value, y_value)
        results.append(averages)

# Create a DataFrame from the results
results_df = pd.DataFrame(results)

# Save the results to a new CSV file
results_df.to_csv('averages_results.csv', index=False)
print("Averages have been calculated and saved to 'averages_results.csv'")


