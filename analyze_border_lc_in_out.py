import pandas as pd
import matplotlib.pyplot as plt

# Read the data from Excel file
file_path = 'C:/Users/Usuario/Desktop/Incidence-Position-On-Phoswich-Detector/output_pha_6_sample_softl1.xlsx' 
df = pd.read_excel(file_path, usecols=['X', 'Y', 'cL1','cL2','cL3','cL4','cC1','cC2','cC3','cC4'])

# Selecting outer-edge and inner points
filtered_data_out = df[(df['X'].between(7.5, 15)) | (df['Y'].between(-7.5, -15))]
filtered_data_in = df[(df['X'].between(-15, 7.5)) | (df['Y'].between(-7.5, 10))]

# Define the columns you want to plot
columns_to_plot = ['cC1', 'cC2', 'cC3', 'cC4']
lines_to_plot = ['cL1', 'cL2', 'cL3', 'cL4']

# Create subplots with 1 row and 2 columns
fig, axes = plt.subplots(2, 2, figsize=(15, 6))

# Plot for filtered_data_out
for column in columns_to_plot:
    bins = pd.cut(filtered_data_out[column], bins=50)
    # Grouping by the bins and calculating the number of entries in the bin
    grouped = filtered_data_out.groupby(bins)[column].sum()
    bin_midpoints = [(interval.left + interval.right) / 2 for interval in grouped.index]
    axes[0,0].plot(bin_midpoints, grouped, marker='o', label=column)

# Plot for filtered_data_in
for column in columns_to_plot:
    bins = pd.cut(filtered_data_in[column], bins=50)
    # Grouping by the bins and calculating the number of entries in the bin
    grouped = filtered_data_in.groupby(bins)[column].sum()
    bin_midpoints = [(interval.left + interval.right) / 2 for interval in grouped.index]
    axes[0,1].plot(bin_midpoints, grouped, marker='o', label=column)


# Plot for filtered_data_out
for line in lines_to_plot:
    bins = pd.cut(filtered_data_out[line], bins=50)
    # Grouping by the bins and calculating the number of entries in the bin
    grouped = filtered_data_out.groupby(bins)[line].sum()
    bin_midpoints = [(interval.left + interval.right) / 2 for interval in grouped.index]
    axes[1,0].plot(bin_midpoints, grouped, marker='o', label=line)

# Plot for filtered_data_in
for line in lines_to_plot:
    bins = pd.cut(filtered_data_in[line], bins=50)
    # Grouping by the bins and calculating the number of entries in the bin
    grouped = filtered_data_in.groupby(bins)[line].sum()
    bin_midpoints = [(interval.left + interval.right) / 2 for interval in grouped.index]
    axes[1,1].plot(bin_midpoints, grouped, marker='o', label=line)

# Add labels and legend for the first subplot
axes[0,0].set_xlabel('Normalized Light Yield')
axes[0,0].set_ylabel('Number of entries')
axes[0,0].set_title('Column: Outer-Edge Points')
axes[0,0].legend()
axes[0,0].grid(True)

# Add labels and legend for the second subplot
axes[0,1].set_xlabel('Normalized Light Yield')
axes[0,1].set_ylabel('Number of entries')
axes[0,1].set_title('Column: Inner Points')
axes[0,1].legend()
axes[0,1].grid(True)

# Add labels and legend for the third subplot
axes[1,0].set_xlabel('Normalized Light Yield')
axes[1,0].set_ylabel('Number of entries')
axes[1,0].set_title('Row: Outer-Edge Points')
axes[1,0].legend()
axes[1,0].grid(True)

# Add labels and legend for the fouth subplot
axes[1,1].set_xlabel('Normalized Light Yield')
axes[1,1].set_ylabel('Number of entries')
axes[1,1].set_title('Row: Inner Points')
axes[1,1].legend()
axes[1,1].grid(True)

# Adjust layout to prevent overlap
plt.tight_layout()

# Show the plot
plt.show()
