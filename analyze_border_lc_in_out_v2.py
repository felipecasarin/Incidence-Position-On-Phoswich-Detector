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
    c_bins_out = pd.cut(filtered_data_out[column], bins=50)
    c_grouped_out = filtered_data_out.groupby(c_bins_out)[column].sum()
    c_bin_midpoints_out = [(interval.left + interval.right) / 2 for interval in c_grouped_out.index]
    axes[0,0].plot(c_bin_midpoints_out, c_grouped_out, marker='o', label=column)

    c_bins_in = pd.cut(filtered_data_in[column], bins=50)
    c_grouped_in = filtered_data_in.groupby(c_bins_in)[column].sum()
    c_bin_midpoints_in = [(interval.left + interval.right) / 2 for interval in c_grouped_in.index]
    axes[1,0].plot(c_bin_midpoints_in, c_grouped_in, marker='o', label=column)    


for line in lines_to_plot:
    l_bins_out = pd.cut(filtered_data_out[line], bins=50)
    l_grouped_out = filtered_data_out.groupby(l_bins_out)[line].sum()
    l_bin_midpoints_out = [(interval.left + interval.right) / 2 for interval in l_grouped_out.index]
    axes[0,1].plot(l_bin_midpoints_out, l_grouped_out, marker='o', label=line)

    l_bins_in = pd.cut(filtered_data_in[line], bins=50)
    l_grouped_in = filtered_data_in.groupby(l_bins_in)[line].sum()
    l_bin_midpoints_in = [(interval.left + interval.right) / 2 for interval in l_grouped_in.index]
    axes[1,1].plot(l_bin_midpoints_in, l_grouped_in, marker='o', label=line)  



# Add labels and legend for the first subplot
axes[0,0].set_xlabel('Normalized Light Yield')
axes[0,0].set_ylabel('Number of entries')
axes[0,0].set_title('Column: Outer-Edge Points')
axes[0,0].legend()
axes[0,0].grid(True)
axes[0,0].locator_params(axis='x', nbins=10)


# Add labels and legend for the second subplot
axes[0,1].set_xlabel('Normalized Light Yield')
axes[0,1].set_ylabel('Number of entries')
axes[0,1].set_title('Row: Inner Points')
axes[0,1].legend()
axes[0,1].grid(True)
axes[0,1].locator_params(axis='x', nbins=10)

# Add labels and legend for the third subplot
axes[1,0].set_xlabel('Normalized Light Yield')
axes[1,0].set_ylabel('Number of entries')
axes[1,0].set_title('Column: Outer-Edge Points')
axes[1,0].legend()
axes[1,0].grid(True)
axes[1,0].locator_params(axis='x', nbins=10)

# Add labels and legend for the fouth subplot
axes[1,1].set_xlabel('Normalized Light Yield')
axes[1,1].set_ylabel('Number of entries')
axes[1,1].set_title('Row: Inner Points')
axes[1,1].legend()
axes[1,1].grid(True)
axes[1,0].locator_params(axis='x', nbins=10)

# Adjust layout to prevent overlap
plt.tight_layout()

# Show the plot
plt.show()
