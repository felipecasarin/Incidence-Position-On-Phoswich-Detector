import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the data from Excel file
file_path = 'C:/Users/Usuario/Desktop/Incidence-Position-On-Phoswich-Detector/diff_values.csv' 

df = pd.read_csv(file_path, usecols=['X', 'Y', 'X_exp', 'Y_exp'])

# Selecting only the range of -10 to 10 for both X and Y axes
#filtered_data = df[(df['X'].between(-8, 8)) & (df['Y'].between(-8, 8))]
filtered_data=df[df['X_exp']==0]
filtered_data=filtered_data[filtered_data['Y_exp']==-6]

# Plotting the heatmap using plt.hist2d()
plt.figure(figsize=(8, 6))
plt.title('Distribution of estimated positions in the Y-axis for the data point (0,-6)')
plt.xlabel('Y (mm)')
plt.ylabel('Number of events')

plt.hist(filtered_data['Y'], bins=30)
plt.grid(True)
plt.show()
