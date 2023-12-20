import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the data from Excel file
file_path = 'C:/Users/Usuario/Desktop/Incidence-Position-On-Phoswich-Detector/output_gp23.xlsx' 
df = pd.read_excel(file_path, usecols=['X', 'Y'])

# Selecting only the range of -10 to 10 for both X and Y axes
filtered_data = df[(df['X'].between(-7.5, 7.5)) & (df['Y'].between(-7.5, 7.5))]

# Plotting the heatmap using plt.hist2d()
plt.figure(figsize=(8, 6))
plt.title('Heatmap of X and Y within the range of -10 to 10')
plt.xlabel('X')
plt.ylabel('Y')

plt.hist2d(filtered_data['X'], filtered_data['Y'], bins=40, cmap='YlOrRd')
plt.colorbar(label='Frequency')
plt.grid(True)
plt.show()
