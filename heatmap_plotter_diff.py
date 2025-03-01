import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the data from csv
file_path = 'C:/Users/Usuario/Desktop/Incidence-Position-On-Phoswich-Detector/diff_values.csv' 


df = pd.read_csv(file_path, usecols=['X_exp', 'Y_exp', 'X_diff', 'Y_diff'])
df = df[df['X_exp'] == 0]

# Plotting the heatmap using plt.hist2d()
plt.figure(figsize=(8, 6))
plt.title('Position deviation in the Y axis with X=0')
plt.xlabel('Y expected (mm)')
plt.ylabel('Y expected - Y calculated (mm)')

plt.hist2d(df['Y_exp'], df['Y_diff'], bins=20, cmap='YlOrRd')
plt.colorbar(label='Frequency')
plt.grid(True)
plt.show()
