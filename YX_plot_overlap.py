import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches


# Read the data from Excel file
file_path = 'C:/Users/Usuario/Desktop/Incidence-Position-On-Phoswich-Detector/output_pha_6_sample_reverse_c.xlsx' 
#file_path = 'C:/Users/Usuario/Desktop/Incidence-Position-On-Phoswich-Detector/row_line_inversion_pha_6/outfiles_0807/output_pha_6_C4_L4.xlsx' 



df = pd.read_excel(file_path, usecols=['X', 'Y'])

# Selecting only the range of -10 to 10 for both X and Y axes
filtered_data=df

fig, ax = plt.subplots()
circle = patches.Circle((0, -0.8), 2.2/2, linewidth=1, edgecolor='green', facecolor='none')
circle1 = patches.Circle((5.4, 5.1), 2.2/2, linewidth=1, edgecolor='green', facecolor='none')

# Plotting the heatmap using plt.hist2d()
ax.set_xlabel('X (mm)')
ax.set_ylabel('Y (mm)')

plt.hist2d(filtered_data['X'], filtered_data['Y'], bins=20, cmap='YlOrRd')
plt.colorbar(label='Frequency')
ax.grid(True)
ax.add_patch(circle)
ax.add_patch(circle1)
plt.show()
