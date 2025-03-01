import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the data from Excel file
#file_path = 'C:/Users/Usuario/Desktop/Incidence-Position-On-Phoswich-Detector/output_pha_6_sample_reverse_c.xlsx' 
#file_path = 'C:/Users/Usuario/Desktop/Incidence-Position-On-Phoswich-Detector/row_line_inversion_pha_6/outfiles_0807/output_pha_6_C4_L4.xlsx' 
file_path = 'C:/Users/Usuario/Desktop/Incidence-Position-On-Phoswich-Detector/2_2_20_07_newdata.xlsx' 

#file_path = 'C:/Users/Usuario/Desktop/Incidence-Position-On-Phoswich-Detector/outfiles_pos_calc/-2_2_pos.xlsx' 


df = pd.read_excel(file_path, usecols=['X', 'Y'])

# Selecting only the range of -10 to 10 for both X and Y axes
#filtered_data = df[(df['X'].between(-8, 8)) & (df['Y'].between(-8, 8))]
filtered_data=df
#filtered_data['X'] = filtered_data['X'].multiply(-1)
#filtered_data['Y'] = filtered_data['Y'].multiply(-1)
# Plotting the heatmap using plt.hist2d()
plt.figure(figsize=(8, 6))
plt.title('(-2,2) arrange')
plt.xlabel('X (mm)')
plt.ylabel('Y (mm)')

plt.hist2d(filtered_data['X'], filtered_data['Y'], bins=20, cmap='YlOrRd')
plt.colorbar(label='Frequency')
plt.grid(True)
plt.show()
