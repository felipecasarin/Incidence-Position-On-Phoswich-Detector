import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the data from Excel file
file_path_x = 'C:/Users/Usuario/Desktop/Incidence-Position-On-Phoswich-Detector/FWHM_X.xlsx' 
file_path_y = 'C:/Users/Usuario/Desktop/Incidence-Position-On-Phoswich-Detector/FWHM_Y.xlsx' 

dfx = pd.read_excel(file_path_x, usecols=['X', 'Y', 'FWHM'])
dfy = pd.read_excel(file_path_y, usecols=['X', 'Y', 'FWHM'])

# Plotting the curves using plt.plot()
plt.figure(figsize=(8, 6))
plt.title('FWHM along both axis at center row and column of measurment points')
plt.xlabel('FWHM (mm) in the given direction')
plt.ylabel('Position of measurement (mm) in the given direction')

plt.plot(dfy['Y'],dfy['FWHM'], marker='o', linestyle='-', color='blue', label='Y direction')
plt.plot(dfx['X'],dfx['FWHM'], marker='o', linestyle='-', color='red', label='X direction')
plt.ylim(0, 5)
plt.legend()
plt.grid(True)
plt.show()
