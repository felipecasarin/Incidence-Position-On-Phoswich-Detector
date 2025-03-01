import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches

# Initialize xtab and ytab
xtab = np.array([-5.3, -3.3, -0.3, 1.7, 4.7])
ytab = np.array([-5.3, -3.3, -0.3, 1.7, 4.7, 6.7])

# Initialize all_points as an empty list to collect points
all_points = []

# Collect points from xtab and ytab
for x in xtab:
    for y in ytab:
        point = np.array([x, y])
        all_points.append(point)

# Convert all_points to a numpy array
all_points = np.array(all_points)

# Create a figure and an axes
fig, ax = plt.subplots()

# Plot the points using matplotlib and set the color to red
ax.scatter(all_points[:, 0], all_points[:, 1], color='red')

# Customize the x and y intervals
x_ticks = np.arange(-10, 10, 2)  # Set x ticks from -10 to 10 with a step of 2
y_ticks = np.arange(-10, 10, 2)  # Set y ticks from -10 to 10 with a step of 2
ax.set_xticks(x_ticks)
ax.set_yticks(y_ticks)

# Set the limits of the axes
ax.set_xlim(-10, 10)
ax.set_ylim(-10, 10)

# Ensure the aspect ratio is equal so the plot is a square
ax.set_aspect('equal')

# Create a square centered at the origin with size 2
square = patches.Rectangle((-8.4, -8.4), 16.8, 16.8, edgecolor='blue', facecolor='none')

# Add the square to the plot
ax.add_patch(square)

ax.set_xlabel('X (mm)')  # Label for x-axis
ax.set_ylabel('Y (mm)')  # Label for y-axis
ax.grid(True)  # Display grid
plt.show()  # Show the plot
