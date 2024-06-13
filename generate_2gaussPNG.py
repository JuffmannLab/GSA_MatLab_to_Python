"""
Martino Zanetti

This code generates a double square Npx*Npx bmp image featuring a double gaussian with inverted intensity.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

#Npx=int(input("Write the number of pixels defining the size of your (square) image:"))
Npx=100

# Define the image size
image_size = (Npx, Npx)

# Create a coordinate grid
x = np.linspace(-1, 1, image_size[1])
y = np.linspace(-1, 1, image_size[0])
x, y = np.meshgrid(x, y)

# Define Gaussian parameters
#fwhm 4.3
#dist 6
zoomfactor = 60
fwhm = 4.3
ratio_fwhm_sigma=2.355
sigma = fwhm/ratio_fwhm_sigma/zoomfactor
mean1 = (-3/zoomfactor, -0.1)
mean2 = (3/zoomfactor, -0.1)

# Define the two Gaussian functions
gaussian1 = np.exp(-((x - mean1[0])**2 + (y - mean1[1])**2) / (2 * sigma**2))
gaussian2 = np.exp(-((x - mean2[0])**2 + (y - mean2[1])**2) / (2 * sigma**2))

# Combine the Gaussians
combined_gaussians = gaussian1 + gaussian2

# Normalize the combined gaussians to the range [0, 1]
norm = Normalize(vmin=0, vmax=combined_gaussians.max())
normalized_image = norm(-combined_gaussians+1)

# Create a figure and axis without borders and padding
figure_size = 1
mydpi = Npx/figure_size
fig, ax = plt.subplots(figsize=(figure_size, figure_size), dpi=mydpi)  # figure_size inches * mydpi DPI = Npx pixels
ax.imshow(normalized_image, cmap='hot', origin='lower')

# Remove axes and extra space
ax.axis('off')
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

# Save the image without extra padding
plt.savefig('2gauss.png', dpi=mydpi, bbox_inches='tight', pad_inches=0, transparent=True)

plt.close(fig)  # Close the figure to free memory
