
# Author: Carlos Arevalo(caeareva)
# Email: carevalo0170@gmail.com

# Required modules
import matplotlib.pyplot as plt 
import matplotlib.patches as mplpatches
import numpy as np 
import argparse

# Make instances of an argument parser and define arguments
parser = argparse.ArgumentParser()
parser.add_argument("-s", "--style_sheet") 
parser.add_argument("-o", "--output_file")

# Read input style sheet
args = parser.parse_args()
style_sheet = args.style_sheet
plt.style.use(args.style_sheet)

# Define figure sizes
figureHeight = 2.0
figureWidth = 3.42

# Generate a wide canvas of arbitrary size
plt.figure(figsize = (figureWidth, figureHeight)) 

# Panel dimensions
panelWidth = 1.0
panelHeight = 1.0

# Normalize axis units:
relativePanelWidth = panelWidth/figureWidth
relaticePanelHeight = panelHeight/figureHeight

# Make sure panels do not overlap: Left, bottom, width, height
panel1 = plt.axes([0.1, 0.2, relativePanelWidth, relaticePanelHeight]) 
panel2 = plt.axes([0.55, 0.2, relativePanelWidth, relaticePanelHeight]) 

# Turn axis labels off:
panel1.tick_params(# axis ="both", which = "both",
                   bottom = False, labelbottom = False,
                   left = False, labelleft = False,
                   right = False, labelright = False,
                   top = False, labeltop = False)

panel2.tick_params(bottom = False, labelbottom = False,
                    left = False, labelleft = False,
                    right = False, labelright = False,
                    top = False, labeltop = False)

# Add points to panel 1:
points = np.arange(0, np.pi/2, 0.065)
for point in points:
    sin = np.sin(point)
    cos = np.cos(point)
    panel1.plot(cos, sin, 
           marker = "o", 
           markerfacecolor = (cos, cos, cos),  
           markersize = 2, 
           markeredgewidth = 0, 
           linewidth = 0)
 
# Add rectangles to panel 2: 10 by 10
rect_coords = np.linspace(0, 0.9, 10, endpoint=True) # np.arange(0, 0.9*np.pi, 10)
for zeta in rect_coords:
    for tau in rect_coords:
        rectangle = mplpatches.Rectangle((zeta, tau), 0.1, 0.1, 
                                        facecolor=(zeta, tau, 1),
                                        edgecolor = "black", 
                                        linewidth = 1)

        panel2.add_patch(rectangle)

# Save figure
plt.savefig("asg1_figure.png", dpi=600)
