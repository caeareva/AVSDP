
################################################################################
#
#   Author: Carlos Arevalo
#   Email: carevalo0170@gmail.com
#
#   Program description:
#   Program takes an input text file in which the first column are genes, second
#   column are the LFC values, and third column, p-values.
#   The program was written to be executed in the command line, but could be 
#   slightly modified to Jupyter Notebook. 
#
#   Program execution:
#   python3 /Users/carlosarevalo/Desktop/volcano_program.py \
#       -i /Users/carlosarevalo/Downloads/test_input_data_2.txt \
#       -o /Users/carlosarevalo/Downloads/volcano_test.png \
#       -s /Users/carlosarevalo/Downloads/stylesheet.mplstyle 
#
################################################################################

import matplotlib.pyplot as plt 
import matplotlib.patches as mplpatches
import numpy as np 
import math
import argparse

# Make instances of an argument parser and define arguments
parser=argparse.ArgumentParser()
parser.add_argument("-i", "--input_file")
parser.add_argument("-o", "--output_file")
parser.add_argument("-s", "--style_sheet") 

# Read input file and style sheet
args=parser.parse_args()
style_sheet = args.style_sheet
plt.style.use(args.style_sheet)
inFile = args.input_file
outFile = args.output_file

# Make lists of input values:
x_values = []
y_values = []
r_xvalues = []
r_yvalues = []
x_labels = []
y_labels = []
labels = []

# Read data and except NA values: 
with open(inFile, "r") as inFile:
    for line in inFile.readlines():
        try:
            splitline = line.strip().split("\t")
            genes = splitline[0]
            x = float(splitline[1])
            y = -1*np.log10(float(splitline[2]))
            # Take absolute value of LFC
            fold_change = 2**abs(x)
            # Label all metabolites with p-value < 0.01
            if x < 6 and fold_change > 0 and y > 2:
                x_labels.append(x)
                y_labels.append(y)
                labels.append(genes)
            elif fold_change > 5 and y > 2:
                r_xvalues.append(x)
                r_yvalues.append(y)
            else:
                x_values.append(x)
                y_values.append(y)
        except ValueError:
            continue

# add +0.3 and -0.3 to values in same y-value to avoid overlapping
#noise_up = 0.3
#noise_down = -0.3
#y_label = []
#for value in y_labels:
#    value = value + noise_up
#    y_label.append(value)
#    if value in y_label:
#        value = value+0.3
#        y_label.append(value)
        #print(value)

# Define figure dimensions:
figureHeight=4
figureWidth=4
plt.figure(figsize=(figureWidth, figureHeight)) 

# Panels size parameters
panelWidth=2.0
panelHeight=2.0
# Normalize axis units:
relativePanelWidth=panelWidth/figureWidth
relativePanelHeight=panelHeight/figureHeight
panel1=plt.axes([0.2, 0.2, relativePanelWidth, relativePanelHeight])

# Axis labels:
panel1.tick_params(bottom=True, labelbottom=True,
                   left=True, labelleft=True,
                   right=False, labelright=False,
                   top=False, labeltop=False)

# Panel axis and labels:
panel1.set_xlim(-8, 8)
panel1.set_ylim(0, 8) #max(y_values)+0.5)
panel1.set_xticks(np.arange(-8, 8.5, 2))
panel1.set_yticks(np.arange(0, 8.5, 1)) #max(y_values)+0.5, 1))
panel1.set_xlabel("$\mathregular{log_{2}}$(fold change)")
panel1.set_ylabel("-$\mathregular{log_{10}}$(p-value)")
panel1.set_title("EA (4:45 hrs)", fontsize=8)

# Add scatters to panel:
panel1.plot(x_values, y_values,
            marker='o',
            markerfacecolor="#CCCCCC", #(0,0,0),
            markeredgecolor="#CCCCCC", #(0,0,0),
            markersize=2,
            markeredgewidth=0,
            linewidth=0)

panel1.plot(r_xvalues, r_yvalues,
            marker='o',
            markerfacecolor="red",
            markeredgecolor="red",
            markersize=2, 
            markeredgewidth=0,
            linewidth=0)  

panel1.plot(x_labels, y_labels,
            marker='o',
            markerfacecolor="red",
            markeredgecolor="red",
            markersize=2,
            markeredgewidth=0,
            linewidth=0) 

# Add labels to scatter:
for i in range(len(labels)):
    x=x_labels[i]
    y=y_label[i]
    panel1.annotate(labels[i],(x, y),
                    xytext=(x-0.25, y),
                    fontsize=4,
                    verticalalignment="bottom",
                    horizontalalignment="center") 
# Save figure
plt.savefig(outFile, dpi=600)
