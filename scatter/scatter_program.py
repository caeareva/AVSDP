
################################################################################
#
#   Author: Carlos Arevalo
#   Email: carevalo0170@gmail.com
#
#   Program description:
#   Program takes an input text file in which the first column are genes, second
#   column are the replicate 1 counts, and third column, replicate 2 counts.
#   The program was written to be run in the command line, but could be slightly 
#   modify to be run in Jupyter notebook. 
#
#   Program execution:
#   python3 /Users/carlosarevalo/Desktop/scatter_program.py \
#       -i /Users/carlosarevalo/Downloads/test_input_data_1.txt \
#       -o /Users/carlosarevalo/Downloads/scatter_test.png \
#       -s /Users/carlosarevalo/Downloads/stylesheet.mplstyle 
#
################################################################################

# Use args to input data
import matplotlib.pyplot as plt 
import matplotlib.patches as mplpatches
import numpy as np 
import argparse

# Make instances of an argument parser and define arguments
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file")
parser.add_argument("-o", "--output_file")
parser.add_argument("-s", "--style_sheet") 

# Read input file and style sheet
args = parser.parse_args()
style_sheet = args.style_sheet
plt.style.use(args.style_sheet)
inFile = args.input_file
# inFile=open("/Users/carlosarevalo/Downloads/test_input_data_1.txt", "r")

# Make lists of input values:
x_values=[]
y_values=[]

# Read input data
with open(inFile, "r") as inFile:
    for line in inFile.readlines():
        splitline=line.strip().split("\t")
        name=splitline[0]
        x_values.append(int(splitline[1]))
        y_values.append(int(splitline[2]))
    pass

x_vals=np.log2([x+1 for x in np.array(x_values)])
y_vals=np.log2([y+1 for y in np.array(y_values)])

# Define figure dimensions:
figureHeight=2
figureWidth=5

# Generate a wide canvas of arbitrary size:
plt.figure(figsize=(figureWidth, figureHeight)) 

# Normalized panel size parameters:
panelWidth=1/figureWidth
panelHeight=1/figureHeight
x_hist_width=panelWidth
x_hist_height=(1/4)/figureHeight
y_hist_width=(1/4)/figureWidth
y_hist_height=panelHeight

# Make sure panels do not overlap: Left, bottom, width, height
panel1 = plt.axes([0.14, 0.15, panelWidth, panelHeight]) 
panel2 = plt.axes([0.54, 0.15, panelWidth, panelHeight])

# Top panels:
x_panel1 = plt.axes([0.14, 0.685, x_hist_width, x_hist_height]) 
x_panel2 = plt.axes([0.54, 0.685, x_hist_width, x_hist_height])

# Side panels:
y_panel1 = plt.axes([0.076, 0.15, y_hist_width, y_hist_height])  
y_panel2 = plt.axes([0.476, 0.15, y_hist_width, y_hist_height]) 

# Panel 1 labels:
panel1.tick_params(bottom=True, labelbottom=True,
                   left=False, labelleft=False,
                   right=False, labelright=False,
                   top=False, labeltop=False)

x_panel1.tick_params(bottom=False, labelbottom=False,
                    left=True, labelleft=True,
                    right=False, labelright=False,
                    top=False, labeltop=False)

y_panel1.tick_params(bottom=True, labelbottom=True,
                    left=True, labelleft=True,
                    right=False, labelright=False,
                    top=False, labeltop=False)
# Panel 2 labels:
panel2.tick_params(bottom=True, labelbottom=True,
                   left=False, labelleft=False,
                   right=False, labelright=False,
                   top=False, labeltop=False)

x_panel2.tick_params(bottom = False, labelbottom = False,
                    left = True, labelleft = True,
                    right = False, labelright = False,
                    top = False, labeltop = False)

y_panel2.tick_params(bottom=True, labelbottom=True,
                    left=True, labelleft=True,
                    right=False, labelright=False,
                    top=False, labeltop=False)

# Panel 1 scatter plot:
panel1.plot(x_vals,y_vals,
    marker='o',
    markerfacecolor=(0,0,0),
    markeredgecolor=(0,0,0),
    markersize=1.5,
    markeredgewidth=0,
    linewidth=0,
    alpha=0.1)    
# panel1.text(0,0,'r='+str(round(stats.spearmanr(xList,yList)[0],2)))
# print(round(stats.spearmanr(xList,yList)[0],2))
panel1.set_xlim(0,max(x_vals))
panel1.set_ylim(0,max(y_vals))


# Add histograms to panel group 1:
x_hist_vals, x_hist_bins = np.histogram(x_vals, np.arange(0,15,0.5))
y_hist_vals, y_hist_bins = np.histogram(y_vals, np.arange(0,15,0.5))

# Plot histograms
for i in np.arange(0,len(x_hist_vals),1):
    x_left=i/2 # len(x_hist_bins)
    x_bottom=0
    x_width=1/2 # len(x_hist_bins)
    x_height=np.log2(x_hist_vals[i]+1)
    y_left=0
    y_bottom=i/2
    y_width=np.log2(y_hist_vals[i]+1)
    y_height=0.5
    rect1 = mplpatches.Rectangle((x_left,x_bottom),x_width,x_height,
                                  facecolor="grey", 
                                  edgecolor=(0,0,0),
                                  linewidth=0.1)
    rect2 = mplpatches.Rectangle((y_left,y_bottom),y_width,y_height,
                                  facecolor="grey", 
                                  edgecolor=(0,0,0),
                                  linewidth=0.1)
    x_panel1.add_patch(rect1)
    y_panel1.add_patch(rect2)

# Panel 1 group labels:
panel1.set_xlim(0,15)
panel1.set_ylim(0,15)
panel1.set_xticks(np.arange(0,16,5))
x_panel1.set_ylim(0,20) 
x_panel1.set_xlim(0,15)
# x_panel1.set_yticks(np.arange(20,-1,20))
y_panel1.set_ylim(0,15)
y_panel1.set_xlim(20,0)
y_panel1.set_yticks(np.arange(0,16,5))

# Panel 2 scatter-heatmap plot:
#for zeta in rect_coords:
#    for tau in rect_coords:
#        rectangle = mplpatches.Rectangle((zeta, tau), 0.1, 0.1, 
#                                        facecolor=(zeta, tau, 1),
#                                        edgecolor = "black", 
#                                        linewidth = 1)
#        panel2.add_patch(rectangle)

panel2.plot(x_vals,y_vals,
    marker='o',
    markerfacecolor=(0,0,0),
    markeredgecolor=(0,0,0),
    markersize=1.5,
    markeredgewidth=0,
    linewidth=0,
    alpha=0.1)    
panel2.set_xlim(0,max(x_vals))
panel2.set_ylim(0,max(y_vals))

for i in np.arange(0,len(x_hist_vals),1):
    x_left=i/2 # len(x_hist_bins)
    x_bottom=0
    x_width=1/2 # len(x_hist_bins)
    x_height=np.log2(x_hist_vals[i]+1)
    y_left=0
    y_bottom=i/2
    y_width=np.log2(y_hist_vals[i]+1)
    y_height=0.5
    rect1 = mplpatches.Rectangle((x_left,x_bottom),x_width,x_height,
                                  facecolor="grey", 
                                  edgecolor=(0,0,0),
                                  linewidth=0.1)
    rect2 = mplpatches.Rectangle((y_left,y_bottom),y_width,y_height,
                                  facecolor="grey", 
                                  edgecolor=(0,0,0),
                                  linewidth=0.1)
    x_panel2.add_patch(rect1)
    y_panel2.add_patch(rect2)

# Panel 2 group labels:
panel2.set_xlim(0,15)
panel2.set_ylim(0,15)
x_panel2.set_ylim(0,20)
x_panel2.set_xlim(0,15)
y_panel2.set_ylim(0,15)
y_panel2.set_xlim(20,0)
panel2.set_xticks(np.arange(0,16,5))

# Save figure
plt.savefig('scatter_test.png', dpi=600)

