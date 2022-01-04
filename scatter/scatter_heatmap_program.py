
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
#   Instead of plotting the scatter points, program plots a heatmap.
# 
#   Program execution:
#   python3 /Users/carlosarevalo/Desktop/scatter_heatmap_program.py \
#       -i /Users/carlosarevalo/Downloads/test_input_data_1.txt \
#       -o /Users/carlosarevalo/Downloads/scatter_heatmap_test.png \
#       -s /Users/carlosarevalo/Downloads/stylesheet.mplstyle 
#
################################################################################

# Use args to input data
import matplotlib.pyplot as plt 
import matplotlib.patches as mplpatches
import numpy as np 
import warnings
import argparse

def general_flooring(number, floor):
  return int(number/floor)*floor

# heatmap function
def heatmap(xList,yList,panel,xmin,xmax,ymin,ymax,binsize):
  white=(1,1,1)
  black=(0,0,0)
  R=np.linspace(white[0],black[0],21)
  B=np.linspace(white[1],black[1],21)
  G=np.linspace(white[2],black[2],21)
  dataDict = {}
  for xBin in np.arange(xmin,xmax,binsize):
    dataDict[xBin]={}
    for yBin in np.arange(ymin,ymax,binsize):
      dataDict[xBin][yBin]=0

  for i in np.arange(0,len(xList),1):
    x=xList[i]
    y=yList[i]
    xFloored=general_flooring(x,binsize) # 7.9 becomes to 7
    yFloored=general_flooring(y,binsize)
    if xFloored in dataDict:
      if yFloored in dataDict[xFloored]:
        dataDict[xFloored][yFloored]+=1 # iterates up by 1

  for xBin in np.arange(xmin,xmax,binsize):
    for yBin in np.arange(ymin,ymax,binsize):
      value=min(dataDict[xBin][yBin],20)
      color=(R[value],G[value],B[value])
      rect=mplpatches.Rectangle((xBin,yBin),binsize,binsize,
                                  facecolor=color, 
                                  edgecolor=(0,0,0),
                                  linewidth=0)
      panel2.add_patch(rect)


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
outFile = args.output_file

# Make lists of input values:
x_values=[]
y_values=[]

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
figureHeight=2.5 # 2
figureWidth=5.5 # 5
# Generates a wide canvas of arbitrary size:
plt.figure(figsize=(figureWidth, figureHeight)) 

# Normalized panel size parameters:
panelWidth=1/figureWidth #*1.5
panelHeight=1/figureHeight #*1.5
x_hist_width=panelWidth
x_hist_height=(1/4)/figureHeight # *1.5
y_hist_width=(1/4)/figureWidth# *1.5
y_hist_height=panelHeight

# Makes sure panels do not overlap: Left, bottom, width, height   
panel1=plt.axes([0.16, 0.2, panelWidth, panelHeight]) 
panel2=plt.axes([0.55, 0.2, panelWidth, panelHeight])
# Top panels:
x_panel1=plt.axes([0.16, 0.63, x_hist_width, x_hist_height]) 
x_panel2=plt.axes([0.55, 0.63, x_hist_width, x_hist_height])
# Side panels:
y_panel1=plt.axes([0.1, 0.2, y_hist_width, y_hist_height])  
y_panel2=plt.axes([0.487, 0.2, y_hist_width, y_hist_height]) 
# scale panel
panel=plt.axes([0.78, 0.2, y_hist_width/2, y_hist_height]) 


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
panel2.tick_params(bottom=True, labelbottom=True,
                   left=False, labelleft=False,
                   right=False, labelright=False,
                   top=False, labeltop=False)
x_panel2.tick_params(bottom=False, labelbottom=False,
                    left=True, labelleft=True,
                    right=False, labelright=False,
                    top=False, labeltop=False)
y_panel2.tick_params(bottom=True, labelbottom=True,
                    left=True, labelleft=True,
                    right=False, labelright=False,
                    top=False, labeltop=False)
panel.tick_params(bottom=False, labelbottom=True,
                  left=True, labelleft=True,
                  right=False, labelright=False,
                  top=False, labeltop=False)

# Panel 1 scatter plot:
panel1.plot(x_vals,y_vals,
    marker='o',
    markerfacecolor=(0,0,0),
    markeredgecolor=(0,0,0),
    markersize=2,
    markeredgewidth=0,
    linewidth=0,
    alpha=0.1) 

panel1.set_xlim(0,max(x_vals))
panel1.set_ylim(0,max(y_vals))

# Add histograms to panel group 1:
x_hist_vals, x_hist_bins = np.histogram(x_vals, np.arange(0,15,0.5))
y_hist_vals, y_hist_bins = np.histogram(y_vals, np.arange(0,15,0.5))

# plot histograms
for i in np.arange(0,len(x_hist_vals),1):
    x_left=i/2 
    x_bottom=0
    x_width=1/2 
    x_height=np.log2(x_hist_vals[i]+1)
    y_left=0
    y_bottom=i/2
    y_width=np.log2(y_hist_vals[i]+1)
    y_height=0.5
    rect1=mplpatches.Rectangle((x_left,x_bottom),x_width,x_height,
                                  facecolor="grey", 
                                  edgecolor=(0,0,0),
                                  linewidth=0.1)
    rect2=mplpatches.Rectangle((y_left,y_bottom),y_width,y_height,
                                  facecolor="grey", 
                                  edgecolor=(0,0,0),
                                  linewidth=0.1)
    x_panel1.add_patch(rect1)
    y_panel1.add_patch(rect2)

for i in np.arange(0,len(x_hist_vals),1):
    x_left=i/2 
    x_bottom=0
    x_width=1/2 
    x_height=np.log2(x_hist_vals[i]+1)
    y_left=0
    y_bottom=i/2
    y_width=np.log2(y_hist_vals[i]+1)
    y_height=0.5
    rect1=mplpatches.Rectangle((x_left,x_bottom),x_width,x_height,
                                  facecolor="grey", 
                                  edgecolor=(0,0,0),
                                  linewidth=0.1)
    rect2=mplpatches.Rectangle((y_left,y_bottom),y_width,y_height,
                                  facecolor="grey", 
                                  edgecolor=(0,0,0),
                                  linewidth=0.1)
    x_panel2.add_patch(rect1)
    y_panel2.add_patch(rect2)


# heatmap scale
colorRange=np.linspace(0,0.8,21)
colors=[(1-i, 1-i, i) for i in colorRange]
for color in colors:
    patch=mplpatches.Rectangle((0, colors.index(color)+0), 1, 21, linewidth=0, facecolor=color)
    panel.add_patch(patch)

# call heatmap
heatmap(x_vals,y_vals,panel2,0,15,0,15,0.33)

panel2.set_xlim(0,max(x_vals))
panel2.set_ylim(0,max(y_vals))

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
# Panel 2 group labels:
panel2.set_xlim(0,15)
panel2.set_ylim(0,15)
panel2.set_xticks(np.arange(0,16,5))
x_panel2.set_ylim(0,20) 
x_panel2.set_xlim(0,15)
y_panel2.set_ylim(0,15)
y_panel2.set_xlim(20,0)
y_panel2.set_yticks(np.arange(0,16,5))

# Panel 1 axis titles
panel1.set_xlabel("Replicate 1(RPM)")
y_panel1.set_ylabel("Replicate 2(RPM)")
y_panel1.set_xlabel("# of \ngenes")
x_panel1.set_ylabel("# of \ngenes")
# Panel 2 axis titles
panel2.set_xlabel("Replicate 1(RPM)")
y_panel2.set_ylabel("Replicate 2(RPM)")
y_panel2.set_xlabel("# of \ngenes")
x_panel2.set_ylabel("# of \ngenes")
# scale panel
warnings.filterwarnings("ignore")
ticks = [x for x in range(0, 20, 10)]
ticks.append(">20")
# print(ticks) 
panel.set_ylim(0,20.00000001,10)
panel.set_yticklabels(ticks) #np.arange(0,21,10))
panel.set_xticks([])

# Save figure
plt.savefig(outFile, dpi=600) 


