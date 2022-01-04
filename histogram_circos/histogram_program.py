
################################################################################
#
#   Author: Carlos Arevalo
#   Email: carevalo0170@gmail.com
#
#   Program description:
#   Program takes an input text file. The program was written to be run in the 
#   command line, but could be slightly modified to be run in Jupyter notebook. 
#
#   Program execution:
#	python3 /Users/carevalo/Desktop/histogram_circos/histogram_program.py \
#        -i /Users/carevalo/Desktop/histogram_circos/test_input_data_4.txt \
#        -o /Users/carevalo/Desktop/histogram_circos/test_figure.png \
#        -s /Users/carevalo/Downloads/stylesheet.mplstyle 
#
################################################################################

# Required modules
import matplotlib.pyplot as plt 
import matplotlib.patches as mplpatches
import numpy as np 
import warnings
import argparse

# Argument parser and definitions
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file")
parser.add_argument("-s", "--style_sheet") 
parser.add_argument("-o", "--output_file")

# Read input file and style sheet
args = parser.parse_args()
style_sheet = args.style_sheet
plt.style.use(args.style_sheet)
inFile = args.input_file
outFile = args.output_file

# Define figure dimentions
figureWidth = 5.5
figureHeight = 3.5
plt.figure(figsize = (figureWidth, figureHeight)) 

# Define panel dimentions
panel_width = 1/figureWidth
panel_heigth = 1/figureHeight
panelHeight = 2.5
panel1_Width = 0.75
panel2_Width = 2.5
relativePanelWidth1 = panel1_Width/figureWidth
relativePanelWidth2 = panel2_Width/figureWidth
relativePanelHeight = panelHeight/figureHeight

# Make figure panels
#panel1 = plt.axes([0.5/figureWidth, 0.1, relativePanelWidth1, relativePanelHeight])
#panel2 = plt.axes([1.75/figureWidth, 0.1, relativePanelWidth2, relativePanelHeight])
panel1 = plt.axes([0.7/figureWidth, 0.1, relativePanelWidth1, relativePanelHeight])
panel2 = plt.axes([2.0/figureWidth, 0.1, relativePanelWidth2, relativePanelHeight])

data_count = []

# Read data and except NA values
with open(inFile, "r") as inFile: 
	for line in inFile.readlines():
		try:
			split_line = line.strip().split("\t")
			random_list = [int(i) for i in split_line[4:12]]
			normalized_y_data = [(j-min(random_list))/(max(random_list)-min(random_list))*100 for j in random_list]
			normalized_y_data.append(float(split_line[13]))
			data_count.append(normalized_y_data)
		except ValueError:
			continue

# sort data based on the Peak_phase(CT) column
data_count = sorted(data_count, key=lambda x:x[-1], reverse=True)

# Build RGB and color palette
yellow = (255,220,0)
white = (255,255,255)
blue = (56,66,156)
R = np.linspace(yellow[0]/255, blue[0]/255, 101)
G = np.linspace(yellow[1]/255, blue[1]/255, 101)
B = np.linspace(yellow[2]/255, blue[2]/255, 101)

# Heatmap panel
y = 0
for value in data_count:
	x = 0
	for position in value[:-1]:
		color = (R[int(position)], G[int(position)], B[int(position)])
		rectangle = mplpatches.Rectangle([x,y],1,1, facecolor=color, linewidth=0)
		panel1.add_patch(rectangle)
		x += 1
	y += 1
	pass


# Circular histogram
for r in np.arange(0.82, 0.98, 0.01):
	xList = []
	yList = []
	for radian in np.linspace(-np.pi/2, np.pi/2, 1000):
		x = np.cos(radian)*r 
		y = np.sin(radian)*r 
		xList.append(-x)
		yList.append(y)
		pass

	panel2.plot(xList, yList, 
				marker = "o",
				markersize = 0.05,
				markeredgewidth = 0, 
				color = "black")

# Histogram 
counts = [i[-1] for i in data_count]
phaseHist = np.histogram(counts, np.arange(0, 26, 2))
arrayCounts = list(phaseHist[0])
tempArray = list(reversed(arrayCounts))

# Add bins to panel2
for i in range(len(phaseHist[0])):
	xList = []
	yList = []
	limit = np.arange(phaseHist[1][i], (phaseHist[1][i] + 2), 0.01)
	lineLimit = np.arange(phaseHist[1][i], phaseHist[1][i] + 2.1, 2) 
	r = (tempArray[i])/100 + 1
	for j in lineLimit:
		scaledVal = (j/24)*(np.pi*2)
		xVal = np.cos(scaledVal+(np.pi/2))*r 
		yVal = np.sin(scaledVal+(np.pi/2))*r
		xIn = [xVal/r, xVal]
		yIn = [yVal/r, yVal]
		panel2.plot(xIn, yIn, 
					markersize = 0,
					markeredgewidth = 0.5,
					color = "black", 
					linewidth = 0.5) 
		pass

	# Fills bins	
	for point in limit:
		scaledVal = (point/24)*(np.pi*2)
		xVal = np.cos(scaledVal+(np.pi/2))*r 
		yVal = np.sin(scaledVal+(np.pi/2))*r
		xList.append(xVal)
		yList.append(yVal)
		xIn = [xVal/r, xVal]
		yIn = [yVal/r, yVal]
		panel2.plot(xIn, yIn, 
					markersize = 0,
					markeredgewidth = 0,
					color = "grey", 
					linewidth = 0.3) 
		pass
	
	# Add histogram top line 
	panel2.plot(xList, yList, 
				markersize = 0,
				markeredgewidth = 0,
				color = "black", 
				linewidth = 0.5)

# Add outer circular lines 
for r in np.arange(2, 5, 1):
	xList = []
	yList = []
	for radian in np.linspace(-np.pi, np.pi, 1000):
		x_pos = np.cos(radian)*r 
		y_pos = np.sin(radian)*r 
		xList.append(x_pos)
		yList.append(y_pos)
		pass

	# Add circular lines
	line, = panel2.plot(xList, yList, 
						markersize = 0,
						markeredgewidth = 0,
						color = "black",
						linewidth = 0.3, 
						linestyle = "--")
	line.set_dashes([6.5, 7, 13, 7])

# Add center circular lines 
for r in [0.8, 1]:
	xList = []
	yList = []
	for radian in np.linspace(-np.pi, np.pi, 1000):
		x_pos = np.cos(radian)*r 
		y_pos = np.sin(radian)*r 
		xList.append(x_pos)
		yList.append(y_pos)
		pass

	panel2.plot(xList, yList,
				markersize = 0,
				markeredgewidth = 0,
				color = "black",
				linewidth = 0.5)

# Add bin labels
panel2.text(-2.219, 0, "100", verticalalignment="center", horizontalalignment="center", fontsize=6)
panel2.text(-3.219, 0, "200", verticalalignment="center", horizontalalignment="center", fontsize=6)
panel2.text(-4.219, 0, "300", verticalalignment="center", horizontalalignment="center", fontsize=6)
# Add centered labels
panel2.text(0, 0, "CT", verticalalignment="center", horizontalalignment="center", fontsize=6)
panel2.text(0, 0.5, "0", verticalalignment="center", horizontalalignment="center", fontsize=6)
panel2.text(0.45, 0.25, "4", verticalalignment="center", horizontalalignment="center", fontsize=6)
panel2.text(0.45, -0.25, "8", verticalalignment="center", horizontalalignment="center", fontsize=6)
panel2.text(0, -0.5, "12", verticalalignment="center", horizontalalignment="center", fontsize=6)
panel2.text(-0.45, -0.25, "16", verticalalignment="center", horizontalalignment="center", fontsize=6)
panel2.text(-0.45, 0.25, "20", verticalalignment="center", horizontalalignment="center", fontsize=6)

# Panel1 axis and labels
panel1.set_xlim([0, 8])
panel1.set_xticks([i + 0.5 for i in range(8)])
warnings.filterwarnings("ignore")
panel1.set_xticklabels(["0", "", "6", "", "12", "", "18", ""])
panel1.set_ylim([0, len(data_count)])
panel1.set_yticks([i for i in range(0, 1201, 200)])
panel1.set_ylabel("Number of genes")
panel1.set_xlabel("CT")

# Panel2 axis
panel2.set_ylim([-4, 4])
panel2.set_xlim([-4, 4])
panel2.set_axis_off()

# Save figure
plt.savefig(outFile, dpi=600)



