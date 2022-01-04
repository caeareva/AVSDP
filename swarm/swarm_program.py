
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
#   python3 /Users/carlosarevalo/Desktop/swarm_program.py \
#       -i /Users/carlosarevalo/Downloads/test_input_data_3.txt \
#       -o /Users/carlosarevalo/Downloads/swarm_test.png \
#       -s /Users/carlosarevalo/Downloads/stylesheet.mplstyle 
#
################################################################################

# required modules
import matplotlib.pyplot as plt 
import matplotlib.patches as mplpatches
import numpy as np 
import sys,random,math
import argparse

# argument parser and definitions
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file")
parser.add_argument("-o", "--output_file")
parser.add_argument("-s", "--style_sheet") 

# read input file and style sheet
args = parser.parse_args()
style_sheet = args.style_sheet
plt.style.use(args.style_sheet)
inFile = args.input_file
outFile = args.output_file

# define figure dimentions
figureHeight=3
figureWidth=7
plt.figure(figsize=(figureWidth, figureHeight)) 

# define panel dimensions
panelHeight=2 
panelWidth=5
relativePanelWidth = panelWidth/figureWidth
relativePanelHeight = panelHeight/figureHeight

# colors:
colors = ["#9EC9E2", "#D12959", "#CDE5D2", "#CF597E", "#FEB24C", \
    "#E9E29C", "#F2ACCA", "#6CB0D6", "dimgrey", "#F0746E", "grey"]

# swarm_panel
swarm_panel = plt.axes([0.1, 0.2, relativePanelWidth, relativePanelHeight])

# dictionaries and lists with conditions, identity values, samples and quality scores
sub_element = {"1":[], "2":[], "3":[], "4":[], "5":[], "6":[], "7":[], "8":[], "9":[], "10":[], ">10":[]}
colors_dict = {"1":["#9EC9E2"], "2":["#D12959"], "3":["#CDE5D2"], "4":["#CF597E"], \
                "5":["#FEB24C"], "6":["#E9E29C"], "7":["#F2ACCA"], "8":["#6CB0D6"], \
                "9":["dimgrey"], "10":["#F0746E"], ">10":["grey"]}

quality_scores = []
element_qualDict = {}
sample_dict = {}

# read data and except NA values
with open(inFile, "r") as inFile: 
    for line in inFile: 
        try:
            splitline = line.strip().split("\t")
            lineElements = splitline[0].split("_") 
            condition = lineElements[3] 
            identity = float(splitline[1])
            if int(condition) <= 10:
                sub_element[condition].append(float(identity))
            else: 
                condition = ">10"
                sub_element[condition].append(float(identity))
            qualScore = float(splitline[0].split("_")[1]) # quality score
            quality_scores.append(qualScore)

            try:
                element_qualDict[condition].append((identity, qualScore))
            except KeyError:
                element_qualDict[condition] = [(identity, qualScore)]
        except ValueError:
            continue

# swarm function
def swarm_plot(identPoints, bin, pointSize):
    placedPoints = [] 
    d = pointSize/72 # point diameter 
    for yPos in identPoints:
        direction = 0
        xPos = float(bin)
        closePoints = []
        placed = False
        for point in placedPoints:
            newPoint = point[1]*(panelHeight/25)
            new_yPoint= yPos*(panelHeight/25)
            if np.abs(newPoint-new_yPoint) <= d:
                closePoints.append(point)
        if len(closePoints)==0:
            placedPoints.append((xPos, yPos))
            continue
        while placed == False:
            dist_list = []
            direction += 1
            shift = d/5
            for point in closePoints:
                print(point)
                dist = abs(xPos*(panelWidth/11.5)-point[0]*(panelWidth/11.5)) # distance
                dist_list.append(dist)
                # print(min(dist_list))
            min_dist = min(dist_list) # defines minimal distance
            if min_dist < d and direction % 2 == 0:
                xPos = shift*(direction/2) + xPos
            elif min_dist < d and direction % 2 == 1:
                xPos = -(shift*(direction/2)) + xPos
            elif min_dist > d:
                placedPoints.append((xPos, yPos))
                placed = True
                xPos = float(bin)
            pass

    # y and x values
    xValue =[point[0] for point in placedPoints]
    yValue =[point[1] for point in placedPoints]
    # color by bins
    color_pos = int(bin)-1
    swarm_panel.plot(xValue, yValue, 
                    marker = "o",
                    markerfacecolor = colors[color_pos], 
                    markeredgecolor = colors[color_pos], 
                    markersize = 0.8,
                    markeredgewidth = 0, 
                    linewidth = 0)

    # add scatters
    swarm_panel.plot([xPos - panelWidth/12.5, xPos + panelWidth/12.5],
        [np.median(yValue), np.median(yValue)], 
        linewidth=1, color ="red")

# add cutoff line at 95% identity
swarm_panel.plot([0,12], [95,95], lw = 0.5, ls = "--", color = "black")

# labels and ticks
swarm_panel.tick_params(bottom=True, labelbottom=True,
                        left=True, labelleft=True,
                        right=False, labelright=False,
                        top=False, labeltop=False)

for key in sub_element.keys(): # keys are elements 1-11
    subset = np.random.choice(sub_element[key], 1000)
    if key == ">10": 
        swarm_plot(subset, "11", 0.6) 
    else: 
        swarm_plot(subset, key, 0.6) 

# define ticks and axis labels:
ticks = [x for x in range(1,11)]
ticks.append('>10') 
swarm_panel.set_xlim(0.5,11)
swarm_panel.set_xticks(np.arange(1,12))
swarm_panel.set_xticklabels(ticks)
swarm_panel.set_xlabel("Subread Coverage")
swarm_panel.set_ylim(75,100)
swarm_panel.set_yticks(np.arange(75,101,5))
swarm_panel.set_ylabel("Identity %")

# save figure
plt.savefig(outFile, dpi=600)

