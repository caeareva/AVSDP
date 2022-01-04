
################################################################################
#
#   Author: Carlos Arevalo
#   Email: carevalo0170@gmail.com
#
#   Program description:
#   Program takes as input a GTF file, two PLS files text files. The program was 
#	written to be run in the command line, but could be slightly modified to be 
#	run in Jupyter notebook. 
#
#	python3 /Users/carevalo/Desktop/genome_browser/genome_browser_program.py \
#        -i1 /Users/carevalo/Desktop/genome_browser/test_input_data_6.psl \
#        -i2 /Users/carevalo/Desktop/genome_browser/test_input_data_6.psl \
#        -g /Users/carevalo/Desktop/genome_browser/refdata-gex-mm10-2020-A-genes_wtMUC1_fsMUC1.gtf \
#        -o /Users/carevalo/Desktop/genome_browser/genome_browser_figure.png \
#        -s /Users/carevalo/Desktop/genome_browser/stylesheet.mplstyle
#
################################################################################


# Required modules
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import matplotlib.image as mpimg
from operator import itemgetter
import numpy as np
import sys
import argparse

# Argument parser and definitions
parser = argparse.ArgumentParser()
parser.add_argument("-i1", "--input_file1")
parser.add_argument("-i2", "--input_file2")
parser.add_argument("-g", "--genome")
parser.add_argument("-o", "--output_file")
parser.add_argument("-s", "--style_sheet") 

# Read style sheet
args = parser.parse_args()
style_sheet = args.style_sheet
plt.style.use(args.style_sheet)

# Read input files
readsFile1 = args.input_file1
readsFile2 = args.input_file2
genomeFile = args.genome
outFile = args.output_file

# Set panels parameters
def panel_params(x_pos, y_pos, figureWidth, figureHeight):
	# defines panel dimensions
	panelHeight=1.25
	panelWidth=10
	relativePanelWidth = panelWidth/figureWidth
	relativePanelHeight = panelHeight/figureHeight
	# sets panel
	panel = plt.axes([x_pos, y_pos, relativePanelWidth, relativePanelHeight], frameon=True)
	panel.set_xticks([])
	# turns all labels and thicks off
	panel.tick_params(bottom=False, labelbottom=False,
                        left=False, labelleft=False,
                        right=False, labelright=False,
                        top=False, labeltop=False)

	return panel

# Read and parse though pls files
def readData(inFile):
	readList = []
	openFile = open(inFile,'r')
	for line in openFile:
		a = line.strip().split('\t')
		chromosome = a[13]
		start = int(a[15])
		end = int(a[16])
		blockstarts = np.array(a[20].split(',')[:-1],dtype=int)
		blockwidths = np.array(a[18].split(',')[:-1],dtype=int)
		read = [chromosome,start,end,blockstarts,blockwidths] 
		readList.append(read)

	return readList

# Read input GTF file
def readGtf(inFile):
    transcriptList=[]
    gtfDict={}
    for line in open(inFile):
        if line[0]!='#':
            a=line.strip().split('\t')
            chromosome=a[0]
            type1=a[2]
            if type1 in  ['exon','CDS']:
                start=int(a[3])
                end=int(a[4])
                transcript=a[8].split(' transcript_id "')[1].split('"')[0]
                if transcript not in gtfDict:
                    gtfDict[transcript]=[]
                gtfDict[transcript].append([chromosome, start, end, type1])

    for transcript, parts in gtfDict.items():
        starts=[]
        ends=[]
        blockstarts=[]
        blockwidths=[]
        types=[]
        for part in parts:
            starts.append(part[1])
            ends.append(part[2])
            blockstarts.append(part[1])
            blockwidths.append(part[2]-part[1])
            chromosome=part[0]
            types.append(part[3])
        transcriptList.append([chromosome, min(starts), max(ends), blockstarts, blockwidths, types])

    return transcriptList

# Plot transcripts
def plotTranscripts(panel, readList, genomicCoord, line_thin, line_thick, width):
    genome_chromosome, genome_start, genome_end = genomicCoord[0], genomicCoord[1], genomicCoord[2]
    plottedReads = [] # Keep list of plotted reads
    for read in readList:
        chromosome,start,end,blockstarts,blockwidths,type1=read[0],read[1],read[2],read[3],read[4], read[5]
        if chromosome == genome_chromosome:
            if genome_start<start<genome_end or genome_start<end<genome_end:
            	# Keep track of plotted reads
            	new_read = [read[0],read[1],read[2],read[3],read[4], False, read[5]]
            	plottedReads.append(new_read)
    
    # Sort plotted reads based on start position
    sortedPlottedReads = sorted(plottedReads, key=itemgetter(2))

    # Make y position go up after every iteration
    for y_pos in range(1, len(sortedPlottedReads), 1):
    	lastPlottedEnd = genome_start-100000
    	for read in sortedPlottedReads:
    		chromosome,start,end,blockstarts,blockwidths,plotted,type1=read[0],read[1],read[2],read[3], \
    																	read[4],read[5],read[6]
    		if plotted is False:
    			if start > lastPlottedEnd:
    				rect = mplpatches.Rectangle((start, y_pos+width),
    											end-start,
    											line_thin,
    											facecolor='black',
    											edgecolor='black',
    											linewidth=0)
    				panel.add_patch(rect)
    				for index in np.arange(0, len(blockstarts), 1):
    					blockstart = blockstarts[index]
    					blockwidth = blockwidths[index]
    					element =  type1[index]
    					if element == "exon":
    						rect1 = mplpatches.Rectangle((blockstart, y_pos+0.12),
    												blockwidth,
    												line_thick,
    												facecolor='black',
    												edgecolor='black',
    												linewidth=0)
    						panel.add_patch(rect1)
    					if element == "CDS":
    						line_thickness = 0.5
    						rect2 = mplpatches.Rectangle((blockstart, y_pos-0.01),
    												blockwidth,
    												line_thickness,
    												facecolor='black',
    												edgecolor='black',
    												linewidth=0)
    						panel.add_patch(rect2)
    				read[5] = True
    				lastPlottedEnd = end

    return y_pos

# Plot reads
def plotReads(panel, readList, genomicCoord, line_thin, line_thick, width):
	genome_chromosome, genome_start, genome_end = genomicCoord[0], genomicCoord[1], genomicCoord[2]
	plottedReads = [] # Keep list of plotted reads
	for read in readList: 
		chromosome,start,end,blockstarts,blockwidths,=read[0],read[1],read[2],read[3],read[4]
		if chromosome == genome_chromosome:
			if genome_start < start < genome_end or genome_start < end < genome_end:
				# Keep track of plotted reads
				new_read = [read[0],read[1],read[2],read[3],read[4], False]
				plottedReads.append(new_read)
		
	# Sort plotted reads based on start position
	sortedPlottedReads = sorted(plottedReads, key=itemgetter(2))

	# Make y position go up after every iteration
	for y_pos in range(1, len(sortedPlottedReads), 1):
		lastPlottedEnd = genome_start
		for read in sortedPlottedReads:
			chromosome,start,end,blockstarts,blockwidths,plotted=read[0],read[1],read[2],read[3],read[4],read[5]
			if plotted is False:
				if start > lastPlottedEnd:
					rect = mplpatches.Rectangle((start, y_pos+width), 
												end-start, 
												line_thin, 
												facecolor='black', 
												edgecolor='black',
												linewidth=0)
					panel.add_patch(rect)
					for index in np.arange(0, len(blockstarts), 1):
						blockstart = blockstarts[index]
						blockwidth = blockwidths[index]
						rect = mplpatches.Rectangle((blockstart, y_pos),
						 							blockwidth,
						 							line_thick, 
													facecolor='black', 
													edgecolor='black', 
													linewidth=0)
						panel.add_patch(rect)
					read[5] = True
					lastPlottedEnd = end

	return y_pos

# Define figure dimentions
figureHeight=5
figureWidth=10
plt.figure(figsize=(figureWidth, figureHeight))

# Set panels: left,bottom, width,height
panel1 = panel_params(0 ,0.65, figureWidth, figureHeight)
panel2 = panel_params(0, 0.35, figureWidth, figureHeight)
panel3 = panel_params(0, 0.05, figureWidth, figureHeight)

# Define genomic region of interest
#genomicCoord = ['chr1', 155184054, 155194688] # MUC1
#genomicCoord = ['chr12', 25205246, 25250929] # KRAS
genomicCoord = ['chr7', 45232945, 45240000]
#genomicCoord = ['chr7', 452329, 4524000]

# Plot transcripts in panel 1
transcriptList = readGtf(genomeFile)
panel1_thin = 0.1
panel1_thick = 0.24
panel1_width = 0.19
bottom = plotTranscripts(panel1, transcriptList, genomicCoord, panel1_thin, panel1_thick, panel1_width)
panel1.set_xlim(genomicCoord[1], genomicCoord[2])
panel1.set_ylim(0.24, 10.05)

# Plot reads for panel 3
readList2 = readData(readsFile1)
panel3_thin = 0.1
panel3_thick = 0.5
panel3_width = 0.2
bottom2 = plotReads(panel3, readList2, genomicCoord, panel3_thin, panel3_thick, panel3_width)
panel3.set_xlim(genomicCoord[1], genomicCoord[2])
panel3.set_ylim(0, 436)

# Plot reads for panel 2
readList1 = readData(readsFile2)
panel2_thin = 0.08
panel2_thick = 0.5
panel2_width = 0.21
bottom1 = plotReads(panel2, readList1, genomicCoord, panel2_thin, panel2_thick, panel2_width)
panel2.set_xlim(genomicCoord[1], genomicCoord[2])
panel2.set_ylim(0.2, 69.5)

# Save figure
plt.savefig(outFile, dpi=1200)

