
################################################################################
#
#   Author: Carlos Arevalo
#   Email: carevalo0170@gmail.com
#
#   Program description:
#   Program takes a fasta file. The program was written to be run in the 
#   command line, but could be slightly modified to be run in Jupyter notebook. 
#
#   Program execution:
#   python3 /Users/carlosarevalo/Desktop/seq_logos_program.py \
#       -i  /Users/carlosarevalo/Downloads/Splice_Sequences.fasta \
#       -p  /Users/carlosarevalo/Downloads/bases \
#       -s  /Users/carlosarevalo/Downloads/stylesheet.mplstyle \
#       -o  /Users/carlosarevalo/Desktop/seq_logos_test.png 
#
################################################################################

# Required modules
import matplotlib.pyplot as plt 
import matplotlib.patches as mplpatches
import matplotlib.image as mpimg
from collections import Counter
import numpy as np 
import warnings
import random
import os
import argparse

# Argument parser and definitions
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file")
parser.add_argument("-p", "--pngs")
parser.add_argument("-s", "--style_sheet") 
parser.add_argument("-o", "--output_file")

# Read input file and style sheet
args = parser.parse_args()
style_sheet = args.style_sheet
plt.style.use(args.style_sheet)
inFile = args.input_file
inLogos = args.pngs
outFile = args.output_file

# Define figure dimentions
figureHeight=3
figureWidth=6
plt.figure(figsize=(figureWidth, figureHeight)) 

# Define panel dimentions
panelHeight=1
panelWidth=2.4

relativePanelWidth = panelWidth/figureWidth
relativePanelHeight = panelHeight/figureHeight

# Make figure panels
panel1 = plt.axes([0.5/6, 0.3, relativePanelWidth, relativePanelHeight])
panel2 = plt.axes([(2.4+1)/6, 0.3, relativePanelWidth, relativePanelHeight])

# Add mid lines
panel1.plot([10 for i in range(10)], list(range(0,10)), color="black", linewidth=0.5)
panel2.plot([10 for i in range(10)], list(range(0,10)), color="black", linewidth=0.5)

# class FastAreader was reused and implemented from a BME160 assignment
class FastAreader:
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname'''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header, sequence = '',''
        with self.doOpen() as fileH:
            header, sequence = '',''
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()
            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence

# Logos dictionary
logos_dict = {"A":[], "T":[], "C":[], "G":[]}

# Read logos and append them to dictionary
for path, dirs, img, in os.walk(inLogos):
    for logosName in img:
        logoPath = os.path.join(path, logosName)
        if logoPath.endswith("A_small.png"):
            A = mpimg.imread(logoPath)
            logos_dict["A"].append(A)
        elif logoPath.endswith("T_small.png"):
            T = mpimg.imread(logoPath)
            logos_dict["T"].append(T)
        elif logoPath.endswith("C_small.png"):
            C = mpimg.imread(logoPath)
            logos_dict["C"].append(C)
        elif logoPath.endswith("G_small.png"):
            G = mpimg.imread(logoPath)
            logos_dict["G"].append(G)
        else:
            continue

# reads fasta sequences
read_obj = FastAreader(inFile)
read_obj.doOpen()
sequenceDict = {}

spliceSite5 = []
spliceSite3 = []

# Separates sequences by splice site
for header, sequence in read_obj.readFasta():
    headerSplit = header.strip().split('_') 
    spliceDirection = headerSplit[0]
    position = int()
    if spliceDirection == "5'":
        spliceSite5.append(sequence)
    else:
        spliceSite3.append(sequence) 
    pass

# Sequences range distionaries
ss5_Dict = {i:[] for i in range(20)}  
ss3_Dict = {i:[] for i in range(20)}

# Assigns sequence range
for sequence in spliceSite5:
    for pos in range(20):
        ss5_Dict[pos].append(sequence[pos])

for sequence in spliceSite3:
    for pos in range(20):
        ss3_Dict[pos].append(sequence[pos])

# Calculates 5' splice site frequency, ratios, and plots SS5
for base in list(ss5_Dict.keys()):
    height_ss5 = 0
    freq = Counter(ss5_Dict[base])
    for i, t in freq.items():
        freq[i]=float(t/len(ss5_Dict[base]))
        height_ss5 += -freq[i]*np.log2(freq[i])
    ratio5 = 2.0 - height_ss5
    pass

    height = 0
    for key, value in sorted(freq.items(), key=lambda item: item[1]):
        imgArray = np.array(logos_dict[key])
        imgLogo = np.squeeze(imgArray) # remove extra dimensions
        panel1.imshow(imgLogo, extent=[base, base+1, height, height+(freq[key]*ratio5)], aspect="auto") 
        height = height + (freq[key]*ratio5)
    pass

# Calculates 3' splice site frequency, ratios, and plots SS3
for base in list(ss3_Dict.keys()):
    height_ss3 = 0
    freq = Counter(ss3_Dict[base])
    for i, t in freq.items():
        freq[i]=float(t/len(ss3_Dict[base]))
        height_ss3 += -freq[i]*np.log2(freq[i])
    ratio3 = 2.0 - height_ss3
    pass

    height = 0
    for key, value in sorted(freq.items(), key=lambda item: item[1]):
        imgArray = np.array(logos_dict[key])
        imgLogo = np.squeeze(imgArray) # remove extra dimensions
        panel2.imshow(imgLogo, extent=[base, base+1, height, height+(freq[key]*ratio3)], aspect="auto")
        height = height + (freq[key]*ratio3)
    pass

# Figure limits
panel1.set_xlim([0, 20])
panel2.set_xlim([0, 20])
panel1.set_ylim([0, 2])
panel2.set_ylim([0, 2])

# Figure tick labels
warnings.filterwarnings("ignore")
panel1.set_xticklabels(np.arange(-10, 11, 5))
panel2.set_xticklabels(np.arange(-10, 11, 5))
panel2.set_yticks([])

# Figure axis titles
panel1.set_title("5'SS")
panel2.set_title("3'SS")
panel1.set_xlabel("Distance to\nSplice Site")
panel2.set_xlabel("Distance to\nSplice Site")
panel1.set_ylabel("Bits")

# Save figure
plt.savefig(outFile, dpi=600)

