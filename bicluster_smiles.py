import sys
import os
import numpy as np
import rdkit
from rdkit import Chem, DataStructs, RDConfig
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole #Needed to show molecules
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions #Only needed if mod
import matplotlib.pyplot as plt
import matplotlib.colors as mcol
import matplotlib.cm as cm
import seaborn as sns
import pandas as pd
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-x", "--xcomp", dest="infile1", help="The x-axis compounds", metavar="FILE")
parser.add_option("-y", "--ycomp", dest="infile2", help="The y-axis compounds", metavar="FILE")
parser.add_option("-o", "--out", dest="outfile", help="The png outputfile", metavar="FILE")
parser.add_option("-s", "--subi", default=100, type=int, dest="subImgSize", help="The compound subimage size")
parser.add_option("-f", "--fsiz", default=1.0, type=float, dest="fontsize", help="The font size")
parser.add_option("-a", action="store_true", default=False, dest="annot")
(options, args) = parser.parse_args()

sns.set(font_scale=options.fontsize)

#calculate 2048bit morgan fingerprints, radius 2
def calcFingerprints(m):
    return AllChem.GetMorganFingerprintAsBitVect(m,2, nBits=2048)

#calculate mols and fps for each infile
def process_infile(in_file):
    smiles = open(in_file).read().splitlines()
    mols= [Chem.MolFromSmiles(smile) for smile in smiles]
    fps = [calcFingerprints(mol) for mol in mols]
    return smiles, mols, fps

#import infiles
if os.path.exists(options.infile1):
    f1_smiles, f1_mols, f1_fps = process_infile(options.infile1)
else:
    print 'Could not find first input file.'
    quit()
if os.path.exists(options.infile2):
    f2_smiles, f2_mols, f2_fps = process_infile(options.infile2)
else:
    print 'Could not find second input file.'
    quit()
if os.path.splitext(options.outfile)[1][1:] != 'png':
    print "outfile  doesn't end with png"
    quit()

#calculate similarity matrix
sim_matrix = [DataStructs.BulkTanimotoSimilarity(fp,f1_fps) for fp in f2_fps]

#generate clustermap
cmap = mcol.LinearSegmentedColormap.from_list("n",[sns.xkcd_rgb["denim blue"],sns.xkcd_rgb["pale red"]])
norm = plt.Normalize(0,1)
g = sns.clustermap(sim_matrix, annot=options.annot, figsize=(10, 10), cmap=cmap, norm=norm, method='complete', cbar_kws={"label": "Tanimoto Coefficient\n(Tc) Similairty\n"})

#reorder inputs into x & y cluster orders
y_order=[int(tick.get_text()) for tick in reversed(g.ax_heatmap.get_ymajorticklabels())]
g.ax_heatmap.set_yticklabels('')
x_order=[int(tick.get_text()) for tick in g.ax_heatmap.get_xmajorticklabels()]
g.ax_heatmap.set_xticklabels('')

#calculate position of axes to draw mols
xl, yl, xh, yh=np.array(g.ax_heatmap.get_position()).ravel()
w=xh-xl
h=yh-yl
size=0.5
siz = 100
subImgSize = (options.subImgSize,options.subImgSize)

#draw mols on x
ax1=g.fig.add_axes([xl, yh-xh-0.05, w, h])
ax1.axison = False
imgplot = ax1.imshow(Draw.MolsToGridImage([f1_mols[idx] for idx in x_order],molsPerRow=len(f1_mols),subImgSize=subImgSize))

#draw mols on y
ax2=g.fig.add_axes([yh-0.075, yl, w, h])
ax2.axison = False
imgplot = ax2.imshow(Draw.MolsToGridImage([f2_mols[idx] for idx in y_order],molsPerRow=1,subImgSize=subImgSize))

#save image
g.savefig(options.outfile,dpi=300,format='png')
