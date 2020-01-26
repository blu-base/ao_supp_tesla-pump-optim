#!/usr/bin/python3
"""
Extended Scatter Matrix plotter for OPAL++ 

This script creates an extended scatterplot matrix for a FullArchive's Gen.csv which
has been created during an Optimization performed with OPAL++.

It autmatically extracts the number objectives and parameters and draws them completely.
In default state, this script does filter by Validity-Value and Contrain Violation.
Further filters on the data can be applied in the Custom Filter section

Main libraries for this script is seaborn, pandas, matplotlib, and numpy.

This work is licensed under CC BY-SA 3.0, https://creativecommons.org/licenses/by-sa/3.0/
Originally posted on Stack Overflow https://stackoverflow.com/a/50690729 by joelostblom https://stackoverflow.com/users/2166823/joelostblom
Adapted and extended to use with OPAL++ data by Sebastian Engel, 2018

"""
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np



import matplotlib as mpl
mpl.use("pgf")

plt.rcParams.update({
    "pgf.texsystem": "lualatex",
    "pgf.preamble": [
#         r"\usepackage[utf8x]{inputenc}",
#         r"\usepackage[T1]{fontenc}",
         r"\usepackage{cmbright}",
         ]
})

genCSVLocation = './Gen.csv'
outputfileName = 'extScatterPlot'

## Switches for validity and contrain violation filter
filterOutInvalids           = True
filterOutContrainViolations = True

## optional Variable. Sets customs headers for the plot indead of o1,o2,o3...,r1,r2,r3,...
## However, this header is only used to draw the plot. Within the script, the default notation (o1,o2,...) is kept.
customHeader = ['ETA','HEA','GAP','RIN','ROU','VOL','RTO']
#customHeader = ['$\eta_{hydr}$','H','gap','VOL','RTO']

## Sanitize Gen.csv Header. There is an one whitespace too much in the header by default...
def inplace_change(fileName,strOld,strNew):
    with open(fileName) as f:
        s = f.read()
        if strOld not in s:
            return
        with open(fileName, 'w') as f:
            s = s.replace(strOld,strNew)
            f.write(s)
inplace_change(genCSVLocation,"Constrain violation", "constrainViolation")


## extracting the amount of objectives and parameters
numObj = 0
numPar = 0
with open(genCSVLocation, 'r') as f:
    s = f.readline()
    cols = s.split(' ')
    for col in cols:
        if col.startswith('o'):
            numObj+=1
        if col.startswith('r'):
            numPar+=1
usecolsInGen = [1,2]
usecolsInGen.extend([ 3 + i for i in range(numObj+numPar)])

## Reads complete CSV
genDF = pd.read_csv(genCSVLocation, usecols=usecolsInGen,delim_whitespace=True,dtype={'Validity':int})


## Filter out invalid data 
if filterOutInvalids == True:
        genDF = genDF[genDF['Validity'] >0]
if filterOutContrainViolations == True:
        genDF = genDF[genDF['constrainViolation'] == 0.0]
## Drop Meta columns for good
genDF = genDF.drop(['Validity'], axis=1)
genDF = genDF.drop(['constrainViolation'], axis=1)


################################################################################
## Custom Filters, see pandas DataFrame documentation for more information
################################################################################
## Setting custom column headers before custom filters. 
## comment these lines out, only if you want to use custom header instead of
## o1,o2, r1,r2,... and so on.
#if 'customHeader' in locals():
#    genDF.columns = customHeader

## Custom Filters, see pandas documentation for syntax
#genDF = genDF[genDF['o1'] > 0 ]
#genDF = genDF[genDF['o1'] < 1300000 ]
#genDF = genDF[genDF['o2'] > 0 ]
#genDF = genDF[genDF['o3'] > 1 ]


## Outlier Detection. Comment out if unwanted
## keeps only the ones that are within +3 to -3 standard deviations in the respective column.
#genDF = genDF[genDF.apply(lambda x: np.abs(x - x.mean()) / x.std() < 3).all(axis=1)]

## Setting custom column headers.
if 'customHeader' in locals():
    genDF.columns = customHeader

## Printing out some information to STDOUT
genDF.info()
genDF_MinMax =[]
print("Minimum and Maximums")
print("Key\tMin\tMax")
for col in genDF.keys():
    colmin = genDF[col].min()
    colmax = genDF[col].max() 
    print(col +'\t'+ colmin.astype(str) + '\t' + colmax.astype(str))
    genDF_MinMax.append([colmin,colmax])
################################################################################
## Helper functions
## function calculating the correlations in upper right triangle. 
## Draws correlation value in circles which size and color is dependent on the correlation. 
## High (absolute) correlation means large circle and font; and the opposite.
################################################################################
def corrdot(*args, **kwargs):
    corr_r = args[0].corr(args[1], 'pearson')
    corr_text = round(corr_r, 2)
    ax = plt.gca()
 #   ax.set_axis_off()
    font_size = abs(corr_r) * 40 + 5
    ax.annotate(corr_text, [.5, .5],xycoords="axes fraction",
            ha='center', va='center', fontsize=font_size)
    ## Optional colorful circle around correlation values
    marker_size = abs(corr_r) * 10000
    ax.scatter([.5], [.5], marker_size, [corr_r], alpha=0.6, cmap="coolwarm",vmin=-1, vmax=1, transform=ax.transAxes)

## Function to add Data Series name. Used in diagonal of the plot matrix
def annotate_colname(x, **kws):
    ax = plt.gca()
    ax.annotate(x.name, xy=(0.5,0.5), xycoords=ax.transAxes, fhorizontalalignment='center',verticalalignment='top', ontsize=15)
    #ax.text(0.5,0.9,, horizontalalignment='center',verticalalignment='top', transform=ax.transAxes,fontsize=15)



################################################################################
## Main Properties of the Plot
################################################################################
## Setting background style
sns.set(style='white', font_scale=1)

## Drawing the plot
g = sns.PairGrid(genDF, aspect=1, diag_sharey=False, despine=False)
g.map_lower(sns.regplot, order=1, truncate=True, ci=95, line_kws={'color': 'red', 'lw': 2},scatter_kws={'color': 'black', 's': 10})
g.map_diag(sns.distplot, color='black', bins = 20, kde_kws={'bw':'scott','color': 'red', 'cut': 0.7, 'lw': 2}, hist_kws={'histtype': 'bar', 'lw': 1, 'edgecolor': 'black', 'facecolor':'green'})
g.map_upper(corrdot)
#g.map_diag(annotate_colname)
g.fig.subplots_adjust(wspace=0.1, hspace=0.1)



# Add titles to the diagonal axes/subplots
for ax, col in zip(np.diag(g.axes), genDF.columns):
    #ax.set_title(col, y=0.5, fontsize=26)
    ax.annotate(col, xy=(0.5,0.8), xycoords=ax.transAxes, fontsize=26)


from math import log10, floor
def round_to_1(x):
    return round(x, -int(floor(log10(abs(x)))))

################################################################################
## Improving axes formats
## specific plots can be accessed by g.axis[i,j]
################################################################################
for i in range(len(genDF_MinMax)):
    n = len(genDF_MinMax) -1  # max index of plot matrix
    axmin = genDF_MinMax[i][0]
    axmax = genDF_MinMax[i][1]
    axmodifier = (axmax - axmin) * 0.05
    steps = 1
    axsteps = (axmax -axmin) / steps
    axticks = [round_to_1(axmin + j*axsteps) for j in range(steps+1)]
    # axes in most left column of plots
    g.axes[i,0].set_ybound(lower=axmin-axmodifier,upper=axmax+axmodifier)
    g.axes[i,0].set_ymargin(0.05)
    g.axes[i,0].set_yticks(axticks)
    g.axes[i,0].minorticks_off()
    g.axes[i,0].ticklabel_format(style='sci', axis='both', scilimits=(-3,-7),useMathText=True)
    g.axes[i,0].tick_params(axis='y',direction='out', length=5,width=1,bottom=True,left=True)
    # axis in bottom row of columns
    g.axes[n,i].set_xbound(lower=axmin-axmodifier,upper=axmax+axmodifier)
    g.axes[n,i].set_xmargin(0.05)
    g.axes[n,i].set_xticks(axticks)
    g.axes[n,i].minorticks_off()
    g.axes[n,i].tick_params(axis='x',direction='out', length=5,width=1,bottom=True,left=True)
    g.axes[n,i].ticklabel_format(style='sci', axis='both', scilimits=(-3,-7),useMathText=True)
    g.axes[n,i].xaxis.set_tick_params(rotation=90)
    print(axticks)


################################################################################
##  Saving the File
################################################################################
## as png
g.savefig(outputfileName + '.png')
## as latex' pgf
#g.savefig(outputfileName + '.pgf')


plt.rcParams['svg.fonttype'] = 'none'
g.savefig(outputfileName + '.svg')
#plt.rcParams['svg.fonttype'] = 'svgfont'
#g.savefig(outputfileName + '2.svg')
#g.savefig(outputfileName + '.pdf')
## in live view
#plt.show()

#plt.style.use("ggplot")

#import tikzplotlib

#tikzplotlib.save("extScatterPlot.tex",figureheight = '\\figureheight', figurewidth = '\\figurewidth')
