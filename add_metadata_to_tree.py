#!/usr/bin/env python

# Use modern print function from python 3.x
from __future__ import print_function

# Usage
import argparse
from argparse import RawTextHelpFormatter
from ete3 import Tree
import pandas as pd
import sys


parser = argparse.ArgumentParser(
    formatter_class=RawTextHelpFormatter,
    description='Add metadata to a one or multiple Newick tree files for editing with FigTree\n'
		'The metadata should be located in a tsv file with the species labels as the first column\n',
    usage='\n  %(prog)s [OPTIONS] -mapfile MAPPING_TABLE -i TREE(s) -o OUTPUT_TREE')
parser.add_argument('-mapfile', metavar='MAPPING_TABLE', default='/home/qiulab/data/CRF_project/data/gene_alignments/meta_144species_map.txt',
                    help='Specify the path to the file that has the species name mapping information. The default file is meta_144species_map.txt')
parser.add_argument('-i', '--input-trees',nargs='+', metavar='TREES', action='store', 
                    help='Path to the original tree file(s) in Newick format')
parser.add_argument('-o','--output-tree', action='store', 
                   help='Specify the path to the output tree file.')
parser.add_argument('-root','--root-species',
                   help='Name of the species that serves as the root')
parser.add_argument('-collapse','--collapse-clade-by',
                   help='Name of the variable [column header in the mapping file] you want to collapse clades by')

parser.add_argument('--version', action='version', version=
'=====================================\n'
'%(prog)s v0.2\n'
'Updated May-2024 by Maeva Perez\n'
'Dependencies: ete3, pandas\n'
'=====================================')

args = parser.parse_args()

import pandas as pd
import ete3



#### Set variables
meta=args.mapfile
input_trees=args.input_trees
outfile=args.output_tree

#### Load metadata
df=pd.read_csv(meta,header=0,sep='\t',index_col=0)
dic=df.to_dict('index')

##### init nexus file

mytree = ete3.PhyloTree(input_trees[0],format=0)

new_nex='#NEXUS\nbegin taxa;\n\tdimensions ntax='+str(len(mytree.get_leaves()))+';\n\ttaxlabels\n\t'
new_nex+='\n\t'.join([node.name for node in mytree.get_leaves()])  
new_nex+='\n;\nend;\n\nbegin trees;'

###### loop through input trees and edit them

for input_tree in input_trees:
    
    new_nex+='\n	tree '+input_tree+' =  [&R] '
    mytree = ete3.PhyloTree(input_tree,format=0)

    ###### reroot-tree
    if args.root_species:
        # annelids=[str(k) for k in dic.keys() if dic[k]['clade_index'] not in [38,39,40,41]]
        # ancestor = mytree.get_common_ancestor(annelids)
        mytree.set_outgroup(str(args.root_species))

    ##### modify tree
    ## Add attributes to leaves
    for node in mytree.iter_leaves():
        for col in df.columns:
            node.add_feature(col,True)
            new_attr=str(dic[int(node.name)][col])
            setattr(node, col, new_attr)
    for node in mytree.traverse():
        if node not in mytree.get_leaves():
        
                    
            elif len({dic[int(leaf.name)][col] for leaf in node.get_leaves()})==1:
                        
                node.add_feature(col,True)
                node.add_feature('clade_name',True)
                node.clade_index=dic[int(node.get_leaves()[0].name)][col]
                
            else:
                node.add_feature('clade_index',False)
                node.add_feature('clade_name',False)


    ##### collapse tree !!! THIS PART IS HARDCODED TO FIT THE ANNELIDA DATASET 
    
    for i in range(1,41):
        val=df[df.clade_index==i]['clade_name'].values[0]
        try:
            node=list(mytree.iter_search_nodes(clade_index=i))[0]
            node.add_feature('clade_genuspec',val)
            node.add_feature('species_name',val)
            node.add_feature('!collapse','{"collapsed"\,0.0}')

        except IndexError:
            continue
    
    #######
    
    for node in mytree.get_leaves():
        node.name=dic[int(node.name)]['IDXgenuspec']
        
    ##### append new_tree
    my_features=[feat for feat in df.columns]+['!collapse']
    new_tree=mytree.write(features=my_features, format=0)
    new_tree=new_tree.replace('&&NHX:','&')
    for col in my_features:
        new_tree=new_tree.replace(':'+col,','+col)
    
    new_nex+=new_tree+'\n'

#### Append options
# with open('FigTree.options.txt', 'r') as f:
#     options=f.read()

options='begin figtree;\n\tset appearance.backgroundColorAttribute="Default";\n\tset appearance.backgroundColour=#ffffff;\n\tset appearance.branchColorAttribute="User selection";\n\tset appearance.branchColorGradient=false;\n\tset appearance.branchLineWidth=1.0;\n\tset appearance.branchMinLineWidth=0.0;\n\tset appearance.branchWidthAttribute="Fixed";\n\tset appearance.foregroundColour=#000000;\n\tset appearance.hilightingGradient=false;\n\tset appearance.selectionColour=#2d3680;\n\tset branchLabels.colorAttribute="User selection";\n\tset branchLabels.displayAttribute="Branch times";\n\tset branchLabels.fontName="sansserif";\n\tset branchLabels.fontSize=8;\n\tset branchLabels.fontStyle=0;\n\tset branchLabels.isShown=false;\n\tset branchLabels.significantDigits=4;\n\tset layout.expansion=0;\n\tset layout.layoutType="RECTILINEAR";\n\tset layout.zoom=0;\n\tset legend.attribute="label";\n\tset legend.fontSize=10.0;\n\tset legend.isShown=false;\n\tset legend.significantDigits=4;\n\tset nodeBars.barWidth=4.0;\n\tset nodeBars.displayAttribute=null;\n\tset nodeBars.isShown=false;\n\tset nodeLabels.colorAttribute="User selection";\n\tset nodeLabels.displayAttribute="Node ages";\n\tset nodeLabels.fontName="sansserif";\n\tset nodeLabels.fontSize=8;\n\tset nodeLabels.fontStyle=0;\n\tset nodeLabels.isShown=false;\n\tset nodeLabels.significantDigits=4;\n\tset nodeShapeExternal.colourAttribute="User selection";\n\tset nodeShapeExternal.isShown=false;\n\tset nodeShapeExternal.minSize=10.0;\n\tset nodeShapeExternal.scaleType=Width;\n\tset nodeShapeExternal.shapeType=Circle;\n\tset nodeShapeExternal.size=4.0;\n\tset nodeShapeExternal.sizeAttribute="Fixed";\n\tset nodeShapeInternal.colourAttribute="User selection";\n\tset nodeShapeInternal.isShown=false;\n\tset nodeShapeInternal.minSize=10.0;\n\tset nodeShapeInternal.scaleType=Width;\n\tset nodeShapeInternal.shapeType=Circle;\n\tset nodeShapeInternal.size=4.0;\n\tset nodeShapeInternal.sizeAttribute="Fixed";\n\tset polarLayout.alignTipLabels=false;\n\tset polarLayout.angularRange=0;\n\tset polarLayout.rootAngle=0;\n\tset polarLayout.rootLength=100;\n\tset polarLayout.showRoot=true;\n\tset radialLayout.spread=0.0;\n\tset rectilinearLayout.alignTipLabels=false;\n\tset rectilinearLayout.curvature=0;\n\tset rectilinearLayout.rootLength=100;\n\tset scale.offsetAge=0.0;\n\tset scale.rootAge=1.0;\n\tset scale.scaleFactor=1.0;\n\tset scale.scaleRoot=false;\n\tset scaleAxis.automaticScale=true;\n\tset scaleAxis.fontSize=8.0;\n\tset scaleAxis.isShown=false;\n\tset scaleAxis.lineWidth=1.0;\n\tset scaleAxis.majorTicks=1.0;\n\tset scaleAxis.minorTicks=0.5;\n\tset scaleAxis.origin=0.0;\n\tset scaleAxis.reverseAxis=false;\n\tset scaleAxis.showGrid=true;\n\tset scaleBar.automaticScale=true;\n\tset scaleBar.fontSize=10.0;\n\tset scaleBar.isShown=true;\n\tset scaleBar.lineWidth=1.0;\n\tset scaleBar.scaleRange=0.0;\n\tset tipLabels.colorAttribute="User selection";\n\tset tipLabels.displayAttribute="Names";\n\tset tipLabels.fontName="sansserif";\n\tset tipLabels.fontSize=8;\n\tset tipLabels.fontStyle=0;\n\tset tipLabels.isShown=true;\n\tset tipLabels.significantDigits=4;\n\tset trees.order=false;\n\tset trees.orderType="increasing";\n\tset trees.rooting=false;\n\tset trees.rootingType="User Selection";\n\tset trees.transform=false;\n\tset trees.transformType="cladogram";\nend;\n\n'

new_nex+='\nend;\n\n'+options

#### write nexus FigTree file

with open(outfile,'w') as fa:
    fa.write(new_nex.replace('\_',','))
    