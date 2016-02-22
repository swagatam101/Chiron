"""
Example script for analysis of HiCCUPS detected peaks for 
specific loop domains in HiC data. HiCCUPS is part of 
Juicebox (http://www.aidenlab.org/software.html). 
The script shows how to use Chiron for topological 
analysis: 
    1. compute genus aat multiple length scales
    2. visualize genus computation 
    3. visualize the genus zero (planar, "non-crossing") loop structure of 
    chromatin where offending links (that contribute to genus, "crossing")
    are displayed in color over the backbone loop domains



Usage: 
    this is an example script. cd into the current bin directory to run it
"""

import os
import sys
import numpy as np
path_to_chiron = os.path.abspath('..')
if  path_to_chiron not in sys.path:
    sys.path.append(path_to_chiron)

import chiron                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
import pandas as pd                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
from collections import defaultdict 
from collections import Counter 
import matplotlib.pyplot as plt
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
#load the example file provided                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
filename = os.path.join(path_to_chiron, 'bin', 'GSE63525_GM12878_primary_HiCCUPS_looplist.txt')
Data = pd.read_csv(filename, sep = '\t', header = 0) 

#which cromosome to analyze? 
ch = 20 
chrname = str(ch)  
sel_data = Data.loc[((Data['chr1']== chrname) & (Data['chr2'] == chrname)), :]
LinkDataAll = sel_data.as_matrix(columns = ['centroid1', 'centroid2']) #the data for this chrm. 


#what length scale to consider? 
window = 10**10

xleft = np.min(LinkDataAll) #or choose any other genome start location
xright = np.round(xleft + window) #the end location of the genomic region considered 
select_ind_left = np.logical_and(LinkDataAll[:,0] >= xleft, LinkDataAll[:,0] < xright)
select_ind_right = np.logical_and(LinkDataAll[:,1] >= xleft, LinkDataAll[:,1] < xright)
select_ind = np.logical_and(select_ind_left, select_ind_right)
sel_LinkData = LinkDataAll[select_ind,:] #the portion of the data within this genomic window 

############################ Genus computation #################################
#compute genus of the link data 
genus, _,_,_ = chiron.computeGenus(sel_LinkData, gap = 2**16) 
#compute which links should be trimmed  (see manual for trimming order) and which 
#can be kept for reduction to genus zero backbone loop structure  
clusters, kept_clusters, kept_links, trimmed_links = chiron.trimToGenusZero(sel_LinkData) 

############################ Loop domain visualization #########################

# class that does all the eesential computation for loop diagrams    
LD = chiron.loopDiagramCreator(sel_LinkData)
#display the loop diagram 
LD.loopDiagramFig(gap_angle = np.pi/10)

######################## Length scale depependent visualization ################

powers = (np.asarray([[6.0, 6.5, 7,7.5], [-np.inf,-np.inf,-np.inf,-np.inf]])).T
length_scales =np.round(10**(powers)).astype(int)
genus_data, count_loops, xbinsAll, _ = chiron. computeGenusLengthScales(LinkDataAll, powers = powers)
chiron.visualizeGenusLengthScales(genus_data, count_loops, xbinsAll, powers, show_counts = True) 


###################### Analyize all autosomal chromosomes at length-scales,
###################### (takes a bit of time) ###################################

powers = (np.asarray([[5.5,5.75,6,6.25,6.5,6.75,7,7.25,7.5], [0,0,0,0,0,0,0,0,0]])).T
length_scales = np.round(10**powers)

genus_accumulator = defaultdict(list) 
count_accumulator = defaultdict(list) 

random_genus_accumulator = defaultdict(Counter)
for ch in np.arange(1,2): #select fewer chromosomes if you are exploring 
    chrname = str(ch)  
    print "working on chrm ", chrname
    sel_data = Data.loc[((Data['chr1']== chrname) & (Data['chr2'] == chrname)), :]
    LinkDataAll = sel_data.as_matrix(columns = ['centroid1', 'centroid2'])
    
    genus_data, count_loops, xbinsOverlap, genus_data_randomize = chiron.computeGenusLengthScales(LinkDataAll, powers = powers, overlapping = True, compute_randomize = True, num_randomize = 50)
    
    
    for i in range(len(powers)): 
        locdata = np.asarray(genus_data[i])
        loc_counts = np.asarray(count_loops[i])
        locdata = locdata[ ~np.isnan(locdata)]
        loc_counts = loc_counts[~np.isnan(loc_counts)]
        genus_accumulator[i].extend(locdata)
        count_accumulator[i].extend(loc_counts)
        random_genus_accumulator[i] += genus_data_randomize[i] 
        

chiron.plotGenusStats(genus_accumulator, count_accumulator, powers, prob_max = 1.1, random_genus_accumulator = random_genus_accumulator)
chiron.plotGenusLengthScale(genus_accumulator, random_genus_accumulator, powers)
plt.show() 