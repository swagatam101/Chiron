"""
Module containing the topological analysis functions in Chiron for genus computation 
of chromatin links (specific interactions)  as input data (Nx2 each row is [p1, p2]
for link origin and destination in genomic positions. 
"""
from __future__ import division 
import numpy as np
import networkx as nx 
import warnings 
from collections import defaultdict 
from collections import Counter 
import copy 



__title__ = 'Chiron Topology'
__author__ = "Swagatam Mukhopadhyay"  
__copyright__ = "Copyright (C) 2015 Swagatam Mukhopadhyay" 
__license__ = "GPL 3.0" 
__version__ = '0.01' 
__credits__ = ['Anirvan M. Sengupta']
__maintainer__ = 'Swagatam Mukhopadhyay'
__status__ = 'Development' 
__copyright__ = 'Copyright 2016 Swagatam Mukhopadhyay'
__docformat__ = 'pdf'


try: 
    from intervaltree import IntervalTree
except: 
    warnings.warn("Cannot use function removeGenusZeroLinks()", ImportWarning)


######### General utilities ####################################################

def uniquerows(Mat):
    """ 
    Find the unique rows of a matrix 
    
    *Args:*
        Mat: 
            matrix input 
    *Returns:*
        Mat:
            Matrix with unique rows
        idx:
            Index of rows that are unique   
    
    """
    temp = np.ascontiguousarray(Mat).view(np.dtype((np.void, Mat.dtype.itemsize * Mat.shape[1])))
    _, idx = np.unique(temp, return_index=True)
    return Mat[idx], idx



######### Utilities used in topological analysis ###############################
def splitAndModify(left, right, G, address, counter, LoopData, gap = 2**16, s1 = 4, s2 = 1): 
    """
    This "splits" the coincident links (at either origin or termination end) and 
    modifies the graph for genus computation appropriately. See algorithm details
    and computeGenus(). Coincident ends are resolved by "splitting" such that new
    crossings (genus-increasing) of links is not introduced as a result.  
    
    
    *Args:*
        left:
            origin of link
        right:
            termination of link, left < right in genomic positions
        G:
            networkx graph for the link 
        address: 
            address map for genome positions onto the odd and even numbers, see
            computeGenus()
        counter: 
            counter for number of coincident link ends
        LoopData:
            The list of lists of edges that correspond to links established between
            "real" and "ghost" points, as a result of splitting etc. 
        gap:
            The gap on the even number line for mapping and splitting link ends,
            see computeGenus()
        s1: 
            gap/s1 is the window within which all splitting is done 
        s2:
            s2 is the odd step from even "real" points to add "ghost" points, 
            see computeGenus()    
            
    
    *Returns:*
       G: 
           modified networkx graph, result of splitting etc. 
       counter:
           modified... 
       LoopData:
           modified by adding new edges to the graph 
        
    """
    assert left < right, 'left point has to be lesser than right point'
    
    #Data ordergin is such that the following is guranteed: 
    #for all coincident left points the nearest right neighbor is encountered first 
    #for all coincident right point the *farthest* neighbor is encountered first       
                      
    #find all nodes that originated from left and right point, that is all children created by splitting 
    # only points that fall between left and right needs to considered for non crossing. All other are splits are such that 
    #they fall outside this range by construction, so no crossing is possible with these loops       
    if counter[left] == 0: #first time encounter this point so node doesn't exist in graph yet
        insert_pos_L = address[left] + gap/s1
        num_CL = 0 
    else:
        num_CL = len([n for n in G if G.node[n]['parent'] == left and G.node[n]['kind'] == 'real' ])        
        right_children_of_left = [n for n in G if G.node[n]['parent'] == left and G.node[n]['kind'] == 'real' and G.node[n]['pos'] > address[left]]
        
        if len(right_children_of_left): #then there are right children to worry about 
            node_pos_RCL = [G.node[n]['pos'] for n in right_children_of_left]  #all node positions of children  
            sorted_pos = np.sort(node_pos_RCL) 
            #insertion/splitting is always between the pos of Left and the first Right children of left 
            insert_pos_L = sorted_pos[0] - s1 
        else: 
            insert_pos_L = address[left] + gap/s1 #this would be the first splitting (bounding range, every other splitting is within this one)
            #arbitrarily chosen to be 1/s1 fraction
       
    
    assert insert_pos_L >= address[left], "increase gap! or decrease s1!" 
        
        
        
    if counter[right] ==0: 
        insert_pos_R = address[right] - s1 
        num_CR = 0 
    else:
        num_CR = len([n for n in G if G.node[n]['parent'] == right and G.node[n]['kind'] == 'real' ])      
        left_children_of_right = [n for n in G if G.node[n]['parent'] == right and G.node[n]['kind'] == 'real' and G.node[n]['pos'] < address[right]]
        
        if len(left_children_of_right): 
            node_pos_LCR = [G.node[n]['pos'] for n in left_children_of_right]  #all node positions of children  
            sorted_pos = np.sort(node_pos_LCR) 
            #insertion/splitting is always between the pos of right and the last left children of Right 
            insert_pos_R = sorted_pos[0] - s1 
        
        else: 
            insert_pos_R = address[right] - s1 
      
    assert np.abs(insert_pos_R - address[right]) < gap, "increase gap! or decrease s1!" 
          
    left_real = str(left) + '_r_' + str(num_CL)
    left_ghost = str(left) + '_g_' + str(num_CL)
    
    right_real = str(right) + '_r_' + str(num_CR)
    right_ghost = str(right) + '_g_' + str(num_CR)
    
     
    if (insert_pos_L%2 or insert_pos_R %2): 
        raise(AssertionError, "Change gap, or s1 and s2")
    
    G.add_node(left_real, pos = insert_pos_L, kind = 'real', parent = left)
    G.add_node(left_ghost, pos = insert_pos_L + s2, kind = 'ghost', parent = left) 
    
    G.add_node(right_real, pos = insert_pos_R, kind = 'real', parent = right) 
    G.add_node(right_ghost, pos = insert_pos_R + s2, kind = 'ghost', parent = right) 
    
    G.add_edge(left_real, right_ghost) #see rule 1 and 2 
    G.add_edge(left_ghost, right_real) 
    
    #stores the points that are connected by loops, as oppoed to backbones 
    LoopData.append([left_real, right_ghost]) 
    LoopData.append([left_ghost, right_real]) 

    counter[left] += 1 
    counter[right] += 1 
    
    return G, counter, LoopData 
    
    
    
################################################################################



def computeGenus(Data, gap = 2**16): 
    """    
    Main function to compute the genus given Data (list of links, N x 2 cotacting 
    pairs of genomic positions).    
    
    Creates a graph G to compute the number of connected components and genus. 
    First consider the case of no coincident ends for loop origins and terminations:  
     
    Then each end in the link list would split into two, "real" (r) and "ghost" (g) 
    where address of ghost on the real line is greater than address of the "real".
    
    Again, in the absence of coincident ends for each link: 
        
    1. The left ends "real" node shares an edge to the right end's "ghost" node
    
    2. The left ends "ghost" node shares an edge to the right end's "real" node, 
       exhausting edges correspoding to links 
       
    3. Along the real line, only "ghost" nodes connect by edge to "real" nodes, 
       in linear order, and in consecutive pairing along the real line (backbone)
        
    4. Count the number of original loops = P (before creating ghosts). Call it P 
    
    5. Count the number of loops (connected components) in the real + ghost graph, 
       call it L 
       
    6. genus :math:`g = (P - L)/2`  
    
    
    Now coming to resolving coincident ends in a manner that introduces no new 
    crossings and doesn't increase genus: 
         
    1. Coincident ends (with n link originating or terminating) will have to be 
       split into n real and n ghost nodes
       
    2. This splitting has to be done in an order such that the splitting itself 
       does not create new link crossings. 
    
    Need to have a strategy for creating nodes such that points are easily ordered. 
    
    Strategy: 
        
    1. Index all original link ends (nodes of G) by large even integers 
    
    2. Create ghosts on large odd numbers 
    
    3. Introduce new real nodes for coincident points in between these large even 
       numbers
       
    4. Ghosts' addresses are always s2 (here 1) greater than reals 
    
    5. gap/s1 (s1 is an option in splitAndModify() function) is the region within 
       which all coincident ends are resolved, increase it if there
       are too many coincident ends
       
    *Args:*
        
        Data: 
            Nx2 link data
        gap:
            Gap between addresses of nodes corresponding to the ends of links 
           
    *Returns:*
        genus:
            computed genus 
        G: 
            networkx graph for computing genus
        LoopData:
            The list of edges corresponding to mapping of links 
        backboneData:
            The list of edges corresponding to mapping of connectivity along the 
            genome 
    """
    
    # cleam up operations 
    # step 1: Order, left < right point along rows, 
    #order rows by left and then by right, so that coincident points have increasing right link 
    Data = np.sort(np.asarray(Data), axis = 1) #sorted for each row so the left point is lesser than right 
    Data, indx = uniquerows(Data) #clean up operation, will use length of Data
    # in genus computation, better have all loops to be unique 
    #print Data         
    Data = Data[np.lexsort((Data[:,1], Data[:,0]))] #this soorts data by first column and then by second column

    G = nx.Graph() 

    points = np.sort(np.unique(Data)) #unique points 
    counter = dict.fromkeys(points, 0) #this is the counter of number of coinident points 
    address = dict(zip(points, np.arange(0, len(points)*gap, gap))) #initialize dict of address for orignial points, 
    #with the gapped index along the line
    #print address 
    LoopData = [] #stores the data to plot in the p1 p2  format  
    for p1, p2 in Data: #order the loop data by left, right in chormosome position order 
        if p1 == p2: 
            raise ValueError('Loop cannot be zero length')     
     
        G, counter, LoopData = splitAndModify(p1, p2, G, address, counter, LoopData, gap = gap) 

    #now run through the points in order of position and introduce the real to ghost backbone edges 

    sorted_graph = sorted(G.nodes(data = True), key = lambda (a, dct): dct['pos']) 
    sorted_names = np.asarray([n for n,dct in sorted_graph])
    #recall, sorted names are real-ghost pairs  
    backbone_ghosts = sorted_names[1:-1:2]  #the first and the last points always 
    #create a cluster together, sorted_names[1] is the first ghost
    backbone_reals = sorted_names[2:-1:2] 
    #error handling, check that backbone_reals are all reals, surely then all ghosts are ghosts 
    test_kinds = np.asarray([G.node[n]['kind'] for n in backbone_reals], dtype = str) 
    #print backbone_reals
    #print backbone_ghosts 
    #print sorted_names 
    assert np.all(test_kinds == 'real'), 'Fatal Error, graph construction is wrong, change gap?' 
    
    #these are guranteed to be of equal length, but for sanity, throw error otherwise
    assert len(backbone_ghosts) == len(backbone_reals),  "Creation of ghosts was wrong" 
    
    backboneData = []         
    for p1, p2 in zip(backbone_ghosts, backbone_reals): 
        G.add_edge(p1, p2)       
        backboneData.append([p1, p2])

    genus = (len(Data) - (nx.number_connected_components(G) -1))/2 
   
    return genus, G, LoopData, backboneData  
    
    
################################################################################

def randomizeLinks(LinkDataAll, chrmRange = None): 
    """
    Randomize LinkData, where randomization is done by maintaining the number of
    links and their lengths, but just scrambling the ends randomly 
    
    *Args:*
        LinkDataAll:
            Nx2 links, rows (p1,p2) 
            
    *Returns:*    
        LinkDataRandom:
            Nx2 random links, rows (q1,q2)
    """
    LinkDataRandom = np.zeros_like(LinkDataAll)
    if chrmRange is None: 
        chrmRange = [np.min(np.ravel(LinkDataAll)), np.max(np.ravel(LinkDataAll)) + 1] #this is the range in which 
        #the random int is drawn for the randomized links, right boundary is exclusive 
        
    loop_len = np.abs(LinkDataAll[:,1] - LinkDataAll[:,0])
    assert np.all(loop_len > 0), "Loop lengths must be > 0"
    
    p1s = np.random.randint(chrmRange[0], high = chrmRange[1], size = len(LinkDataAll))
       
    p2s = p1s  + loop_len
    bad_ind = p2s >= chrmRange[1] 
    p2s[bad_ind] = p1s[bad_ind] - 2*loop_len[bad_ind] #factor of 2 because you had added loop length before...
    #we are trying to fold it back on the other side 
    LinkDataRandom[:,0] = p1s
    LinkDataRandom[:,1] = p2s
    loop_len = np.abs(LinkDataRandom[:,1] - LinkDataRandom[:,0])
    assert np.all(loop_len > 0), "Random Loop lengths must be > 0, fatal error" #just sanity check 
    return LinkDataRandom  
    
################################################################################
    
def computeGenusLengthScales(LinkDataAll, powers = (np.asarray([[6.0, 6.5, 7,7.5], [-np.inf,0,0,0]])).T,
    overlapping = False, compute_randomize = False, num_randomize = 10):
    """
    Computes genus at different lengthscales, returns two dicts, of genus and of 
    counts of loops found in windows of lengthscale steps. 
    
    *Args*:
        LinkDataAll: 
            Nx2 links, rows (p1,p2) 
        powers:
            log10(length scales), first row is long distance scale of link sizes, 
            second column is short distance scale, analysis is for pairs (rows) 
        overlapping: 
            Binary, if true creates windows that overlap one-half window size
        compute_randomize: 
            Binary, whether to randomize links and compute genus for each window 
        num_randomize: 
            Number of times to randomize the links 
             
    *Returns*:
        genus_data: dict of lists of genus computation, keys are length scales, 
            List is for each window along the genome 
        count_loops:
            Dict of lists, counts of links in each window, same keys as above 
        xbinsAll: 
            Dict of dict of list. Highest dict is 'left' and 'right', for left
            and right edges of windows, the next level dict is identical to above
        genus_data_random:
            Dict of list, where the list is a counter for each genus value (integers)  
    """
    xRange = [np.min(LinkDataAll[:,0]), np.max(LinkDataAll[:,1])] 
    genus_data = defaultdict(list) 
    count_loops = defaultdict(list) 

    xbinsAllL = defaultdict(list)
    xbinsAllR = defaultdict(list)

    genus_data_random = dict()  
    length_scales  = np.round(10**(powers)).astype(int)     
    for i,(upper_ls,lower_ls) in enumerate(length_scales): #run over the long distance cutoff  
        
        loc_counter = Counter()
        
        if not overlapping:  
            xbins = np.arange(xRange[0],xRange[1],upper_ls)   
            xbins = np.append(xbins, xRange[1]+1) #all the bins within which the 
                #contacts are considered
            xbinsL = xbins[:-1]
            xbinsR = xbins[1:]
            
        else: 
            xbinsL = np.arange(xRange[0],xRange[1],np.round(0.5*upper_ls)) 
            xbinsR = xbinsL + upper_ls  
            

        loc_data = [np.nan]*(len(xbinsL)) 
        loc_count = [np.nan]*(len(xbinsL))
        
        for j,x in enumerate(xbinsL): 
            #find all contact point pairs that fall within this range
            select_ind_left = np.logical_and(LinkDataAll[:,0] >= x, LinkDataAll[:,0] < xbinsR[j])
            select_ind_right = np.logical_and(LinkDataAll[:,1] >= x, LinkDataAll[:,1] < xbinsR[j])
            select_ind = np.logical_and(select_ind_left, select_ind_right)
            
            if lower_ls > 0: #short distance cutoff  
                pass_ind = np.abs(LinkDataAll[:,0]- LinkDataAll[:,1]) > lower_ls
                select_ind = np.logical_and(select_ind, pass_ind)
                
            if np.count_nonzero(select_ind): 
                loc_LinkData = LinkDataAll[select_ind,:] - x 
                genus, G, LoopData, backboneData = computeGenus(loc_LinkData, gap = 2**16) 
                loc_data[j] = genus
                loc_count[j] = len(loc_LinkData)
                if compute_randomize and np.count_nonzero(select_ind) >1: 
                    for k in range(num_randomize):  
                        loc_link_random = randomizeLinks(loc_LinkData, chrmRange = [0, upper_ls])
                        genus, G, LoopData, backboneData = computeGenus(loc_link_random, gap = 2**16) 
                        loc_counter[genus] += 1 
                        
        genus_data[i] = loc_data
        count_loops[i] = loc_count 
        genus_data_random[i] = loc_counter 
        xbinsAllL[i] = xbinsL
        xbinsAllR[i] = xbinsR
    
   
    xbinsAll = dict()
    xbinsAll['left']  = xbinsAllL
    xbinsAll['right'] = xbinsAllR
            
    return genus_data, count_loops, xbinsAll, genus_data_random 




################################################################################

def trimToGenusZero(LinkData): 
    """
    Culls link such that the new set of links is guranteed to be genus zero, 
    maximally crossing links are culled first, with preference to retaining long 
    range links.
    
    *Args:*
        LinkData: 
            Nx2 links, rows (p1, p2) 
            
    *Returns;*
         clusters: 
             Clusters of links that isolate along the genone 
         kept_clusters:
             Clusters after trimming
        kept_links: 
            Links retained after culling, genus zero links
        trimmed_links:
            Links culled, genus changing 
    """ 
    # need to find the set of all overlapping regions and of disjoint regions   
    # order the bonds first by the left orgin and then by the right origin 
    # then find footporint of all overlapping regions (clusters) 
    clusters= defaultdict(list)
    Data = np.sort(LinkData, axis = 1) #sorted for each row 
    Data, indx = uniquerows(Data) #clean up operation, will use length of Data
    Data = Data[np.lexsort((Data[:,1], Data[:,0]))] #this sorts data by first column and then by second column
    counter = 0         
    clusters[0].append([Data[0,0], Data[0,1]])
    
    footprint = [Data[0,0], Data[0,1]] #inital footprint is first bond  
    
    for i, (p1, p2) in enumerate(Data[1:]): 
    
        if p1 <= footprint[1]:  #i.e., the region overlaps the previous region
            clusters[counter].append([p1, p2])
            footprint[1] = np.max((p2, footprint[1]))  # thelexsort above ensures that p1 > footprint[0]
        else: 
            counter +=1
            footprint = [p1, p2] #new singleton cluster to grown, hence footprint is the new link 
            clusters[counter].append([p1,p2])
            
                        
    genus_all = np.zeros(len(clusters))
    
    kept_clusters =  defaultdict(list)                                              
    trimmed_links = []
    kept_links = []                                                                                               
                                                                                                                                                                                                
    for c in clusters: 
        genus, _, _,_ = computeGenus(np.asarray(clusters[c]))  
        genus_all[c] = genus      
        
        loc_data = clusters[c] 
        num_entries = len(loc_data)
        
        if num_entries > 1 and genus > 0:
            scores = np.zeros((num_entries, 4), dtype = np.int)     
            
            loc_data = np.asarray(loc_data)
            for i, (p1, p2) in enumerate(loc_data):   
                #find the number of crossings!!
                right_overlap = np.logical_and(np.logical_and(loc_data[:,0] >= p1, loc_data[:,0] <= p2), loc_data[:,1] > p2) 
                left_overlap = np.logical_and(np.logical_and(loc_data[:,1] >= p1, loc_data[:,1] <= p2), loc_data[:,0] < p1)
                scores[i,:] = (p1, p2, np.count_nonzero(np.logical_or(left_overlap, right_overlap)), np.abs(p1-p2))
                
            #now sort scores, but first select only those with >0 crossings 
            scores_S = scores[scores[:,3]  > 0]  
            sort_order = np.lexsort((scores_S[:,3], -scores_S[:,2]))  #this sorts first by crossing in descending and then by range of link in ascending
            
            ind_to_remove = []
            genus_running = copy.copy(genus)
            pointer = 0
            
            #print '---------------------------'
            #print
            #print sort_order
            #print
            #print scores_S
            #print scores_S[sort_order]
            
            
            while genus_running > 0: 
                ind_to_remove.append(sort_order[pointer])
                rem_loc = np.delete(loc_data, ind_to_remove, axis = 0)
                genus_running ,_ ,_ ,_ = computeGenus(rem_loc)
                pointer += 1
                    
            ind_to_keep = [x for x in range(len(loc_data)) if x not in ind_to_remove]        
            
            
            to_keep = [[i,j] for (i,j) in loc_data[ind_to_keep]]
            to_remove = [[i,j] for (i,j) in loc_data[ind_to_remove]]
            kept_clusters[c] = to_keep
            kept_links.extend(to_keep)
            trimmed_links.extend(to_remove)
        else: 
            kept_clusters[c] = loc_data
            kept_links.extend(loc_data)    
                
    return clusters, kept_clusters, kept_links, trimmed_links
    
    
    
################################################################################    
def removeGenusZeroLinks(LinkData): 
    """
    Removes all links that are are equivalent toplogically to other links and
    do not contribute to genus! Uses IntervalTree to find all crossing links to a 
    given link (in O(N log N) vs. O(N^2) for N links). 
    
    Required IntervalTree package! 
    
    *Args*:
        LinkData: 
            Nx2 links, rows (p1,p2) 
    *Returns*:
        removal_linkData:
            The list of links retained 
    """
    
    org_tree = IntervalTree.from_tuples(list(map(tuple, LinkData)) )  
    removal_tree = IntervalTree.from_tuples(list(map(tuple, LinkData)) )
    
    for m,t in enumerate(org_tree): 
        
        loc_set = org_tree[t[0]:t[1]]
    
        for l in loc_set: 
            if t[0] >= l[0] and t[1]<= l[1] and t!=l: #find if i is contained in l
                
                #now find the intervals overlapping Interval(i[0], l[0]) and i[1], l[1]
                left_set = org_tree[l[0]:t[0]]
                right_set = org_tree[t[1]:l[1]]
                if len(right_set) == 1 and len(left_set) ==1: # the right and left set has one overlap: 
                    # the interval l 
                    removal_tree.remove(t)
                
    removal_linkData = [[i[0], i[1]] for i in removal_tree]
    return removal_linkData    
            
################################################################################