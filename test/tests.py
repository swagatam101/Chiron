from chiron import *
import numpy as np
from itertools import cycle

def testGenusComputation(): 
    #should be converted to unittest  
    #this is the form of p1, p2, in contact
    ## CHANGE THIS TO GIVE ANY DATA YOU LIKE 
    #genus 0 
    Data = np.asarray([[0,10], [2,8], [4,6]])
    print 'Known genus is 0 ', 
    genus, G, LoopData, backboneData = computeGenus(Data) 
    
    print "computed genus is: ", genus 
    visualize(genus, G, LoopData, backboneData, Data)

    Data = np.asarray([[0,10], [2,6], [4,8]])
    print 'Known genus is 1 ', 
    genus, G, LoopData, backboneData = computeGenus(Data) 
    
    print "computed genus is: ", genus 
    visualize(genus, G, LoopData, backboneData, Data)

    Data = np.asarray([[0,8], [2,6], [4,10]])
    print 'Known genus is 1 ', 
    genus, G, LoopData, backboneData = computeGenus(Data) 
    
    print "computed genus is: ", genus 
    visualize(genus, G, LoopData, backboneData, Data)
    
    Data = np.asarray([[0,8], [0,2], [2,10], [0, 3], [0, 5], [5, 10]])
    print 'Known genus is 2 ', 
    genus, G, LoopData, backboneData = computeGenus(Data) 
    
    print "computed genus is: ", genus 
    visualize(genus, G, LoopData, backboneData, Data)
    
################################################################################

def testEllipseConstruction(ratio=0.5, a=10, center = [0,0], angles = [0,45, 90], arc_fracs = [1/10.0, 1/4.0, 1/3.0], gridsize = 1000): 
    """
    A test fucntion to understand and check accuracy of finding points on 
    ellpse, clockwise from base, and rotations 
    """
    col_gen = cycle('bgrcmk')
    figsize = (20, 20) 
    fig = plt.figure(figsize=figsize) 
    ax = fig.add_subplot(1,1,1)
    
    all_phi, all_path_length = loopDiagramCreator.precomputeEllipticArc(gridsize, ratio)
    arc_lengths = np.asarray(arc_fracs)*all_path_length[-1]*a
    
    for angle in angles: 
        for i, ar in enumerate(arc_lengths):  
        
            pac = Ellipse(center, 2*a, 2*ratio*a, angle= angle, fill = False, linewidth = 2) #IMPORTANT: this angle is in degrees! 
            ax.add_patch(pac)
            
            x1, y1 = loopDiagramCreator.findxy(all_phi, all_path_length, ar, a, ratio, theta = np.radians(angle)) #this has to bein radians 
            x,y = x1 + center[0], y1 + center[0]
            ax.annotate('x' + str(i), xy=(x, y), xytext=(x+1, y-1),
                        arrowprops=dict(facecolor='black', shrink=0.001, alpha = 0.5),
                        )
            ax.plot([x, center[0]], [y, center[1]], '--', linewidth = 2, color =col_gen.next() )
            
        
    ax.set_xlim([center[0]-2*a, center[0]+ 2*a]) 
    ax.set_ylim([center[1]-2*a, center[1]+ 2*a])
    ax.set_aspect('equal') 
    ax.set_title("test of placement of point along ellipse")
    plt.grid()
    plt.show()
    
    
    