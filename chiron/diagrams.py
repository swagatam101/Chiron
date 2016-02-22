"""
Module containing the visualization function in Chiron for genus computation 
visualizations of chromatin loop domains. Input data is list of links identified 
to be in two-point specific contacts (by other programs)  
"""
from __future__ import division 
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx 
from collections import defaultdict 
import scipy.special as scsp 
import copy 
import bisect
import matplotlib.patches as patches
from matplotlib.patches import Ellipse
import matplotlib.cm as cmx
from pylab import get_cmap
from matplotlib import colors
from matplotlib.path import Path
import matplotlib.gridspec as gridspec
from matplotlib import colors

from .topology import * 

__title__ = 'Chiron Diagrams for Topology'
__author__ = "Swagatam Mukhopadhyay"  
__copyright__ = "Copyright (C) 2015 Swagatam Mukhopadhyay" 
__license__ = "GPL 3.0" 
__version__ = '0.01' 
__credits__ = ['Anirvan M. Sengupta']
__maintainer__ = 'Swagatam Mukhopadhyay'
__status__ = 'Development' 
__copyright__ = 'Copyright 2016 Swagatam Mukhopadhyay'
__docformat__ = 'pdf'



################################################################################
def computeMajor(ratio, circum): 
    """
    Use the complete elliptic integrals of second kind to figure out what value
    of major axis corresponds to the circumference, given b/a ratio 
    
    :math:`circum = 4 a E(k)` where :math:`k = \\sqrt{ 1- b^2/a^2 } = \\sqrt{1 - ratio^2}` 
    therefore, :math:`a = circum/(4E(k))`
    
    *Args:*
        ratio: 
            ratio of semi-minor to semi-major axis of Ellipse, ratio = b/a <= 1 
        circum: 
            circumference/perimeter of ellipse 
    
    *Returns:*
        a: 
            semi-major axis of ellipse 
    """
    k = np.sqrt(1- ratio**2)
    a = circum/(4*scsp.ellipe(k))
    
    return a 

################################################################################

                                                                        
def createRainbow(LinkDataAll, backbone = None, xRange = None, shiftRatio = 0.01, yMaxFrac = 2, ax = None, unit = 1, title = '', cmap = 'spectral'):   
    """
    Creates custom "rainbow plots" displaying specific interaction pairs, has 
    the option of diplaying weights of interactions. Can display genus 
    computation connected-components, grey, or color.  
    
    *Args:*
        LinkDataAll: 
            List of lists [p1, p2] or [p1, p2, w12] where p1 and p2 
            are genomic locations of interacting points. LinkData all is Nx2 or Nx3
            array, weights w12 are optional. For connected-components w12
            the component index. Color figure is plotted if weights provided.  
        
        backbone: 
            Used only when displaying genus computation where "vertices" form of 
            interaction links and the intervening polymer (backbone), see 
            computeGenus, ignore if visualizing multi-C data
        
        xRange: 
            Range of genomic region; defaults to, 
            xRange = [np.min(LinkDataAll[:,:2]), np.max(LinkDataAll[:,:2])]
            
        shiftRatio:
            Is the fraction of the xRange (the length of genomic region viewed) 
            that the baseline is shifted by,
              
            shift = shiftRatio*(xRange[1] - xRange[0]) 
        
        yMaxFrac:
            Determines the y-range as a fraction of the largest loop, useful to zoom
            in if there is one very large loop; yMax = max_loop_size*yMaxFrac,
            ax.set_ylim([0, yMax]) 
         
        ax: 
            Can pass axis of a figure to plot
            
        unit: 
            Is the genomic unit for xticks, set to 10^3 to get xticks in Kbs 
            for example 
        
        title: 
            Title of the plot
        
        cmap: 
            Color map for display of weights 
        
    *Returns:*
        figure axis 
        
    *Examples:* 
        >>> createRainbow(np.asarray([[1,7][5,10][12, 15]]))
    
    """
    
    
    colored = False #indicator whether need to draw colored diagram, i.e. LinkDataAll is 3 column    
    LinkDataAll = np.asarray(LinkDataAll) 
    assert LinkDataAll.ndim == 2, "Data should be Nx2 or Nx3 array"
    
    dict_colors = dict() 
    
    
    grey_color = '#301c0c'  #default color 
    if LinkDataAll.shape[1] == 2: 
        #if LinkDataAll is two columns (i.e., no color) then add thrid column with zeros  
        LinkDataAll =np.hstack((LinkDataAll, np.zeros((len(LinkDataAll),1))))
        dict_colors[0] = grey_color 
        
    elif LinkDataAll.shape[1]  ==3 : 
        colored = True
        my_map = get_cmap(cmap) 
        cNorm  = colors.Normalize(vmin=np.min(LinkDataAll[:,2]), vmax= np.max(LinkDataAll[:,2]))
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=my_map)
        sm = plt.cm.ScalarMappable(cmap=my_map, norm = cNorm)
    
    
        cluster_index = np.unique(LinkDataAll[:,2]) 
    
        for c in cluster_index:
            dict_colors[c] = scalarMap.to_rgba(c)   
    else: 
        raise TypeError('LinkData can only be N x2 or Nx3 matrix (list of list or numpy array), each row being [p1, p2, weight/value]')   
    
    
    
    if ax is None: 
        figsize = (20, 5) 
        fig = plt.figure(figsize=figsize) 
        ax = fig.add_subplot(1,1,1)

    
    codes = [Path.MOVETO,
            Path.CURVE4,
            Path.CURVE4,
            Path.CURVE4,
            ]
    
    
    if xRnage is None: 
        xRange = [np.min(LinkDataAll[:,:2]), np.max(LinkDataAll[:,:2])] 
        
    shift = shiftRatio*(xRange[1] - xRange[0]) 
    
    for p1, p2, w in LinkDataAll:                
        
        dist = np.abs(p1-p2)
        
        verts = [
                (p1, shift),  # P0
                (p1, dist+shift), # P1
                (p2, dist+shift), # P2
                (p2, shift), # P3
                ]



        path = Path(verts, codes)
        patch = patches.PathPatch(path, facecolor = 'none', edgecolor = dict_colors[w], lw=2)
        ax.add_patch(patch)


    max_loop_size =  np.max(np.abs(LinkDataAll[:,0]- LinkDataAll[:,1]))
    yMax = max_loop_size*yMaxFrac 
    ax.set_ylim([0, yMax]) 
    ax.set_xlim(xRange) 
    
    if backbone is None: 
        ax.axhline(shift, color = grey_color, linewidth = 1)  
    else: 
        for b1, b2, w in backbone: 
            #ax.axhline(y = shift, xmin = float(np.min((b1,b2))-xRange[0])/(xRange[1] - xRange[0]), xmax = float(np.max((b1, b2))- xRange[0])/(xRange[1] - xRange[0]), color = colors[w], linewidth = 2)  
            ax.plot([np.min((b1,b2)), np.max((b1, b2))], [shift, shift], color = dict_colors[w], linewidth = 2)  
           
    ax.get_yaxis().set_visible(False) 
    

    numLabels = np.round(np.min((10, (xRange[1] - xRange[0]))))  
    skipper   = int((xRange[1] - xRange[0])/numLabels)
    xlabels   = ['%i' % (float(x)/unit) for x in np.arange(xRange[0], xRange[1], skipper)]
    xlocation = np.arange(xRange[0], xRange[1], skipper)
    
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(xlocation)
    ax.set_xticklabels(xlabels, fontsize = 15)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.set_title(title, fontsize = 15) 
    
    
    if colored and ax is None:
        c_bar = fig.add_axes([0.90, 0.70, 0.01, 0.20])
        sm._A = [] 
        myCB = plt.colorbar(mappable = sm, cax = c_bar)
        myCB.ax.tick_params(labelsize = 15) 
        #myCB.ax.set_ylabel('weight', fontsize = 15)   
            
    plt.show() 
    return ax 
        


################################################################################

# Rules:
# matplotlibs ellipse patch has width height and center of ellipse, 
# I decided on the convention that width is the 2*semi_major, and height is
# 2*semi_minor, therefore at zero angle, this looks like a sitting egg  
# the angle of the ellipse is the rotaiton of the horizontal axis (semi-major)
# clockwise 

class loopDiagramCreator:
    """
    Does all the computation to create loop diagrams (see paper). 
    Loop diagrams are constructed by minimally reducing the list of links 
    (specific interactions) called "Data" here, such that the reduced set is
    genus zero (planar, non-crossing). The outermost loop for a "cluster" of
    loop-within-loop is the "footprint" loop. The "footprint" loops become
    ellipses (with realistic perimeter), and the planar loops are drawn as
    links in one color(blue) and the offending loops (genus-increasing) drawn in
    another(orange). The backbone (non-looped section) are displayed on a (semi)
    circle.  
     
    """
    
    
    
    def __init__(self, Data):
        
        """
        *Args:*
            Data: 
                Nx2 numpy array with N links between points [p1, p2]
                
        *Attributes:*
            kept_links: 
                Result of trimtoZero(), the list of links that maintain genus
                zero 
            trimmed_links: 
                The list of links that have to trimmed to make genus zero, 
                order of trimming is in trimToGenusZero()
        """
        self.Data = Data #this is the 2 column matrix with link points 
        _,_, self.kept_links, self.trimmed_links = trimToGenusZero(Data)    
        
        
    def findFootPrints(self): 
        """
        Determines the "footprint" of planar links; isolated loops along the 
        chromatin which are independent loop domains.
        
        *Attributes:*
            FPs: 
                numpy array of "footprints"
            xRange: 
                The range of the backbone excluding the looped out domains  
          
        """
        loc_Data = np.sort(self.kept_links, axis = 1) #sorted for each row 
        loc_Data, indx = uniquerows(loc_Data) #clean up operation, will use length of Data
        loc_Data = loc_Data[np.lexsort((loc_Data[:,1], loc_Data[:,0]))] #this sorts data by first column and then by second column
        
                    
        footprints = defaultdict(list)
        #self.clusters = defaultdict(list)
        counter = 0 
        #self.clusters[0].append([loc_Data[0,0], loc_Data[0,1]])
        footprints[0] = [loc_Data[0,0], loc_Data[0,1]] #inital footprint is first bond  
        
        
        for i, (p1, p2) in enumerate(loc_Data[1:]): 
        
            if p1 <= footprints[counter][1]:  #i.e., the region overlaps the previous region
            
                footprints[counter][1] = np.max((p2, footprints[counter][1]))  # thelexsort above ensures that p1 > footprint[0]
                #self.clusters[counter].append([p1, p2])
            else: 
                counter +=1
                footprints[counter] = [p1, p2] #new singleton cluster to grown, hence footprint is the new link 
                #self.clusters[counter].append([p1,p2])
        
        self.FPs = np.asarray([footprints[c] for c in footprints])
        self.xRange = [np.min(self.Data), np.ceil(np.max(self.Data) - np.sum(np.abs(self.FPs[:,1] - self.FPs[:,0])))] 
        
    #---------------------------------------------------------------------------
    @staticmethod 
    def precomputeEllipticArc(grid_size, ratio):
        """
        Precomputes the incomplete elliptic integrals on a grid of theta so that 
        one can find the adress on a ellipse (x, y coordinates) for a given arc 
        length, :math:`k = \\sqrt{1 - b^2/a^2}` where b is semi-minor and a is 
        semi-major axis of the ellipse 
        
      
        *Args:*  
            grid_size: 
                The number of points (on a=1 ellipse) where the 
                incomplete elliptic integrals is computed 
            ratio: 
                b/a (semi-minor/semi-major) of Ellipse 
                
        *Returns:* 
            all_phi: 
                the grid of angles (ellipse parameter) on which the arc 
                lengths are evaluated 
            all_path_lengths: 
                the evaluated arc lengths 
    
        """
        k = np.sqrt(1- ratio**2)
        phi = np.linspace(0, np.pi/2, grid_size+1)
        ellip = scsp.ellipeinc(phi, k)    
        # now the phi runs from 0 - pi/2 owing to symmetry of the ellipse
        # compute the arc length of the ellipse along the perimeter for phi from
        # 0 to 2pi, weave together 0-pi/2, and then pi/2-0 etc. 
        
        max_val = ellip[-1]
        ellip = np.delete(ellip, -1)
        all_path_length = np.zeros(4*(grid_size) + 1 )  # four quarters of 0 - (pi/2 - grid) 
        all_phi = np.zeros(4*(grid_size) + 1) # the last value is 2*pi
        
        all_path_length[0:grid_size]  = ellip #the last element is for pi/2
        all_path_length[(grid_size):(2*grid_size)] = (max_val - ellip[::-1])  + max_val #weaving requres reversing,
        #owing to the nature of arc length on an ellipse
        all_path_length[(2*grid_size):(3*grid_size)] = ellip  + 2*max_val
        all_path_length[(3*grid_size):(4*grid_size)] = (max_val - ellip[::-1])  + 3*max_val 
        
        all_path_length[-1] = 4*max_val 
        for i in range(4):   
            all_phi[i*(grid_size):(i+1)*(grid_size)] = phi[:-1] + i*phi[-1]                    
        
        all_phi[-1] = 2*np.pi
        
        return all_phi, all_path_length 
    #-----------------------------------------------------------------------

    @staticmethod
    def findxy(all_phi, all_path_length, arc_length, semi_major, ratio, theta = -0.5*np.pi):
        """
    
        Finds location (x,y) on an ellipse (center at (0,0) with major axis along 
        x-axis at zero rotation) which is rotated anti-clockwise by angle theta; 
        given the the grid of ellipse paramter phi, the evaluation of perimeter 
        along the ellipse for phi, the semi_major axis of the ellipse, and the 
        arc_length to match.
        
        *Args:* 
            all_phi: 
                see precomputeEllipticArc() method
            all_path_length: 
                see precomputeEllipticArc() method
            arc_length: 
                the arc length upto the point for which location is queried
            semi_major: 
                a of Ellipse 
            ratio: 
                b/a
            theta: 
                angle of anti-clockwise rotation of ellipse 
                
        *Returns:*
            np.array([x,y]): 
                location
        
        """
        
        #first find the (x,y) that would best fit, in conventional ellipse 
        #parametrization, center at (0,0) major axis along x-axis 
        scaled_path_length = semi_major*all_path_length
        phi = all_phi[bisect.bisect(scaled_path_length, arc_length)-1] # all_path_length 
    
        #is ordered array, bisect finds the best place to insert arc_length
        
        # x = a cos(phi)
        # y = b sin(phi) 
        
        #now, first reverse the angle 
        x = semi_major*np.cos(-phi)
        y = ratio*semi_major*np.sin(-phi)
        
        #The rotation matrix counter clockwise is 
        # R = [ cos(theta), -sin(theta)
        #       sin(theta), cos(theta)]
    
        R = np.asarray([[np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]])                                                                                                                                                                                                                                                                                                                                     
        x, y = np.dot(R, [x, y])                                                                                                                                                   
        
        return np.array([x,y])
        
    #-----------------------------------------------------------------------
            
    def _findPositionOfLinkEnds(self, links):
        """
        Private method to find address of link points (p1,p2) on the ellipse or 
        the backbone arc, as the case may be 
        
        *Args:*
            links: 
                numpy array Nx2 of N links (p1,p2)
                
        *Returns:*
            numpy array of Nx4 of positions, (x1,y1,x2,y2) where p1 is at (x1,y1) etc. 
        """
        positions  = np.zeros((len(links), 4)) # [x1, y1, x2, y2] is the convention for
        #adress of p1 and p2 
        for i, (p1, p2) in enumerate(links): 
            #find the x,y point on the polymer 
            #any point is either entirely in the loop
            #or crosses the loop domains 
            indA = bisect.bisect(self.FPintervals, p1)
            indB = bisect.bisect(self.FPintervals, p2)
            #if both indA and indB are odd then 
            #the bond is entirely within one of the
            #FPs
            closest_loop_point = self.FPintervals[indA-1]
            arclength = p1 - closest_loop_point
            
            
            
            if indA%2: 
            #left point is in interval, need to compute both x and y
            
                a1 = self.majorAxis[closest_loop_point]
                x1, y1 = self.findxy(self.all_phi, self.all_path_length, arclength, a1, self.ratio, theta = np.pi + self.loop_at_angles[closest_loop_point])
                #needed to add np.pi to the rotation because the the horizontantal left facing ellipse is at pi angle (going counterclockwise)
                #and the horizontal right facing is at zero. the alignment of the ellipse is pi off from the radial direction 
                
                Ecenter_x, Ecenter_y = self.Ecenters[closest_loop_point]
                x1 += Ecenter_x
                y1 += Ecenter_y
            
            else: 
            #need to compute position of point along backbone arc 
                arc_span = arclength/self.radius 
                place_at_angle = self.loop_at_angles[closest_loop_point] - arc_span # we are going clockise, therefore, -ve arc_span
                x1, y1 = (self.radius)*np.cos(place_at_angle), (self.radius)*np.sin(place_at_angle)
        
                
            closest_loop_point = self.FPintervals[indB-1] 
            arclength = p2 - closest_loop_point 
            
            if indB%2: 
                #left point is in interval, need to compute both x and y
                
                a2 = self.majorAxis[closest_loop_point]
                x2, y2 = self.findxy(self.all_phi, self.all_path_length, arclength, a2, self.ratio, theta = np.pi + self.loop_at_angles[closest_loop_point])
                Ecenter_x, Ecenter_y = self.Ecenters[closest_loop_point]
                x2 += Ecenter_x
                y2 += Ecenter_y
                
            else: 
                #need to compute position of point along backbone arc 
                arc_span = arclength/self.radius 
                place_at_angle = self.loop_at_angles[closest_loop_point] - arc_span # we are going clockise, therefore, -ve arc_span
                x2, y2 = (self.radius)*np.cos(place_at_angle), (self.radius)*np.sin(place_at_angle)
                
            positions[i] = [x1, y1, x2, y2]
                
        return positions   
    
    #-----------------------------------------------------------------------        
            
    def _computeEllipseAddresses(self):
        """
        private method to compute the addresses of the loop domain ellipses 
        
        *attributes:*
            radius: 
                Radius of backbone (semi-)circle
            Ecenters:
                Center locations of the ellipses 
            loop_at_angles:
                Angles at which the ellipses have to oriented (radially along
                backbone semi-circle) 
            majorAxis: 
                major axis of all the ellipses 
            max_loop_size: 
                Maximum perimeter of ellipses  
            FPintervals: 
                An array of "footprints" in increasing order, [f1, g1, f2, g2, .. ] 
                where [f1,g1] is the first footprint, etc. 
            all_phi:
                see preComputeEllipticArc()
            all_path_length:
                see preComputeEllipticArc()
            pos_kept_links: 
                computed positions of kept_links 
            pos_trimmed_links: 
                computed positions of trimmed_links 
                    
        """
        
        
        print "working on circular config..."
        #the ellipses are drawn clockwise 
        self.findFootPrints() #find the enclosing loops that "canopy" isolated clusters 
        self.radius = (self.xRange[1] - self.xRange[0])/(2*np.pi - 2*self.gap_angle)
        self.Ecenters  = dict()
        self.loop_at_angles = dict()
        self.majorAxis = dict() 
        
    
        running_sum = 0
        for i, (p1, p2) in enumerate(self.FPs): #running over the footprints                 
        
            xp1 = p1 - running_sum
            a = computeMajor(self.ratio, np.abs(p1-p2))
            running_sum += np.abs(p2 - p1)

            arc_angle = self.start_angle - (xp1-self.xRange[0])/self.radius  
            Ecenter_x, Ecenter_y = (self.radius+a)*np.cos(arc_angle), (self.radius+a)*np.sin(arc_angle)
        
            self.Ecenters[p1] = (Ecenter_x, Ecenter_y)
            self.Ecenters[p2] = (Ecenter_x, Ecenter_y) 
            self.loop_at_angles[p1] = arc_angle  
            self.loop_at_angles[p2] = arc_angle
            self.majorAxis[p1] = a 
            self.majorAxis[p2] = a
            
            
        self.max_loop_size = np.max(np.abs(self.FPs[:,0] - self.FPs[:,1]))         
        self.FPintervals = np.ravel(self.FPs[:,0:2], order = 'C') #these are pairs of beginning and
        #end of the FPs in linear order, [a1,b1, a2,b2, ...] where a and b a re FP start and end
                            
        # compute the elliptic arc length 
        self.all_phi, self.all_path_length = self.precomputeEllipticArc(self.grid_size, self.ratio)
        
        
        self.pos_kept_links = self._findPositionOfLinkEnds(self.kept_links)
        self.pos_trimmed_links = self._findPositionOfLinkEnds(self.trimmed_links)
            
            
    #-----------------------------------------------------------------------        
                
    def loopDiagramFig(self, figsize = (40,20), circle_center = [0,0], ratio = 0.33, grid_size = 1000, gap_angle = np.pi/2, unit = 10**6, annotation_angle_gap = 3): 
        """
        The plotting function for the loop diagrams, change the colors in the 
        source file if needed 
        
        
        *Args:*
            figsize: 
                size of figure 
            circle-center: 
                center of the (semi-) circle for the backbone 
            ratio:
                b/a for all eliipses to be display loop domains 
            grid_size:
                grid size for precomuting the elliptic arc lengths against ellipse
                parametric angle 
            gap_angle: 
                The angle of gap left from bottom of circle to left clockwise
                for the backbone semi-circle drawing 
            unit:
                Genomic unit, 10^6 is 1 Mb for example
            
            annotation_angle_gap:
                Minimum angle left between succesive annotation of genomic position
                to avoid over crowding   
            
        
        *Attributes:*
            ratio: 
            grid_size:
            gap_angle:
            circle_center:
            start_angle: 
                3/2 pi - gap_angle
            angle_span_deg:
                [start angle, end_angle] in degrees 
                
        """
        
        color_arc = '#726255'
        color_text = '#372e29'
        color_trim_links = '#eb6841'
        color_kept_links = '#00a0b0'
        
        self.ratio = ratio #thsi is b/a: semi-minor/semi-major axis 
        self.grid_size = grid_size #grid size for computing the incomplete elliptical integrals for perimeter along ellipse 
        self.gap_angle = gap_angle 
        self.circle_center = circle_center
        self.start_angle = 1.5*np.pi - self.gap_angle   #3/2 pi if gap is zero 
        self.angle_span_deg =np.degrees(np.asarray([self.start_angle + 2*self.gap_angle, self.start_angle]))  #this is anti-clockwise, 
        self._computeEllipseAddresses()
        
        

        fig = plt.figure(figsize=figsize) 
        ax = fig.add_subplot(1,1,1)
        ax.set_aspect('equal') 
        patch = patches.Arc(self.circle_center, 2*self.radius, 2*self.radius, angle = 0, theta1 = self.angle_span_deg[0], theta2 = self.angle_span_deg[1], linewidth = 3, edgecolor = color_arc)
        ax.add_patch(patch)
        for k,v in self.Ecenters.iteritems(): 
            a = self.majorAxis[k]
            arc_angle = self.loop_at_angles[k]
            patch = Ellipse((v[0], v[1]), 2*a, 2*self.ratio*a, angle = np.degrees(arc_angle), linewidth =2, fill = False, edgecolor = color_arc)  
            ax.add_patch(patch)
            
            
        angle_old = -np.inf    
        for p1, p2 in self.FPs: 
            ps1 = '%.2f' % (float(p1)/unit)
            ps2 = '%.2f' % (float(p2)/unit)
            arc_angle = self.loop_at_angles[p1]     
            base_x, base_y = (0.99*self.radius)*np.cos(arc_angle), (0.99*self.radius)*np.sin(arc_angle)
            base_xt, base_yt = (0.90*self.radius)*np.cos(arc_angle), (0.90*self.radius)*np.sin(arc_angle)
            text_rot_angle = np.degrees(arc_angle) 
            if arc_angle > np.pi/2: 
                text_rot_angle = np.degrees(-(np.pi - arc_angle))
                
            if np.abs(text_rot_angle - angle_old) > annotation_angle_gap: 
                ax.annotate(ps1 + '-' + ps2, xy=(base_x, base_y), xytext=(base_xt, base_yt), rotation = text_rot_angle,  horizontalalignment='center', size = 10, color = color_text, 
                verticalalignment='center', arrowprops=dict(facecolor=color_arc, arrowstyle="-"))
            angle_old = copy.copy(text_rot_angle)
 
        for (x1, y1, x2, y2) in self.pos_kept_links: 
            #print 'kept', x1,y1, x2, y2
            ax.plot([x1, x2], [y1, y2], color = color_kept_links, linewidth = 2, linestyle = '-')
            
        for (x1, y1, x2, y2) in self.pos_trimmed_links: 
            #print 'trimmed', x1,y1, x1, y2                                                                        
            ax.plot([x1, x2], [y1, y2], color = color_trim_links, linewidth = 2, linestyle = '-')                                                         
                                                                                    
        ax.get_yaxis().set_visible(False) 
        ax.get_xaxis().set_visible(False)
        

        max_ellipse_major = np.max(self.majorAxis.values())
        
        ax.set_xlim([-2*self.radius - 2*max_ellipse_major , 2*self.radius + 2*max_ellipse_major]) 
        ax.set_ylim([-2*self.radius - 2*max_ellipse_major , 2*self.radius + 2*max_ellipse_major])
        
        ax.set_aspect('equal') 
        plt.show()
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
################################################################################   

    
def visualize(genus, G, LoopData, backboneData, Data, **kwargs): 
    
    """
    Visualize computed genus using custom rainbow diagrams
    
    *Args:*
      genus:
          computed genus
      G:
          networkx graph returned by genus computation, see topology module 
      LoopData:
          Which of the edges of the graph are links correspoding to loops
      backboneData: 
          which of the edges of the graph contribute to the "backbone" along the
          genome
      Data: 
          Orginal data of links (Nx2 of each row [p1,p2]) 
      **kwargs:
          Args. passed on to createRainbow()  
    
    """
    CC = list(nx.connected_component_subgraphs(G)) 

    NodeColor= {} 
    count = 0 
    for g in CC: 
        for n in g.nodes(): 
            NodeColor[n] = count
        count +=1 
    
    
    LinkDataColor = [] 
    
    for p1, p2 in LoopData: 
        LinkDataColor.append([G.node[p1]['pos'], G.node[p2]['pos'], NodeColor[p1]]) 
        
    backboneColor = [] 
    for p1, p2 in backboneData: 
        backboneColor.append([G.node[p1]['pos'], G.node[p2]['pos'], NodeColor[p1]]) 
        
    figsize = (20, 15) 
    fig, ax = plt.subplots(2, 1, figsize=figsize) 
    createRainbow(Data, ax = ax[0], title = 'orginal', **kwargs)
    createRainbow(LinkDataColor, backbone = backboneColor, ax = ax[1], title = 'split and looped: genus: ' +str(genus), **kwargs )
  
    
################################################################################
        
def visualizeGenusLengthScales(genus_data, count_loops, xBinsAll, powers, units = 10**7, prob_max = 1, show_counts = True, display = 'upper'): 
    """
    Visualize the genus computation at various lengthscales and the summary 
    histograms of genus (and counts) of loops at those lengthscales; 
    plots a heatmap along the genome with local genus variaiton form median
    value at that length-scale. One row for each length-scale.
    
    Use computeGenusLengthScales() to compute the args. of this.   
    
    *Args:*
        genus_data: 
            The genus computed locally, dict. of lists with keys for length scales
        count_loops: 
            The number of links that featured in the above locla genus computation,
            dict of list as above 
        xBinsAll:
            The bin edge positions along the genome (dict of list) for the above
        powers: 
            Lx2 array of powers (of 10) for upper (column 0) and lower length scale
            (column 1) 
        units:
            Units in genome length 
        prob_max:
            The y limit of probability (y-axis) in histograms 
        show_counts: 
            Binary, whether to plot the distribution of count_loops
        display: 
            "upper" or "lower" for upper or lower length scale to display as titles
             
   
    """
    

    numrows = len(powers) 
    numlabels = 5 
    
    if show_counts: 
        gs = gridspec.GridSpec(numrows,3, width_ratios = [10,1,1], hspace = 0.2, wspace = 0.02)
    else: 
        gs = gridspec.GridSpec(numrows,2, width_ratios = [10,1], hspace = 0.2, wspace = 0.02)
        
    fig = plt.figure(figsize=(20,8)) 
    
    
    my_cmap = cmx.get_cmap('RdBu')
    bad_cl = [205./255., 201./255., 201./255.,1]
    cNorm  = colors.Normalize(vmin=-5, vmax=5)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=my_cmap)
    sm = plt.cm.ScalarMappable(cmap=my_cmap, norm = cNorm)
    
    
    heatmap_ax = [None]*numrows
    heatmap_ax[0] = plt.subplot(gs[0,0])
    bar_ax1 = [None]*numrows
    bar_ax1[0] = plt.subplot(gs[0,1])
    
    if show_counts: 
        bar_ax2 = [None]*numrows
        bar_ax2[0] = plt.subplot(gs[0,2])
    
    for i in range(1, numrows): 
        heatmap_ax[i] = plt.subplot(gs[i,0], sharex = heatmap_ax[0], sharey = heatmap_ax[0])
        bar_ax1[i] = plt.subplot(gs[i,1], sharey = bar_ax1[0])
        if show_counts: 
            bar_ax2[i] = plt.subplot(gs[i,2], sharey = bar_ax1[0]) #all bar plots share y axis 
    
    for i in range(numrows): 
        loc_data = genus_data[i]
        loc_data = np.asarray(loc_data).reshape([1, len(loc_data)])
        med = np.median(loc_data)
        xbinsL = xBinsAll['left'][i]
        xbinsR = xBinsAll['right'][i]
        for j,x in enumerate(xbinsL): 
            if np.isnan(loc_data[0,j]):
                color = bad_cl
            else:
                color  = scalarMap.to_rgba((loc_data[0,j] - med))
                
            heatmap_ax[i].axvspan(x, xbinsR[j], color = color)
            heatmap_ax[i].axvline(x, color = 'w', linewidth = 0.3)
        
        heatmap_ax[i].set_yticks([])
        heatmap_ax[i].set_yticklabels([])
        
        heatmap_ax[i].set_xticks(range(xbinsL[0], xbinsR[-1], units))
        
        if i== numrows-1: 
            heatmap_ax[i].set_xticklabels(np.round(np.arange(xbinsL[0], xbinsR[-1], units)/units).astype(int), fontsize = 15)
        else:    
            heatmap_ax[i].set_xticklabels(np.round(np.arange(xbinsL[0], xbinsR[-1], units)/units).astype(int), visible = False)
        
        if display == "upper": 
            heatmap_ax[i].set_ylabel("$10^{" + str(powers[i,0]) + "}$", fontsize = 15 )
        elif display == "lower":
            heatmap_ax[i].set_ylabel("$10^{" + str(powers[i,1]) + "}$", fontsize = 15 )
        
        ####### now do the bar plots ##############
        
        loc_data =  loc_data[ ~np.isnan(loc_data)]
        bins = np.arange(np.min(loc_data), np.max(loc_data)+1)  
        hist, bins = np.histogram(loc_data, bins) 
        hist = hist/float(np.sum(hist))
        
        
        bar_ax1[i].bar(bins[:-1], hist, width =1, align = 'center', color = 'k', alpha = 0.5, ec = 'w')
        
        if len(bins) <= numlabels:   
            bar_ax1[i].set_xticks(bins[:-1])
        else: 
            skipper = np.round((bins[-1] - bins[0])/float(numlabels))
            bar_ax1[i].set_xticks(bins[:-1:skipper]) 
            bar_ax1[i].set_xticklabels(bins[:-1:skipper].astype(int), fontsize = 12)
        
       
        bar_ax1[i].yaxis.tick_right()
        bar_ax1[i].set_yticks(np.arange(0, prob_max, 0.2))    
        bar_ax1[i].set_yticklabels(np.arange(0, prob_max, 0.2), fontsize = 12, visible = False) 
        bar_ax1[i].set_xlim(left = np.max((0,np.min(bins)))-0.5)
        bar_ax1[i].set_ylim([0,prob_max])
        
        if (not show_counts) and i == numrows-1: #last barplot has a y label 
            bar_ax1[i].set_yticklabels(np.arange(0, prob_max,0.2), fontsize = 12, visible = True) 
        
        if i == 0: 
                bar_ax1[i].set_title('Genus', fontsize = 12) 
            
        ###########################
        
        if show_counts: 
            loc_data = count_loops[i]
            loc_data = np.asarray(loc_data).reshape([1, len(loc_data)])
            loc_data =  loc_data[ ~np.isnan(loc_data)]
            
            bins = np.arange(np.min(loc_data), np.max(loc_data) +1)   
            hist, bins = np.histogram(loc_data, bins) 
            
            hist = hist/float(np.sum(hist))
            bar_ax2[i].bar(bins[:-1], hist, width =1, align = 'center', color = 'k', alpha = 0.5, ec = 'w')
            
            if len(bins) <= numlabels:   
                bar_ax2[i].set_xticks(bins[:-1])
            else: 
                skipper = np.round((bins[-1] - bins[0])/float(numlabels))
                bar_ax2[i].set_xticks(bins[:-1:skipper]) 
                bar_ax2[i].set_xticklabels(bins[:-1:skipper].astype(int), fontsize = 12)
            
            
            bar_ax2[i].yaxis.tick_right()
            bar_ax2[i].set_yticks(np.arange(0, prob_max,0.2))    
            bar_ax2[i].set_yticklabels(np.arange(0, prob_max,0.2), fontsize = 12, visible = False) 
            bar_ax2[i].set_xlim(left = np.max((0,np.min(bins)))-0.5)
            bar_ax2[i].set_ylim([0,prob_max])
            
            if i == numrows -1: 
                bar_ax2[i].set_yticklabels(np.arange(0, prob_max,0.2), fontsize = 12, visible = True) 
        
            if i == 0: 
                bar_ax2[i].set_title('Counts', fontsize = 12) 
        
    c_bar = fig.add_axes([0.91, 0.70, 0.01, 0.20])
    sm._A = [] 
    myCB = plt.colorbar(mappable = sm, cax = c_bar)
    myCB.ax.tick_params(labelsize = 12) 
    
 


################################################################################

def plotGenusStats(genus_data, count_data, powers, prob_max = 1, random_genus_accumulator = None, display = 'upper'): 
    """
    Plots a grid of histograms of genus, genus of randomized links and counts of 
    observed links at length scales along rows, columns are different length-scales 
    
    *Args:*
        genus_data: 
            The genus computed locally, dict. of lists with keys for length scales
        count_loops: 
            The number of links that featured in the above locla genus computation,
            dict of list as above 
        powers: 
            Lx2 array of powers (of 10) for upper (column 0) and lower length scale
            (column 1) 
        prob_max:
            The y limit of probability (y-axis) in histograms 
        randon_genus_accumulator: 
            If not none, then plots the random data provided here, corresponding 
            to the randomization of the genus data 
        display: 
            "upper" or "lower" for upper or lower length scale to display as titles
             
   
    """
    color1 = [250/255.0,199/255.0,120/255.0]
    color2 = [153/255.0,203/255.0,186/255.0]
    color3 = [228/255.0,108/255.0,89/255.0]
    import matplotlib.gridspec as gridspec
    switch_to_plot = 30 #num bins at which the plot switches to stem plot instead of bar plot
    numcols = len(powers) 
    
    numlabels = 5 
    
    if random_genus_accumulator is None: 
        gs = gridspec.GridSpec(2, numcols, hspace = 0.15, wspace = 0.04)
    else: 
        gs = gridspec.GridSpec(3, numcols, hspace = 0.15, wspace = 0.04)    
        
    fig = plt.figure(figsize=(20,9)) 
    bar_ax1 = [None]*numcols
    bar_ax1[0] = plt.subplot(gs[0,0])
    
    bar_ax1[0].set_yticks(np.arange(0, prob_max, 0.2)) 
    bar_ax1[0].set_yticklabels(np.arange(0, prob_max,0.2), fontsize = 12, visible = True) 
    
    bar_ax2 = [None]*numcols
    
    for i in range(0, numcols): 
    
        bar_ax1[i] = plt.subplot(gs[0,i], sharey = bar_ax1[0])
        bar_ax2[i] = plt.subplot(gs[-1,i], sharey = bar_ax1[0]) #all bar plots share y axis 
    
    if random_genus_accumulator is not None:
        bar_ax3 = [None]*numcols 
        for i in range(0, numcols): 
            bar_ax3[i] = plt.subplot(gs[1,i], sharey = bar_ax1[0]) 
    
    for i in range(len(powers)): 
        loc_data = genus_data[i]
        bins = np.arange(np.min(loc_data), np.max(loc_data))   
        hist, bins = np.histogram(loc_data, bins) 
    
        hist = hist/float(np.sum(hist))
        if len(bins) > switch_to_plot:
            good_ind = hist>0
            markerline, stemlines, baseline = bar_ax1[i].stem(bins[:-1][good_ind], hist[good_ind], linestyle = '--', markerfmt  = " ")
            #plt.setp(markerline, 'markerfacecolor',  color1, 'markersize', 3, 'markeredgecolor', color1) 
            plt.setp(stemlines, 'color', color1, 'linewidth', 3)
            plt.setp(baseline, 'linewidth', 0.0) 
        else: 
            bar_ax1[i].bar(bins[:-1], hist, width =1, align = 'center', color = color1, ec = 'w')
        
        if len(bins) <= numlabels:   
            bar_ax1[i].set_xticks(bins[:-1])
        else: 
            skipper = np.round((bins[-1] - bins[0])/float(numlabels))
            bar_ax1[i].set_xticks(bins[:-1:skipper]) 
            bar_ax1[i].set_xticklabels(bins[:-1:skipper].astype(int), fontsize = 12)
        
        
        #bar_ax1[i].yaxis.tick_right()
         
        if i == 0: 
            bar_ax1[i].set_ylabel('Observed', fontsize = 12) 
        else:  
            plt.setp(bar_ax1[i].get_yticklabels(), visible=False) 
            
        bar_ax1[i].set_xlim([np.max((0,np.min(bins)))-1, np.max(bins)+1])
        if display == 'upper': 
            bar_ax1[i].set_title("$10^{" + str(powers[i,0]) + "}$", fontsize = 16 )
        elif display == 'lower': 
            bar_ax1[i].set_title("$10^{" + str(powers[i,1]) + "}$", fontsize = 16 )
        else: 
            bar_ax1[i].set_title("$10^{" + str(powers[i,0]) + "}$" + ','  + "$10^{" + str(powers[i,1]) + "}$", fontsize = 16 )
        
        ###########################
        
        
        loc_data = count_data[i]
        loc_data = np.asarray(loc_data).reshape([1, len(loc_data)])
        loc_data =  loc_data[ ~np.isnan(loc_data)]
        
        bins = np.arange(np.min(loc_data), np.max(loc_data))   
        hist, bins = np.histogram(loc_data, bins) 
        
        hist = hist/float(np.sum(hist))
        
        if len(bins) > switch_to_plot:
            good_ind = hist >0
            markerline, stemlines, baseline = bar_ax2[i].stem(bins[:-1][good_ind], hist[good_ind], linestyle = '--', markerfmt = " ")
            #plt.setp(markerline, 'markerfacecolor',  color2, 'markersize', 3, 'markeredgecolor', color2) 
            plt.setp(stemlines, 'color', color2, 'linewidth', 3)
            plt.setp(baseline, 'linewidth', 0.0) 
        else: 
            bar_ax2[i].bar(bins[:-1], hist, width =1, align = 'center', color = color2, ec = 'w', hatch = '/')
        
        if len(bins) <= numlabels:   
            bar_ax2[i].set_xticks(bins[:-1])
        else: 
            skipper = np.round((bins[-1] - bins[0])/float(numlabels))
            bar_ax2[i].set_xticks(bins[:-1:skipper]) 
            bar_ax2[i].set_xticklabels(bins[:-1:skipper].astype(int), fontsize = 12)
        
        
        #bar_ax2[i].yaxis.tick_right()
        if i == 0: 
            bar_ax2[i].set_ylabel('Count dist.', fontsize = 12)
        else:   
            plt.setp(bar_ax2[i].get_yticklabels(), visible=False)
            
        bar_ax2[i].set_xlim([np.max((0,np.min(bins)))-1, np.max(bins)+1])
    
        
        if random_genus_accumulator is not None: 
            loc_data = random_genus_accumulator[i]
            loc_hist = np.asarray([[x,y] for x,y in loc_data.items()])
            bins = loc_hist[:,0]
            hist = loc_hist[:,1]
            hist = hist/float(np.sum(hist)) 
            if len(bins) > switch_to_plot:
                good_ind = hist >0
                markerline, stemlines, baseline = bar_ax3[i].stem(bins[good_ind], hist[good_ind], linestyle = '--', markerfmt = " ")
                #plt.setp(markerline, 'markerfacecolor',  color2, 'markersize', 3, 'markeredgecolor', color2) 
                plt.setp(stemlines, 'color', color3, 'linewidth', 3)
                plt.setp(baseline, 'linewidth', 0.0) 
            else: 
                bar_ax3[i].bar(bins, hist, width =1, align = 'center', color = color3, ec = 'w', hatch = '\\')
            
            if len(bins) <= numlabels:   
                bar_ax3[i].set_xticks(bins)
            else: 
                skipper = np.round((bins[-1] - bins[0])/float(numlabels))
                bar_ax3[i].set_xticks(bins[::skipper]) 
                bar_ax3[i].set_xticklabels(bins[::skipper].astype(int), fontsize = 12)
            
            
            #bar_ax2[i].yaxis.tick_right()
            if i == 0: 
                bar_ax3[i].set_ylabel('Randomized', fontsize = 12)
            else:   
                plt.setp(bar_ax3[i].get_yticklabels(), visible=False)
                 
            bar_ax3[i].set_xlim([np.max((0,np.min(bins)))-1, np.max(bins)+1])
            
            
            
    return fig 


    
################################################################################    
    
def plotGenusLengthScale(genus_accumulator, random_genus_accumulator, powers, display = 'upper'): 
    """
    Summary plot of length scale vs genus with error bars for observed and randomized
    results
    
    *Args:*
        genus_accumulator: 
            dict of list of genus computation in bins at various length-scales (keys
            of dict) 
        random_genus_accumulator: 
            Randomized links genus
        powers:
            Lx2 powers of 10 for L (i.e. log10(length-scales)) with rows for
            upper and lower length-scale cutoffs. 
            
    
    """
    obs_data = []   
    random_data = [] 
    max_genus = 50 
    fig, ax = plt.subplots(1, 1, figsize=(15,10))
    
    color1 = [1,0,0]
    color2 = [0,0,1]
    
    if display == 'upper':
        powers_loc = powers[:,0]
    elif display == 'lower': 
        powers_loc = powers[:,1]
    for i in range(len(powers)): 
        
        loc_data = genus_accumulator[i]
        bins = np.arange(np.min(loc_data), np.max(loc_data))   
        hist, bins = np.histogram(loc_data, bins)
        hist = hist/float(np.sum(hist))
        obs_data.append([powers_loc[i], np.median(loc_data)])
        
        loc_random = random_genus_accumulator[i]
        loc_hist = np.asarray([[x,y] for x,y in loc_random.items()])
        rand_bins = loc_hist[:,0]
        rand_hist = loc_hist[:,1]
        rand_hist = rand_hist/float(np.sum(rand_hist))
        med_val = rand_bins[np.argmax(rand_hist)]
        random_data.append([powers_loc[i], med_val])
        colora = [color1 + [c] for c in hist]
        ax.scatter([powers_loc[i]+0.025]*len(bins), bins, marker = 's', s = 100, color = colora, edgecolor = colora)
        colorb = [color2 + [c] for c in rand_hist]
        ax.scatter([powers_loc[i]-0.025]*len(rand_bins), rand_bins, marker = 's', s = 100, color = colorb, edgecolor = colorb)
        
    obs_data = np.asarray(obs_data)
    random_data = np.asarray(random_data)    
    ax.plot(obs_data[:,0], obs_data[:,1], color = color1, marker = '', linewidth = 4, label = 'Obs.')
    ax.plot(random_data[:,0], random_data[:,1], color = color2, marker = '', linewidth = 4, linestyle = '--', label = "Rand.")
    ax.set_xticks(powers_loc)
    xticks = ['$10^{' + str(p) + '}$' for p in powers_loc]
    ax.set_xticklabels(xticks, fontsize = 18)      
    ax.set_ylim(bottom = -1, top = max_genus)
    plt.tick_params(axis='both', labelsize=18)
    ax.set_xlabel("Length Scale", fontsize = 18)
    ax.set_ylabel("Genus", fontsize = 18)
    plt.legend()



################################################################################