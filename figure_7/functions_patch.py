import numpy as np
from numpy import linalg as LA

def sum_squared_displacements(path, maxstep):
    #get subpaths such that succesive steps have step lengt less than maxstep
    #skip timepoints where there are not patch. 
    subpaths=[]
   
    #get first position that contains data
    for i in range(0,len(path)):
        if ~np.isnan(path[i,0]):
            break
    i_ini=i
    
    
    for i in range(i_ini+1,len(path)):
        #if current and previous point have data and the step lenght is long, a new subpath, starts: store the old one, and start a new one. 
        if ~np.isnan(path[i,0]) and ~np.isnan(path[i-1,0]) and LA.norm(path[i,:]-path[i-1,:]) > maxstep :
            subpaths.append([i_ini,i-1])
            i_ini=i
            
        #if current point is nan and the previous point has data close subpath 
        if np.isnan(path[i,0]) and ~np.isnan(path[i-1,0]):
            subpaths.append([i_ini,i-1])
            
        #if the current point has data but the previous point was nan, set i_ini      
        if ~np.isnan(path[i,0]) and np.isnan(path[i-1,0]):
            i_ini=i
    
    #add last subpath 
    if ~np.isnan(path[len(path)-1,0]):       
        subpaths.append([i_ini,len(path)-1])
        
    sum_norm_disp = np.zeros(len(path)-1)
    sum_squared_disp = np.zeros(len(path)-1)
    sum_squared_disp_squared = np.zeros(len(path)-1)    
    #number of each intervals computed
    n_squared_disp = np.zeros(len(path)-1)
    
    for subp in range(len(subpaths)):
        #get msds in subpath
        sublen = subpaths[subp][1] - subpaths[subp][0]
        for interval in range(1,sublen+1):
            #compute along the path
            #for pos in range( first, frist + sublen - interval + 1)
            for pos in range( subpaths[subp][0], subpaths[subp][0] + sublen - interval + 1):
                norm_disp = LA.norm( path[pos,:] - path[pos+interval,:])
                sum_norm_disp[interval-1] += norm_disp
                sum_squared_disp[interval-1] += norm_disp**2
                sum_squared_disp_squared[interval-1] += norm_disp**4
                n_squared_disp[interval-1] += 1
    return sum_norm_disp, sum_squared_disp, sum_squared_disp_squared, n_squared_disp

def geodesic_dist(p1,p2,R):
    #compute large circle distance (geodetic distance )
    dist = LA.norm(p1-p2)
    phi=np.arcsin(dist/2/R)
    gcdist=2*phi*R
    return gcdist

def sum_squared_displacements_sphere(path, maxstep, R):
    #get subpaths such that succesive steps have step lengt less than maxstep
    #skip timepoints where there are not patch. 
    subpaths=[] 
   
    #get first position that contains data
    for i in range(0,len(path)):
        if ~np.isnan(path[i,0]):
            break
    i_ini=i
    
    
    for i in range(i_ini+1,len(path)):
        #if current and previous point have data and the step lenght is long, a new subpath, starts: store the old one, and start a new one. 
        if ~np.isnan(path[i,0]) and ~np.isnan(path[i-1,0]) and geodesic_dist(path[i,:],path[i-1,:],R) > maxstep :
            subpaths.append([i_ini,i-1])
            i_ini=i
            
        #if current point is nan and the previous point has data close subpath 
        if np.isnan(path[i,0]) and ~np.isnan(path[i-1,0]):
            subpaths.append([i_ini,i-1])
            
        #if the current point has data but the previous point was nan, set i_ini      
        if ~np.isnan(path[i,0]) and np.isnan(path[i-1,0]):
            i_ini=i
    
    #add last subpath 
    if ~np.isnan(path[len(path)-1,0]):       
        subpaths.append([i_ini,len(path)-1])
        
    sum_norm_disp = np.zeros(len(path)-1)
    sum_squared_disp = np.zeros(len(path)-1)
    sum_squared_disp_squared = np.zeros(len(path)-1)    
    #number of each intervals computed
    n_squared_disp = np.zeros(len(path)-1)
    
    for subp in range(len(subpaths)):
        #get msds in subpath
        sublen = subpaths[subp][1] - subpaths[subp][0]
        for interval in range(1,sublen+1):
            #compute along the path
            #for pos in range( first, frist + sublen - interval + 1)
            for pos in range( subpaths[subp][0], subpaths[subp][0] + sublen - interval + 1):
                norm_disp = geodesic_dist( path[pos,:],path[pos+interval,:],R)
                sum_norm_disp[interval-1] += norm_disp
                sum_squared_disp[interval-1] += norm_disp**2
                sum_squared_disp_squared[interval-1] += norm_disp**4
                n_squared_disp[interval-1] += 1
    return sum_norm_disp, sum_squared_disp, sum_squared_disp_squared, n_squared_disp
        


def get_centroid(mat):
    N = np.shape(mat)[0]
    #coordinated of the distribution to make the calculation
    [xr,yr]=np.where(mat>0.5*np.max(mat))
    inds = np.ravel_multi_index(np.vstack((xr,yr)), np.shape(mat))
    #COMPUTATION OF THE X COMPONENT OF THE CENTROID
    if (np.max(xr) - np.min(xr)) > N/2: #space between found patches is larger than half of the box (assuming the diameter of the patch is smaller than N/2)->the patch is going over the boundary
        [binsleft]=np.where(xr < N/2) #bins found on the left side
        [binsright]=np.where(xr >= N/2) #bins found on the right side
        ind_bins_left = np.ravel_multi_index(np.vstack((xr[binsleft],yr[binsleft])), np.shape(mat))
        ind_bins_right = np.ravel_multi_index(np.vstack((xr[binsright],yr[binsright])), np.shape(mat))
        
        if len(binsright) > len(binsleft): #%majority of bins are on the right, easier to shift to the right            
            xsum_left = np.sum(mat.flat[ind_bins_left]*(xr[binsleft]+ N) )        
            xsum_right = np.sum(mat.flat[ind_bins_right]*(xr[binsright]))
            xmean= (xsum_right + xsum_left)/np.sum(mat[xr,yr])
            if xmean > N - 1:
                xmean=xmean-N
        else:
            xsum_left = np.sum(mat.flat[ind_bins_left]*(xr[binsleft]) )        
            xsum_right = np.sum(mat.flat[ind_bins_right]*(xr[binsright] - N))
            xmean= (xsum_right + xsum_left)/np.sum(mat[xr,yr])
            if xmean < 0 :
                xmean=xmean + N        
    else:
        xmean= np.sum(mat.flat[inds]*(xr))/np.sum(mat[xr,yr])            

    #SAME COMPUTATION FOR THE Y COMPONENT OF THE CENTROID, READ DOWN AND UP INSTEAD OF LEFT AND RIGHT    
    if (np.max(yr) - np.min(yr)) > N/2: #space between found patches is larger than half of the box (assuming the diameter of the patch is smaller than N/2)->the patch is going over the boundary
        [binsleft]=np.where(yr < N/2) #bins found on the left side
        [binsright]=np.where(yr >= N/2) #bins found on the right side
        ind_bins_left = np.ravel_multi_index(np.vstack((xr[binsleft],yr[binsleft])), np.shape(mat))
        ind_bins_right = np.ravel_multi_index(np.vstack((xr[binsright],yr[binsright])), np.shape(mat))
        
        if len(binsright) > len(binsleft): #%majority of bins are on the right, easier to shift to the right            
            ysum_left = np.sum(mat.flat[ind_bins_left]*(yr[binsleft]+ N) )        
            ysum_right = np.sum(mat.flat[ind_bins_right]*(yr[binsright]))
            ymean= (ysum_right + ysum_left)/np.sum(mat[xr,yr])
            if ymean > N - 1:
                ymean=ymean-N
        else:
            ysum_left = np.sum(mat.flat[ind_bins_left]*(yr[binsleft]) )        
            ysum_right = np.sum(mat.flat[ind_bins_right]*(yr[binsright] - N))
            ymean= (ysum_right + ysum_left)/np.sum(mat[xr,yr])
            if ymean < 0 :
                ymean=ymean + N        
    else:
        ymean= np.sum(mat.flat[inds]*(yr))/np.sum(mat[xr,yr])            
    return xmean,ymean