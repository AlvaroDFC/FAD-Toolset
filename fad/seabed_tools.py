"""A set of functions for processing seabed information for a Project."""


import os
import matplotlib.pyplot as plt
import numpy as np
import yaml
from .helpers import getFromDict

def loadBoundary(filename, lat0, lon0):
        '''
        Load a lease area boundary for the project from an input file.
        
        Parameters
        ----------
        filename : path
            path/name of file containing bathymetry data (format TBD)
        '''
        
        # load data from file
        Xs, Ys = processBoundary(filename, lat0, lon0)
        
        boundary = setBoundary(Xs, Ys)
    
        return boundary


def setBoundary(Xs, Ys):
    '''Set the boundaries of the project based on x-y polygon vertices.'''
    
    # check compatibility with project grid size
    
    # save as project boundaries
    boundary = np.vstack([[Xs[i],Ys[i]] for i in range(len(Xs))])
    # self.boundary = np.vstack([Xs, Ys])
    
    # if the boundary doesn't repeat the first vertex at the end, add it
    if not all(boundary[0,:] == boundary[-1,:]):
        boundary = np.vstack([boundary, boundary[0,:]])
    
    return boundary
    
    # figure out masking to exclude grid data outside the project boundary

def setGrid(xs, ys, grid_x, grid_y, grid_depth):
    '''
    Set up the rectangular grid over which site or seabed
    data will be saved and worked with. Directions x and y are 
    generally assumed to be aligned with the East and North 
    directions, respectively, at the array reference point.
    
    Parameters
    ----------        
    xs : float array
        x coordinates relative to array reference point [m]
    ys : float array
        y coordinates relative to array reference point [m]
    '''
    
    # Create a new depth matrix with interpolated values from the original
    depths = np.zeros([len(ys), len(xs)])  # note: indices are iy, ix
    for i in range(len(ys)):
        for j in range(len(xs)):
            depths[i,j], nvec = getDepthFromBathymetry(xs[j], ys[i], 
                                grid_x, grid_y, grid_depth)
    
    # Replace the grid data with the updated values
    grid_x = np.array(xs)
    grid_y = np.array(ys)
    grid_depth = depths

    return grid_x, grid_y, grid_depth

def loadBathymetry(filename, interpolate=False):
        '''
        Load bathymetry information from an input file (format TBD), convert to
        a rectangular grid, and save the grid to the floating array object (TBD).
        
        Paramaters
        ----------
        filename : path
            path/name of file containing bathymetry data (format TBD)
        '''
        
        # load data from file
        Xs, Ys, Zs = readBathymetryFile(filename)  # read MoorDyn-style file
        # Xs, Ys, Zs = sbt.processASC(filename, self.lat0, self.lon0)
        
        # ----- map to existing grid -----
        # if no grid, just use the bathymetry grid
        if not interpolate: #len(self.grid_x) == 0: 
            grid_x = np.array(Xs)
            grid_y = np.array(Ys)
            grid_depth = np.array(Zs)
            
        else:
        # interpolate onto grid defined by grid_x, grid_y
            for i, x in enumerate(grid_x):
                for j, y in enumerate(grid_y):
                    grid_depth[i,j], _ = getDepthFromBathymetry(x, y, Xs, Ys, Zs)

        return grid_x, grid_y, grid_depth



def readBathymetryFile(filename, dtype=float):

    with open(filename, 'r') as f:
        # skip the header
        line = next(f)
        # collect the number of grid values in the x and y directions from the second and third lines
        line = next(f)
        nGridX = int(line.split()[1])
        line = next(f)
        nGridY = int(line.split()[1])
        # allocate the Xs, Ys, and main bathymetry grid arrays
        bathGrid_Xs = np.zeros(nGridX)
        bathGrid_Ys = np.zeros(nGridY)
        bathGrid = np.zeros([nGridY, nGridX], dtype=dtype)  # MH swapped order June 30
        # read in the fourth line to the Xs array
        line = next(f)
        bathGrid_Xs = [float(line.split()[i]) for i in range(nGridX)]
        strlist = []
        # read in the remaining lines in the file into the Ys array (first entry) and the main bathymetry grid
        for i in range(nGridY):
            line = next(f)
            entries = line.split()
            bathGrid_Ys[i] = entries[0]
            if dtype==float:
                bathGrid[i,:] = entries[1:]
            if dtype==str:
                strlist.append(entries[1:])
        if dtype==str:
            bathGrid = np.array(strlist)
    
    return bathGrid_Xs, bathGrid_Ys, bathGrid

            
def getSoilTypes(filename, soil_mode='layered', profile_source=None):
    '''
    Load soil properties or layered profiles depending on soil_mode.

    Parameters
    ----------
    filename : str
        Path to .txt file containing grid and profile/soil label definitions
    soil_mode : str
        'uniform' or 'layered'
    profile_source : str or None
        Path to YAML file with layered soil profiles (used only for 'layered')

    Returns
    -------
    soilProps : dict
        Dictionary of soil type properties (uniform) or layered profiles (layered)
    '''
    soilProps = {}
    used_labels = []

    # NOTE: why passing file name instead of the already read soil grid?
    with open(filename, 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if line.strip().startswith('---') and 'SOIL TYPES' in line.upper():
            break

    # Extract used labels from the SOIL TYPES section
    for line in lines[i+3:]:
        if '---' in line:
            break
        entries = line.strip().split()
        label = entries[0]
        used_labels.append(label); print(label)

    if soil_mode == 'uniform':
        var_names = lines[i+1].split()
        for line in lines[i+3:]:
            if '---' in line:
                break
            entries = line.strip().split()
            label = entries[0]
            soilProps[label] = {}
            for iv, var in enumerate(var_names[1:]):
                val = entries[iv+1]
                soilProps[label][var] = [float(val)] if val != '-' else [0.0]

    elif soil_mode == 'layered':
        if profile_source is None:
            raise ValueError("profile_source (path to YAML) is required for layered mode.")

        # Load the full YAML file of profiles
        with open(profile_source, 'r') as f:
            all_profiles = yaml.safe_load(f)

        # Reassign each label to the actual layer list directly
        for label in used_labels:
            if label not in all_profiles:
                raise KeyError(f'Profile ID {label} not found in YAML: {profile_source}')
            soilProps[label] = all_profiles[label]['layers']  # now a list of layer dicts

        print(f"[DEBUG] Loaded profiles from YAML: {list(soilProps.keys())}")
        if used_labels:
            print(f"[DEBUG] Example layers for {used_labels[0]}: {soilProps[used_labels[0]]}")
        else:
            print("[WARNING] No profile labels were found in the soil grid.")


    else:
        raise ValueError(f"Unrecognized soil_mode '{soil_mode}'")

    return soilProps


def processBoundary(filename, lat, lon,meters=True):
    '''Reads boundary information from a CSV file and stores the boundary 
    coordinate list in a set of arrays. This function can be extended to
    deal with multiple boundary sets.
        
    Parameters
    ----------
    filename : string
        Filename containing columns of x and y coordinates of boundary.
    lat : float
        lattitude of reference point to use for array y grid
    long : float
        lattitude of reference point to use for array x grid

    Returns
    -------
    Xs : array
        x values of grid points [m]
    Ys : array
        y values of grid points [m]
    '''
    
    import pandas as pd
    
    zerozero = (lat, lon)  # lattitude and longitude of reference point (grid origin)
    
    delin = pd.read_csv(filename)
    longs = np.array(delin['X_UTM10'])
    lats = np.array(delin['Y_UTM10'])
    
    if meters:
        Xs = longs
        Ys = lats
    #else:
        #Xs, Ys = convertLatLong2Meters(zerozero, lats, longs)
    
    return Xs, Ys









def resampleGrid(x_new, y_new, x_old, y_old, grid_values):
    '''Interpolate an existing array of values on a rectangular grid to a new
    rectangular grid.
    
    Parameters
    ----------
    x_new : list
        x values of the new grid to interpolate to
    y_new : list
        y values of the new grid to interpolate to
    x_old : list
        x values of the original grid
    y_old : list
        x values of the original grid
    grid_values : 2D array
        The values on the old grid to be interpolated from (dimensions must
        match the length of y_old and x_old, in that order).
    
    Returns
    -------
    grid_values_new : 2D array
        Interpolated grid values on y_new and x_new grid lines.
    '''
    
    grid_values_new = np.zeros([len(y_new), len(x_new)])
    
    for i in range(len(y_new)):
        for j in range(len(x_new)):
            grid_values_new[i,j], _ = getDepthFromBathymetry(x_new[j], y_new[i],
                                                   x_old, y_old, grid_values)
    
    return grid_values_new


def getInterpNums(xlist, xin, istart=0):  # should turn into function in helpers
    '''
    Paramaters
    ----------
    xlist : array
        list of x values
    xin : float
        x value to be interpolated
    istart : int
        first lower index to try
    
    Returns
    -------
    i : int
        lower index to interpolate from
    fout : float
        fraction to return   such that y* = y[i] + fout*(y[i+1]-y[i])
    '''
    
    if np.isnan(xin):
        raise Exception('xin value is NaN.')
    
    nx = len(xlist)
  
    if xin <= xlist[0]:  #  below lowest data point
        i = 0
        fout = 0.0
  
    elif xlist[-1] <= xin:  # above highest data point
        i = nx-1
        fout = 0.0
  
    else:  # within the data range
 
        # if istart is below the actual value, start with it instead of 
        # starting at 0 to save time, but make sure it doesn't overstep the array
        if xlist[min(istart,nx)] < xin:
            i1 = istart
        else:
            i1 = 0

        for i in range(i1, nx-1):
            if xlist[i+1] > xin:
                fout = (xin - xlist[i] )/( xlist[i+1] - xlist[i] )
                break
    
    return i, fout


def interpFromGrid(x, y, grid_x, grid_y, values):
    '''Interpolate from a rectangular grid of values.'''

    # get interpolation indices and fractions for the relevant grid panel
    ix0, fx = getInterpNums(grid_x, x)
    iy0, fy = getInterpNums(grid_y, y)

    # handle end case conditions
    if fx == 0:
        ix1 = ix0
    else:
        ix1 = min(ix0+1, values.shape[1])  # don't overstep bounds
    
    if fy == 0:
        iy1 = iy0
    else:
        iy1 = min(iy0+1, values.shape[0])  # don't overstep bounds
    
    # get corner points of the panel
    c00 = values[iy0, ix0]
    c01 = values[iy1, ix0]
    c10 = values[iy0, ix1]
    c11 = values[iy1, ix1]

    # get interpolated points and local value
    cx0    = c00 *(1.0-fx) + c10 *fx
    cx1    = c01 *(1.0-fx) + c11 *fx
    c0y    = c00 *(1.0-fy) + c01 *fy
    c1y    = c10 *(1.0-fy) + c11 *fy
    value  = cx0 *(1.0-fy) + cx1 *fy

    # get local slope
    dx = grid_x[ix1] - grid_x[ix0]
    dy = grid_y[iy1] - grid_y[iy0]
    
    # deal with being on an edge or a zero-width grid increment
    if dx > 0.0:
        dc_dx = (c1y-c0y)/dx
    else:
        dc_dx = c0y*0  # maybe this should raise an error
    
    if dy > 0.0:
        dc_dy = (cx1-cx0)/dy
    else:
        dc_dy = cx0*0  # maybe this should raise an error
    
    # return the interpolated value, the derivatives, and the grid indices
    return value, dc_dx, dc_dy, ix0, iy0



def getDepthFromBathymetry(x, y, grid_x, grid_y, grid_depth, index=False):
    ''' interpolates local seabed depth and normal vector
    
    Parameters
    ----------
    x, y : float
        x and y coordinates to find depth and slope at [m]
    
    Returns
    -------        
    depth : float
        local seabed depth (positive down) [m]
    nvec : array of size 3
        local seabed surface normal vector (positive out) 
    index : bool, optional
        If True, will also retun ix and iy - the indices of the intersected
        grid panel.
    '''
    
    # Call general function for 2d interpolation
    depth, dc_dx, dc_dy, ix0, iy0 = interpFromGrid(x, y, grid_x, grid_y, grid_depth)
    
    # Compute unit vector of the seabed panel
    nvec = np.array([dc_dx, dc_dy, 1.0])/np.linalg.norm([dc_dx, dc_dy, 1.0])
    
    if index:
        return depth, nvec, ix0, iy0
    else:
        return depth, nvec
    
def loadSoil(filename=None, yaml=None, soil_mode='uniform', profile_source=None, dir=None):
        '''
        Load geotechnical information from input file or YAML.
        Supports two soil modes: 'uniform' and 'layered'.

        Parameters
        ----------
        filename : str, optional
            Path to .txt/.dat file with soil labels/profile IDs and coordinates
        yaml : dict, optional
            Dictionary containing soil data and properties (used when filename is None)
        soil_mode : str
            Either 'uniform' or 'layered'
        profile_source : str, optional
            Path to YAML file with layered profile definitions (only used if soil_mode='layered')
        dir : str, optional
            Directory to use for relative file paths
        '''
        xs = None
        ys = None
        soil_names = None
        soilProps = None

        if profile_source is not None and not os.path.isabs(profile_source):
            profile_source = os.path.join(dir, profile_source)

        # Case 1: File input (grid + properties)
        if filename is not None:
            if filename.endswith('.shp'):
                raise ValueError("Shapefiles not supported in Project class")

            elif filename.endswith('.txt') or filename.endswith('.dat'):
                # Load label/profile_id grid
                xs, ys, soil_names = readBathymetryFile(filename, dtype=str)

                # Load soil properties
                soilProps = getSoilTypes(filename, soil_mode=soil_mode, profile_source=profile_source)

            if yaml:
                soilProps = yaml.get('soil_types', soilProps)  # allow overwriting via YAML

        # Case 2: YAML only (no filename)
        elif filename is None:
            if yaml:
                xs = yaml['x']
                ys = yaml['y']
                soil_names = yaml['type_array']
                raw_soil_types = yaml['soil_types']
        
                # Ensure all soil types have a 'layers' field
                soilProps = {}
                for key, entry in raw_soil_types.items():
                    if 'layers' in entry:
                        soilProps[key] = entry 
                    else:
                        # Wrap old flat format into single-layer profile (optional fallback)
                        layer = dict(entry)
                        layer.setdefault('top', 0)
                        layer.setdefault('bottom', 50)
                        layer.setdefault('soil_type', key)
                        soilProps[key] = {'layers': [layer]}
            else:
                print('[Warning] No soil input provided — using default values')
                xs = [0]
                ys = [0]
                soil_names = [['mud']]  # note: should be 2D to match grid structure
                soilProps = {
                    'mud': {'layers': [{
                        'soil_type': 'clay',
                        'top': 0, 'bottom': 50,
                        'gamma_top': 10, 'gamma_bot': 10,
                        'Su_top': 2.39, 'Su_bot': 59.39
                    }]},
                    'rock': {'layers': [{
                        'soil_type': 'rock',
                        'top': 0, 'bottom': 50,
                        'UCS_top': 5, 'UCS_bot': 5,
                        'Em_top': 7, 'Em_bot': 7
                    }]}
                }


        else:
            raise ValueError("Invalid combination of filename/yaml inputs")

        # --- Set defaults only for uniform mode (when values are missing) ---
        if soil_mode == 'uniform':
            for key, props in soilProps.items():
                props['Su0']   = getFromDict(props, 'Su0',   shape=-1, dtype=list, default=[2.39])
                props['k']     = getFromDict(props, 'k',     shape=-1, dtype=list, default=[1.41])
                props['alpha'] = getFromDict(props, 'alpha', shape=-1, dtype=list, default=[0.7])
                props['gamma'] = getFromDict(props, 'gamma', shape=-1, dtype=list, default=[8.7])
                props['phi']   = getFromDict(props, 'phi',   shape=-1, dtype=list, default=[0.0])
                props['UCS']   = getFromDict(props, 'UCS',   shape=-1, dtype=list, default=[7.0])
                props['Em']    = getFromDict(props, 'Em',    shape=-1, dtype=list, default=[50.0])

                # ensure no array-like leftovers
                for k, prop in props.items():
                    if hasattr(prop, '__array__'):
                        props[k] = np.array(prop)

        if xs is not None:
            soil_x = np.array(xs)
            soil_y = np.array(ys)
            soil_names = np.array(soil_names)

        print(f"Loaded soilProps keys: {list(soilProps.keys())}")

        return soilProps, soil_x, soil_y, soil_names, soil_mode








if __name__ == '__main__':
    
    centroid = (40.928, -124.708)  #humboldt    
    xs = np.arange(-30000,30001,400)
    ys = np.arange(-40000,40001,400)
    
    xs, ys, depths = processGeotiff('humboldt.tif', centroid[0], centroid[1], xs=xs, ys=ys, outfilename='test output.txt')
    
    import moorpy as mp
    ms = mp.System(depth=np.max(depths), bathymetry='test output.txt')
    ms.initialize()
    ms.plot(hidebox=True, args_bath={'cmap':'viridis'})
    '''
    # try converting to a different grid
    x_new = np.arange(-20000, 20001, 800)
    y_new = np.arange(-20000, 20001, 800)
    depths_new = resampleGrid(x_new, y_new, xs, ys, depths)
    '''
    plt.show()
