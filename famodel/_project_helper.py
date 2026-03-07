'''
Helper methods for project management.
'''

import numpy as np
from moorpy.helpers import loadLineProps
from famodel.famodel_base import rotationMatrix
from famodel.famodel_base import Node


def _parse_array(array=None, ua=None):
    '''build array info dictionary from uniform array'''
    arrayInfo = []

    if array is not None:
        arrayInfo   = [dict(zip(array['keys'], row)) for row in array['data']]
    
    elif ua is not None:
        WestStart   = ua['west_start'] 
        NorthStart  = ua['north_start']
        xSpacing    = ua['spacing_x']
        ySpacing    = ua['spacing_y']
        topID       = ua['topsideID'] 
        pfID        = ua['platformID']
        moorID      = ua['mooringID']
        pfhead      = ua['heading_adjust']
    
        # get locations of platforms
        xs = WestStart + np.arange(0, ua['n_cols']) * xSpacing
        ys = NorthStart + np.arange(0, ua['n_rows']) * ySpacing
    
        xlocs,ylocs = np.meshgrid(xs,ys)

        outx = np.hstack(xlocs)
        outy = np.hstack(ylocs)

        # Parse the array to table format
        arrayInfo = []

        for i in range(ua['n_rows']*ua['n_cols']):
            arrayInfo.append({'ID':'fowt'+str(i), 'topsideID':topID, 'platformID':pfID,
                                'mooringID':moorID, 'x_location':outx[i], 'y_location':outy[i],
                                'heading_adjust':pfhead})
        
    return arrayInfo

def _parse_all_cable_info(array_cables_data=None, dyn_cable_configs_data=None, cable_types_data=None, cable_appendages_data=None):
    # Variables: cables, dyn_cable_configs, cable_types, cable_appendages as dictionaries.

    # Load dictionary for variable arrayCableInfo or return empty list if not present
    arrayCableInfo = []
    if array_cables_data is not None and 'data' in array_cables_data:
        arrayCableInfo = [dict(zip( array_cables_data['keys'], row))
                         for row in array_cables_data['data']]
        
    # Load dictionary for variable dyn_cable_configs or return empty dict if not present
    dyn_cable_configs = {}
    if dyn_cable_configs_data is not None:
        for k, v in dyn_cable_configs_data.items():
            dyn_cable_configs[k] = v

    # Load dictionary for variable cable_types or return empty dict if not present
    cable_types = {}
    if cable_types_data is not None:
        for k, v in cable_types_data.items():
            cable_types[k] = v
        
    # Load dictionary for variable cable_appendages or return empty dict if not present
    cable_appendages = {}
    if cable_appendages_data is not None:
        for k, v in cable_appendages_data.items():
            cable_appendages[k] = v

    return arrayCableInfo, dyn_cable_configs, cable_types, cable_appendages

def _parse_all_mooring_info(lineProps, array_mooring_data=None, mooring_systems_data=None, mooring_line_types_data=None, mooring_connector_types_data=None, mooring_line_configs_data=None, arrayInfo=None):
    # Variables: arrayMooring, mooring_systems, mooring_line_types, mooring_connector_types, mooring_line_configs as dictionaries.

    # Load dictionary for variable arrayMooring or return empty dict if not present
    arrayMooring = {}
    # for mooring lines: save a list of dictionaries from each row in the data section
    if array_mooring_data and array_mooring_data.get("line_data") is not None:
        arrayMooring = [dict(zip(array_mooring_data['line_keys'], row)) for row in array_mooring_data['line_data']]

    # Load dictionary for variable mooring_line_types or return empty dict if not present
    lineTypes = {}
    if mooring_line_types_data is not None and 'mooring_line_properties_file' in mooring_line_types_data:
        if 'mooring_line_properties_file' in mooring_line_types_data:
            mp_file = mooring_line_types_data['mooring_line_properties_file']
            lineProps = loadLineProps(mp_file)
            # remove this entry to keep everything below working properly
            mooring_line_types_data.pop('mooring_line_properties_file')

        # check if table format was used at all
        if 'keys' and 'data' in mooring_line_types_data: # table-based
            dt = mooring_line_types_data # save location for code clarity
            # save a list of dictionaries from each row in the data section
            ms_info = [dict(zip(dt['keys'], row)) for row in dt['data']]
            # save the list into lineTypes dictionary and rename the index as the linetype name
            for k in range(0,len(ms_info)):
                lineTypes[ms_info[k]['name']] = ms_info[k]
                
        # read in line types from list format as well(will overwrite any repeats from table)
        for k, v in mooring_line_types_data.items():
            # set up line types dictionary
            lineTypes[k] = v
        
    # Load dictionary for variable connectorTypes or return empty dict if not present
    connectorTypes = {}
    if mooring_connector_types_data is not None:
        for k, v in mooring_connector_types_data.items():
            connectorTypes[k] = v

    # Load dictionary for variable lineConfigs or return empty dict if not present
    lineConfigs = {}
    if mooring_line_configs_data is not None:
        for k, v in mooring_line_configs_data.items():
                # set up mooring config
                lineConfigs[k] = v
                # check line types listed in line configs matches those in linetypes section
                if lineTypes: # if linetypes section is included in dictionary
                    for j in range(0,len(v['sections'])): # loop through each line config section
                        if 'type' in v['sections'][j]: # check if it is a connector or line config
                            if not v['sections'][j]['type'] in lineTypes: # check if they match
                                raise Exception(f"Mooring line type '{v['sections'][j]['type']}' listed in mooring_line_configs is not found in mooring_line_types")
        # check line configurations listed in mooring systems matches those in line configs list
        # NOTE: I think other than a check this is dead code. I think mSystems is not the intedend output but rather msys
        mSystems = {}
        if mooring_systems_data: # if mooring_systems section is included in dictionary
            for j,m_s in enumerate(mooring_systems_data): # loop through each mooring system
                for i in range(0, len(arrayInfo)): # loop through each entry in array
                    if m_s == arrayInfo[i]['mooringID']:
                        mSystems[m_s] = [dict(zip(mooring_systems_data[m_s]['keys'], row)) for row in mooring_systems_data[m_s]['data']]
                        for i in range(0,len(mSystems[m_s])): #len(mSystems[m_s]['data'])): # loop through each line listed in the system
                                if not mSystems[m_s][i]['MooringConfigID'] in lineConfigs: # check if they match
                                
                                    raise Exception(f"Mooring line configuration '{mSystems[m_s][i]['MooringConfigID']}' listed in mooring_systems is not found in mooring_line_configs")
                    
    return arrayMooring, mSystems, lineTypes, lineProps, connectorTypes, lineConfigs

def _parse_anchor_data(anchor_data=None, anchor_keys=None, anchor_type_data=None):
    # for anchors: save a list of dictionaries from each row in the data section
    arrayAnchor = {}
    if anchor_data and anchor_keys:
        arrayAnchor = [dict(zip(anchor_keys, row)) for row in anchor_data]
        
    # ----- anchor types -----
    anchorTypes = {}
    if anchor_type_data:
        for k, v in anchor_type_data.items():
            anchorTypes[k] = v

    return arrayAnchor, anchorTypes

def _parse_platform_data(platform_data=None, platforms_data=None):
    RAFTDict    = {} # dictionary for raft platform information
    platforms   = [] # dictionary of platform information

    assert any(x is None for x in (platform_data, platforms_data)), "Cannot read in items for both 'platforms' and 'platform' keywords. Use either 'platform' keyword for one platform or 'platforms' keyword for a list of platforms."
    
    if platform_data:
        # checks
        if type(platform_data) is list and len(platform_data)>1:
            raise Exception("'platform' section keyword must be changed to 'platforms' if multiple platforms are listed")
        
        else:
            if isinstance(platform_data,list):
                platforms.append(platform_data[0])
            else:
                platforms.append(platform_data)
            RAFTDict['platform'] = platform_data

    # load list of platform dictionaries into RAFT dictionary
    elif platforms_data:
        Warning('I believe this method is not correct. You extend dicts not lists')
        platforms.extend(platforms_data)
        RAFTDict['platforms'] = platforms_data

    return platforms, RAFTDict
        
def _parse_topsides(topside_data=[]):
    turbines    = []
    substations = []

    for ts in topside_data:
        if 'TURBINE' in ts['type'].upper():
            turbines.append(ts)
        elif 'SUBSTATION' in ts['type'].upper():
            substations.append(ts)

    return turbines, substations


def _build_platforms_instance(arrayInfo_row, platformTypes, platformList):
    """
    Build Platform instances from arrayInfo.

    Parameters
    ----------
    arrayInfo_row : dict
        Parsed array row (includes platformID, ID, x_location, y_location, etc.).
    platformTypes : list[dict]
        Platform type definitions (RAFT-style). This is the same object as Project.platformTypes.
    platformList : dict[str, Platform]
        Dict to be populated/updated in-place. This is the same object as Project.platformList.

    Returns
    -------
    platforms : list[Platform]
        Platforms created in the same order as arrayInfo.
    """
    from famodel.platform.platform import Platform

    pfID = int(arrayInfo_row["platformID"] - 1)
    pf_type = platformTypes[pfID]

    # Position (z priority: row -> platformType -> 0)
    if "z_location" in arrayInfo_row:
        r = [arrayInfo_row["x_location"], arrayInfo_row["y_location"], arrayInfo_row["z_location"]]
    elif "z_location" in pf_type:
        r = [arrayInfo_row["x_location"], arrayInfo_row["y_location"], pf_type["z_location"]]
    else:
        r = [arrayInfo_row["x_location"], arrayInfo_row["y_location"], 0]

    hydrostatics = pf_type.get("hydrostatics", {}) or {}

    # ID: required by Platform.addPlatform. If missing, make one here.
    pid = arrayInfo_row.get("ID", None)
    if pid is None:
        pid = f"fowt{len(platformList)}"

    # Optional: allow per-row RAFT platform dict override (keeps your old behavior if needed)
    # If you don't want platformTypes mutated here, delete this block.
    raft_dict = arrayInfo_row.get("raft_platform_dict", {}) or {}
    platform_type_index = pfID
    if raft_dict:
        raft_dict = dict(raft_dict)     # defensive copy
        raft_dict["type"] = pf_type.get("type", "")
        platformTypes.append(raft_dict)
        platform_type_index = len(platformTypes) - 1

    platform = Platform.addPlatform(
        r=r,
        id=pid,
        phi=arrayInfo_row.get("heading_adjust", 0),
        entity=pf_type.get("type", ""),
        rFair=pf_type.get("rFair", 0),
        zFair=pf_type.get("zFair", 0),
        platform_type=platform_type_index,
        hydrostatics=hydrostatics,
    )

    platformList[platform.id] = platform

    return platform

def _build_fairleads_list(platformTypes, pfID, platform):
    '''Build list of Fairlead instances for a platform. 
    Calling the class method Fairlead.addFairlead'''

    from famodel.platform.fairlead import Fairlead

    pf_fairs = []
    fct = 0 # fairlead count for id numbering

    for fl in platformTypes[pfID].get("fairleads", []):
        # if headings provided, adjust r_rel with headings
        if 'headings' in fl:
            for head in fl['headings']:
                R = rotationMatrix(0,0,np.radians(90-head))
                # apply to unrotated r_rel
                r_rel = np.matmul(R, fl['r_rel'])
                pf_fairs.append(Fairlead.addFairlead(id=platform.id+'_F'+str(fct+1), 
                                                    platform=platform, 
                                                    r_rel=r_rel))
                fct += 1
        # otherwise, just use r_rel as-is
        elif 'r_rel' in fl:
            pf_fairs.append(Fairlead.addFairlead(id=platform.id+'_F'+str(fct+1), 
                                                platform=platform, 
                                                r_rel=fl['r_rel']))
            fct += 1

    return pf_fairs

def _build_jtubes_list(platformTypes, pfID, platform):
    '''Build list of Jtube instances for a platform. 
    Calling the class method Jtube.addJtube'''

    from famodel.cables.components import Jtube

    pf_jtubes = []

    for jct, jt in enumerate(platformTypes[pfID].get("JTubes", [])):
        # if headings provided, adjust r_rel with headings
        if 'headings' in jt:
            for head in jt['headings']:
                R = rotationMatrix(0,0,np.radians(90-head))
                # apply to unrotated r_rel
                r_rel = np.matmul(R, jt['r_rel'])
                pf_jtubes.append(Jtube.addJtube(id=platform.id+'_JT'+str(jct+1), 
                                                platform=platform, 
                                                r_rel=r_rel))
                jct += 1
        # otherwise, just use r_rel as-is
        elif 'r_rel' in jt:
            pf_jtubes.append(Jtube.addJtube(id=platform.id+'_JT'+str(jct+1), 
                                            platform=platform, 
                                            r_rel=jt['r_rel']))
            jct += 1

    return pf_jtubes


def _build_turbine_and_substation_instances(
    arrayInfo_row, row_index, platform, topsides, turbineList, turbineTypes, substationList):

    from famodel.turbine.turbine import Turbine
    from famodel.substation.substation import Substation
    from famodel.famodel_base import Node  # wherever Node is

    entity      = (platform.entity or "").upper()
    topside_id  = arrayInfo_row.get("topsideID", 0)

    if not topside_id or topside_id <= 0:
        return

    # resolve topside_dd
    if isinstance(topsides, list):
        topside_dd = topsides[topside_id - 1]
    else:
        topside_dd = topsides

    if entity == "FOWT":
        # Define name and ID for turbine/substations 
        name    = f"T{topside_id}_{row_index}"
        typeID  = topside_id

        # rotor diameter + turbineTypes append (matches original)
        rotor_diameter = 0

        if topside_dd and "blade" in topside_dd:
            # Update project.turbineTypes
            turbineTypes.append(topside_dd)

            # Define the rotor diameter for the turbine initialization
            blade = topside_dd["blade"]

            if isinstance(blade, list):
                rotor_diameter = blade[0]["Rtip"] * 2
            else:
                rotor_diameter = blade["Rtip"] * 2

        # Update project.turbineList
        turbine_instance    = Turbine.addTurbine(dd=topside_dd, name=name, D=rotor_diameter, typeID=typeID)
        turbineList[name]   = turbine_instance

        return turbine_instance

    if entity == "SUBSTATION":
        name                    = f"S{topside_id}_{row_index}"
        
        # Update project.substationList
        substationInstance      = Substation.addSubstation(dd=topside_dd, name=name)
        substationList[name]    = substationInstance
        
        return substationInstance

    # fallback node (When we add WECs or other topside types, we can expand this if/else block)
    name    = f"N{topside_id}_{row_index}"
    node    = Node(name)
    node.dd = topside_dd
    platform.attach(node)
