'''
Helper methods for project management.
'''
import os
import numpy as np
from moorpy.helpers import loadLineProps
from fad.famodel_base import rotationMatrix
from fad.famodel_base import Node
import fad.seabed_tools as sbt
from .helpers import getAnchors, getMoorings, attachFairleads, calc_heading, calc_midpoint, \
                    getDynamicCables, getStaticCables, getFromDict

from fad.cables.cable import Cable

from copy import deepcopy
from matplotlib import pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap


def _parse_array(array=None, ua=None):
    '''build array info dictionary from uniform array'''
    arrayInfo = []

    if array is not None and array.get('data') is not None:
        arrayInfo   = [dict(zip(array['keys'], row)) for row in array['data']]
    
    elif array is not None and array.get('data') is None:
        return arrayInfo
    
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
    if mooring_line_types_data is not None:
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
    mSystems = {}
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
    from fad.platform.platform import Platform

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

    from fad.platform.fairlead import Fairlead

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

    from fad.cables.components import Jtube

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

    from fad.turbine.turbine import Turbine
    from fad.substation.substation import Substation
    from fad.famodel_base import Node  # wherever Node is

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

def _build_mooring_anchor_instance_nonshared(arrayInfo_row, 
                            mSystems,
                            lineConfigs,
                            alph, 
                            connectorTypes,
                            lineProps,
                            lineTypes,
                            platform,
                            platFormTypes_pfID,
                            pf_fairs,
                            anchorTypes,
                            **kwargs):
    """
    Creates one Mooring + one Anchor instance.
    - Does NOT attach ends (leave attachTo in project.py as requested)
    - DOES register the mooring in project.mooringList (passed by reference on project)
    """

    from fad.mooring.mooring import Mooring  # adjust import if needed
    from fad.anchors.anchor import Anchor

    mooringList = {}
    anchorList = {}

    if lineConfigs and mSystems and arrayInfo_row['mooringID'] != 0: #if not fully shared mooring on this platform
        m_s = arrayInfo_row['mooringID'] # get mooring system ID

        # get the mooring system dictionary for this platform
        mySys = mSystems[m_s] # get the mooring system for this platform from the mooring_systems section of the input dictionary
        # get mooring headings (need this for platform class)
        headings = []
        for ii in range(0,len(mySys)):
            headings.append(np.radians(mySys[ii]['heading']))
        
        # get the mooring line information 
        for j in range(0,len(mySys)): # loop through each line in the mooring system
            # - - -  create anchor first
            # set anchor info
            lineAnch = mySys[j]['anchorType'] # get the anchor type for the line
            
            # add anchor class instance to anchorList in project class
            array_name = str(arrayInfo_row['ID'])+alph[j]
        
            # get the configuration for that line in the mooring system
            mooring_name = mySys[j]['MooringConfigID']
        
            # create mooring and connector dictionary
            mdd = getMoorings(mooring_name, lineConfigs, 
                                connectorTypes, array_name,
                                lineProps=lineProps, 
                                lineTypes=lineTypes, 
                                rho=kwargs['rho_water'], g=kwargs['g'])
            
            # Add entries to mooring dictionary that depend on project class and platform info
            mdd['zAnchor']  = -kwargs['depth'] 
            mdd['rad_fair'] = platform.rFair
            mdd['z_fair']   = platform.zFair
            
            moor = Mooring.addMooring(id=array_name,
                                        rel_heading=headings[j],
                                        dd=mdd, 
                                        lineProps=lineProps)
            
            # NOTE: the following methods are supposed to return something which is not used
            if 'fairlead' in mySys[j]:
                attachFairleads(moor,
                                1,
                                platform,
                                fair_ID_start=platform.id+'_F',
                                fair_inds=mySys[j]['fairlead'])
                
            elif 'rFair' in platFormTypes_pfID and 'zFair' in platFormTypes_pfID:
                moor.attachTo(platform, 
                                r_rel=[platform.rFair,
                                        0,
                                        platform.zFair], 
                                end='b')
                
            elif pf_fairs:
                attachFairleads(moor,
                                1,
                                platform,
                                fair_ID = pf_fairs[j].id)

            else:
                print('Warning: platform definition did not include either rFair & zFair or fairlead definitions. Assuming 0 fairlead radius and depth.')
            
            
            # Position the subcomponents along the Mooring
            moor.positionSubcomponents()

            # Anchor stuff
            ad, mass = getAnchors(lineAnch, anchorTypes) # call method to create anchor dictionary
            anch = Anchor.addAnchor(id=array_name, dd=ad, mass=mass)

            # attach ends
            moor.attachTo(anch, end='A')

            # Add mooring instance to mooringList in project class
            mooringList[moor.id] = moor
            anchorList[anch.id] = anch

    return mooringList, anchorList

def _build_array_cable_instances(
    arrayCableInfo,
    dyn_cable_configs,
    cable_types,
    cable_appendages,
    platformList,
    jtube_by_platform,
    *,
    depth,
    rho_water,
    g,
    start_index=0,
):
    """
    Build Cable instances from the array_cables table.

    Notes
    -----
    - Does NOT require the Project instance.
    - Does NOT reposition cables.
    - Returns:
        1) cable_dict: {cable_id: Cable}
        2) reposition_args: {cable_id: {'rad_fair': [rJTubeA, rJTubeB]}}
    """

    cable_dict = {}
    reposition_args = {}

    if not arrayCableInfo:
        return cable_dict, reposition_args

    for i, cab in enumerate(arrayCableInfo):
        A = None
        dd = {"cables": [], "joints": []}

        dyn_cabA = cab["DynCableA"] if "NONE" not in cab["DynCableA"].upper() else None
        dyn_cabB = cab["DynCableB"] if "NONE" not in cab["DynCableB"].upper() else None
        stat_cab = cab["cableType"] if "NONE" not in cab["cableType"].upper() else None

        JtubeA = cab["JTubeA"] if "JTubeA" in cab else None
        JtubeB = cab["JTubeB"] if "JTubeB" in cab else None

        rJTubeA = None
        rJTubeB = None

        attachA = cab["AttachA"]
        attachB = cab["AttachB"]

        # Keep the same assumption as the current code:
        # headings are taken from platform phi, so these attach IDs must exist in platformList.
        if attachA not in platformList:
            raise Exception(
                f"AttachA {attachA} for array cable {i} must match a platform "
                f"to compute heading/phi in this builder."
            )
        if attachB not in platformList:
            raise Exception(
                f"AttachB {attachB} for array cable {i} must match a platform "
                f"to compute heading/phi in this builder."
            )

        A_phi = platformList[attachA].phi
        B_phi = platformList[attachB].phi

        if dyn_cabA:
            Acondd, jAcondd = getDynamicCables(
                dyn_cable_configs[dyn_cabA],
                cable_types,
                cable_appendages,
                depth,
                rho_water=rho_water,
                g=g,
            )

            if stat_cab or dyn_cabB:
                dd["joints"].append(jAcondd)
            else:
                # suspended cable
                Acondd["headingB"] = np.radians(cab["headingB"]) + B_phi

            Acondd["headingA"] = np.radians(cab["headingA"]) + A_phi

            if "rJTube" in dyn_cable_configs[dyn_cabA]:
                rJTubeA = dyn_cable_configs[dyn_cabA]["rJTube"]
                Acondd["rJTube"] = rJTubeA

            dd["cables"].append(Acondd)
            A = Acondd["A"]

        if stat_cab:
            dd["cables"].append(
                getStaticCables(
                    stat_cab,
                    cable_types,
                    rho_water=rho_water,
                    g=g,
                    A=A,
                )
            )

        if dyn_cabB:
            Bcondd, jBcondd = getDynamicCables(
                dyn_cable_configs[dyn_cabB],
                cable_types,
                cable_appendages,
                depth,
                rho_water=rho_water,
                g=g,
            )

            if "rJTube" in dyn_cable_configs[dyn_cabB]:
                rJTubeB = dyn_cable_configs[dyn_cabB]["rJTube"]
                Bcondd["rJTube"] = rJTubeB

            Bcondd["headingB"] = np.radians(cab["headingB"]) + B_phi
            dd["cables"].append(Bcondd)
            dd["joints"].append(jBcondd)

        cable_id = f"cable{start_index + i}"
        cable = Cable(cable_id, d=dd)

        # Attach end A
        if attachA in platformList:
            if jtube_by_platform.get(attachA) and JtubeA:
                cable.subcomponents[0].attachTo(
                    jtube_by_platform[attachA][JtubeA - 1],
                    end="A",
                )
            else:
                cable.attachTo(platformList[attachA], end="A")
        elif attachA in cable_appendages:
            pass
        else:
            raise Exception(
                f"AttachA {attachA} for array cable {i} does not match any "
                f"platforms or appendages."
            )

        # Attach end B
        if attachB in platformList:
            if jtube_by_platform.get(attachB) and JtubeB:
                cable.subcomponents[-1].attachTo(
                    jtube_by_platform[attachB][JtubeB - 1],
                    end="B",
                )
            else:
                cable.attachTo(platformList[attachB], end="B")
        elif attachB in cable_appendages:
            pass
        else:
            raise Exception(
                f"AttachB {attachB} for array cable {i} does not match any "
                f"platforms or appendages."
            )

        cable_dict[cable_id] = cable
        reposition_args[cable_id] = {"rad_fair": [rJTubeA, rJTubeB]}

    return cable_dict, reposition_args

def _build_cable_instances_from_cables_section(
    cableInfo,
    dyn_cable_configs,
    cable_types,
    cable_appendages,
    platformList,
    jtube_by_platform,
    *,
    depth,
    rho_water,
    g,
    start_index=0,
):
    """
    Build Cable instances from the higher-level `cables` section.

    Notes
    -----
    - Does NOT require the Project instance.
    - Does NOT call reposition.
    - Returns:
        1) cable_dict: {cable_id: Cable}
        2) reposition_args: {cable_id: {'rad_fair': [rJTubeA, rJTubeB]}}
    """

    cable_dict = {}
    reposition_args = {}

    if not cableInfo:
        return cable_dict, reposition_args

    for i, cab in enumerate(cableInfo):
        rJTubeA = None
        rJTubeB = None
        A = None

        JtubeA = cab["endA"]["JTube"] if "JTube" in cab["endA"] else None
        JtubeB = cab["endB"]["JTube"] if "JTube" in cab["endB"] else None

        dd = {"cables": [], "joints": []}

        dyn_cabA = (
            cab["endA"]["dynamicID"]
            if "NONE" not in cab["endA"]["dynamicID"].upper()
            else None
        )
        dyn_cabB = (
            cab["endB"]["dynamicID"]
            if "NONE" not in cab["endB"]["dynamicID"].upper()
            else None
        )
        stat_cab = cab["type"] if "NONE" not in cab["type"].upper() else None

        attachA = cab["endA"]["attachID"]
        attachB = cab["endB"]["attachID"]

        if attachA not in platformList:
            raise Exception(
                f"AttachA {attachA} for cable {cab['name']} must match a platform "
                f"to compute heading/phi in this builder."
            )
        if attachB not in platformList:
            raise Exception(
                f"AttachB {attachB} for cable {cab['name']} must match a platform "
                f"to compute heading/phi in this builder."
            )

        A_phi = platformList[attachA].phi
        B_phi = platformList[attachB].phi

        # end A dynamic section
        if dyn_cabA:
            Acondd, jAcondd = getDynamicCables(
                dyn_cable_configs[dyn_cabA],
                cable_types,
                cable_appendages,
                depth,
                rho_water=rho_water,
                g=g,
            )

            if stat_cab or dyn_cabB:
                dd["joints"].append(jAcondd)
            else:
                # suspended cable
                Acondd["headingB"] = np.radians(cab["endB"]["heading"]) + B_phi

            Acondd["headingA"] = np.radians(cab["endA"]["heading"]) + A_phi

            if "rJTube" in dyn_cable_configs[dyn_cabA]:
                rJTubeA = dyn_cable_configs[dyn_cabA]["rJTube"]
                Acondd["rJTube"] = rJTubeA

            dd["cables"].append(Acondd)
            # A = dyn_cable_configs[dyn_cabA].get("A", Acondd.get("A")) #old version
            A = dyn_cable_configs[dyn_cabA]['A']

        # static section
        if stat_cab:
            statcondd = getStaticCables(
                stat_cab,
                cable_types,
                rho_water=rho_water,
                g=g,
                A=A,
            )

            if "routing_x_y_r" in cab and cab["routing_x_y_r"]:
                statcondd["routing"] = cab["routing_x_y_r"]

            if "burial" in cab and cab["burial"]:
                statcondd["burial"] = cab["burial"]

            dd["cables"].append(statcondd)

        # end B dynamic section
        if dyn_cabB:
            Bcondd, jBcondd = getDynamicCables(
                dyn_cable_configs[dyn_cabB],
                cable_types,
                cable_appendages,
                depth,
                rho_water=rho_water,
                g=g,
            )

            Bcondd["headingB"] = np.radians(cab["endB"]["heading"]) + B_phi

            if "rJTube" in dyn_cable_configs[dyn_cabB]:
                rJTubeB = dyn_cable_configs[dyn_cabB]["rJTube"]
                Bcondd["rJTube"] = rJTubeB

            dd["cables"].append(Bcondd)
            dd["joints"].append(jBcondd)

        cable_id = 'cable' + f"{start_index + i}"
        cable = Cable(cable_id, d=dd)

        # attach end A
        if attachA in platformList:
            if jtube_by_platform.get(attachA) and JtubeA is not None:
                cable.subcomponents[0].attachTo(
                    jtube_by_platform[attachA][JtubeA],
                    end="A",
                )
            else:
                cable.attachTo(platformList[attachA], end="A")
        elif attachA in cable_appendages:
            pass
        else:
            raise Exception(
                f"AttachA {attachA} for cable {cab['name']} does not match any "
                f"platforms or appendages."
            )

        # attach end B
        if attachB in platformList:
            if jtube_by_platform.get(attachB) and JtubeB is not None:
                cable.subcomponents[-1].attachTo(
                    jtube_by_platform[attachB][JtubeB],
                    end="B",
                )
            else:
                cable.attachTo(platformList[attachB], end="B")
        elif attachB in cable_appendages:
            pass
        else:
            raise Exception(
                f"AttachB {attachB} for cable {cab['name']} does not match any "
                f"platforms or appendages."
            )

        cable_dict[cable_id] = cable
        reposition_args[cable_id] = {"rad_fair": [rJTubeA, rJTubeB]}

    return cable_dict, reposition_args

def _build_raft_dict(
    raft,
    RAFTDict,
    turbineTypes,
    arrayInfo,
    site_data,
    *,
    depth,
    rho_water,
    rho_air,
    mu_air,
    name="Project_Array",
):
    """
    Build and return the RAFT input dictionary.

    Notes
    -----
    - Does NOT require the Project instance.
    - Does NOT call self.getRAFT(...).
    - Returns a new dictionary; does not mutate the input RAFTDict in place.
    """

    raft_dict = deepcopy(RAFTDict) if RAFTDict is not None else {}

    if not raft:
        return raft_dict

    # turbine(s)
    nt = len(turbineTypes) if turbineTypes is not None else 0
    if nt == 1:
        raft_dict["turbine"] = turbineTypes[0]
    elif nt > 1:
        raft_dict["turbines"] = turbineTypes

    # settings / cases
    if site_data.get("RAFT_settings"):
        raft_dict["settings"] = site_data["RAFT_settings"]

    if site_data.get("RAFT_cases"):
        raft_dict["cases"] = site_data["RAFT_cases"]

    # array table
    if arrayInfo:
        raftinfo = deepcopy(arrayInfo)

        for row in raftinfo:
            if "heading_adjust" in row:
                row["heading_adjust"] = -row["heading_adjust"]

        raft_dict["array"] = {
            "keys": list(raftinfo[0].keys()),
            "data": [list(row.values()) for row in raftinfo],
        }

    # site block
    raft_dict["site"] = {
        "water_depth": depth,
        "rho_water": rho_water,
        "rho_air": rho_air,
        "mu_air": mu_air,
    }
    raft_dict["site"]["shearExp"] = getFromDict(
        site_data.get("general", {}),
        "shearExp",
        default=0.12,
    )

    raft_dict["name"] = name
    raft_dict["type"] = "input file for RAFT"

    return raft_dict

#### Plotting helper function ####
# NOTE: later move to a separate plotting helper file
def _nearest_soil_names_on_bathymetry_grid(grid_x, grid_y, soil_x, soil_y, soil_names):
    """
    Map categorical soil names onto the bathymetry grid using nearest-neighbor lookup.

    Parameters
    ----------
    grid_x, grid_y : 1D array-like
        Bathymetry grid coordinates.
    soil_x, soil_y : 1D array-like
        Soil grid coordinates.
    soil_names : 2D array-like of str
        Soil labels on the soil grid, shape (len(soil_y), len(soil_x)).

    Returns
    -------
    soil_on_bathy : np.ndarray[str]
        Soil labels resampled onto the bathymetry grid,
        shape (len(grid_y), len(grid_x)).
    """
    soil_on_bathy = np.empty((len(grid_y), len(grid_x)), dtype=object)

    for iy, y in enumerate(grid_y):
        jy = int(np.argmin(np.abs(np.asarray(soil_y) - y)))
        for ix, x in enumerate(grid_x):
            jx = int(np.argmin(np.abs(np.asarray(soil_x) - x)))
            soil_on_bathy[iy, ix] = soil_names[jy, jx]

    return soil_on_bathy

# TODO: move to a separate script eventually
FOLK7_SOIL_COLOR_MAP = {
    '-':    '#f2f2f2',  # empty / outside
    '11.0': '#6f78d8',  # Mud
    '12.0': '#8edcc4',  # Sandy Mud
    '13.0': '#d7ee72',  # Muddy Sand
    '2.0':  '#f1e36b',  # Sand
    '3.0':  '#b7ad2f',  # Coarse substrate
    '4.0':  '#d9b6a3',  # Mixed sediment
    '5.0':  '#a31717',  # Rock & Boulders
    '6.0':  '#ffffff',  # No data
    '9.0':  '#d9d2e9',  # Restricted data
    '10.0': '#eeeeee'   # Unpublic data
}


FOLK7_SOIL_LABELS_PRETTY = {
    '11.0': '1.1 Mud',
    '12.0': '1.2 Sandy Mud',
    '13.0': '1.3 Muddy Sand',
    '2.0':  '2 Sand',
    '3.0':  '3 Coarse substrate',
    '4.0':  '4 Mixed sediment',
    '5.0':  '5 Rock & Boulders',
    '6.0':  '6 No data',
    '9.0':  '9 Restricted data',
    '10.0': '10 Unpublic data',
    '-':    'Outside'
}

def _build_soil_facecolors(soil_names_2d, cmap_soil=None, soil_alpha=0.85):
    """
    Convert a 2D soil label array into RGBA facecolors plus categorical legend data.

    Colors are taken from FOLK7_SOIL_COLOR_MAP and tick labels from
    FOLK7_SOIL_LABELS_PRETTY. The ``cmap_soil`` argument is ignored.

    Parameters
    ----------
    soil_names_2d : 2D array-like of str
        Soil labels on the target plotting grid.
    cmap_soil : ignored
        Kept for backward compatibility; colors come from FOLK7_SOIL_COLOR_MAP.
    soil_alpha : float, optional
        Alpha channel to apply to all soil colors.

    Returns
    -------
    facecolors : np.ndarray
        RGBA array with shape soil_names_2d.shape + (4,).
    soil_types_pretty : list[str]
        Pretty labels for each unique soil category (suitable for colorbar ticks).
    soil_type_to_int : dict
        Mapping from raw soil key to integer category index.
    cmap_obj : ListedColormap
        Matplotlib colormap built from FOLK7_SOIL_COLOR_MAP.
    norm : BoundaryNorm
        Norm suitable for discrete categorical colorbars.
    """
    soil_names_2d = np.asarray(soil_names_2d)

    # Preserve first-appearance order instead of np.unique sorting
    soil_types_raw = list(dict.fromkeys(soil_names_2d.ravel().tolist()))

    # Build a ListedColormap using the fixed FOLK7 hex colors
    hex_colors = [FOLK7_SOIL_COLOR_MAP.get(s, '#cccccc') for s in soil_types_raw]
    cmap_obj = ListedColormap(hex_colors)

    soil_type_to_int = {name: i for i, name in enumerate(soil_types_raw)}
    soil_int = np.vectorize(soil_type_to_int.get)(soil_names_2d)

    bounds = np.arange(len(soil_types_raw) + 1)
    norm = BoundaryNorm(bounds, cmap_obj.N)

    facecolors = cmap_obj(norm(soil_int))
    facecolors[..., 3] = soil_alpha

    # Resolve pretty labels for colorbar tick display
    soil_types_pretty = [FOLK7_SOIL_LABELS_PRETTY.get(s, s) for s in soil_types_raw]

    return facecolors, soil_types_pretty, soil_type_to_int, cmap_obj, norm

def _build_bathy_grid_and_depth(site, depth, dir, defaults_dict, interpolate=False):
    grid_x = defaults_dict.get('grid_x', np.array([0]))
    grid_y = defaults_dict.get('grid_y', np.array([0]))
    grid_depth = defaults_dict.get('grid_depth', np.array([[depth]]))

    bathy_dict = site.get('bathymetry', None)

    if not bathy_dict:
        return grid_x, grid_y, grid_depth
    
    # load bathymetry information, if provided
    if 'file' in bathy_dict and bathy_dict['file']: # make sure there was a file provided even if the key is there
        filename = bathy_dict['file']
        
        # if it's a relative file location, specify the root directory
        if not os.path.isabs(filename): 
            filename = os.path.join(dir, filename)
        grid_x, grid_y, grid_depth = sbt.loadBathymetry(filename, interpolate=interpolate)
        
    elif 'x' in bathy_dict and 'y' in bathy_dict:
        grid_x = np.array(bathy_dict['x'])
        grid_y = np.array(bathy_dict['y'])
        grid_depth = np.array(bathy_dict['depths'])
        # setGrid(xs,ys)
    else:
        # assume a flat bathymetry
        grid_depth  = np.array([[depth]])
        sbt.setGrid([0],[0], grid_x, grid_y, grid_depth)

    return grid_x, grid_y, grid_depth

def _build_boundary(site, lat0, lon0, default_boundary):
    boundary_dict = site.get('boundaries', None)

    if boundary_dict is None:
        return default_boundary

    if 'file' in boundary_dict and boundary_dict['file']:  # load boundary data from file if filename provided
            # TODO: move to project helper. By passing lon and lat, the method is project class independent.
            boundary = sbt.loadBoundary(boundary_dict['file'], lat0, lon0)

    elif 'x_y' in boundary_dict and boundary_dict['x_y']:  # process list of boundary x,y vertices
            xy = boundary_dict['x_y']
            boundary = np.zeros([len(xy),2])
            for i in range(len(xy)):
                boundary[i,0] = float(xy[i][0])
                boundary[i,1] = float(xy[i][1])

    else:
        raise ValueError("\nBoundary information provided in site block, but no valid 'file' or 'x_y' key found. Check your input data.")
    
    return boundary

def _build_soil(site, dir, defaults_dict):
    soilProps   = defaults_dict.get('soilProps', None)
    soil_x      = defaults_dict.get('soil_x', None)
    soil_y      = defaults_dict.get('soil_y', None)
    soil_names  = defaults_dict.get('soil_names', None)
    soil_mode   = defaults_dict.get('soil_mode', None)

    seabed_dict = site.get('seabed', None)

    if seabed_dict is None:
        return soilProps, soil_x, soil_y, soil_names, soil_mode

    # if there's a file listed in the seabed dictionary
    if 'file' in seabed_dict and seabed_dict['file']:
        filename = seabed_dict['file']

        # In case relative file location fails, use the root directory
        if not os.path.isabs(filename): 
            filename = os.path.join(dir, filename)

        # without reading the file to tell whether it has soil property information listed, check to see if soil property information is given
        if 'soil_types' in seabed_dict and seabed_dict['soil_types']:     # if the yaml does have soil property information
            soilProps, soil_x, soil_y, soil_names, soil_mode = sbt.loadSoil(filename=filename, yaml=seabed_dict, soil_mode='uniform', dir=dir)

        elif 'profile_source' in seabed_dict and seabed_dict['profile_source']:
            soilProps, soil_x, soil_y, soil_names, soil_mode = sbt.loadSoil(filename=filename, soil_mode='layered', profile_source=str(seabed_dict['profile_source']), dir=dir)
            
        else:       # if the yaml doesn't have soil property information, read in just the filename to get all the information out of that
            soilProps, soil_x, soil_y, soil_names, soil_mode = sbt.loadSoil(filename=filename)
    # if there's no file listed in the seabed dictionary, load in just the yaml information (assuming the ['x', 'y', and 'type_array'] information are there)
    else:
        # NOTE: this is missing specifying information such as soil_mode at least
        soilProps, soil_x, soil_y, soil_names, soil_mode = sbt.loadSoil(yaml=seabed_dict)

    return soilProps, soil_x, soil_y, soil_names, soil_mode

def _build_marine_growth(site):
    marine_growth_dict  = site.get('marine_growth', None)
    marine_growth       = None
    marine_growth_buoys = None

    if marine_growth_dict is None or 'data' not in marine_growth_dict:
        return None, None

    marine_growth = [dict(zip(marine_growth_dict['keys'], row))
                            for row in marine_growth_dict['data']]
    
    if 'buoys' in marine_growth_dict:
        marine_growth_buoys = marine_growth_dict['buoys']

    return marine_growth, marine_growth_buoys





