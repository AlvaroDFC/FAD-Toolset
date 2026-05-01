# class for cable components

import numpy as np
from fad.famodel_base import Node, Edge


class Joint(Node, dict):
    '''Subsea joint for power cables, such as might join dynamic and static cables together.
    '''
    def __init__(self,id, r=None,**kwargs):
        '''
        Connectors inherit from dict, so properties can be passed in as arguments
        and they will be assigned like dictionary entries. 
        
        Parameters
        ----------
        dd : dictionary, optional
            Design dictionary The default is None.
        r : list, optional
            x,y,z location of the connector. The default is [0,0,0].
        kwargs
            Additional optional parameters, such as m
        '''
        from fad.project import getFromDict
        dict.__init__(self, **kwargs)  # initialize dict base class (will put kwargs into self dict)
        Node.__init__(self, id)  # initialize Node base class
        # Joint position and orientation

        # set defaults if they weren't provided
        self['m'  ] = getFromDict(self, 'm'  , default=0)
        
        # MoorPy Point Object for Joint
        self.mpConn = None
        
        # Dictionary for failure probability
        self.failure_probability = {}
        
    #this might be useful for the ends of dynamic cables
    def makeMoorPyConnector(self, ms):
        '''Create a MoorPy connector object in a MoorPy system
        Parameters
        ----------
        ms : class instance
            MoorPy system
        
        Returns
        -------
        ms : class instance
            MoorPy system 

        '''
        from fad.project import getFromDict
        # create connector as a point in MoorPy system
        ms.addPoint(1,self['r'])
        # assign this point as mpConn in the anchor class instance
        self.mpConn = ms.pointList[-1]

        self.mpConn.m = getFromDict(self,'m',default=10000)
        self.mpConn.v = getFromDict(self,'v',default=0)
        self.mpConn.CdA = getFromDict(self,'CdA',default=0)

        return(ms)

id="xytg13"
class Jtube(Node, dict):
    def __init__(self, id, r=None, **kwargs):
        dict.__init__(self, **kwargs)
        Node.__init__(self, id)

    @classmethod
    def addJtube(
        cls,
        id=None,
        platform=None,
        r_rel=(0, 0, 0),
        cable=None,
        end="b",
        **kwargs,
    ):
        """
        Create a Jtube object and optionally attach it to a platform and/or a cable.

        Parameters
        ----------
        id : str | None
            Jtube id. If None and platform is provided, an id is generated.
        platform : Node-like | None
            Parent platform to attach under (expects .attach(), .id, .attachments).
        r_rel : array-like length 3
            Location relative to platform.
        cable : object | None
            Cable object with .subcomponents and attachTo() on end connectors.
        end : {'a','b',0,1,'A','B'}
            Which end of the cable to connect to.
        kwargs : dict
            Extra fields to store in the Jtube dict (since Jtube inherits dict).
        """
        # create an id if needed
        if id is None:
            if platform is None:
                raise ValueError("Jtube.addJtube: id is None and platform is None, cannot auto-generate id.")
            id = platform.id + str(len(platform.attachments))

        # create J-tube object (store kwargs into dict)
        jt = cls(id=id, **kwargs)

        # attach subordinately to platform and provide relative location
        if platform is not None:
            platform.attach(jt, r_rel=r_rel)

        # attach equally to cable end connector
        if cable is not None:
            if end in ["a", "A", 0]:
                cable.subcomponents[0].attachTo(jt)
            elif end in ["b", "B", 1]:
                cable.subcomponents[-1].attachTo(jt)
            else:
                raise ValueError(f"Jtube.addJtube: invalid end={end!r}")

        return jt