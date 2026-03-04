from famodel.famodel_base import Node

class Fairlead(Node):
    def __init__(self, id):
        super().__init__(id)

    @classmethod
    def addFairlead(
        cls,
        id=None,
        platform=None,
        r_rel=(0, 0, 0),
        mooring=None,
        end="b",
    ):
        """
        Create a Fairlead object and optionally attach it to a platform and/or a mooring.

        Parameters
        ----------
        id : str | None
            Fairlead id. If None and platform is provided, an id is generated.
        platform : Node-like | None
            Parent platform to attach under (expects .attach(), .id, .r, .attachments).
        r_rel : array-like length 3
            Location relative to platform.
        mooring : object | None
            Mooring object with .subcomponents and join() on end connectors.
        end : {'a','b',0,1,'A','B'}
            Which end of the mooring to connect to.
        """
        # create an id if needed
        if id is None:
            if platform is None:
                raise ValueError("Fairlead.addFairlead: id is None and platform is None, cannot auto-generate id.")
            id = platform.id + str(len(platform.attachments))

        fl = cls(id=id)

        # attach subordinately to platform and provide relative location
        if platform is not None:
            platform.attach(fl, r_rel=r_rel)
            # ensure vector math works (platform.r might be numpy)
            fl.r = platform.r + r_rel  # absolute location

        # attach equally to mooring end connector
        if mooring is not None:
            if end in ["a", "A", 0]:
                mooring.subcomponents[0].join(fl)
            elif end in ["b", "B", 1]:
                mooring.subcomponents[-1].join(fl)
            else:
                raise ValueError(f"Fairlead.addFairlead: invalid end={end!r}")

        return fl