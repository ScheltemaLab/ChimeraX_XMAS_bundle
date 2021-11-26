from chimerax.atomic.pbgroup import Pseudobonds

class Pseudobonds:

    @property
    def dashes(self, value):

        groups = self.groups
        for group in pbs.groups:
            group.dashes = value 
