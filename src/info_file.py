# vim: set expandtab shiftwidth=4 softtabstop=4:
# === UCSF ChimeraX Copyright ===
# Copyright 2016 Regents of the University of California.
# All rights reserved. This software provided pursuant to a license 
# agreement containing restrictions on its disclosure, duplication and 
# use. For details see:
# http://www.rbvi.ucsf.edu/chimerax/docs/licensing.html
# This notice must be embedded in or attached to all copies, including 
# partial copies, of the software or any revisions or derivations 
# thereof.
# === UCSF ChimeraX Copyright ===

import pandas as pd

class InfoFile:
    
    
    def __init__(self, path):
        
        self.columns = ["Row in evidence file", "Pseudobond", "Overlap category", "Distance"]
        self.df = pd.DataFrame(columns=self.columns)
        self.path = path
    
    
    def add(self, row_number, value, category="", distance=""):
        
        data = [[row_number, value, category, distance]]
        df_add = pd.DataFrame(data, columns=self.columns)
        self.df = self.df.append(df_add, ignore_index=True)
        
        return len(self.df.index) - 1
        
    
    def create_file(self):
               
        self.df.sort_values(["Row in evidence file", "Pseudobond"], inplace=True)
        self.df.to_csv(self.path, sep="\t", index=False)
        

