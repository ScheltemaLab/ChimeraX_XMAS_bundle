# Copyright 2022 Scheltema LAB
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import pandas as pd

# The name and content of the first (reference) column depends on the search
# engine used. This dictionary contains the reference column names.
ref_columns = {"Proteome Discoverer": "Row in evidence file",
               "pLink": "Peptide_Order",
              "Xi": "PeptidePairID",
              "Xi_alternative": "PSMID",
              "mzIdentML": "Peptide id"}


# Each mapping information file is created and maintainted with an instance if 
# this class.
class InfoFile:
    
    
    def __init__(self, path, engine):
        # Create a dataframe for the evidence file
        self.ref_column = ref_columns[engine]
        self.columns = [self.ref_column, "Pseudobond", "Overlap category", 
                        "Distance (A)"]
        self.df = pd.DataFrame(columns=self.columns)
        self.path = path
    
    
    def add(self, row_number, value, category="", distance=""):
        # Add a row to the dataframe
        data = [[row_number, value, category, distance]]
        df_add = pd.DataFrame(data, columns=self.columns)
        self.df = pd.concat([self.df, df_add], ignore_index=True)
        
        return len(self.df.index) - 1
        
    
    def create_file(self):
        # Create the tsv file from the dataframe
        self.df.sort_values([self.ref_column, "Pseudobond"], inplace=True)
        self.df.to_csv(self.path, sep="\t", index=False)