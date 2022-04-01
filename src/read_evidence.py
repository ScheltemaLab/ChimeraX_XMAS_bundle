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


from csv import Sniffer
from operator import attrgetter
import pandas as pd

double_attributes = {"Sequences": ["SequenceA", "SequenceB"],
                  "Positions": ["XLinkPositionA", "XLinkPositionB"]}


# Enable evidence files from multiple results files and input formats to be
# read.
class Evidence:

    
    def __init__(self, evidence_file):
        
        # Use a specific class to parse mzIdentML files
        if evidence_file.endswith(".mzid"):
            cls = MzIdentML
        else:
            cls = Tabular
        
        evidence = cls(evidence_file)
        self.peptide_pairs = evidence.peptide_pairs
        self.engine = evidence.engine
        
        
# Read evidence from mzIdentML files            
class MzIdentML:    

    
    def __init__(self, evidence_file):
        
        if not "parse_xl_peptides" in dir():       
            from .mzidentml import parse_xl_peptides
        
        self.peptide_pairs = parse_xl_peptides(evidence_file)
        self.engine = "mzIdentML"
        for peptide_pair in self.peptide_pairs:
            sort_peptides(peptide_pair)
                
    
# Read evidence from files with tabular format (*.csv, *.tsv, *.txt, *.xls
# *.xlsx)        
class Tabular:

    
    def __init__(self, evidence_file):
        
        import os
        
        # Dictionary dictating, for each search engine, which column names to be 
        # used from the header (index=0), which function to be used to obtain 
        # the relevant information (index=1), and which function to further 
        # parse this information (index=2)
        # Dictionary is slightly different for pLink, because this engine 
        # requires a different approach for parsing
        self.engines = {"Proteome Discoverer": [["Sequence A", "Sequence B", 
                                                 "Max. XlinkX Score", 
                                                 "Is Decoy"],
                                                self.parse_pd_xi_seqs_scores, 
                                                self.parse_pd_pos_ids],
                        "pLink": [["Peptide"], self.parse_plink],
                        "Xi": [["Peptide1", "Peptide2", "Score", "IsDecoy"],                                
                               self.parse_pd_xi_seqs_scores,                               
                               self.parse_xi_pos_ids],
                        # Two different layouts are supported for Xi
                        "Xi_alternative": [["PepSeq1", "PepSeq2", 
                                            "Score", "IsDecoy"],
                                           self.parse_pd_xi_seqs_scores, 
                                           self.parse_xi_pos_ids]}
        
        # Get header and delimiter
        extension = os.path.splitext(evidence_file)[1][1:]
        self.is_excel = (extension == "xls" or extension == "xlsx") 
        header, delimiter = self.parse_headers(evidence_file)
        
        # Get engine with header
        engine = self.parse_engine(header)
        
        if engine == "":
            print("Unsupported evidence file format")
            return
        
        self.engine = engine
        
        # Make the dataframe
        df = self.make_df(evidence_file, delimiter, header)
        
        # Obtain the relevant information by calling the search engine-specific
        # function
        self.engines[engine][1](df)
        
        
    def parse_pd_xi_seqs_scores(self, df):
        
        # PD and Xi evidence files have a similar structure. This method parses
        # the peptide sequences and scores of both of them.

        xi_alternative = self.engine == "Xi_alternative" 
        parameters = self.engines[self.engine]
        col_names, function = parameters[0], parameters[2]
             
        attributes = (double_attributes["Sequences"] + ["Score", "IsDecoy"])
        range_xlinks = range(len(df.index))
        # Each peptide pair is stored as an instance of class PeptidePair
        peptide_pairs = [PeptidePair() for i in range_xlinks]
            
        for i in range(len(col_names)):
            try:
                values = df[col_names[i]].tolist()
            except:
                continue
            for j in range_xlinks:
                peptide_pair = peptide_pairs[j]
                setattr(peptide_pair, attributes[i], values[j])
            
        function(peptide_pairs, df, xi_alternative)
        self.peptide_pairs = peptide_pairs
        
        
    def parse_pd_pos_ids(self, peptide_pairs, *args):
        
        # Parse the crosslink positions and peptide pair references of PD
        # evidence files
        
        for i, peptide_pair in enumerate(peptide_pairs):
            peptide_pair.Ref = i + 2
            for j, seq_attr in enumerate(double_attributes["Sequences"]):
                seq = getattr(peptide_pair, seq_attr)
                pos_attr = double_attributes["Positions"][j]
                if peptide_pair.invalid(seq_attr):
                    pos = ""
                    seq_new = ""
                elif seq.count("[") == 0: 
                    pos = 0
                    seq_new = seq
                else:
                    pos = seq.index("[")
                    seq_new = seq.replace("[", "").replace("]", "")
                setattr(peptide_pair, pos_attr, pos)
                setattr(peptide_pair, seq_attr, seq_new)
            sort_peptides(peptide_pair)
            
            
    def parse_plink(self, df):
        
        # Parse pLink evidence files
        
        from numpy import isnan
        import re
        
        # pLink evidence files have a two-dimensional header. The upper header
        # corresponds to unique peptide pairs. The lower header corresponds to 
        # spectra belonging to a peptide pair. The lower (sub-) header is read 
        # and, subsequently, removed from the dataframe for convenience
        subheader = df.iloc[0,:].values.tolist()
        df.drop(index=df.index[0], axis=0, inplace=True)
        peptide_pairs = []  
        
        # Helper function to extract sequences and crosslink positions
        def parse_seq_pos(value):
            sequences = re.findall("[A-Z]+", value)
            positions = re.findall("[0-9]+", value)
            return sequences, positions

        for row in df.itertuples():
            ref = row.Peptide_Order
            # If the Peptide_Order cell contains a value, the row contains
            # peptide pair data. If not, it contains spectra data
            if not isnan(ref):
                peptide_pair = PeptidePair()
                peptide_pairs.append(peptide_pair)
                peptide_pair.Ref = int(ref)
                sequences, positions = parse_seq_pos(row.Peptide)
                for i, seq_attr in enumerate(double_attributes["Sequences"]):
                    setattr(peptide_pair, seq_attr, sequences[i])
                    setattr(peptide_pair, double_attributes["Positions"][i], 
                            int(positions[i]) - 1)
                sort_peptides(peptide_pair)
                continue
            data = list(row)[1:]
            sub_pair = dict(zip(subheader, data))
            score = float(sub_pair["Score"])
            if (peptide_pair.Score == "" or peptide_pair.Score < score):
                peptide_pair.Score = score
                
            self.peptide_pairs = peptide_pairs
        
            
    def parse_xi_pos_ids(self, peptide_pairs, df, alternative):
        
        # Parse crosslink positions and peptide references for Xi evidence 
        # files
        
        # Dictionary containing the relevant column names for Xi (False) and
        # Xi_alternative (True)
        col_names = {False: ["FromSite", "ToSite", "PeptidePairID"],
                     True: ["LinkPos1", "LinkPos2", "PSMID"]}
        keys = double_attributes["Positions"] + ["Ref"]
        
        params = {}
        for i, key in enumerate(keys):
            col_name = col_names[alternative][i]
            params[key] = df[col_name].tolist()
            
        for i, peptide_pair in enumerate(peptide_pairs):
            peptide_pair.Ref = params["Ref"][i]
            for j, key in enumerate(keys[:2]):
                seq_attr = double_attributes["Sequences"][j]
                pos_attr = double_attributes["Positions"][j]
                if peptide_pair.invalid(key):
                    seq = ""
                    pos = ""
                else:
                    seq = getattr(peptide_pair, seq_attr)
                    residues = [res for res in seq if res.isupper()]
                    seq = "".join(residues)
                    pos = params[key][i] - 1
                setattr(peptide_pair, seq_attr, seq)
                setattr(peptide_pair, pos_attr, pos)
                
        
    def parse_headers(self, evidence_file):
        
        if self.is_excel:
            df = pd.read_excel(evidence_file, nrows=0)
            header = df.columns.values.tolist()
            return header, None
        
        with open(evidence_file) as f:
            line = f.readline().rstrip()
        dialect = Sniffer().sniff(line)
        delimiter = dialect.delimiter
        header = line.split(delimiter)
            
        return header, delimiter
                
                
    def parse_engine(self, header):
        
        engine = ""
        for i, eng in enumerate(self.engines):
            col_name = self.engines[eng][0][0]
            if col_name not in header:
                continue
            engine = eng
            break
        
        return engine
                
                
    def make_df(self, evidence_file, delimiter, header):
        
        # Write evidence file to a Pandas dataframe. Approach differs between
        # Excel (xlsx) files, non-Excel non-pLink files, and non-Excel pLink
        # files
        if self.is_excel:
            df = pd.read_excel(evidence_file)
            return df
        
        if self.engine != "pLink":        
            df = pd.read_csv(evidence_file, sep=delimiter, 
                             decimal=".")
            score_col = self.engines[self.engine][0][2]
            if (df[score_col].dtypes != "float64" and delimiter != ","):
                df = pd.read_csv(evidence_file, sep=delimiter, decimal=",")
            return df
        
        # For non-Excel pLink files, header might need to be adjusted, since  
        # the first row (the sub-header) can contain more cells than the 
        # header, which causes problems when converting to dataframe.
        with open(evidence_file) as f:
            f.readline()
            subheader = f.readline().rstrip().split(delimiter)
        missing_cells = len(subheader) - len(header)
        if missing_cells > 0:
            addition = ["Unnamed: %s" % i for i in range(missing_cells)]
            col_names = header + addition
        else:
            col_names = header   
        df = pd.read_csv(evidence_file, sep=delimiter, header=0, 
                         names=col_names, decimal=".")
            
        return df
    
    
class PeptidePair:
    
    
    def __init__(self):
        
        self.Ref = ""
        self.Score = ""
        self.CrossLinker = ""
        self.CrosslinkType = ""
        self.SequenceA = ""
        self.ModificationsA = ""
        self.AccessionA = ""
        self.XLinkPositionA = ""
        self.PositionA = ""
        self.SequenceB = ""
        self.ModificationsB = ""
        self.XLinkPositionB = ""
        self.AccessionB = ""
        self.PositionB = ""
        self.ProteinDescriptionsA = ""
        self.ProteinDescriptionsB = ""
        self.IsDecoy = False
        self.QValue = 0
        self.NumCSMs = 1
        self.AlignmentsA = []
        self.AlignmentsB = []
        
        
    def get_info(self, attr="Sequence"):
        
        # Return multiple attributes
        f = attrgetter(attr + "A", attr + "B")
        
        return f(self)  
    
    
    def invalid(self, seq_attr):
        
        # If sequence is from a decoy peptide or is not a string, it cannot be
        # used for mapping(/alignment)
        seq = getattr(self, seq_attr)
        return (self.IsDecoy or not isinstance(seq, str))
            
                    
def sort_peptides(peptide_pair):
    
    # Sort peptides of a peptide pair based on sequence
    
    seq_attributes = double_attributes["Sequences"]
    pos_attributes = double_attributes["Positions"]
    
    peptides_unsorted = [None] * 2
    pos_unsorted = [None] * 2
    
    for i, seq_attr in enumerate(seq_attributes):
        peptides_unsorted[i] = getattr(peptide_pair, seq_attr)
        pos_unsorted[i] = getattr(peptide_pair, pos_attributes[i])
    peptides_sorted = sorted(peptides_unsorted)
    
    if peptides_sorted == peptides_unsorted:
        return
    
    pos_sorted = pos_unsorted[::-1]
    
    for i, seq_attr in enumerate(seq_attributes):
        setattr(peptide_pair, seq_attr, peptides_sorted[i])
        setattr(peptide_pair, pos_attributes[i], pos_sorted[i])      
        
        
