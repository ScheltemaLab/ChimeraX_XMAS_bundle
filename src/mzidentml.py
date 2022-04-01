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


from .read_evidence import PeptidePair
import re
import xml.dom.minidom as minidom


class XlPeptide(PeptidePair):
    
    def __init__(self):

        super().__init__()
        
    
    @staticmethod
    def merge_XL_DON_RCV_peptide(don, rcv, key):
        don.Id = key
        don.SequenceB = rcv.SequenceB
        don.ModificationsB = rcv.ModificationsB
        don.XLinkPositionB = rcv.XLinkPositionB
        don.AccessionB = rcv.AccessionB
        don.PositionB = rcv.PositionB
        don.ProteinDescriptionsB = rcv.ProteinDescriptionsB
        if (rcv.Score > don.Score):
            don.Score = rcv.Score
        if (rcv.NumCSMs > don.NumCSMs):
            don.NumCSMs = rcv.NumCSMs
        return don
    

    def as_dictionary(self):
        return {"Ref" : self.Ref,
            "Score" : self.Score,
            "CrossLinker" : self.CrossLinker,
            "CrosslinkType" : self.CrosslinkType,
            "SequenceA" : self.SequenceA,
            "ModificationsA" : self.ModificationsA,
            "AccessionA" : self.AccessionA,
            "XLinkPositionA" : self.XLinkPositionA,
            "PositionA" : self.PositionA,
            "SequenceB" : self.SequenceB,
            "ModificationsB" : self.ModificationsB,
            "AccessionB" : self.AccessionB,
            "XLinkPositionB" : self.XLinkPositionB,
            "PositionB" : self.PositionB,
            "ProteinDescriptionsA" : self.ProteinDescriptionsA,
            "ProteinDescriptionsB" : self.ProteinDescriptionsB,
            "IsDecoy" : self.IsDecoy,
            "QValue" : self.QValue,
            "NumCSMs" : self.NumCSMs,
            "AlignmentsA": self.AlignmentsA,
            "AlignmentsB": self.AlignmentsB}
    

class DBSequence:
    def __init__(self):
       self.Accession = ""
       self.Id = ""
       self.ProteinDescription = ""
       

class XlPeptideEvidence:
    def __init__(self):
        self.PeptideRef = ""
        self.DBSequenceRef = ""
        self.Id = ""
        self.End = 0
        self.Start = 0
        

class XlSpectrumIdentification:
    def __init__(self):
        self.PeptideRef = ""
        self.Score = ""
        self.SpectraCount = 1
        

# For some reason mzIdentML uses their unique sequence id insead of input 
# accession id
# This function returns dictionary with ID as key, from which we can find 
# protein accession and description if needed
def parse_sequence_dbs(mzidentfile):
    doc = minidom.parse(mzidentfile)
    db_sequences = {}
    protein_sequence_data = doc.getElementsByTagName("DBSequence")
    for sequence in protein_sequence_data:
        sq = DBSequence()
        sq.Accession = sequence.getAttribute("accession")
        sq.Id = sequence.getAttribute("id")
        cv_param = sequence.getElementsByTagName("cvParam")
        sq.ProteinDescription = cv_param[0].getAttribute("value")
        db_sequences[sq.Id] = sq
    return db_sequences


# Only real purpose to read this data blob is to find to which protein the 
# peptide belongs and its location in sequence
def parse_xl_peptides_evidence(mzidentfile):
    doc = minidom.parse(mzidentfile)
    evidences = {}
    peptide_evidence = doc.getElementsByTagName("PeptideEvidence")
    for evidence in peptide_evidence:
        if ("XL" in evidence.getAttribute("peptide_ref")):
            ev = XlPeptideEvidence()
            ev.DBSequenceRef  = evidence.getAttribute("dBSequence_ref")
            ev.Id = evidence.getAttribute("id")
            ev.PeptideRef= evidence.getAttribute("peptide_ref")
            ev.Start = evidence.getAttribute("start")
            ev.End = evidence.getAttribute("end")
            evidences[ev.PeptideRef] = ev
    return evidences


# Store the best crosslink score and CSM count from spectra hits
def parse_spectrum_identification_result(mzidentfile):

    # Helper function to extract info for both x-link peptides. DON has index 0 
    # and RCV index 1.
    def extract_info (specindex, best_spectra_match):
        xlspecid = XlSpectrumIdentification()
        xlspecid.PeptideRef = spec_id[specindex].getAttribute("peptide_ref")
        cvparam = spec_id[specindex].getElementsByTagName("cvParam")
        for cv in cvparam:
            if (cv.getAttribute("name") == "xi:score"):
                xlspecid.Score = cv.getAttribute("value")
                xlspecid
                # Not sure if this is needed, because it looks like at least PD 
                # exports only the best hit
                if xlspecid.PeptideRef in best_spectra_match:
                    if (xlspecid.Score > best_spectra_match[xlspecid.PeptideRef].Score):
                        xlspecid.SpectraCount = xlspecid.SpectraCount + 1
                        best_spectra_match[xlspecid.PeptideRef] = xlspecid
                else:
                    best_spectra_match[xlspecid.PeptideRef] = xlspecid
        return best_spectra_match

    doc = minidom.parse(mzidentfile)
    best_spectra_match = {}
    spectra_matched = doc.getElementsByTagName("SpectrumIdentificationResult")
    for spectramatch in spectra_matched:
        spec_id = spectramatch.getElementsByTagName("SpectrumIdentificationItem")
        
        pep_ref = spec_id[0].getAttribute("peptide_ref")
        if ("XL_" in pep_ref):
            extract_info (0, best_spectra_match)
            extract_info (1, best_spectra_match)

    return best_spectra_match


def parse_xl_peptides(mzidentfile):
    doc = minidom.parse(mzidentfile)
    dbs = parse_sequence_dbs(mzidentfile)
    evidences = parse_xl_peptides_evidence(mzidentfile)
    spectra_scores = parse_spectrum_identification_result(mzidentfile)

    peptides = doc.getElementsByTagName("Peptide")

    don_peptides = {}
    rcv_peptides = {}

    for peptide in peptides:
        if ("XL_" in peptide.getAttribute("id")):
            pep = XlPeptide()
            
            key = peptide.getAttribute("id")
            peptidesequence = peptide.getElementsByTagName("PeptideSequence")[0].firstChild.nodeValue
            
            #Mod name [location]
            modifications = []
            pepmods = peptide.getElementsByTagName("Modification")
            for mod in pepmods:
                loc = mod.getAttribute("location")
                cvparams = mod.getElementsByTagName("cvParam")
                name = cvparams[0].getAttribute("name")
                modifications.append("".join((name," [", loc, "]")))
                if (cvparams[0].getAttribute("cvRef") == "XLMOD"):
                    xlinkposition = int(loc)
            
            if (key in evidences):
                evidence = evidences[key]
                accessions = dbs[evidence.DBSequenceRef].Accession
                position = evidence.Start
                proteindescription = dbs[evidence.DBSequenceRef].ProteinDescription
            else:
                accessions = ""
                position = ""
                proteindescription = ""
            
            if (key in spectra_scores):
                pep.Score = float(spectra_scores[key].Score)
                pep.NumCSMs = int(spectra_scores[key].SpectraCount)

            if ("_DON" in key):
                key = re.sub("_DON", "", key)
                pep.Ref = key
                pep.AccessionA = accessions
                pep.ModificationsA = ";".join(modifications)
                pep.PositionA = position
                pep.ProteinDescriptionsA = proteindescription
                pep.SequenceA = peptidesequence
                pep.XLinkPositionA = xlinkposition - 1
                don_peptides[key] = pep
            else:
                key = re.sub("_RCV", "", key)
                pep.Ref = key
                pep.AccessionB = accessions
                pep.ModificationsB = ";".join(modifications)
                pep.PositionB = position
                pep.ProteinDescriptionsB = proteindescription
                pep.SequenceB = peptidesequence
                pep.XLinkPositionB = xlinkposition - 1
                rcv_peptides[key] = pep

    xlinks = []
    for key in don_peptides:
        xlinks.append(XlPeptide.merge_XL_DON_RCV_peptide(don_peptides[key], 
                                                         rcv_peptides[key], 
                                                         key))
       
    return xlinks