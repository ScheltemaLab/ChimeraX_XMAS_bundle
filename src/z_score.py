from PyQt5.QtWidgets import QDialog, QLabel, QLineEdit, QPushButton, QDialogButtonBox, QHBoxLayout, QVBoxLayout, QFileDialog
from pathlib import Path
from chimerax.core.commands import run
from chimerax.core.models import MODEL_POSITION_CHANGED
import re
import os

class ZScoreSelector:


    def __init__(self, session):

        self.session = session
        
        self.main_dialog = QDialog()
        self.main_dialog.setWindowTitle("Create HADDOCK input from DisVis output")

        label = QLabel("Select a DisVis output folder:")
        line_edit = QLineEdit()
        
        # Start for testing
        line_edit.setText("C:/Users/Nathan/Documents/zscore/disvis_out")
        self.ok_clicked(line_edit)
        return
    # End for testing
    
        line_edit.setMinimumWidth(200)
        folder_button = QPushButton("Select folder")
        folder_button.clicked.connect(lambda: self.file_dialog(line_edit))
        ok_cancel = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        ok_cancel.accepted.connect(lambda: self.ok_clicked(line_edit))
        ok_cancel.rejected.connect(self.main_dialog.close)

        file_layout = QHBoxLayout()
        file_layout.addWidget(line_edit)
        file_layout.addWidget(folder_button)
        
        main_layout = QVBoxLayout()
        main_layout.addWidget(label)
        main_layout.addLayout(file_layout)
        main_layout.addWidget(ok_cancel)

        self.main_dialog.setLayout(main_layout)
        self.main_dialog.show()
        

    def file_dialog(self, line_edit):

        folder = QFileDialog.getExistingDirectory(self.main_dialog, 
                                                  "Select a DisVis output folder", 
                                                  "", QFileDialog.ShowDirsOnly
                                                  | QFileDialog.DontResolveSymlinks)
        line_edit.setText(folder)


    def ok_clicked(self, line_edit):
        
        self.folder = line_edit.text() + "/"
        
        if self.folder == "/":
            return
        
        extensions = ["*.pdb", "*.pb", "*.mrc"]
        for extension in extensions:
            for file in self.get_files(extension):
                command = "open %s" % file
                run(self.session, command)
                
        log_path = self.folder + "disvis.log"
        log_file = open(log_path, "r")
        keys = "fix", "scan"
        chains = {}
        for line in log_file:
            line_lower = line.lower()
            if not line_lower.count("pdb") > 0:
                continue
            for key in keys:
                if (line_lower.count(key)) > 0:
                    chains[key] = self.get_name(line, "pdb")
                    break

        pb_path = os.path.normpath(list(self.get_files("*.pb"))[0])
        pb_model = self.get_name(pb_path, "pb")      
        for model in self.session.models:
            if model.name == pb_model:
                pbs = model.pseudobonds
                self.pb_model = model
        number_of_pbs = len(pbs)
        
        input_path = list(self.get_files("*.txt"))[0]
        input_lines = self.list_from_file(input_path,
                                          "\S+\s+" * 5 + "\S+",
                                          number_of_pbs)
        
        self.get_input_distances(input_path, pbs)
        self.trigger_handler(chains.values(), pbs)
        
        zscore_path = self.folder + "z-score.out"
        zscores = self.list_from_file(zscore_path, "\S+(?=\n)", number_of_pbs)
        
        self.find_zscores(pbs, chains, input_lines, zscores)
        
        
    def trigger_handler(self, chain_names, pbs):
        
        triggerset = self.session.triggers
        self.movement_handler = triggerset.add_handler(
            MODEL_POSITION_CHANGED,
            lambda trigger, trigger_data, cn=chain_names, p=pbs: 
                self.handle_movement(trigger, trigger_data, cn, p))
        
        
    def handle_movement(self, trigger, trigger_data, chain_names, pbs):
        
        # Maybe it only works if fixed or scanning are moved?
        if not trigger_data.name in chain_names:
            return
        
        self.color_pbs(pbs)
        
        
    def get_input_distances(self, path, pbs):
        
        for line in open(path, "r"):
            first_line = line
            break
        
        self.minimum = float(self.find_string("\S+(?=\s+\S+\n)", first_line))
        self.maximum = float(self.find_string("\S+(?=\n)", first_line))
        
        self.color_pbs(pbs)
        
        
    def color_pbs(self, pbs):
        
        for pb in pbs:
            if pb.length >= self.minimum and pb.length <= self.maximum:
                color = [255, 255, 0, 255]
            else:
                color = [170, 0, 0, 255]
            pb.color = color
        
        
    def find_string(self, search_string, whole_string):
        
        return re.search(search_string, whole_string).group(0)
    
        
    def find_zscores(self, pbs, chains, input_lines, zscores):
        
        for pb in pbs:
            strings = [None] * len(pb.atoms)
            for atom in pb.atoms:
                model = atom.structure.name
                if model == chains["fix"]:
                    i = 0
                elif model == chains["scan"]:
                    i = 1
                residue = atom.residue
                chain_id = residue.chain_id
                number = str(residue.number)
                name = atom.name
                string = chain_id + " " + number + " " + name
                strings[i] = string
            pb_line = " ".join(strings)
            j = input_lines.index(pb_line)
            pb.zscore = float(zscores[j])
            
    
    def list_from_file(self, path, search_string, number):
                
        file = open(path, "r")
        lst = [None] * number
        i = 0
        for line in file:
            item = self.find_string(search_string, line)
            lst[i] = item
            i += 1
            
        return lst
        
        
    def get_name(self, line, extension):
        
        search_string = "\w+\.%s" % extension
        
        return self.find_string(search_string, line)      
    
        
    def get_files(self, extension):
        
        files = Path(self.folder).glob(extension)
        
        return files
