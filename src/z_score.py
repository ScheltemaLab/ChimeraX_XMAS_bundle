from PyQt5.QtWidgets import QLabel, QLineEdit, QPushButton, QDialogButtonBox, QHBoxLayout, QVBoxLayout, QFileDialog, QSizePolicy
from PyQt5.QtCore import Qt
from pathlib import Path
from chimerax.core.commands import run
from chimerax.core.models import MODEL_POSITION_CHANGED
import re
import os
from .tool import Slider, ExportSlider

class ZScoreSelector:


    def __init__(self, session, tool_window):

        self.session = session
        
        self.main_dialog = tool_window.create_child_window("Create HADDOCK input from DisVis output")

        label = QLabel("Select a DisVis output folder:")
        line_edit = QLineEdit()
    
        line_edit.setMinimumWidth(200)
        folder_button = QPushButton("Select folder")
        folder_button.clicked.connect(lambda: self.file_dialog(line_edit))
        use_button = QPushButton("Use folder")
        use_button.setSizePolicy(QSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred))
        use_button.clicked.connect(lambda: self.ok_clicked(line_edit))
        ok_cancel = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        ok_cancel.accepted.connect(self.create_haddock_output)
        ok_cancel.rejected.connect(self.main_dialog.destroy)

        file_layout = QHBoxLayout()
        file_layout.addWidget(line_edit)
        file_layout.addWidget(folder_button)
        
        main_layout = self.main_dialog.layout = QVBoxLayout()
        main_layout.addWidget(label)
        main_layout.addLayout(file_layout)
        main_layout.addWidget(use_button)
        main_layout.setAlignment(use_button, Qt.AlignRight)
        main_layout.addWidget(ok_cancel)
        
        self.main_dialog.ui_area.setLayout(main_layout)
        self.main_dialog.cleanup = self.cleanup
        self.main_dialog.manage("side")
        
        # Start for testing
        # line_edit.setText("C:/Users/ilsel/OneDrive/Documenten/MCLS/Bioinformatics_profile/zscore/disvis_out")
        # self.ok_clicked(line_edit)
        # return
        # End for testing
        
        
    def cleanup(self):
        
        if not hasattr(self, "pbs"):
            return
        
        self.pbs.displays = True
        self.pbs.colors = self.colors["Main"]
        self.triggerset.remove_handler(self.movement_handler)
        
        
    def create_haddock_output(self):
        
        if not hasattr(self, "pbs"):
            return
        
        setters = self.setters
        values = [None] * len(setters)
        for i, setter in enumerate(setters):
            values[i] = float(setter.text())
            
        chosen_restraints = []
        
        for pb in self.pbs:
            if (pb.zscore < values[0] or pb.zscore > values[1]):
                continue
            chosen_restraints.append(pb.restraint_number)
        
        chosen_restraints.sort()
        for restraint in chosen_restraints:
            print(restraint)
            

    def file_dialog(self, line_edit):

        folder = QFileDialog.getExistingDirectory(None, 
                                                  "Select a DisVis output folder", 
                                                  "", QFileDialog.ShowDirsOnly
                                                  | QFileDialog.DontResolveSymlinks)
        line_edit.setText(folder)


    def ok_clicked(self, line_edit):
        
        self.folder = line_edit.text() + "/"
        
        if self.folder == "/":
            return
        
        number_of_pbs, disvis_chains = self.read_log()
                          
        extensions = [disvis_chains["fix"], disvis_chains["scan"], ".mrc"]
        for extension in extensions:
            for file in self.get_files(extension):
                command = "open %s" % file
                run(self.session, command)
                
        disvis_chains = self.get_chains(disvis_chains)
                
        zscore_path = self.folder + "z-score.out"
        zscores = self.get_zscores(zscore_path, "\S+(?=\n)")
        
        input_path = list(self.get_files(".txt"))[0]
        self.pbs = self.make_pbs(input_path, disvis_chains, zscores)
        
        self.get_input_distances(input_path)
        self.colors = {"Main": [255, 255, 0, 255], "Cutoff": [255, 0, 0, 255]}
        self.color_pbs()
        self.trigger_handler(disvis_chains.values())
        self.prepare_slider(zscores)
        
        
    def prepare_slider(self, zscores):
        
        minimum = min(zscores)
        maximum = max(zscores)
        slider = ZScoreSlider("zscore", True, minimum, maximum, self.pbs)
        main_layout = self.main_dialog.layout
        rows = main_layout.count()
        main_layout.insertLayout(rows - 1, slider.layout)
        
        self.setters = slider.setters
        
        
    def get_chains(self, disvis_chains):
        
        for model in self.session.models:
            if not model.name in disvis_chains.values():
                continue
            index = list(disvis_chains.values()).index(model.name)
            key = list(disvis_chains.keys())[index]
            disvis_chains[key] = model
            
        return disvis_chains
        
    
    def make_pbs(self, path, disvis_chains, zscores):      
            
        file = open(path, "r") 
        pb_manager = self.session.pb_manager
        name = os.path.basename(path).replace(".txt", "")
        group = pb_manager.get_group(name)
        group.radius = 0.5
        group.color = [255, 255, 0, 255]
        for i, line in enumerate(file):
            split_line = line.split(" ")
            atoms = [None] * len(disvis_chains)
            j = 0
            for k, chain in enumerate(disvis_chains):
                chain_id = split_line[j]
                res_number = int(split_line[j + 1])
                atom_name = split_line[j + 2]
                model = disvis_chains[chain]
                residue = model.find_residue(chain_id, res_number)
                for atom in residue.atoms:
                    if atom.name != atom_name:
                        continue
                    atoms[k] = atom
                j += 3
            pb = group.new_pseudobond(atoms[0], atoms[1])
            pb.zscore = zscores[i]
            pb.restraint_number = i + 1
        self.session.models.add([group])
        
        file.close()
        return group.pseudobonds
        
    
    def read_log(self):      
        
        log_path = self.folder + "disvis.log"
        log_file = open(log_path, "r")
        keys = "fix", "scan"
        disvis_chains = {}
        for line in log_file:
            line_lower = line.lower()
            if not line_lower.count("pdb") > 0:
                if line_lower.count("distance") > 0:
                    number_of_pbs = self.find_string("\d+(?=\n)", line)
            for key in keys:
                if (line_lower.count(key)) > 0:
                    disvis_chains[key] = self.find_string("\w+\.pdb", line)
                    break
                
        log_file.close()
        return number_of_pbs, disvis_chains
        
        
    def trigger_handler(self, chains):
        
        self.triggerset = self.session.triggers
        self.movement_handler = self.triggerset.add_handler(
            MODEL_POSITION_CHANGED,
            lambda trigger, trigger_data, c=chains: 
                self.handle_movement(trigger, trigger_data, c))
        
        
    def handle_movement(self, trigger, trigger_data, chains):
        
        # Maybe it only works if fixed or scanning are moved?
        if not trigger_data in chains:
            return
        
        self.color_pbs()
        
        
    def get_input_distances(self, path):
        
        file = open(path, "r")
        for line in file:
            first_line = line
            break
        
        file.close()
        
        self.minimum = float(self.find_string("\S+(?=\s+\S+\n)", first_line))
        self.maximum = float(self.find_string("\S+(?=\n)", first_line))
        
        
    def color_pbs(self):
        
        for pb in self.pbs:
            if pb.length >= self.minimum and pb.length <= self.maximum:
                color = self.colors["Main"]
            else:
                color = self.colors["Cutoff"]
            pb.color = color
        
        
    def find_string(self, search_string, whole_string):
        
        return re.search(search_string, whole_string).group(0)
            
    
    def get_zscores(self, path, search_string):
        
        # Doesn't work if you use file instead of open(path, "r)
        zscores = [None] * sum(1 for line in open(path, "r"))
        file = open(path, "r")
        for i, line in enumerate(file):
            zscores[i] = float(self.find_string(search_string, line))
        
        file.close()            
        return zscores    
    
        
    def get_files(self, extension):
        
        files = Path(self.folder).glob("*" + extension)
        
        return files
    
    
class ZScoreSlider(Slider):
    
    
    def __init__(self, value_type="zscore", enabled=True, minimum=None, maximum=None, pbs=None):
        
        self.add = -minimum
        super().__init__(value_type, enabled, minimum, maximum, pbs)
        self.function = ExportSlider.display_pseudobonds
        
    
    def within_range(self, value_type, pbs, function=None):

        score_range = self.slider.value()
        real_range = self.get_real_values(score_range)
        within_range = {}
        for pb in pbs:
            zscore = pb.zscore
            is_outside_range = (zscore < real_range[0] 
                                or zscore > real_range[1])
            within_range[pb] = is_outside_range
            
        self.function(self, within_range)
        
    
    def get_real_values(self, values):
        
        real_values = [None] * len(values)
        
        for i, value in enumerate(values):
            real_value = value / 1000 - self.add
            real_values[i] = round(real_value, 3)
            
        return real_values
    
        
    def get_slider_values(self, minimum, maximum):

        minimum = (minimum + self.add) * 1000
        maximum = (maximum + self.add) * 1000
        
        return minimum, maximum