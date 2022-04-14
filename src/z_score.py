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


from chimerax.core.commands import run
from chimerax.core.models import MODEL_POSITION_CHANGED
import os
from pathlib import Path
from Qt.QtCore import Qt
from Qt.QtWidgets import (QLabel, QLineEdit, QPushButton, QDialogButtonBox,
                             QHBoxLayout, QVBoxLayout, QFileDialog,
                             QSizePolicy)
import re


class ZScoreSelector:

    # A window is opened in which the user can select a DisVis output folder.
    # From this folder, the .pdb and .mrc files are opened, and the pseudobonds
    # are drawn on the structures. The user-specified distance restraints are a
    # cut-off for pseudobond color. Colors are adjusted upon model movement. A
    # slider allows selection of pseudobonds based on z-score. A .pb file is
    # generated with the selected pseudobonds.

    def __init__(self, xmas):
        
        # Create a child window of the main tool window
        self.xmas = xmas
        self.session = xmas.session
        tool_window = xmas.tool_window        
        self.main_dialog = tool_window.create_child_window("Create HADDOCK "
                                                           "input from DisVis "
                                                           "output")

        label = QLabel("Select a DisVis output folder:")
        line_edit = QLineEdit()
    
        line_edit.setMinimumWidth(200)
        folder_button = QPushButton("Select folder")
        folder_button.clicked.connect(lambda: self.file_dialog(line_edit))
        use_button = QPushButton("Use folder")
        use_button.setSizePolicy(QSizePolicy(QSizePolicy.Maximum, 
                                             QSizePolicy.Preferred))
        use_button.clicked.connect(lambda: self.use_folder_clicked(line_edit))
        ok_cancel = QDialogButtonBox(QDialogButtonBox.Ok 
                                     | QDialogButtonBox.Cancel)
        ok_cancel.accepted.connect(self.create_pb_file)
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
            

    def file_dialog(self, line_edit):
        
        # Get the DisVis output folder
        
        get_directory = QFileDialog.getExistingDirectory
        folder = get_directory(None, "Select a DisVis output folder", "",
                               QFileDialog.ShowDirsOnly 
                               | QFileDialog.DontResolveSymlinks)
        line_edit.setText(folder)


    def use_folder_clicked(self, line_edit):
        
        # Open the .pdb and .mrc files, and create the pseudobonds with 
        # corresponding z-scores
        
        self.folder = line_edit.text() + "/"
        
        if self.folder == "/":
            return
        
        disvis_chains = self.read_log()
                          
        extensions = [disvis_chains["fix"], disvis_chains["scan"], ".mrc"]
        for extension in extensions:
            for file in self.get_files(extension):
                command = "open %s" % file
                run(self.session, command)
                
        disvis_chains = self.get_chains(disvis_chains)
               
        zscores = self.get_zscores()
        
        input_path = list(self.get_files(".txt"))[0]
        self.pbs = self.make_pbs(input_path, disvis_chains, zscores)
        
        # Obtain input distances and color pseudobonds accordingly
        self.get_input_distances(input_path)
        self.colors = {"Main": [255, 255, 0, 255], "Cutoff": [255, 0, 0, 255]}
        self.color_pbs()
        # Re-color when a chain is moved
        self.trigger_handler(disvis_chains.values())
        
        # Make the z-score slider
        self.prepare_slider(zscores) 
        
    
    def read_log(self):   
        
        # Obtain the number of restraints (= pseudobonds), and the fixed and
        # scanning chain from the DisVis log file
        
        log_path = self.folder + "disvis.log"
        # If the folder does contains a log file, the DisVis job was run on the
        # grid version
        try:
            log_file = open(log_path, "r")
            self.grid = True
        except:
            self.grid = False
            return {"fix": "fixed_chain.pdb",
                    "scan": "scanning_chain.pdb"}
        keys = "fix", "scan"
        disvis_chains = {}
        for line in log_file:
            line_lower = line.lower()
            for key in keys:
                if (line_lower.count(key)) > 0:
                    # Store the name of the chain in a dictionary
                    disvis_chains[key] = self.find_string("\S+\.pdb", line)
                    break
                
        log_file.close()
        return disvis_chains
        
        
    def find_string(self, search_string, whole_string):
        
        return re.search(search_string, whole_string).group(0)
    
        
    def get_files(self, extension):
        
        # Get the files with the specified extension from the DisVis output
        # folder
        
        files = Path(self.folder).glob("*" + extension)
        
        return files
        
        
    def get_chains(self, disvis_chains):
        
        # Store the ChimeraX structures that correspond to the fixed and 
        # scanning chains in the dictionary (overwrite the names)
        
        for model in self.session.models:
            if not model.name in disvis_chains.values():
                continue
            index = list(disvis_chains.values()).index(model.name)
            key = list(disvis_chains.keys())[index]
            disvis_chains[key] = model
        
        return disvis_chains
            
    
    def get_zscores(self):
        
        # Get the zscores from their file (depends on whether the job was run
        # on the grid or local DisVis version)
        
        if self.grid:
            zscore_path = self.folder + "z-score.out"
            zscores = [None] * sum(1 for line in open(zscore_path, "r"))
            file = open(zscore_path, "r")
            for i, line in enumerate(file):
                zscores[i] = float(self.find_string("\S+(?=\n)", line))
                
        else:
            zscore_path = self.folder + "results.html"
            zscores = []
            file = open(zscore_path, "r")
            for line in file:
                strings = re.findall("<td>[^<]+</td>\s+<td>[^<]+</td>\s+<td>[^"
                                     "<]+</td>\s+<td>[^<]+</td>\s+<td>[^<]+</td"
                                     ">", line)
                if len(strings) != 1:
                    continue
                string = strings[0]
                zscore = float(self.find_string("[-\.0-9]+(?=</td>\Z)", 
                               string))
                zscores.append(zscore)
        
        file.close()            
        return zscores   
        
    
    def make_pbs(self, path, disvis_chains, zscores):     
        
        # Make the pseudobonds model
        
        # Use the DisVis input .txt file containing the restraints to find the
        # atoms
        file = open(path, "r") 
        pb_manager = self.session.pb_manager
        self.name = os.path.basename(path).replace(".txt", ".pb")
        group = pb_manager.get_group(self.name)
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
        self.session.models.add([group])
        
        file.close()
        return group.pseudobonds
        
        
    def get_input_distances(self, path):
        
        # Obtain the distance restraints specified in the DisVis search
        
        file = open(path, "r")
        for line in file:
            first_line = line
            break
        
        file.close()
        
        self.minimum = float(self.find_string("\S+(?=\s+\S+\n)", first_line))
        self.maximum = float(self.find_string("\S+(?=\n)", first_line))
        
        
    def color_pbs(self):
        
        # If distance restraint is violated, the pseudobond is colored red. If
        # not, it is colored yellow.
        
        for pb in self.pbs:
            if pb.length >= self.minimum and pb.length <= self.maximum:
                color = self.colors["Main"]
            else:
                color = self.colors["Cutoff"]
            pb.color = color
        
        
    def trigger_handler(self, chains):
        
        # Add a trigger handler for moved models
        
        self.triggerset = self.session.triggers
        self.movement_handler = self.triggerset.add_handler(
            MODEL_POSITION_CHANGED,
            lambda trigger, trigger_data, c=chains: 
                self.handle_movement(trigger, trigger_data, c))
        
        
    def handle_movement(self, trigger, trigger_data, chains):
        
        # If one of the DisVis chains is moved, pseudobonds are re-colored
        # according to their distances

        if not trigger_data in chains:
            return

        self.color_pbs()
        
        
    def prepare_slider(self, zscores):
        
        # Make a slider for the z-scores. Its range is adapted to the minimum
        # and maximum z-scores in the DisVis output.
        
        minimum = min(zscores)
        maximum = max(zscores)
        slider = ZScoreSlider("zscore", True, minimum, maximum, self.pbs)
        main_layout = self.main_dialog.layout
        rows = main_layout.count()
        main_layout.insertLayout(rows - 1, slider.layout)
        
        self.setters = slider.setters
        
        
    def create_pb_file(self):
        
        # The user has clicked OK; pseudobonds that are within the specified
        # z-score range are written in a new .pb file.
        
        if not hasattr(self, "pbs"):
            return
        
        setters = self.setters
        values = [None] * len(setters)
        for i, setter in enumerate(setters):
            values[i] = float(setter.text())
            
        chosen_restraints = []
        
        for pb in self.pbs:
            if pb.outside_range:
                continue
            pb_line = self.xmas.create_pb_line(pb)
            chosen_restraints.append(pb_line)
            
        if len(chosen_restraints) == 0:
            print("No restraints found for this z-score range")
            return
            
        path = self.folder + "Selected_restraints/"
        
        if not os.path.exists(path):
            os.makedirs(path)
            
        file_path = path + self.name.replace(".pb", "_selected.pb")
        
        self.xmas.write_file(file_path, chosen_restraints, "export")
        print("Selected restraints saved in %s" % file_path)
        
        
    def cleanup(self):
        
        # Code runs when the window is closed
        
        if not hasattr(self, "pbs"):
            return
        
        self.pbs.displays = True
        self.pbs.colors = self.colors["Main"]
        self.triggerset.remove_handler(self.movement_handler)
    

from .tool import Slider, ExportSlider
   
class ZScoreSlider(Slider):
    
    # Interactive slider to specify a range in z-scores. If a pseudobond's 
    # z-score is outside the specified range, it is not displayed.    
    
    def __init__(self, value_type="zscore", enabled=True, minimum=None,
                 maximum=None, pbs=None):
        
        # ZScoreSlider inherits from Slider. 
        # The slider behaves weirdly when the minimum values is below zero. 
        # However, z-scores can be negative values. Therefore, we will set the
        # slider minimum to zero and perform a correction using the minimum 
        # z-score multiplied by -1 (self.add)
        
        self.add = -minimum
        super().__init__(value_type, enabled, minimum, maximum, pbs)
        self.function = ExportSlider.display_pseudobonds
        
        
    def get_initial_texts(self, minimum, maximum):
        
        return self.scale(minimum, maximum)
        
        
    def scale(self, minimum, maximum):
        
        return minimum, round(maximum + 0.001, 3)
    
        
    def get_slider_values(self, minimum, maximum):
        
        # Convert real z-scores to slider values

        minimum = (minimum + self.add) * 1000
        maximum = (maximum + self.add) * 1000
        
        return minimum, maximum
        
    
    def within_range(self, value_type, pbs, function):
        
        # Determine whether the pseudobonds' z-scores are within the specified
        # range
        
        invert = self.invert.isChecked()
        score_range = self.slider.value()
        real_range = self.get_real_values(score_range)
        
        for pb in pbs:
            zscore = pb.zscore
            is_outside_range = self.check_value(invert, zscore, real_range)
            pb.outside_range = is_outside_range
            
        self.function(self, pbs)
        
    
    def get_real_values(self, values):
        
        # Convert slider values to real z-scores
        
        real_values = [None] * len(values)
        
        for i, value in enumerate(values):
            real_value = value / 1000 - self.add
            real_values[i] = round(real_value, 3)
            
        return real_values