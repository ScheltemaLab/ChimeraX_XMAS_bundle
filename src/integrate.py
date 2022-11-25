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


from chimerax.pdb.pdb import save_pdb
import os
from Qt.QtGui import QDoubleValidator
from Qt.QtWidgets import (QVBoxLayout, QPushButton, QRadioButton, QButtonGroup, 
                          QDialogButtonBox, QLabel, QTreeWidgetItemIterator,
                          QHBoxLayout, QCheckBox, QFileDialog, QGridLayout,
                          QLineEdit)


class Integrate:
    
    # A window is opened in which the user can select which Integrate step to
    # perform; creating DisVis input, analyzing DisVis output, and creating
    # HADDOCK input
    
    def __init__(self, xmas_instance):
        
        # Create the window with buttons for each step
        
        self.xmas_instance = xmas_instance
        
        title = "Select an integration step"
        tool_window = xmas_instance.tool_window
        self.steps_window = xmas_instance.create_child_window(tool_window, 
                                                              title)
        layout = QVBoxLayout()     
        
        texts = ("Create DisVis input", "Analyze DisVis output", 
                 "Create HADDOCK input")
        functions = (self.disvis_input, self.disvis_output, self.haddock_input)

        for i, text in enumerate(texts):
            button = QPushButton(text)
            button.clicked.connect(functions[i])
            layout.addWidget(button)
            
        cancel = QDialogButtonBox(QDialogButtonBox.Cancel)
        cancel.rejected.connect(self.steps_window.destroy)
        layout.addWidget(cancel)

        self.steps_window.ui_area.setLayout(layout)
        self.steps_window.manage(None)
        
        
    def disvis_input(self):
        
        # Clicking "Create DisVis input" directs the user to the Export window,
        # in which the checkbox for DisVis is already checked
        
        xmas = self.xmas_instance
        # First check whether pseudobonds have been selected; selected 
        # pseudobonds are required to open the Export window
        empty = xmas.is_selection_empty(xmas.show_subset_dialog)
        
        if empty:
            return
        
        xmas.checkbox_pb.setChecked(False)
        xmas.checkbox_disvis.setChecked(True)
        
        self.steps_window.destroy()
    
    
    def disvis_output(self):
        
        # The user has clicked "Analyze DisVis output"; show the corresponding
        # dialog (built in ZScoreSelector)
        
        from .z_score import ZScoreSelector
        
        self.z_score_selector = ZScoreSelector(self.xmas_instance)
        self.steps_window.destroy()
    
    
    def haddock_input(self):
        
        # The user had clicked "Create HADDOCK input"; show the corresponding
        # dialog. First, the user can select a pseudobonds model. Only single-
        # body docking is supported by XMAS. Therefore, only pseudobonds models
        # that connect two molecular models can be selected
        
        xmas = self.xmas_instance
        tool_window = xmas.tool_window
        title = "Create HADDOCK input"
        self.haddock_window = xmas.create_child_window(tool_window, title)        
        self.haddock_layout = QVBoxLayout()
        
        self.haddock_layout.addWidget(QLabel("Select pseudobonds model:"))
        # Add QRadioButtons for each valid model to the dialog
        self.add_pb_models(xmas)
        
        self.ok_cancel = QDialogButtonBox(QDialogButtonBox.Ok
                                          | QDialogButtonBox.Cancel)
        # When the user has selected a pseudobonds model, chain A, and the  
        # input type(s), and has clicked "OK", the output is created
        self.ok_cancel.accepted.connect(self.create_output)
        self.ok_cancel.rejected.connect(self.haddock_window.destroy)
        # "OK" button is only be enabled if a pseudobonds model has been
        # selected
        self.ok_cancel.button(QDialogButtonBox.Ok).setEnabled(False)
        self.haddock_layout.addWidget(self.ok_cancel)
        self.haddock_layout.expanded = False
        
        self.steps_window.destroy()
        self.haddock_window.ui_area.setLayout(self.haddock_layout)
        self.haddock_window.manage(None)
        
        
    def add_pb_models(self, xmas):
        
        # The user is allowed to select pseudobond models that connect two
        # molecular models. QRadioButtons are shown for each.
        # The pseudobond models are taken from XMAS's main window, and then
        # checked for suitability for HADDOCK input (i.e., whether they
        # connect two molecular models)
        
        pbonds_menu = xmas.pbonds_menu
        iterator = QTreeWidgetItemIterator(pbonds_menu)
        available_models = False
        self.pb_group = QButtonGroup(self.haddock_layout)
        
        while iterator.value():
            item = iterator.value()
            model = item.model
            unique_structures = list(model.pseudobonds.unique_structures)
            if not len(unique_structures) == 2:
                continue
            available_models = True
            id_string = model.id_string
            text = item.text(0) + " (" + id_string + ")"
            button = QRadioButton(text)
            button.model = model
            self.haddock_layout.addWidget(button)
            self.pb_group.addButton(button)
            iterator += 1
            
        if not available_models:
            label = "No valid pseudobond models available"
            self.haddock_layout.addWidget(QLabel(label))
            return
        
        # When a model is selected, options to select chain A, the output 
        # type(s), and the distances for in the Restraints file become 
        # available for the selected model
        self.pb_group.buttonClicked.connect(self.select_chain_a)      
    
        
    def select_chain_a(self, button):
        
        # HADDOCK requires the chain ID to be A for the first molecular model,
        # and B for the second molecular model. The user can now select which
        # molecular model will be chain A. The other model will be chain B
        
        # Before any pseudobonds model has been selected, the options to select 
        # chain A, the input type(s), and the distances for in the Restraints 
        # file are not inserted in the layout yet ("expanded" attribute set to 
        # False). They are inserted when a model's QRadioButton is clicked.
        # However, when the layout is already expanded, and a different model
        # is selected, different chains for selecting chain A is displayed. 
        # Therefore, the corresponding buttons are deleted
        if self.haddock_layout.expanded:
            for chain_button in self.chain_buttons:
                self.haddock_layout.removeWidget(chain_button)
                subtract = 6
                add_boxes = False
                self.ok_cancel.button(QDialogButtonBox.Ok).setEnabled(False)
        else:
            subtract = 1
            add_boxes = True
            widgets = [QLabel(""), QLabel("Select chain A:")]
            for widget in widgets:
                index = self.haddock_layout.count() - subtract
                self.haddock_layout.insertWidget(index, widget)
                
        self.haddock_layout.expanded = True
        
        # The pseudobonds model's unique structures (i.e. the molecular models
        # that it connects to each other) are extracted. A QRadioButton is 
        # inserted in the layout for each structure, so that one of the models
        # can be selected as chain A
        unique_structures = list(button.model.pseudobonds.unique_structures)
        unique_structures.sort(key=lambda s: s.id_string)
        
        self.chain_group = QButtonGroup(self.haddock_layout)
        self.chain_group.buttonClicked.connect(self.enable_ok)
        
        # Widgets are inserted above the "OK" and "Cancel" buttons
        self.chain_buttons = [None] * len(unique_structures)
        for i, structure in enumerate(unique_structures):
            button_text = structure.name + " (" + structure.id_string + ")"
            button = QRadioButton(button_text)
            index = self.haddock_layout.count() - subtract
            self.haddock_layout.insertWidget(index, button)
            self.chain_group.addButton(button)
            button.id_string = structure.id_string
            self.chain_buttons[i] = button
            button.chain = structure
        
        # The first time that the user has selected a pseudobonds model, the
        # options to select the input type(s) and the distances for in the 
        # restraints file appear
        if add_boxes:
            self.select_output()
        
        
    def enable_ok(self):
        
        # "OK" is only enabled when the user has selected chain A
        
        self.ok_cancel.button(QDialogButtonBox.Ok).setEnabled(True)
        
    
    def select_output(self):
        
        # Widgets to select the input type(s) and the distances for in the 
        # restraints file are inserted in the layout
        
        # The user can choose to print interface residue numbers to the 
        # ChimeraX log, create a Restraints file, or create PDB files for the 
        # molecular models, or any combination thereof
        self.output_layout = QHBoxLayout()
        outputs = ("Residue numbers", "Restraints file", "PDB files")
        self.checkboxes = [None] * len(outputs)
        for i, output in enumerate(outputs):
            box = QCheckBox(output)
            self.checkboxes[i] = box
            self.output_layout.addWidget(box)
            box.setChecked(True)
        
        # Insert QLineEdits to specify the lower, median, and upper distance
        # to be used in the Restraints file
        distance_layout = QGridLayout()
        distance_layout.addWidget(QLabel("Define restraint distances:"), 0, 0, 
                                  1, 2)
        distances = ("Lower:", "Median:", "Upper:")
        presets = ("5.0", "10.0", "25.0")
        self.line_edits = [None] * len(distances)
        for i, distance in enumerate(distances):
            row = i + 1
            distance_layout.addWidget(QLabel(distance), row, 0)
            line_edit = QLineEdit()
            line_edit.setValidator(QDoubleValidator(0.0, float("inf"), 1))
            line_edit.setPlaceholderText(presets[i])
            distance_layout.addWidget(line_edit, row, 1)
            self.line_edits[i] = line_edit
        
        # Widgets are inserted above the "OK" and "Cancel" buttons
        index = self.haddock_layout.count() - 1
        self.haddock_layout.insertWidget(index, QLabel(""))
        self.haddock_layout.insertWidget(index + 1, 
                                         QLabel("Select input type(s):"))
        self.haddock_layout.insertLayout(index + 2, self.output_layout)
        self.haddock_layout.insertWidget(index + 3, QLabel("")) 
        self.haddock_layout.insertLayout(index + 4, distance_layout)
            
                    
    def create_output(self):
        
        # The user has selected all parameters and clicked "OK"
        
        # The input type(s) are extracted
        types = (numbers, restraints_file, pdbs) = self.output_types()
        if not True in types:
            print("Please select input type(s)")
            return
        
        # The pseudobonds model, the molecular model that will be chain A, and 
        # the molecular model that will be chain B are extracted
        pb_model, chain_a, chain_b = self.input_models()
        
        # A copy is made for each molecular model. From these copies, the
        # residues will be renumbered from 1 and the PDB files created
        copy_a = chain_a.copy("Copy A")
        copy_b = chain_b.copy("Copy B")
                
        copied_models = {chain_a.id_string: copy_a, chain_b.id_string: copy_b}
        
        session = self.xmas_instance.session
        
        # A copy of the selected pseudobonds model, connecting the copies of
        # the molecular models, is created
        group = self.create_group(session, pb_model, copied_models)
        
        # The user should select a folder in which to save the restraints and 
        # PDB files
        if restraints_file or pdbs:
            get_folder = QFileDialog.getExistingDirectory
            folder = get_folder(None, "Select a folder for HADDOCK input", "",
                                QFileDialog.ShowDirsOnly 
                                | QFileDialog.DontResolveSymlinks)
            if folder == "":
                group.delete()
                return
            self.path = folder + "/"
        
        # Renumber residues and create PDBs from the copies of the molecular
        # models
        self.renumber_and_create_pdbs(session, copied_models.values(), pdbs)
        
        # Create a dictionary from the pseudobonds model copy, to extract the
        # information for the interface residues and Restraints file
        if numbers or restraints_file:
            atoms_dict = self.create_atoms_dict(group, copy_a)
            if numbers:
                self.print_residue_numbers(atoms_dict, session)
            if restraints_file:
                pb_model_name = os.path.splitext(pb_model.name)[0]
                self.create_restraints_file(atoms_dict, pb_model_name)
        
        for model in [group] + list(copied_models.values()):
            model.delete()
        self.haddock_window.destroy()
        
        
    def output_types(self):
        
        # Find out which output checkboxes are checked
        
        checked = [box.isChecked() for box in self.checkboxes]
        
        return tuple(checked)
        
        
    def input_models(self):
        
        # Check which pseudobonds model is selected and which of its unique
        # structures is selected as chain A
        
        button = self.pb_group.checkedButton()
        pb_model = button.model
        
        chains = [button.chain for button in self.chain_group.buttons()]
        chain_a_id = self.checked_button_id(self.chain_group)
        for chain in chains:
            if chain.id_string == chain_a_id:
                chain_a = chain
            else:
                chain_b = chain
                
        return pb_model, chain_a, chain_b
            
            
    def checked_button_id(self, button_group):
        
        # Check which molecular model the user has selected as chain A
        
        return button_group.checkedButton().id_string
    
    
    def create_group(self, session, pb_model, copied_models):
        
        # Create a copy of the selected pseudobonds model, so that the residue
        # numbers can be extracted from this model after the residues have
        # been renumbered from 1. The copy connects the copies of the molecular
        # model
        
        group = session.pb_manager.get_group("group")
        
        for pb in pb_model.pseudobonds:
            atom1, atom2 = pb.atoms
            # Intralinks are excluded
            if atom1.structure == atom2.structure:
                continue
            atoms = [atom1, atom2]
            new_atoms = [None] * 2
            for i, atom in enumerate(atoms):
                residue = atom.residue
                model = copied_models[residue.structure.id_string]
                new_residue = model.find_residue(residue.chain_id, 
                                                 residue.number)
                new_atom = new_residue.find_atom(atom.name)
                new_atoms[i] = new_atom
            group.new_pseudobond(new_atoms[0], new_atoms[1])
            
        return group
    
    
    def renumber_and_create_pdbs(self, session, models, pdbs):
        
        # HADDOCK requires the residues from both models to be renumbered 
        # from 1. It also requires the chain ID to be A for the first molecular 
        # model, and B for the second molecular model
        
        # Renumber the residues
        for i, model in enumerate(models):
            start = 1
            chains = model.chains
            chain_id = ["A", "B"][i]
            replace_chains = False
            for chain in chains:
                residues = chain.existing_residues
                for residue in residues:
                    residue.number = start
                    start += 1
                # Store whether chain IDs should be changed
                if (chain.chain_id != chain_id and not replace_chains):
                    replace_chains = True
            
            if not pdbs:
                continue
            
            # PDB files are created from the molecular models. If chain IDs
            # should be replaced, a temporary file with the original chain IDs
            # is first created. From this temporary file, the final PDB file
            # with correct chain IDs is created and the temporary file is 
            # deleted
            ext = {True: "_TEMP", False: ""}
            path = (self.path + "Chain" + chain_id + ext[replace_chains] 
                    + ".pdb")
            save_pdb(session, path, models=[model])
            
            if not replace_chains:  
                continue
            
            # Replace chain IDs in the temporary PDB file
            temp = open(path, "r")
            haddock = open("".join(path.rsplit("_TEMP", 1)), "w")      
            for line in temp:
                if line.startswith(("ANISOU", "ATOM", "HETATM", "TER")):
                    chain = line[22]
                    if chain != chain_id:
                        line = line[:21] + chain_id + line[22:]
                haddock.write(line)
            temp.close()
            haddock.close()
            os.remove(path)
            
            
    def create_atoms_dict(self, group, chain_a):
        
        # Create a dictionary from the pseudobonds model copy. The atoms 
        # belonging to chain A are the keys. For each key atom, a list with the
        # connected atoms that belong to chain B is stored
        
        atoms_dict = {}
        
        for pb in group.pseudobonds:
            atom1, atom2 = pb.atoms
            atoms = [atom1, atom2]
            for atom in atoms:
                if atom.structure == chain_a:
                    key = atom
                else:
                    value = atom
            if key not in atoms_dict.keys():
                atoms_dict[key] = []
            atoms_dict[key].append(value)
            
        return atoms_dict
        
        
    def print_residue_numbers(self, atoms_dict, session):
        
        # Print the residues numbers for the interface residues in the 
        # ChimeraX log
        
        numbers_a = [int(atom.residue.number) for atom in atoms_dict.keys()]
        numbers_a = self.get_residues_string(numbers_a)
        
        numbers_b = set()
        for atoms in atoms_dict.values():
            numbers = [atom.residue.number for atom in atoms]
            for number in numbers:
                numbers_b.add(number)
        numbers_b = self.get_residues_string(list(numbers_b))
        session.logger.info("<b>Residue numbers chain A:</b><br>%s<br>"
                            "<b>Residue numbers chain B:</b><br>%s" 
                            % (numbers_a, numbers_b), is_html=True)
                            
                            
    def get_residues_string(self, numbers_list):
        
        # Create a string from a list of residue numbers
        
        numbers_list.sort()
        numbers_list = [str(number) for number in numbers_list]
        
        return ", ".join(numbers_list)
        
        
    def create_restraints_file(self, atoms_dict, pb_model_name):
        
        # Create the Restraints file
        
        distances = self.get_distances()
        
        file_name = self.path + pb_model_name + ".tbl"
        file = open(file_name, "w")
        file.write("! HADDOCK AIR restraints\n!\n")
    
        for atom in atoms_dict:
            file.write("assign ( name %s and resid %s  and segid A)"
                       "\n       (\n" 
                        % (atom.name.lower(), atom.residue.number))
            other_atoms = atoms_dict[atom]
            for i, other_atom in enumerate(other_atoms):
                file.write("        ( name %s and resid %s  and segid B)\n" 
                            % (other_atom.name.lower(), 
                               other_atom.residue.number))
                if i < len(other_atoms) - 1:
                    file.write("      or\n")
                else:
                    file.write("       )  %s\n\n" % distances)
                    
        file.close()
    
    
    def get_distances(self):
        
        # Get the distances that the user has specified in the dedicated 
        # QLineEdits for the Restraints file
        
        distances = [None] * len(self.line_edits)
        
        for i, line_edit in enumerate(self.line_edits):
            distance = line_edit.text()
            if distance == "":
                distance = line_edit.placeholderText()
            distances[i] = str(float(distance))
            
        return distances[1] + " " + distances[0] + " " + distances[2]