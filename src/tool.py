# vim: set expandtab shiftwidth=4 softtabstop=4:

# === UCSF ChimeraX Copyright ===
# Copyright 2016 Regents of the University of California.
# All rights reserved.  This software provided pursuant to a
# license agreement containing restrictions on its disclosure,
# duplication and use.  For details see:
# http://www.rbvi.ucsf.edu/chimerax/docs/licensing.html
# This notice must be embedded in or attached to all copies,
# including partial copies, of the software or any revisions
# or derivations thereof.
# === UCSF ChimeraX Copyright ===

from chimerax.core.tools import ToolInstance

class CrosslinkMapper(ToolInstance):

    # Inheriting from ToolInstance makes us known to the ChimeraX tool mangager,
    # so we can be notified and take appropriate action when sessions are closed,
    # saved, or restored, and we will be listed among running tools and so on.
    
    SESSION_ENDURING = False    # Does this instance persist when session closes
    SESSION_SAVE = True         # We do save/restore in sessions
    help = "help:user/tools/tutorial.html"
                                # Let ChimeraX know about our help page

    def __init__(self, session, tool_name):
        # "session"   - chimerax.core.session.Session instance
        # "tool_name" - string

        # Initialize base class
        super().__init__(session, tool_name)

        # Set name displayed on title bar (defaults to tool_name)
        # Must be after the superclass init, which would override it
        self.display_name = "Crosslink Mapper"

        # Create the main window for our tool. The window object will have
        # a "ui_area" where we place the widgets composing our interface.
        # The window isn't shown until its "manage" method is called.
        from chimerax.ui import MainToolWindow
        self.tool_window = MainToolWindow(self)

        # Method to build and show the main window
        self._build_ui()

        # Call trigger handler to take action when certain triggers fire
        self.trigger_handler()
        
        # Override the "cleanup" method to perform additional actions upon closing the main tool window
        self.tool_window.cleanup = self.cleanup


    def _build_ui(self):

        # Put our widgets in the tool window
 
        from PyQt5.QtWidgets import QVBoxLayout, QGridLayout, QHBoxLayout, QTreeWidget, QAbstractItemView, QPushButton, QLabel
        
        outer_layout = QVBoxLayout()
        top_layout = QGridLayout()
        pbonds_layout = QVBoxLayout()
        buttons_layout = QHBoxLayout()

        minimum_tree_height = 120
        
        # A treewidget that will contain all structural protein models open in the session
        self.model_selector = QTreeWidget()
        self.model_selector.setHeaderLabels(["Name", "ID"])
        self.model_selector.setMinimumHeight(minimum_tree_height)

        # A treewidget that will contain PD evidence files that the user has selected
        self.file_selector = QTreeWidget()
        self.file_selector.setHeaderLabels(["Name", "ID"])
        self.file_selector.setColumnWidth(0, 200)
        self.file_selector.setMinimumHeight(minimum_tree_height)

        file_button = QPushButton()
        file_button.setText("Click to select files")
        # Upon clicking the file button, a file dialog is shown, where the user can select
        # .xlsx files
        file_button.clicked.connect(self.dialog)        

        map_button = QPushButton()
        map_button.setText("Click to map crosslinks")
        # Upon clicking the map button, the map_button_clicked method is called twice,
        # each time with different arguments
        map_button.clicked.connect(lambda: 
                                    self.map_button_clicked(selector=self.model_selector, selector_type="model"))
        map_button.clicked.connect(lambda: 
                                    self.map_button_clicked(selector=self.file_selector, selector_type="file"))
        
        top_layout.addWidget(QLabel("Available models"), 0, 0)
        top_layout.addWidget(QLabel("Available files"), 0, 1)
        top_layout.addWidget(self.model_selector, 1, 0)
        top_layout.addWidget(self.file_selector, 1, 1)
        top_layout.addWidget(file_button, 2, 1)
        top_layout.addWidget(map_button, 3, 1)

        # In this treewidget, pseudond models from .pb files are shown;
        # both models that are created with Crosslink Mapper, as well as models that are opened 
        # independently of Crosslink Mapper.
        self.pbonds_menu = QTreeWidget()
        self.pbonds_menu.setHeaderLabels(["Name", "Code (model IDs-file IDs)"])
        self.pbonds_menu.setColumnWidth(0, 300)
        self.pbonds_menu.setMinimumHeight(minimum_tree_height)
        # When a pseudobond model is (de)selected in the menu, it should also be (de)selected in
        # the ChimeraX session. Call "check_signal" method to arrange this        
        self.pbonds_menu.itemChanged.connect(self.check_signal)

        pbonds_layout.addWidget(QLabel("Crosslink models"))
        pbonds_layout.addWidget(self.pbonds_menu)

        subset_button = QPushButton()
        subset_button.setText("Export data subsets (IN PROGRESS)")

        # buttons_layout.addWidget(subset_button) - This button is currently not in use
        # Upon clicking the subset button, a dialog is opened where the user can select the parameters
        # for the data subset required
        subset_button.clicked.connect(self.show_subset_dialog)

        outer_layout.addLayout(top_layout)
        outer_layout.addLayout(pbonds_layout)
        outer_layout.addLayout(buttons_layout)

        # Set the layout as the contents of our window
        self.tool_window.ui_area.setLayout(outer_layout)

        # Show the tool in the ChimeraX main window. Because of its relatively large size, the tool is not shown
        # on the user-preferred side of the window, which would have been arranged by 'self.tool_window.manage("side")'
        self.tool_window.manage(None)

        # Create an empty dictionary where pb models made with the Crosslink Mapper will be stored
        self.created_models = {}

        # Add models open in session to the window with the "add_models" method 
        self.manage_models(self.session.models)

    
    def dialog(self):

        # The user has clicked the file button: with this dialog, evidence files can be selected

        from PyQt5.QtWidgets import QFileDialog, QTreeWidgetItem
        from PyQt5.QtCore import Qt

        file_dialog = QFileDialog()
        file_dialog.setFileMode(QFileDialog.ExistingFiles)
        selected_files, _ = file_dialog.getOpenFileNames(None, "Select evidence files", "", "Excel (*.xlsx)")

        # Create a dictionary containing the evidence files to avoid showing the same file multiple times
        if not hasattr(self, "evidence_files"):
            self.evidence_files = {}
        for selected_file in selected_files:
            if selected_file in self.evidence_files.values():
                selected_files.remove(selected_file)
        number_of_files = len(selected_files)
        # Store the number of files in the treewidget to give each file a proper ID
        if not hasattr(self, "number_of_items"):
            self.number_of_items = 0
        i = 0
        # Files get a unique letter sequence as ID, which is created by the "file_id_generator" method, and are
        # added to the treewidget
        for ID in self.file_id_generator(self.number_of_items, self.number_of_items + number_of_files):
            path_to_file = selected_files[i]
            short_name = self.get_short_filename(path_to_file)
            item = QTreeWidgetItem(self.file_selector)
            item.setText(0, short_name)
            item.setCheckState(0, Qt.Unchecked)
            item.setText(1, ID)
            i += 1
            self.evidence_files[ID] = path_to_file
            self.number_of_items += 1


    def file_id_generator(self, start, end):

        # Creates file IDs for user-selected evidence files. The start index is specified, so that the generator
        # does not start over every time new files are selected. The end index is specified, so that the number of 
        # IDs created corresponds to the number of newly added files.

        import string
        import itertools

        for i in itertools.count(1):
            for p in itertools.product(string.ascii_uppercase, repeat=i):
                start -= 1
                if (end > 0 and start < 0):
                    yield "".join(p)
                    end -= 1
                elif (end > 0 and start > -1):
                    end -= 1
                else:
                    break
            if end < 1:
                break


    def get_short_filename(self, path_to_file):

        # Method to obtain a short name that can be shown to represent a file or model in the tool

        import re

        short_name = re.search("[^/]+\Z", path_to_file).group(0)

        return short_name


    def map_button_clicked(self, selector, selector_type):

        # The user has clicked the map button: the selected models and files are extracted from their
        # respective selector treewidgets

        from PyQt5.QtWidgets import QTreeWidgetItemIterator

        iterator = QTreeWidgetItemIterator(selector, QTreeWidgetItemIterator.Checked)
        # If no models/files have been selected, print a message
        if not iterator.value():
            print("Please check one or more %ss" % selector_type)
            self.missing_data = True
        else:
            self.checked = []
            while iterator.value():
                item = iterator.value()
                ID = item.text(1)
                iterator += 1
                self.checked.append(ID)
            if selector_type == "model":
                self.checked_models = self.checked
                self.missing_data = False
            elif selector_type == "file":
                number_of_checked_files = len(self.checked)
                self.checked_files = [None]*number_of_checked_files
                for i in range(number_of_checked_files):
                    self.checked_files[i] = self.evidence_files[self.checked[i]]
                # If no models have been selected, mapping should not be performed.
                # Therefore, missing_data is first checked. Mapping is performed by a call to 
                # the "map_crosslinks" method
                if not self.missing_data:
                    self.map_crosslinks()


    def map_crosslinks(self):
        
        # The user has selected one or multiple models and evidence files, and clicked the map button: 
        # map crosslinked peptides

        import pandas as pd
        from PyQt5.QtCore import Qt

        # Each checked file is mapped to all checked models
        for j in range(len(self.checked_files)):

            evidence_file = self.checked_files[j]

            # Display bold log message to signify which file is being mapped
            self.session.logger.info("<b>Peptide pair mapping of PD output file: %s</b>" % evidence_file, is_html=True)
            
            # Read the file and extract the peptide pairs from it
            dataframe = pd.read_excel(evidence_file)
            
            input_peptides_A = dataframe["Sequence A"].tolist()
            input_peptides_B = dataframe["Sequence B"].tolist()

            input_pairs = []
            # Keep track of peptide pairs with lacking sequence information
            sequence_info_lacking = 0

            # Peptides within the peptide pairs are sorted, to be able to remove inversely identical peptide pairs
            for i in range(len(dataframe.index)):
                peptide_A = input_peptides_A[i]
                peptide_B = input_peptides_B[i]
                try:
                    input_pairs.append(sorted([peptide_A, peptide_B]))
                except:
                    sequence_info_lacking += 1

            # Print a message when one or multiple peptide pairs lack sequence information
            if sequence_info_lacking == 1:
                print("1 peptide pair is disregarded due to lacking sequence information")
            elif sequence_info_lacking >= 1:
                print("%s peptide pairs are disregarded due to lacking sequence information" % str(sequence_info_lacking))

            # Function that removes duplicate peptide pairs from peptide pair list
            def deduplicate(peptide_pair_list):
                peptide_pair_list.sort()
                last = object()
                for item in peptide_pair_list:
                    if item == last:
                        continue
                    yield item
                    last = item

            input_pairs_deduplicated = list(deduplicate(input_pairs))

            # Print a message stating how many peptide pairs where unique
            number_of_deduplicated = len(input_pairs_deduplicated)
            print("Unique peptide pairs: %s out of %s" % (number_of_deduplicated, len(input_pairs)))

            # Store all peptide pairs in list of lists
            # Each peptide is stored as a dictionary, that contains the peptide sequence and the position of the crosslinked residue.

            def create_peptide_dictionary(input_peptide):
                if input_peptide.count("[") == 1:
                    peptide_sequence = input_peptide.replace("[", "").replace("]", "")
                    crosslink_position = input_peptide.index("[")
                else:
                    peptide_sequence = input_peptide
                    crosslink_position = 0
                peptide_dictionary = {"Peptide sequence":peptide_sequence, "Crosslink position":crosslink_position, "Alignments":[]}
                return peptide_dictionary

            peptide_pairs = [None] * number_of_deduplicated

            for i in range(number_of_deduplicated):
                peptide_pairs[i] = [create_peptide_dictionary(input_pairs_deduplicated[i][0]), 
                                    create_peptide_dictionary(input_pairs_deduplicated[i][1])]

            # Now we align all peptides to the sequences of all chains open in the ChimeraX session
            #
            # Make a dictionary containing strings for the sequences of all chains of all selected models. The peptide sequences are compared
            # to these sequence strings. The index of a crosslinked residue on the sequence string is always different than its 'number' 
            # attribute of the corresponding Residue object in ChimeraX, since numbering Residue objects commences with the number 1 or higher, 
            # instead of 0. To ensure that the right residue number is used in the .pb file, the first residue number of each chain is stored, 
            # which is then added to the index found on the sequence string. Furthermore, some residues are present in the sequence string, but 
            # their Residue objects are not shown in the Chimerax structure. These Residue objects are of type NoneType, and ChimeraX does not 
            # enable creating pseudobonds between NoneType residues. To prevent adding these crosslinks to the .pb file, the positions of all 
            # NoneType residues is also stored in the dictionary. 

            sequences = {}

            for model in self.session.models:
                model_id = model.id_string
                if model_id not in self.checked_models:
                    continue
                # Create a new dictionary for each model
                sequences[model_id] = {}
                for chain in model.chains:
                    # Keep track of the number of NoneType residues at the start of the sequence. These will influence the first residue
                    # number         
                    preceding_nonetypes = 0
                    first_residue_number_found = False
                    sequence = chain.characters
                    nonetype_positions = []
                    residues = chain.residues
                    for i in range(len(residues)):
                        residue = residues[i]
                        if (residue is None and not first_residue_number_found):
                            preceding_nonetypes += 1
                            nonetype_positions.append(i)
                        elif (residue is None and first_residue_number_found):
                            nonetype_positions.append(i)
                        elif (residue is not None and not first_residue_number_found):
                            first_residue_number = residue.number - preceding_nonetypes
                            first_residue_number_found = True
                    # For each chain, create a dictionary containing its sequence, the position of the first residue, 
                    # and the NoneType positions
                    sequences[model_id][chain.chain_id] = {"Chain sequence":sequence, 
                                                            "First residue number":first_residue_number, 
                                                            "NoneType positions":nonetype_positions}

            # Loop all peptide sequences over all chains to find perfect alignments
            for model in sequences:
                all_chains = sequences[model]
                for chain in all_chains:
                    chain_sequence = all_chains[chain]["Chain sequence"]
                    first_residue_number = all_chains[chain]["First residue number"]
                    chain_reference = "#" + model + "/" + chain
                    for peptide_pair in peptide_pairs:
                        for peptide in peptide_pair:
                            peptide_sequence = peptide["Peptide sequence"]
                            peptide_length = len(peptide_sequence)
                            for start in range(len(chain_sequence) - peptide_length + 1):
                                end = start + peptide_length
                                # If the crosslinked residue is not present in the structure, the pseudobond cannot be 
                                # mapped, and therefore we will disregard this alignment
                                if (start + peptide["Crosslink position"] not in all_chains[chain]["NoneType positions"] 
                                        and chain_sequence[start:end] == peptide_sequence):
                                    peptide["Alignments"].append([start + first_residue_number, 
                                                                    end - 1 + first_residue_number, 
                                                                    chain_reference])
            
            # Continue with creating all valid pseudobonds for all peptide pairs and store them in a list
            # If both peptides of a pair are aligned on the same chain, these alignments are checked for overlap. 
            # Make separate lists for peptide pairs with non-overlapping peptides and those with overlapping peptides

            pbonds_filtered = []
            pbonds_filtered_overlapping = []
            pbonds_lists = [pbonds_filtered, pbonds_filtered_overlapping]

            # Function to write an identifying string for each atom
            def write_pseudobond(alignments, overlapping=False):
                pos_A = alignments[0][0] + peptide_pair[0]["Crosslink position"]
                pos_B = alignments[1][0] + peptide_pair[1]["Crosslink position"]
                pbond_sorted = sorted([pbond[0][2] + ":" + str(pos_A) + "@CA", pbond[1][2] + ":" + str(pos_B) + "@CA"])
                if overlapping == False:
                    pbonds_filtered.append(pbond_sorted)
                else:
                    pbonds_filtered_overlapping.append(pbond_sorted)

            # The number of perfectly aligned peptide pairs is counted
            number_of_aligned_pairs = 0

            for peptide_pair in peptide_pairs:
            # First select the peptide pairs for which both peptides have alignments
                if (peptide_pair[0]["Alignments"] and peptide_pair[1]["Alignments"]):
                    number_of_aligned_pairs += 1
                    pbonds_unfiltered = [[a, b] for a in peptide_pair[0]["Alignments"] for b in peptide_pair[1]["Alignments"]]
                    # Then write overlapping and non-overlapping peptide pairs in separate lists
                    for pbond in pbonds_unfiltered:
                        if (pbond[0][2] == pbond[1][2] 
                                and max(pbond[0][0], pbond[1][0]) < min(pbond[0][1], pbond[1][1]) + 1):
                            write_pseudobond(pbond, overlapping=True)
                        else:
                            write_pseudobond(pbond)
            
            # Deduplicate the lists of pseudobonds
            for i in range(len(pbonds_lists)):
                if len(pbonds_lists[i]) > 0:
                    pbonds_lists[i] = list(deduplicate(pbonds_lists[i]))
            
            # Print a log message stating how many peptide pairs could be aligned
            print("Unique peptide pairs with pseudobonds: %s" % number_of_aligned_pairs)
            
            # Store pseudobonds in .pb file:
            #
            # Pseudobonds from peptide pairs with non-overlapping peptides:

            # Function to create a string for each line in the .pb file
            def create_pb_line(pbond):
                return pbond[0] + " " + pbond[1] + "\n"

            # Create a .pb file containing all filtered pbonds
            if len(pbonds_filtered) > 0:
                # Create a code for the model with all model IDs and the evidence file ID that was used
                model_ids = ",".join([str(model_id) for model_id in self.checked_models])
                file_id = self.checked[j]
                pb_file_code = model_ids + "-" + file_id
                pb_file_path = evidence_file.replace(".xlsx", "_%s.pb" % pb_file_code)
                pb_file_name = self.get_short_filename(pb_file_path)
                # Store the model and its code in the "created_models" dictionary
                if pb_file_name not in self.created_models.keys():
                    self.created_models[pb_file_name] = pb_file_code
                  
                pb_file = open(pb_file_path, "w")
                for pbond in pbonds_filtered:
                    pb_file.write(create_pb_line(pbond))
                pb_file.close()
                # Open this .pb file in ChimeraX. It will be added to the pbonds menu via the "manage_models" method that is called
                # due to the "ADD_MODELS" trigger
                from chimerax.core.commands import run
                run(self.session, "open %s" % pb_file_path)
                # Show the code of this file in the pbonds menu
                item = self.pbonds_menu.findItems(pb_file_name, Qt.MatchExactly, column=0)[0]
                item.setText(1, pb_file_code)

            # Pseudobonds from peptide pairs with overlapping peptides:       
            number_of_overlapping_peptides = len(pbonds_filtered_overlapping)

            if number_of_overlapping_peptides > 0:
                # Print a log message stating how many pseudobonds were not mapped due to overlapping peptides
                if number_of_overlapping_peptides == 1:
                    print("1 pseudobond was not mapped due to overlapping peptides.")
                elif number_of_overlapping_peptides > 1:
                    print("%s pseudobonds were not mapped due to overlapping peptides." % number_of_overlapping_peptides)

                # Self-links and non-self-links should be stored in separate files, 
                # since ChimeraX cannot open self-links
                nonselflinks_file_exists = False
                selflinks_file_exists = False

                # Function to to open and/or write .pb files for overlapping peptides
                def create_pb_file_overlapping(pbond, link_type=None, file_exists=False, pb_file=None):
                    if not file_exists:
                        if link_type == "intralink":
                            end_file_name = ""
                        elif link_type == "selflink":
                            end_file_name = "_selflinks"
                        new_pb_file_path = pb_file_path.replace(".pb", "_overlapping%s.pb" % end_file_name)
                        pb_file = open(new_pb_file_path, "w")
                        file_exists = True
                        pb_file.write(create_pb_line(pbond))
                        return new_pb_file_path, pb_file, file_exists
                    else:
                        pb_file.write(create_pb_line(pbond))


                for pbond in pbonds_filtered_overlapping:
                    if (pbond[0] != pbond[1] and nonselflinks_file_exists == False):
                        # Create the file for non-self-links and add the pseudobond to this file
                        (pb_file_path_nonselflinks, 
                         pb_file_nonselflinks, 
                         nonselflinks_file_exists
                         ) = create_pb_file_overlapping(pbond, "intralink")
                    elif (pbond[0] != pbond[1] and nonselflinks_file_exists == True):
                        # Add pseudobond to the file for non-self-links
                        create_pb_file_overlapping(pbond, file_exists=True, pb_file=pb_file_nonselflinks)
                    elif (pbond[0] == pbond[1] and selflinks_file_exists == False):
                        # Create the file for self-links and add the pseudobond to this file
                        number_of_selflinks = 1
                        (pb_file_path_selflinks, 
                         pb_file_selflinks, 
                         selflinks_file_exists
                         ) = create_pb_file_overlapping(pbond, "selflink")
                    elif (pbond[0] == pbond[1] and selflinks_file_exists == True):
                        # Add pseudobond to the file for self-links
                        number_of_selflinks += 1
                        create_pb_file_overlapping(pbond, file_exists=True, pb_file=pb_file_selflinks)

                # Tell the user in which files the pseudobonds from overlapping peptides can be found and close these files
                if (nonselflinks_file_exists and not selflinks_file_exists):
                    print("These pseudobonds were stored in %s." % pb_file_path_nonselflinks)
                    pb_file_nonselflinks.close()
                elif (nonselflinks_file_exists and selflinks_file_exists):
                    print("%s of these pseudobonds were self-links.\n"
                          "Self-links were stored in %s.\n"
                          "The remaining pseudobonds were stored in %s." 
                          % (number_of_selflinks, pb_file_path_selflinks, pb_file_path_nonselflinks)
                          )
                    pb_file_nonselflinks.close()
                    pb_file_selflinks.close()
                elif (not nonselflinks_file_exists and selflinks_file_exists):
                    print("All of these pseudobonds were self-links, and stored in %s." % pb_file_path_selflinks)
                    pb_file_selflinks.close()


    def check_signal(self, item, column):

        # (De)select models depending on the item (de)selected in the pbonds menu

        from PyQt5.QtCore import Qt

        for model in self.session.models:
            if (model.name == item.text(column) and item.checkState(column) == Qt.Checked):
                model.selected = True
            elif (model.name == item.text(column) and item.checkState(column) == Qt.Unchecked):
                model.selected = False
                

    def show_subset_dialog(self):

        # Through this dialog, the user can export a subset of the models selected in the pbonds menu

        from PyQt5.QtWidgets import QGridLayout, QDialog, QListWidget, QListWidgetItem, QLabel
        from PyQt5.QtCore import Qt

        layout = QGridLayout()

        subset_dialog = QDialog()
        subset_dialog.setWindowTitle("Export a subselection of pseudobonds")

        # Menu to select the models included in the subset    
        self.dialog_model_selector = QListWidget()
        # Create an item that can check all models at once
        self.item_all = QListWidgetItem(self.dialog_model_selector)
        self.item_all.setText("All models")
        self.item_all.setCheckState(Qt.Unchecked)
        # Connect a change of selection to the "check_all_models" method. When "All models" has been checked,
        # all models underneath it should also be checked
        self.dialog_model_selector.itemChanged.connect(self.check_all_models)

        # Function to make items for a given QListWidget instance
        def create_listwidgetitems(listwidget, text_list):
            for text in text_list:
                item = QListWidgetItem(listwidget)
                item.setText(text)
                item.setCheckState(Qt.Unchecked)

        # Get the IDs of all protein models that the selected pseudobond models are connected to
        model_ids = self.get_model_ids()

        # Create an item for each model in the list of model IDs
        if len(model_ids) > 0:
            for model_id in model_ids:
                create_listwidgetitems(self.dialog_model_selector, [model_id])

        # Menu to select whether only intralinks, only interlinks, or both need to be exported
        link_selector = QListWidget()
        create_listwidgetitems(link_selector, ["Intralinks", "Interlinks"])

        layout.addWidget(QLabel("Models:"), 0, 0)
        layout.addWidget(self.dialog_model_selector, 1, 0)
        layout.addWidget(QLabel("Link types:"), 0, 1)
        layout.addWidget(link_selector, 1, 1)

        subset_dialog.setLayout(layout)

        subset_dialog.exec_()


    def check_all_models(self, item):

        # if the "All models" checkbox has been (un)checked, all other checkboxes also need to be (un)checked

        from PyQt5.QtCore import Qt

        if item == self.item_all:
            for i in range(1, self.dialog_model_selector.count()):
                current_item = self.dialog_model_selector.item(i)
                if item.checkState() == Qt.Checked:
                    current_item.setCheckState(Qt.Checked)
                else:
                    current_item.setCheckState(Qt.Unchecked)


    def get_model_ids(self):

        # For each selected pseudobond, the IDs of the models that it is connected to are added to a set.
        # This is then converted to a list, containing the IDs of all models that the selected pseudobonds are connected to.

        from chimerax.atomic.pbgroup import PseudobondGroup, selected_pseudobonds
        from chimerax.atomic.structure import Structure
        
        model_ids = set()

        for model in self.session.models:
            if (isinstance(model, PseudobondGroup) or not isinstance(model, Structure)):
                continue
            for pb in selected_pseudobonds(self.session):
                for atom in pb.atoms:
                    if atom in model.atoms:
                        model_ids.add(model.id_string)

        model_ids = list(model_ids)
        model_ids.sort()

        return model_ids


    def manage_models(self, models, model_selector_only=False):

        # The main tool window contains two treewidgets showing models; the model selector for non-pseudobond structural models,
        # and the pbonds menu for pseudobond models. With this method, both treewidgets of the main tool window can be filled.

        from PyQt5.QtWidgets import QTreeWidgetItemIterator

        for model in models:
            # Call method determining to which treewidget the model should be added (if any)
            treewidget, is_model_selector = self.determine_treewidget(model)
            # Skip this model when:
            if (treewidget is None
                # The model does not belong in any of the treewidgets
                or (not is_model_selector and model_selector_only)):
                # or (in the case that we only want to fill the model selector),
                # the model does not belong in the model selector
                continue

            iterator = QTreeWidgetItemIterator(treewidget)
            already_present = False
            while iterator.value():
                item = iterator.value()
                # Check if the model is already present in the dedicated treewidget
                if (item.text(0) == model.name
                    # The model name is already present 
                    and ((is_model_selector and item.text(1) == model.id_string)
                    # In the model selector, different models can have identical names, so also check for identical IDs 
                    or not is_model_selector)):
                    # Pseudobond models should not have the same name
                    already_present = True
                    # If the model is already present:
                    break
                iterator += 1

            # If the model is not yet present, add it to the dedicated treewidget by caling the "add_item" method.
            if not already_present:
                self.add_item(treewidget, model, is_model_selector)


    def determine_treewidget(self, model):

        # Determine which treewidget the model belongs to

        from chimerax.atomic.structure import Structure
        from chimerax.atomic.pbgroup import PseudobondGroup

        # Structural, non-pseudobond models go in the model selector
        if (isinstance(model, Structure) and not isinstance(model, PseudobondGroup)):
            treewidget = self.model_selector
            is_model_selector = True
        # Models made from .pb files go in the pbonds_menu
        elif model.name[-3:] == ".pb":
            treewidget = self.pbonds_menu
            is_model_selector = False
        else:
            treewidget = None
            is_model_selector = False

        return treewidget, is_model_selector             


    def add_item(self, treewidget, model, is_model_selector):

        # Adds an item for a given model to a given treewidget instance

        from PyQt5.QtWidgets import QTreeWidgetItem
        from PyQt5.QtCore import Qt

        model_name = model.name

        item = QTreeWidgetItem(treewidget)
        item.setText(0, model_name)
        item.setCheckState(0, Qt.Unchecked)
        # For models in the model selector, show the model's "id_string" attribute as ID
        if is_model_selector:
            column_1_text = model.id_string
        # For models in the pbonds menu that have been created with Crosslink Mapper, show the model's code
        # NB: if Crosslink Mapper is closed and reopened in the same session, we will start with a new, empty, dictionary
        elif (not is_model_selector and model_name in self.created_models.keys()):
            column_1_text = self.created_models[model_name]
        # Models in the pbonds menu that have not been created with the current instance of Crosslink Mapper
        # will not get a code
        else:
            column_1_text = "N/A"
        item.setText(1, column_1_text)


    def trigger_handler(self):

        # Create trigger handlers for:
        # - change of selection, to change the CheckState of items in the pbonds menu upon (de)selection of
        #   pseudobond models in the session
        # - models being added to and removed from the session, to be added and removed from the tool window

        from chimerax.core.selection import SELECTION_CHANGED
        from chimerax.core.models import ADD_MODELS, REMOVE_MODELS

        self.triggerset = self.session.triggers             
        self.change_selection_handler = self.triggerset.add_handler(SELECTION_CHANGED, self.selection_handler)
        self.add_model_handler = self.triggerset.add_handler(ADD_MODELS, self.model_handler)
        self.remove_model_handler = self.triggerset.add_handler(REMOVE_MODELS, self.model_handler)


    def selection_handler(self, trigger, trigger_data):

        # Called upon change of model selection, to adjust the CheckState of models in the pbonds menu, if necessary

        from chimerax.atomic.pbgroup import PseudobondGroup
        from PyQt5.QtWidgets import QTreeWidgetItemIterator
        from PyQt5.QtCore import Qt

        # Get names of selected pseudobond models in the session
        selected_pb_models = []
        for model in self.session.models:
            if (isinstance(model, PseudobondGroup) and model.selected == True):
                selected_pb_models.append(model.name)
        # Check the selected models and uncheck the deselected ones in the pbonds menu
        iterator = QTreeWidgetItemIterator(self.pbonds_menu)
        while iterator.value():
            item = iterator.value()
            if item.text(0) in selected_pb_models:
                item.setCheckState(0, Qt.Checked)
            else:
                item.setCheckState(0, Qt.Unchecked)
            iterator += 1


    def model_handler(self, trigger, trigger_data):

        # Called when models are added to or removed from the session. The trigger_data is a list of the models that were added or removed

        # If models are added, call "manage_models" method
        if trigger == "add models":
            self.manage_models(trigger_data)
        # If models are removed, call "remove_item" method
        else:
            self.remove_item(trigger_data)


    def remove_item(self, models):

        # Called when models are removed from the session

        from PyQt5.QtCore import Qt

        # Make sure that the "model_selector" treewidget is only updated once, since multiple times would be unnecessary
        model_selector_updated = False

        for model in models:
            # First determine to which treewidget the removed model belongs (model selector or pbonds menu)
            treewidget, is_model_selector = self.determine_treewidget(model)
            # If the model belongs to the model selector treewidget, clear this treewidget and refill it by calling the "manage_models" method
            if is_model_selector and not model_selector_updated:
                treewidget.clear()
                self.manage_models(self.session.models, model_selector_only=True)
                model_selector_updated = True
            # If the model belongs to the pbonds treewidget, remove this item selectively
            elif not is_model_selector and treewidget is not None:
                root = treewidget.invisibleRootItem()
                item = treewidget.findItems(model.name, Qt.MatchExactly, 0)[0]
                root.removeChild(item)                


    def cleanup(self):
        
        # Called when main tool window is closed. Trigger handlers are removed to prevent errors when the trigger fires

        self.triggerset.remove_handler(self.change_selection_handler)
        self.triggerset.remove_handler(self.add_model_handler)
        self.triggerset.remove_handler(self.remove_model_handler)
