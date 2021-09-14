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

        # Initialize base class.
        super().__init__(session, tool_name)

        # Set name displayed on title bar (defaults to tool_name)
        # Must be after the superclass init, which would override it.
        self.display_name = "Crosslink Mapper"

        # Create the main window for our tool.  The window object will have
        # a "ui_area" where we place the widgets composing our interface.
        # The window isn't shown until we call its "manage" method.

        from chimerax.ui import MainToolWindow
        self.tool_window = MainToolWindow(self)

        # Override the destroy method to perform additional action upon destruction of window
        self.tool_window.cleanup = self.cleanup

        # Our user interface is simple enough that we could probably inline
        # the code right here, but for any kind of even moderately complex
        # interface, it is probably better to put the code in a method so
        # that this __init__ method remains readable.
        self._build_ui()

        # Call trigger handler
        self.trigger_handler()

    def _build_ui(self):
        # Put our widgets in the tool window
 
        from PyQt5.QtWidgets import QVBoxLayout, QGridLayout, QHBoxLayout, QTreeWidget, QAbstractItemView, QPushButton, QLabel
        
        outer_layout = QVBoxLayout()
        top_layout = QGridLayout()
        pbonds_layout = QVBoxLayout()
        buttons_layout = QHBoxLayout()

        minimum_tree_height = 120
        
        self.model_selector = QTreeWidget()
        self.model_selector.setHeaderLabels(["Name", "ID"])
        self.model_selector.setMinimumHeight(minimum_tree_height)

        self.file_selector = QTreeWidget()
        self.file_selector.setHeaderLabels(["Name", "ID"])
        self.file_selector.setColumnWidth(0, 200)
        self.file_selector.setMinimumHeight(minimum_tree_height)

        file_button = QPushButton()
        file_button.setText("Click to select files")

        map_button = QPushButton()
        map_button.setText("Click to map crosslinks")
        
        top_layout.addWidget(QLabel("Available models"), 0, 0)
        top_layout.addWidget(QLabel("Available files"), 0, 1)
        top_layout.addWidget(self.model_selector, 1, 0)
        top_layout.addWidget(self.file_selector, 1, 1)
        top_layout.addWidget(file_button, 2, 1)
        top_layout.addWidget(map_button, 3, 1)

        self.pbonds_menu = QTreeWidget()
        self.pbonds_menu.setHeaderLabels(["Name", "Code (model IDs-file IDs)"])
        self.pbonds_menu.setColumnWidth(0, 300)
        self.pbonds_menu.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.pbonds_menu.setMinimumHeight(minimum_tree_height)        
        self.pbonds_menu.itemChanged.connect(self.check_signal)

        pbonds_layout.addWidget(QLabel("Crosslink models"))
        pbonds_layout.addWidget(self.pbonds_menu)

        subset_button = QPushButton()
        subset_button.setText("Export data subsets (IN PROGRESS)")

        #buttons_layout.addWidget(subset_button)

        outer_layout.addLayout(top_layout)
        outer_layout.addLayout(pbonds_layout)
        outer_layout.addLayout(buttons_layout)

        file_button.clicked.connect(self.dialog)
        map_button.clicked.connect(lambda: self.map_button_clicked(selector=self.model_selector, selector_type="model"))
        map_button.clicked.connect(lambda: self.map_button_clicked(selector=self.file_selector, selector_type="file"))
        subset_button.clicked.connect(self.show_subset_dialog)

        # Set the layout as the contents of our window
        self.tool_window.ui_area.setLayout(outer_layout)

        self.tool_window.manage(None)

        # Create an empty dictionary where pb models generated with the tool will be stored
        self.generated_models = {}

        # Add models open in session to the window with the add_models method 
        self.manage_models(self.session.models)


    def show_subset_dialog(self):
        from PyQt5.QtWidgets import QGridLayout, QDialog, QListWidget, QListWidgetItem, QLabel
        from PyQt5.QtCore import Qt
        layout = QGridLayout()

        subset_dialog = QDialog()
        subset_dialog.setWindowTitle("Export a subselection of pseudobonds")
            
        self.dialog_model_selector = QListWidget()
        # model_selector.itemChanged.connect(self.check_all_models(model_selector, item_all))

        self.item_all = QListWidgetItem(self.dialog_model_selector)
        self.item_all.setText("All models")
        self.item_all.setCheckState(Qt.Unchecked)
        self.dialog_model_selector.itemChanged.connect(self.check_all_models)

        def generate_listwidgetitems(listwidget, text_list):
            for text in text_list:
                item = QListWidgetItem(listwidget)
                item.setText(text)
                item.setCheckState(Qt.Unchecked)

        model_ids = self.get_model_ids()

        if len(model_ids) > 0:
            for model_id in model_ids:
                generate_listwidgetitems(self.dialog_model_selector, [model_id])

        link_selector = QListWidget()

        generate_listwidgetitems(link_selector, ["Intralinks", "Interlinks"])

        layout.addWidget(QLabel("Models:"), 0, 0)
        layout.addWidget(self.dialog_model_selector, 1, 0)
        layout.addWidget(QLabel("Link types:"), 0, 1)
        layout.addWidget(link_selector, 1, 1)

        subset_dialog.setLayout(layout)

        subset_dialog.exec_()


    def get_model_ids(self):
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

    
    def check_all_models(self, item):
        from PyQt5.QtCore import Qt

        if item == self.item_all:
            for i in range(1, self.dialog_model_selector.count()):
                current_item = self.dialog_model_selector.item(i)
                if item.checkState() == Qt.Checked:
                    current_item.setCheckState(Qt.Checked)
                else:
                    current_item.setCheckState(Qt.Unchecked)


    def manage_models(self, models, model_selector_only=False):
        from PyQt5.QtWidgets import QTreeWidgetItemIterator

        for model in models:
            self.determine_treewidget(model)
            if self.treewidget is None:
                continue
            if not self.is_model_selector and model_selector_only:
                continue

            iterator = QTreeWidgetItemIterator(self.treewidget)
            already_present = False
            while iterator.value():
                item = iterator.value()
                if (item.text(0) == model.name and ((self.is_model_selector and item.text(1) == model.id_string) or not self.is_model_selector)):
                    already_present = True
                    break
                iterator += 1

            if not already_present:
                self.add_item(self.treewidget, model)   


    def determine_treewidget(self, model):
        from chimerax.atomic.structure import Structure
        from chimerax.atomic.pbgroup import PseudobondGroup

        # Structural, non-pseudobond models go in the model selector
        if (isinstance(model, Structure) and not isinstance(model, PseudobondGroup)):
            self.treewidget = self.model_selector
            self.is_model_selector = True
        # Models made from .pb files go in the pbonds_menu
        elif model.name[-3:] == ".pb":
            self.treewidget = self.pbonds_menu
            self.is_model_selector = False
        else:
            self.treewidget = None
            self.is_model_selector = False          


    def add_item(self, treewidget, model):
        from PyQt5.QtWidgets import QTreeWidgetItem
        from PyQt5.QtCore import Qt

        model_name = model.name

        item = QTreeWidgetItem(treewidget)
        item.setText(0, model_name)
        item.setCheckState(0, Qt.Unchecked)
        if self.is_model_selector:
            column_1_text = model.id_string
        elif (not self.is_model_selector and model_name in self.generated_models.keys()):
            column_1_text = self.generated_models[model_name]
        else:
            column_1_text = "N/A"
        item.setText(1, column_1_text)


    def remove_item(self, models):
        from PyQt5.QtCore import Qt

        model_selector_updated = False

        for model in models:
            self.determine_treewidget(model)
            if self.is_model_selector and not model_selector_updated:
                self.treewidget.clear()
                self.manage_models(self.session.models, model_selector_only=True)
                model_selector_updated = True
            elif not self.is_model_selector and self.treewidget is not None:
                root = self.treewidget.invisibleRootItem()
                item = self.treewidget.findItems(model.name, Qt.MatchExactly, 0)[0]
                root.removeChild(item)


    def check_signal(self, item, column):
        from PyQt5.QtCore import Qt

        for model in self.session.models:
            if (model.name == item.text(column) and item.checkState(column) == Qt.Checked):
                model.selected = True
            elif (model.name == item.text(column) and item.checkState(column) == Qt.Unchecked):
                model.selected = False


    def file_id_generator(self, start, end):
        import string, itertools
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

        
    def dialog(self):
        from PyQt5.QtWidgets import QFileDialog, QTreeWidgetItem
        from PyQt5.QtCore import Qt
        import re

        if not hasattr(self, "input_files"):
            self.input_files = {}
        file_dialog = QFileDialog()
        file_dialog.setFileMode(QFileDialog.ExistingFiles)
        selected_files, _ = file_dialog.getOpenFileNames(None, "Select input files", "", "Excel (*.xlsx)")

        # Make sure that files aren't added multiple times
        for selected_file in selected_files:
            if selected_file in self.input_files.values():
                selected_files.remove(selected_file)
        number_of_files = len(selected_files)
        if not hasattr(self, "number_of_items"):
            self.number_of_items = 0
        i = 0
        for ID in self.file_id_generator(self.number_of_items, self.number_of_items + number_of_files):
            file_name = selected_files[i]
            file_name_short = re.search("[^/]+\Z", file_name).group(0)
            item = QTreeWidgetItem(self.file_selector)
            item.setText(0, file_name_short)
            item.setCheckState(0, Qt.Unchecked)
            item.setText(1, ID)
            i += 1
            self.input_files[ID] = file_name
            self.number_of_items += 1

    
    # Function for displaying log messages
    def display_log_message(self, msg):
        self.session.logger.info(msg, is_html=True)


    def map_button_clicked(self, selector, selector_type):
        from PyQt5.QtWidgets import QTreeWidgetItemIterator

        iterator = QTreeWidgetItemIterator(selector, QTreeWidgetItemIterator.Checked)
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
                    self.checked_files[i] = self.input_files[self.checked[i]]
                if not self.missing_data:
                    self.map_crosslinks()


    def map_crosslinks(self):
        
        # The user has clicked the "Run" button; map crosslinked peptides
        import pandas as pd
        from chimerax.atomic import Structure
        import re
        from PyQt5.QtWidgets import QTreeWidgetItem
        from PyQt5.QtCore import Qt

        for j in range(len(self.checked_files)):

            input_file = self.checked_files[j]

            # Display first log message
            self.display_log_message("<b>Peptide pair mapping of PD output file: %s</b>" % input_file)
            
            # Read the file that the user has selected and extract the peptide pairs from it
            dataframe = pd.read_excel(input_file)
            
            input_peptides_A = dataframe["Sequence A"].tolist()
            input_peptides_B = dataframe["Sequence B"].tolist()

            input_pairs = []
            # Keep track of peptide pairs with missing sequences
            pairs_with_missing_sequence = 0

            # Peptides within the peptide pairs are sorted, to be able to remove inversely identical peptide pairs
            for i in range(len(dataframe.index)):
                peptide_A = input_peptides_A[i]
                peptide_B = input_peptides_B[i]
                try:
                    input_pairs.append(sorted([peptide_A, peptide_B]))
                # In case of missing sequence:
                except:
                    pairs_with_missing_sequence += 1

            if pairs_with_missing_sequence == 1:
                print("1 peptide pair is disregarded due to missing peptide sequence")
            elif pairs_with_missing_sequence >= 1:
                print("%s peptide pairs are disregarded due to missing peptide sequence" % str(pairs_with_missing_sequence))

            # Duplicate peptide pairs are removed
            def deduplicate(lst):
                lst.sort()
                last = object()
                for item in lst:
                    if item == last:
                        continue
                    yield item
                    last = item

            input_pairs_deduplicated = list(deduplicate(input_pairs))

            # Display a message if duplicates have been removed
            number_of_deduplicated = len(input_pairs_deduplicated)
            self.display_log_message("Unique peptide pairs: %s out of %s" % (number_of_deduplicated, len(input_pairs)))

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
                peptide_pairs[i] = [create_peptide_dictionary(input_pairs_deduplicated[i][0]), create_peptide_dictionary(input_pairs_deduplicated[i][1])]

            # Now we align all peptides to the sequences of all chains open in the ChimeraX session
            #
            # Make a dictionary containing all the sequences of all chains of all structural models and the starting positions, i.e. the
            # position in the sequence of the first residue of the structure. The positions of the residues that are absent from the structure
            # are also stored in the dictionary as unmapped positions.
            sequences = {}

            for model in self.session.models:
                model_id = model.id_string
                if model_id not in self.checked_models:
                    continue
                sequences[model_id] = {}
                for chain in model.chains:          
                    preceding_nonetypes = 0
                    starting_position_found = False
                    sequence = chain.characters
                    unmapped_positions = []
                    current_position = -1
                    for residue in chain.residues:
                        current_position += 1
                        if (residue is None and not starting_position_found):
                            preceding_nonetypes += 1
                            unmapped_positions.append(current_position)
                            continue
                        elif (residue is None and starting_position_found):
                            unmapped_positions.append(current_position)
                        elif (residue is not None and not starting_position_found):
                            starting_position = residue.number - preceding_nonetypes
                            starting_position_found = True
                    sequences[model_id][chain.chain_id] = {"Chain sequence":sequence, "Starting position":starting_position, "Unmapped positions":unmapped_positions}

            # Loop all peptide sequences over all chains for the alignment
            for structure in sequences:
                all_chains = sequences[structure]
                for chain in all_chains:
                    chain_sequence = all_chains[chain]["Chain sequence"]
                    starting_position = all_chains[chain]["Starting position"]
                    chain_reference = "#" + structure + "/" + chain
                    for peptide_pair in peptide_pairs:
                        for peptide in peptide_pair:
                            peptide_sequence = peptide["Peptide sequence"]
                            peptide_length = len(peptide_sequence)
                            for start in range(len(chain_sequence) - peptide_length + 1):
                                end = start + peptide_length
                                # If the crosslinked residue is not present in the structure, the pseudobond cannot be 
                                # mapped, and therefore we will disregard this alignment
                                if (start + peptide["Crosslink position"] not in all_chains[chain]["Unmapped positions"] and chain_sequence[start:end] == peptide_sequence):
                                    peptide["Alignments"].append([start + starting_position, end - 1 + starting_position, chain_reference])
            
            # Continue with creating all valid pseudobonds for all peptide pairs and store them in a list
            # If both peptides of a pair are aligned on the same chain, these alignments are checked for overlap. 
            # Make separate lists for peptide pairs with non-overlapping peptides and those with overlapping peptides

            pbonds_filtered = []
            pbonds_filtered_overlapping = []
            pbonds_lists = [pbonds_filtered, pbonds_filtered_overlapping]

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
            # First select the peptide pairs for which both peptides have alignments.
                if (peptide_pair[0]["Alignments"] and peptide_pair[1]["Alignments"]):
                    number_of_aligned_pairs += 1
                    pbonds_unfiltered = [[a, b] for a in peptide_pair[0]["Alignments"] for b in peptide_pair[1]["Alignments"]]
                    for pbond in pbonds_unfiltered:
                        if (pbond[0][2] == pbond[1][2] and max(pbond[0][0], pbond[1][0]) < min(pbond[0][1], pbond[1][1]) + 1):
                            write_pseudobond(pbond, overlapping=True)
                        else:
                            write_pseudobond(pbond)
            
            # Deduplicate the lists of pseudobonds
            for i in range(len(pbonds_lists)):
                if len(pbonds_lists[i]) > 0:
                    pbonds_lists[i] = list(deduplicate(pbonds_lists[i]))
            
            # Display a log message stating how many peptide pairs could be aligned
            self.display_log_message("Unique peptide pairs with pseudobonds: %s" % number_of_aligned_pairs)
            
            # Store pseudobonds in .pb file

            # From peptide pairs with non-overlapping peptides

            def generate_pb_line(pbond):
                return pbond[0] + " " + pbond[1] + "\n"

            if len(pbonds_filtered) > 0:
                model_ids = ",".join([str(model_id) for model_id in self.checked_models])
                file_id = self.checked[j]
                pb_file_code = model_ids + "-" + file_id
                pb_file_path = input_file.replace(".xlsx", "_%s.pb" % pb_file_code)
                pb_file_name = re.search("[^/]+\Z", pb_file_path).group(0)
                if pb_file_name not in self.generated_models.keys():
                    self.generated_models[pb_file_name] = pb_file_code
                  
                pb_file = open(pb_file_path, "w")
                for pbond in pbonds_filtered:
                    pb_file.write(generate_pb_line(pbond))
                pb_file.close()
                # Open this .pb file in ChimeraX. It will automatically be added to the pbonds menu
                from chimerax.core.commands import run
                run(self.session, "open %s" % pb_file_path)
                # Show the code of this file in the pseudobonds menu
                item = self.pbonds_menu.findItems(pb_file_name, Qt.MatchExactly, column=0)[0]
                item.setText(1, pb_file_code)

            # From peptide pairs with overlapping peptides        
            number_of_overlapping_peptides = len(pbonds_filtered_overlapping)

            if number_of_overlapping_peptides > 0:
                # Display a log message stating how many pseudobonds were not mapped due to overlap
                if number_of_overlapping_peptides == 1:
                    self.display_log_message("1 pseudobond was not mapped due to overlapping peptides.")
                elif number_of_overlapping_peptides > 1:
                    self.display_log_message("%s pseudobonds were not mapped due to overlapping peptides." % number_of_overlapping_peptides)

                # Self-links and non-self-links (here called intralinks) should be stored in separate files, since ChimeraX cannot process self-links
                intralinks_found = False
                selflinks_found = False           

                for pbond in pbonds_filtered_overlapping:
                    if (pbond[0] != pbond[1] and intralinks_found == False):
                        pb_file_path_intralinks = pb_file_path.replace(".pb", "_overlapping.pb")
                        pb_file_intralinks = open(pb_file_path_intralinks, "w")
                        pb_file_intralinks.write(generate_pb_line(pbond))
                        intralinks_found = True
                    elif (pbond[0] != pbond[1] and intralinks_found == True):
                        pb_file_intralinks.write(generate_pb_line(pbond))
                    elif (pbond[0] == pbond[1] and selflinks_found == False):
                        number_of_selflinks = 1
                        pb_file_path_selflinks = pb_file_path.replace(".pb", "_overlapping_selflinks.pb")
                        pb_file_selflinks = open(pb_file_path_selflinks, "w")
                        pb_file_selflinks.write(generate_pb_line(pbond))
                        selflinks_found = True
                    elif (pbond[0] == pbond[1] and selflinks_found == True):
                        number_of_selflinks += 1
                        pb_file_selflinks.write(generate_pb_line(pbond))

                # Tell the user in which files the pseudobonds from overlapping peptides can be found and close these files
                if (intralinks_found == True and selflinks_found == False):
                    self.display_log_message("These pseudobonds were stored in %s." % pb_file_path_intralinks)
                    pb_file_intralinks.close()
                elif (intralinks_found == True and selflinks_found == True):
                    self.display_log_message("%s of these pseudobonds were self-links.<br>Self-links were stored in %s.<br>The remaining pseudobonds were stored in %s." % (number_of_selflinks, pb_file_path_selflinks, pb_file_path_intralinks))
                    pb_file_intralinks.close()
                    pb_file_selflinks.close()
                elif (intralinks_found == False and selflinks_found == True):
                    self.display_log_message("All of these pseudobonds were self-links, and stored in %s." % pb_file_path_selflinks)
                    pb_file_selflinks.close()


    # Create a trigger handler for change of selection, to change selection of pbonds menu upon (de)selection of
    # pseudobond models outside of the bundle
    def trigger_handler(self):
        from chimerax.core.selection import SELECTION_CHANGED
        from chimerax.core.models import ADD_MODELS, REMOVE_MODELS

        self.triggerset = self.session.triggers             
        self.change_selection_handler = self.triggerset.add_handler(SELECTION_CHANGED, self.selection_handler)
        self.add_model_handler = self.triggerset.add_handler(ADD_MODELS, self.model_handler)
        self.remove_model_handler = self.triggerset.add_handler(REMOVE_MODELS, self.model_handler)


    def selection_handler(self, trigger, trigger_data):
            from chimerax.atomic.pbgroup import PseudobondGroup
            from PyQt5.QtWidgets import QTreeWidgetItemIterator
            from PyQt5.QtCore import Qt

            # Get names of selected pseudobond models
            selected_pb_models = []
            for model in self.session.models:
                if (isinstance(model, PseudobondGroup) and model.selected == True):
                    selected_pb_models.append(model.name)
            # Check the selected models and uncheck the deselected ones
            iterator = QTreeWidgetItemIterator(self.pbonds_menu)
            while iterator.value():
                item = iterator.value()
                if item.text(0) in selected_pb_models:
                    item.setCheckState(0, Qt.Checked)
                else:
                    item.setCheckState(0, Qt.Unchecked)
                iterator += 1


    def model_handler(self, trigger, trigger_data):
        if trigger == "add models":
            self.manage_models(trigger_data)
        else:
            self.remove_item(trigger_data)                


    def cleanup(self):
        self.triggerset.remove_handler(self.change_selection_handler)
        self.triggerset.remove_handler(self.add_model_handler)
        self.triggerset.remove_handler(self.remove_model_handler)
