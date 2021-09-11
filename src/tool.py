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

        # We will be adding an item to the tool's context menu, so override
        # the default MainToolWindow fill_context_menu method
        self.tool_window.fill_context_menu = self.fill_context_menu

        # Our user interface is simple enough that we could probably inline
        # the code right here, but for any kind of even moderately complex
        # interface, it is probably better to put the code in a method so
        # that this __init__ method remains readable.
        self._build_ui()

    def _build_ui(self):
        # Put our widgets in the tool window
 
        from PyQt5.QtWidgets import QLabel, QGridLayout, QPushButton
        
        layout = QGridLayout()
        
        self.file_button = QPushButton()
        self.file_button.setText("Click to select")
        self.run_button = QPushButton()
        self.run_button.setText("Run")
        
        layout.addWidget(QLabel("Input file:"), 0, 0)
        layout.addWidget(self.file_button, 0, 1)
        layout.addWidget(self.run_button, 1, 1)
        
        # Arrange for our "dialog" method to be called when the "Click to select"
        # button is clicked
        self.file_button.clicked.connect(self.dialog)
        
        # Arrange for our "run_button_clicked" method to be called when the
        # user presses the Return key
        self.run_button.clicked.connect(self.run_button_clicked)

        # Set the layout as the contents of our window
        self.tool_window.ui_area.setLayout(layout)

        # Show the window on the user-preferred side of the ChimeraX
        # main window
        self.tool_window.manage("side")

        
    def dialog(self):
        from PyQt5.QtWidgets import QFileDialog
        self.input_file_path, check = QFileDialog.getOpenFileName(None, "Select an input file", "", "Excel (*.xlsx)")

    
    # Function for displaying log messages
    def display_log_message(self, msg):
        self.session.logger.info(msg, is_html=True)

    
    def run_button_clicked(self):
        
        # The user has clicked the "Run" button; map crosslinked peptides
        import pandas as pd
        from chimerax.atomic import Structure

        # Display first log message
        self.display_log_message("<b>Peptide pair mapping of PD output file: %s</b>" % self.input_file_path)
        
        # Read the file that the user has selected and extract the peptide pairs from it
        dataframe = pd.read_excel(self.input_file_path)
        
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
            if not isinstance(model, Structure):
                continue
            model_id = model.id_string
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
            input_file_name = self.input_file_path.replace(".xlsx", "")      
            pb_file_path = input_file_name + "_pseudobonds.pb"
            pb_file = open(pb_file_path, "w")
            for pbond in pbonds_filtered:
                pb_file.write(generate_pb_line(pbond))
            pb_file.close()
            # Open this .pb file in ChimeraX
            from chimerax.core.commands import run
            run(self.session, "open %s" % pb_file_path)

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
                    pb_file_path_intralinks = input_file_name + "_pseudobonds_overlapping.pb"
                    pb_file_intralinks = open(pb_file_path_intralinks, "w")
                    pb_file_intralinks.write(generate_pb_line(pbond))
                    intralinks_found = True
                elif (pbond[0] != pbond[1] and intralinks_found == True):
                    pb_file_intralinks.write(generate_pb_line(pbond))
                elif (pbond[0] == pbond[1] and selflinks_found == False):
                    number_of_selflinks = 1
                    pb_file_path_selflinks = input_file_name + "_pseudobonds_overlapping_selflinks.pb"
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


    def fill_context_menu(self, menu, x, y):
        # Add any tool-specific items to the given context menu (a QMenu instance).
        # The menu will then be automatically filled out with generic tool-related actions
        # (e.g. Hide Tool, Help, Dockable Tool, etc.) 
        #
        # The x,y args are the x() and y() values of QContextMenuEvent, in the rare case
        # where the items put in the menu depends on where in the tool interface the menu
        # was raised.
        from PyQt5.QtWidgets import QAction
        clear_action = QAction("Clear", menu)
        clear_action.triggered.connect(lambda *args: self.input_file_path.clear())
        menu.addAction(clear_action)


    def take_snapshot(self, session, flags):
        # Remembers the input file name when a session is saved
        return {
            "version": 1,
            "current input file": self.input_file_path
        }


    @classmethod
    def restore_snapshot(class_obj, session, data):
        # When a saved session is reopened
        inst = class_obj(session, "Test Bundle 5")
        inst.input_file_path = data["current input file"]
        return inst
