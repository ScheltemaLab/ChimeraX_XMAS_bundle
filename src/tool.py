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


class TestBundle5(ToolInstance):

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
        self.display_name = "Test Bundle 5"

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
        self.input_file, check = QFileDialog.getOpenFileName(None, "Select an input file", "", "Excel (*.xlsx)")

    def run_button_clicked(self):
        
        # The user has clicked the "Run" button; map crosslinked peptides
        import pandas as pd
        from chimerax.atomic import Structure
        
        # Read the file that the user has selected and extract the peptide pairs from it
        df = pd.read_excel(self.input_file)

        input_peptides_A = df["Sequence A"].tolist()
        input_peptides_B = df["Sequence B"].tolist()

        input_pairs = []

        for i in range(len(df.index)):
            input_pairs.append(sorted([input_peptides_A[i], input_peptides_B[i]]))

        # Deduplicate the peptide pairs ****
        def deduplicate(lst):
            lst.sort()
            last = object()
            for item in lst:
                if item == last:
                    continue
                yield item
                last = item

        input_pairs_deduplicated = list(deduplicate(input_pairs))
        
        # Store all peptide pairs in list of lists.
        # Each peptide is stored as a dictionary, that contains the peptide sequence and the position of the crosslinked residue.
        def create_peptide_dictionary(input_peptide):
            peptide_sequence = input_peptide.replace("[", "").replace("]", "")
            crosslink_position = input_peptide.index("[")
            peptide_dictionary = {"Peptide sequence":peptide_sequence, "Crosslink position":crosslink_position, "Alignments":[]}
            return peptide_dictionary

        peptide_pairs = []

        for i in range(len(input_pairs_deduplicated)):
            peptide_pairs.append([create_peptide_dictionary(input_pairs_deduplicated[i][0]), create_peptide_dictionary(input_pairs_deduplicated[i][1])])

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
                            if (start + peptide["Crosslink position"] not in all_chains[chain]["Unmapped positions"] and chain_sequence[start:end] == peptide_sequence):
                                peptide["Alignments"].append([start + starting_position, end - 1 + starting_position, chain_reference])
        
        # Continue with creating all pseudobonds for all peptide pairs. 
        # If both peptides of a pair are aligned on the same chain, these alignments are checked for overlap. 
        # If they do overlap, the corresponding pseudobond is filtered out.
        pbonds_filtered = []

        for peptide_pair in peptide_pairs:
        # First select the peptide pairs for which both peptides have alignments.
            if (peptide_pair[0]["Alignments"] and peptide_pair[1]["Alignments"]):
                pbonds_unfiltered = [[a, b] for a in peptide_pair[0]["Alignments"] for b in peptide_pair[1]["Alignments"]]
                for pbond in pbonds_unfiltered:
                # Skip overlapping peptides
                    if (pbond[0][2] == pbond[1][2] and max(pbond[0][0], pbond[1][0]) < min(pbond[0][1], pbond[1][1]) + 1):
                        continue
                    # Append the remaining peptides to the list of filtered pseudobonds, in the format that is required for the .pb file
                    pos_A = pbond[0][0] + peptide_pair[0]["Crosslink position"]
                    pos_B = pbond[1][0] + peptide_pair[1]["Crosslink position"]
                    pbonds_filtered.append(sorted([pbond[0][2] + ":" + str(pos_A) + "@CA", pbond[1][2] + ":" + str(pos_B) + "@CA"]))
                    
        # Deduplicate them
        pbonds_filtered_deduplicated = list(deduplicate(pbonds_filtered))
        
        # Store the pseudobonds in a pb file        
        pb_file_name = self.input_file.replace(".xlsx", "_pseudobonds.pb")
        pb_file = open(pb_file_name, "w")
        
        for pbond in pbonds_filtered_deduplicated:
            pb_file.write(pbond[0] + " " + pbond[1] + "\n")
            
        pb_file.close()
        
        # Open this .pb file in ChimeraX
        from chimerax.core.commands import run
        run(self.session, "open %s" % pb_file_name)

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
        clear_action.triggered.connect(lambda *args: self.input_file.clear())
        menu.addAction(clear_action)

    def take_snapshot(self, session, flags):
        # Remembers the input file name when a session is saved
        return {
            "version": 1,
            "current input file": self.input_file
        }

    @classmethod
    def restore_snapshot(class_obj, session, data):
        # When a saved session is reopened
        inst = class_obj(session, "Test Bundle 5")
        inst.input_file = data["current input file"]
        return inst
