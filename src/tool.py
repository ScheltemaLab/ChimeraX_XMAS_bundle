# vim: set expandtab shiftwidth=4 softtabstop=4:
# === UCSF ChimeraX Copyright ===
# Copyright 2016 Regents of the University of California.
# All rights reserved. This software provided pursuant to a license 
# agreement containing restrictions on its disclosure, duplication and 
# use. For details see:
# http://www.rbvi.ucsf.edu/chimerax/docs/licensing.html
# This notice must be embedded in or attached to all copies, including 
# partial copies, of the software or any revisions or derivations 
# thereof.
# === UCSF ChimeraX Copyright ===


      
from .info_file import InfoFile
from .matplotlib_venn._venn2 import venn2
from .matplotlib_venn._venn3 import venn3
from .read_evidence import Evidence
from chimerax.atomic.molarray import Pseudobonds
from chimerax.atomic.pbgroup import selected_pseudobonds, PseudobondGroup
from chimerax.atomic.structure import Structure
from chimerax.color_key.model import ColorKeyModel
from chimerax.core.colors import Colormap
from chimerax.core.commands import run
from chimerax.core.models import ADD_MODELS, REMOVE_MODELS
from chimerax.core.selection import SELECTION_CHANGED
from chimerax.core.tools import ToolInstance
from chimerax.ui import MainToolWindow
from chimerax.ui.widgets.color_button import MultiColorButton
import itertools
import matplotlib.pyplot as plt
from matplotlib.text import Text
import numpy as np 
import operator
import os
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QDoubleValidator, QIntValidator
from PyQt5.QtWidgets import (QVBoxLayout, QGridLayout, QHBoxLayout,
                             QTreeWidget, QAbstractItemView, QPushButton,
                             QLabel, QCheckBox, QDialogButtonBox,
                             QLineEdit, QSpacerItem, QSizePolicy, 
                             QTreeWidgetItemIterator, QFileDialog,
                             QTreeWidgetItem, QListWidget, QListWidgetItem,
                             QComboBox, QButtonGroup, QRadioButton, 
                             QStyledItemDelegate)     
from qtrangeslider import QRangeSlider
import re
import seaborn as sns
import string
from venn import venn 

class XMAS(ToolInstance):
    # Inheriting from ToolInstance makes us known to the ChimeraX tool 
    # manager, so we can be notified and take appropriate action when 
    # sessions are closed, saved, or restored, and we will be listed 
    # among running tools and so on.
    
    # Does this instance persist when session closes
    SESSION_ENDURING = False  
    # We do save/restore in sessions  
    SESSION_SAVE = True         
                                

    def __init__(self, session, tool_name):
        # "session"   - chimerax.core.session.Session instance
        # "tool_name" - string

        # Initialize base class
        super().__init__(session, tool_name)

        # Set name displayed on title bar (defaults to tool_name)
        # Must be after the superclass init, which would override it
        self.display_name = "XMAS"

        # Create the main window for our tool. The window object will 
        # have a "ui_area" where we place the widgets composing our 
        # interface. The window isn't shown until its "manage" method 
        # is called.
        self.tool_window = MainToolWindow(self)
        
        self.pb_manager = self.session.pb_manager

        # Call trigger handler to take action when certain triggers fire
        self.trigger_handler()
        
        # Override the "cleanup" method to perform additional actions 
        # upon closing the main tool window
        self.tool_window.cleanup = self.cleanup
        
        cmap = Colormap(None, ((1, 0, 0, 1), (1, 1, 0, 1), (0, 1, 0, 1)))
        self.score_cmap = cmap.linear_range(0, 200)
        
        # Method to build and show the main window
        self._build_ui()


    def _build_ui(self):

        # Put our widgets in the tool window
        
        outer_layout = QVBoxLayout()
        top_layout = QGridLayout()
        pbonds_layout = QVBoxLayout()
        buttons_layout = QHBoxLayout()
        
        # A treewidget that will contain all structural protein models 
        # open in the session
        self.model_selector = QTreeWidget()
        self.model_selector.setHeaderLabels(["Name", "ID"])
        self.model_selector.setColumnWidth(0, 200)

        # A treewidget that will contain PD evidence files that the user
        # has selected
        self.file_selector = QTreeWidget()
        self.file_selector.setHeaderLabels(["Name"])
        self.file_selector.setColumnWidth(0, 200)
        self.file_selector.setSelectionMode(QAbstractItemView.MultiSelection)

        file_buttons = {"Import files":self.dialog,
                        "Remove selected files":self.remove_files}
        row_index= 2
        col_index = 2

        for key in file_buttons:
            button = QPushButton()
            button.setText(key)
            action = file_buttons[key]
            button.clicked.connect(action)
            top_layout.addWidget(button, row_index, col_index)
            col_index += 1

        map_button = QPushButton()
        map_button.setText("Map crosslinks")
        # Upon clicking the map button, the map_button_clicked method is
        # called twice, each time with different arguments
        map_dict = {"model":self.model_selector, "file":self.file_selector}
        for key in map_dict:
            selector = map_dict[key]
            function = lambda _, s=selector, t=key: self.map_button_clicked(s,
                                                                            t)
            map_button.clicked.connect(function)
        
        top_layout.addWidget(QLabel("Molecular models"), 0, 0)
        top_layout.addWidget(QLabel("Evidence files"), 0, 2)
        top_layout.addWidget(self.model_selector, 1, 0, 1, 2)
        top_layout.addWidget(self.file_selector, 1, 2, 1, 2)
        top_layout.addWidget(map_button, 3, 2, 1, 2)

        # In this treewidget, pseudond models from .pb files are shown;
        # both models that are created with XMAS, as well as
        # models that are opened independently of XMAS.
        self.pbonds_menu = QTreeWidget()
        self.pbonds_menu.setHeaderLabels(["Name", "Linked models IDs"])
        self.pbonds_menu.setColumnWidth(0, 300)
        # When a pseudobond model is (de)selected in the menu, it should
        # also be (de)selected in the ChimeraX session. Call
        # "check_signal" method to arrange this        
        self.pbonds_menu.itemChanged.connect(self.check_signal)

        pbonds_layout.addWidget(QLabel("Crosslink models"))
        pbonds_layout.addWidget(self.pbonds_menu)

        buttons_dict = {"Analyze":self.show_analyze_dialog,
                        "Export": self.show_subset_dialog,
                        "Visualize": self.show_visualize_dialog}

        for key in buttons_dict:
            button = QPushButton()
            button.setText(key)
            buttons_layout.addWidget(button)
            function = (lambda _, f=buttons_dict[key]:
                        self.is_selection_empty(f))
            button.clicked.connect(function)

        self.integrate_button = QPushButton("Integrate")
        buttons_layout.addWidget(self.integrate_button)
        self.integrate_button.clicked.connect(self.show_integrate_dialog)

        outer_layout.addLayout(top_layout)
        outer_layout.addLayout(pbonds_layout)
        outer_layout.addLayout(buttons_layout)

        # Set the layout as the contents of our window
        self.tool_window.ui_area.setLayout(outer_layout)

        # Show the tool in the ChimeraX main window. Given its 
        # relatively large size, the tool is not shown on the 
        # user-preferred side of the window, which would have been 
        # arranged by 'self.tool_window.manage("side")'
        self.tool_window.manage("left")

        # Create an empty dictionary where pb models made with Crosslink
        # Mapper will be stored
        self.created_models = {}

        # Add models open in session to the window with the "add_models"
        # method 
        self.add_models(self.session.models)


    def show_integrate_dialog(self):
        
        from .z_score import ZScoreSelector
        
        ZScoreSelector(self)
        
        
    def create_child_window(self, title):
        
        window = self.tool_window.create_child_window(title)
        
        return window


    def show_visualize_dialog(self, pbs):
        
        title = "Visualization settings of pseudobonds in the available set"
        self.visualize_dialog = self.create_child_window(title)
        visualize_dialog = self.visualize_dialog
        visualize_dialog.gradient_active = False

        outer_layout = QVBoxLayout()
        settings_layout = QGridLayout()

        pointers = np.array([pb.cpp_pointer for pb in pbs], dtype=np.uintp)
        pbs = VisualizePseudobonds(pointers)
        
        self.distance_policy = "distance"
        self.score_policy = "score"
        policies = [self.distance_policy, self.score_policy]
        checkboxes_text = ["Gradient", "Cut-off"]
        sliders = list(self.make_sliders(pbs, VisualizeSlider,
                                         self.visualize_dialog))
        visualize_dialog.sliders = dict(zip(policies, sliders))
        visualize_dialog.color_keys = dict(zip(policies, 
                                               [None] * len(policies)))
        checkboxes = [None] * len(policies) * len(checkboxes_text)  
        gradient_group = QButtonGroup(outer_layout)
        gradient_group.setExclusive(False)
        for i in range(len(checkboxes)):
            checkboxes[i] = QCheckBox()
        lower_layout, reset_values = self.create_lower_layout(pbs)
        gradient_group.buttonToggled.connect(lambda button, checked: 
                                             self.uncheck(button, checked, 
                                                          gradient_group, pbs, 
                                                          reset_values))

        row = 0
        for i, policy in enumerate(policies):
            label = "Color by " + policy
            settings_layout.addWidget(QLabel(label), row, 0, 1, 2)
            slider = sliders[i]
            settings_layout.addLayout(slider.layout, row, 2, 2, 1)
            for j, text in enumerate(checkboxes_text):
                index = row + j
                checkbox = checkboxes[index]
                checkbox.setText(text)
                settings_layout.addWidget(checkbox, row + 1, j)
                if not slider.enabled:
                    checkbox.setEnabled(False)
                gradient_box = checkboxes[row]
                gradient_box.policy = policy
                gradient_group.addButton(gradient_box)
                cutoff_box = checkboxes[row + 1]
                cutoff_box.clicked.connect(lambda checked, s=slider:
                                           self.show_slider(checked, s))
            row += 2

        outer_layout.addWidget(QLabel("Value-based coloring:"))    
        outer_layout.addLayout(settings_layout)
        outer_layout.addWidget(QLabel(""))
        outer_layout.addWidget(QLabel("Customized styling:"))
        outer_layout.addLayout(lower_layout)
        visualize_dialog.ui_area.setLayout(outer_layout)
        visualize_dialog.manage("side")
        

    def uncheck(self, button, checked, gradient_group, pbs, reset_values):

        for but in gradient_group.buttons():
            if but == button:
                continue
            if checked:
                but.setChecked(False)
                reset_values = None
            elif (not checked and but.isChecked()):
                reset_values = None

        self.color_gradient(checked, pbs, button.policy, reset_values)


    def create_lower_layout(self, pbs):

        color_button = MultiColorButton(has_alpha_channel=True,
                                        max_size=(16,16))
        color_button2 = MultiColorButton(has_alpha_channel=True,
                                         max_size=(16,16))

        visualize_dialog = self.visualize_dialog

        def color_button_function(button, attribute, row_name, function):
            group = self.pb_manager.get_group("Temp")
            if row_name == visualize_dialog.row_names[0]:
                if group.num_pseudobonds > 0:
                    group.clear()
                for pb in pbs:
                    a1, a2 = pb.atoms
                    new_pb = group.new_pseudobond(a1, a2)
                    new_pb.color = pb.color
                color = group.model_color
                group.delete()
            elif row_name == visualize_dialog.row_names[1]:
                color = visualize_dialog.custom_values[row_name][attribute][0]
            button.set_color(color)
            button.color_changed.connect(lambda rgba: function(rgba, attribute,
                                                               pbs, None,
                                                               row_name))
        
        def radii_edit_function(line_edit, attribute, row_name, function):
            line_edit.setValidator(QDoubleValidator())
            line_edit.setMaximumWidth(50)
            line_edit.textChanged.connect(lambda text: function(text,
                                                                attribute, pbs,
                                                                float,
                                                                row_name))

        def dashes_edit_function(line_edit, attribute, row_name, function):
            line_edit.setValidator(QIntValidator())
            line_edit.setMaximumWidth(50)
            line_edit.textChanged.connect(lambda text: function(text,
                                                                attribute,
                                                                pbs, int,
                                                                row_name))

        layout = QGridLayout()
        labels = ["Color", "Radius", "Dashes (changes complete model(s))"]
        for i, label in enumerate(labels):
            layout.addWidget(QLabel(label), 0, i + 1)                    
        row_names = visualize_dialog.row_names = ["Main", "Cut-off"]
        widget_rows = {
            row_names[0]: [color_button, QLineEdit(), QLineEdit()], 
            row_names[1]: [color_button2, QLineEdit(), QLineEdit()]
            }
        functions = [color_button_function, radii_edit_function,
                     dashes_edit_function]
        attributes = visualize_dialog.attributes = {"color":"colors",
                                                    "radius":"radii",
                                                    "dashes":"dashes"}
        cutoff_values = [[170, 0, 0, 255], 0.5]

        # The next part is inelegant, but at least it works...
                
        custom_values = visualize_dialog.custom_values = {}
        for row_name in widget_rows:
            custom_values[row_name] = {}
            for i, attribute in enumerate(attributes):
                if attribute == "dashes":
                    break
                if row_name == row_names[0]:
                    custom_attribute = attributes[attribute]
                    value = getattr(pbs, custom_attribute)
                else:
                    value = [cutoff_values[i]] * len(pbs)
                custom_values[row_name][attribute] = value

        reset_values = {}
        for i, row_name in enumerate(widget_rows):
            row = i + 1
            layout.addWidget(QLabel(row_name + ":"), row, 0)
            widgets = widget_rows[row_name]
            for j, widget in enumerate(widgets):
                attribute = list(attributes.keys())[j]
                if (row_name == row_names[1] and attribute == "dashes"):
                    break
                col = j + 1
                layout.addWidget(widget, row, col)
                if row_name == row_names[0]:
                    reset_attribute = attributes[attribute]
                    value = getattr(pbs, reset_attribute)
                    reset_values[reset_attribute] = value
                functions[j](widget, attribute, row_name,
                             self.change_pseudobonds_style)

        # "apply" attribute to control whether the changes are accepted
        visualize_dialog.apply = False
        save_box = QCheckBox("Save pb file")
        apply_cancel = QDialogButtonBox(QDialogButtonBox.Apply 
                                        | QDialogButtonBox.Cancel, Qt.Vertical)        
        apply_cancel.clicked.connect(lambda button:
                                     self.apply_clicked(button, "Apply", 
                                                        save_box.checkState(),
                                                        pbs))
        closing_function = lambda: self.reset_style(pbs, reset_values)
        visualize_dialog.cleanup = closing_function
        apply_cancel.rejected.connect(visualize_dialog.destroy)
        
        spacer = QSpacerItem(0, 0, QSizePolicy.Expanding)
        no_cols = layout.columnCount()
        layout.addItem(spacer, 1, no_cols, 1, 2)
        layout.addItem(spacer, 2, no_cols - 1, 1, 2)
        layout.addWidget(save_box, 2, layout.columnCount())
        layout.addWidget(apply_cancel, 1, layout.columnCount(), 2, 1)

        return layout, reset_values


    def change_pseudobonds_style(self, value, attribute, pbs,
                                 convert_function=None, row_name=""):
       
        if str(value) == "":
            return

        elif convert_function is not None:
            value = convert_function(value)
            
        dialog = self.visualize_dialog

        if attribute == list(dialog.attributes.keys())[2]:
            setattr(pbs, attribute, value)
        else:
            row_names = dialog.row_names
            is_cutoff = {row_names[0]: False, row_names[1]: True}
            custom_values = dialog.custom_values
            custom_values[row_name][attribute] = [value] * len(pbs)
            for pb in pbs:
                if not pb.outside_range == is_cutoff[row_name]:
                    continue
                if (attribute == "color" and not pb.outside_range 
                    and dialog.gradient_active):
                    pb.color = pb.gradient_color
                    continue
                setattr(pb, attribute, value)


    def color_gradient(self, checked, pbs, policy, reset_values):

        dialog = self.visualize_dialog        
        row_key = dialog.row_names[0]
        attr_key = list(dialog.attributes.keys())[0]

        if not checked:
            dialog.gradient_active = False
            for i, pb in enumerate(pbs):
                if pb.outside_range:
                    continue
                color = dialog.custom_values[row_key][attr_key][i]
                pb.color = color
            if dialog.color_keys[policy] is None:
                return
            color_key = dialog.color_keys[policy]
            color_key.delete()
            dialog.color_keys[policy] = None
            return
        
        dialog.gradient_active = True
        
        if policy == self.distance_policy:
            rgbas = ((1, 1/3 , 1, 1), (5/6, 1/2, 5/6, 1), (2/3, 2/3, 2/3, 1),
                     (1/3, 2/3, 5/6, 1), (0, 2/3, 1, 1))
            slider = dialog.sliders[policy].slider
            maximum = slider.value()[1]
            color_range = 0, maximum
            attribute = "length"
        elif policy == self.score_policy:
            rgbas = ((1, 0, 0, 1), (1, 1/2, 0, 1), (1, 1, 0, 1),
                     (1/2, 1, 0, 1), (0, 1, 0, 1))
            color_range = 0, 200 
            attribute = "xlinkx_score"

        color_key = ColorKeyModel(self.session)
        color_key._rgbas_and_labels = self.get_rgbas_and_labels(rgbas,
                                                                color_range)
        color_key._ticks = True
        color_key._tick_thickness = 2
        color_key._label_offset = -15
        color_key.name = policy.title() + " gradient"
        self.session.models.add([color_key])
        dialog.color_keys[policy] = color_key
        colors = tuple([rgba for i, rgba in enumerate(rgbas) if i % 2 == 0])
        cmap = Colormap(None, colors)
        self.cmap = cmap.linear_range(min(color_range), max(color_range))
        values = [getattr(pb, attribute) for pb in pbs]
        rgba8_list = self.cmap.interpolated_rgba8(values)

        for i, pb in enumerate(pbs):
            color = rgba8_list[i]
            pb.gradient_color = color            
            if pb.outside_range:
                continue
            pb.color = color

    
    def get_rgbas_and_labels(self, rgbas, color_range):

        minimum = min(color_range)
        maximum = max(color_range)
        sum_values = minimum + maximum
        if not sum_values % 2 == 0:
            maximum += 1
        middle = str(int((minimum + maximum)/2))
        labels = [str(minimum), "", middle, "", str(maximum)]
        rgbas_and_labels = [(rgba, labels[i]) for i, rgba in enumerate(rgbas)]

        return rgbas_and_labels


    def apply_clicked(self, button, text, checkstate, pbs):
        
        if button.text() != text:
            return
        
        if checkstate == Qt.Checked:
            title = "Save pseudobond styles"
            extension = "*.pb"
            file_path, _ = QFileDialog.getSaveFileName(None, 
                                                       title, "", extension)
            if file_path != "":
                models = [group[0].id_string for group in pbs.by_group]
                spec = "#" + ",".join(models)
                for model in models:
                    run(self.session, "save \"%s\" %s" % (file_path, spec))
        
        self.visualize_dialog.apply = True
        self.visualize_dialog.destroy()


    def reset_style(self, pbs, reset_values):
        
        dialog = self.visualize_dialog
        
        if dialog.apply:
            return

        for attribute in reset_values:
            value = reset_values[attribute]
            setattr(pbs, attribute, value)

        color_keys = dialog.color_keys

        for policy in color_keys:
            color_key = color_keys[policy]
            if color_key is None:
                continue
            try:
                color_key.delete()
            except:
                pass

    
    def show_slider(self, checked, slider):

        widgets = slider.widgets
        widget = widgets[0]
        enabled = widget.isEnabled()

        if (enabled and not checked):
            set_enabled = False
            if slider.invert.isChecked():
                slider.invert.click()
        elif (not enabled and checked):
            set_enabled = True
        else:
            return
        
        slider.enabled = set_enabled
        slider.within_range(slider.value_type, None, slider.function)
        
        for widget in widgets:
            widget.setEnabled(set_enabled)


    def remove_files(self):

        iterator = QTreeWidgetItemIterator(self.file_selector,
                                           QTreeWidgetItemIterator.Selected)

        if not iterator.value():
            return

        remove = []

        while iterator.value():
            item = iterator.value()
            remove.append(item)
            del self.evidence_files[item.text(0)]      
            iterator += 1

        root = self.file_selector.invisibleRootItem()      

        for item in remove:
            root.removeChild(item) 


    def find_shortest(self, pbs_dict, names):
        
        if hasattr(self.analyze_dialog, "find_shortest"):
            return
        
        self.analyze_dialog.find_shortest = True
        allow_layout = QHBoxLayout()
        label_front = QLabel("Allow difference of")
        line_edit = QLineEdit()
        line_edit.setText("2")
        line_edit.setValidator(QDoubleValidator(0.0, float("inf"), 1000))
        label_end = QLabel("Ångström")
        widgets = [label_front, line_edit, label_end]
        for widget in widgets:
            allow_layout.addWidget(widget)
        ok = QDialogButtonBox(QDialogButtonBox.Ok)
        ok.accepted.connect(lambda: self.show_shortest(pbs_dict, names,
                                                       line_edit.text()))
        main_layout = self.analyze_dialog.layout
        main_layout.insertLayout(3, allow_layout)
        main_layout.insertWidget(4, ok)

    
    def show_shortest(self, pbs_dict, names, allow_value):

        i = 0

        for model in pbs_dict:

            if not hasattr(model, "XMAS_made"):
                print("Find shortest option unavailable for %s: model has not "
                      "been generated in-session with XMAS" % model.name)
                continue

            name = names[i] + "_shortest.pb"
            group = self.get_pseudobonds_model(name)
            peptide_pairs_dict = {}

            for pb in pbs_dict[model]:
                peptide_pairs = pb.peptide_pairs
                for pair in peptide_pairs:
                    if not pair in peptide_pairs_dict:
                        peptide_pairs_dict[pair] = []
                    peptide_pairs_dict[pair].append(pb)

            shortest_pbs = []

            for pair in peptide_pairs_dict:
                minimum = float("inf")
                pbs = peptide_pairs_dict[pair]              
                for pb in pbs:
                    if pb.length < minimum:
                        minimum = pb.length 
                for pb in pbs:
                    if pb.length > minimum + float(allow_value):
                        continue
                    shortest_pbs.append(pb)

            pbs_atoms_dict = self.pbs_atoms(shortest_pbs, 
                                            operation="find shortest")
            
            for atoms in pbs_atoms_dict:
                pb = group.new_pseudobond(atoms[0], atoms[1])
                xlinkx_scores = pbs_atoms_dict[atoms]
                pb.xlinkx_score = max(xlinkx_scores)

            group_item = self.pbonds_menu.findItems(
                    group.name, Qt.MatchExactly, column=0)[0]
            group_item.setText(1, model.item.text(1))

            i += 1


    def create_venn(self, pbs_dict, names):

        number_of_groups = len(pbs_dict)

        if (number_of_groups == 2 or number_of_groups == 3):
            function = self.venn_matplotlib
        elif (number_of_groups > 3 and number_of_groups < 7):
            function = self.venn_venn
        else:
            print("Venn diagrams can only be plotted for two to six groups")
            return

        sets = self.get_plotting_data(pbs_dict, distance=False)
        
        function(sets, names)
        self.create_plot("Venn diagram")
        plt.show()
        
        # To do: make the labels in the Venn diagram draggable
        # self.drag_handler = DragHandler(plt.gcf())
        
    
    def venn_matplotlib(self, sets, names):
        
        if len(sets) == 2:
            function = venn2
        else:
            function = venn3
            
        function(sets, tuple(names))
        
    
    def venn_venn(self, sets, names):
        
        data = {}
        
        for i, name in enumerate(names):
            data[name] = sets[i]
            
        venn(data)


    def create_distance_plot(self, pbs_dict, names):

        distances = self.get_plotting_data(pbs_dict)    
        
        plt.figure()
        self.create_plot("Distance plot")
        sns.set(style="whitegrid")
        sns.boxenplot(data=distances, orient="h")
        sns.stripplot(data=distances, orient="h", size=4, color="gray")
        plt.xlabel("Distance (Å)")
        plt.yticks(plt.yticks()[0], names)
        plt.tight_layout()
        plt.show()
          
        
    def create_plot(self, title):

        plt.ion()
        manager = plt.get_current_fig_manager()
        manager.set_window_title(title)


    def get_plotting_data(self, pbs_dict, distance=True):

        values = [None] * len(pbs_dict)

        for i, model in enumerate(pbs_dict):
            pbs = pbs_dict[model]
            if distance:
                values[i] = [pb.length for pb in pbs]
            else:
                values[i] = set(pb.string() for pb in pbs)

        return values
    
    
    def update_distances(self, pbs_dict, _):
        
        for model in pbs_dict:
            pbs = pbs_dict[model]
            files = set()
            for pb in pbs:
                try:
                    file = pb.info_file
                except:
                    continue
                file = pb.info_file
                files.add(file)
                distance = pb.length
                indices = pb.indices
                for index in indices:
                    file.df.at[index, "Distance (A)"] = distance
            for file in files:
                file.create_file()
                print("Distances updated in %s" % file.path)

    
    def dialog(self):

        # The user has clicked the file button: with this dialog, 
        # evidence files can be selected

        file_dialog = QFileDialog()
        # Set file mode to existing files to enable selection of multiple files
        file_dialog.setFileMode(QFileDialog.ExistingFiles)
        get_file_names = QFileDialog.getOpenFileNames
        selected_files, _ = get_file_names(None, "Select evidence files", "", 
                                           "Evidence File (*.csv *.mzid *.tsv "
                                           "*.txt *.xlsx)")

        # Create a dictionary containing the evidence files to avoid 
        # showing the same file multiple times
        if not hasattr(self, "evidence_files"):
            self.evidence_files = {}
        for selected_file in selected_files:
            if selected_file in self.evidence_files.values():
                selected_files.remove(selected_file)
        number_of_files = len(selected_files)
        # Store the number of files in the treewidget to give each file 
        # a proper ID
        if not hasattr(self, "number_of_items"):
            self.number_of_items = 0
        i = 0
        # Files get a unique letter sequence as ID, which is created by 
        # the "file_id_generator" method, and are added to the 
        # treewidget
        for ID in self.file_id_generator(
            self.number_of_items, self.number_of_items + number_of_files):
            path_to_file = selected_files[i]
            short_name = self.get_short_filename(path_to_file)
            item = QTreeWidgetItem(self.file_selector)
            item.setText(0, short_name)
            item.setCheckState(0, Qt.Checked)
            i += 1
            self.evidence_files[short_name] = path_to_file
            self.number_of_items += 1


    def file_id_generator(self, start, end):

        # Creates file IDs for user-selected evidence files. The start 
        # index is specified, so that the generator does not start over 
        # every time new files are selected. The end index is pecified, 
        # so that the number of IDs created corresponds to the number 
        # of newly added files.

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

        # Method to obtain a short name to represent a file or model in
        # the tool

        short_name = os.path.basename(path_to_file)

        return short_name


    def map_button_clicked(self, selector, selector_type):

        # The user has clicked the map button: the selected models and 
        # files are extracted from their respective selector treewidgets

        iterator = QTreeWidgetItemIterator(
            selector, QTreeWidgetItemIterator.Checked)
        # If no models/files have been selected, print a message
        if not iterator.value():
            print("Please check one or more %ss" % selector_type)
            self.missing_data = True
        else:
            checked_items = []
            while iterator.value():
                item = iterator.value()
                if selector_type == "model":
                    ID = item.text(1)
                elif selector_type == "file":
                    ID = item.text(0)
                iterator += 1
                checked_items.append(ID)
            if selector_type == "model":
                self.checked_models = checked_items
                self.missing_data = False
            elif selector_type == "file":
                checked_files = [None] * len(checked_items)
                for i, item in enumerate(checked_items):
                    checked_files[i] = self.evidence_files[item]
                # If no models have been selected, mapping should not
                # be performed. Therefore, "missing_data" has to be
                # False. Mapping is performed by a call to the
                # "map_crosslinks" method
                if not self.missing_data:
                    self.map_crosslinks(self.checked_models, checked_files)
    

    def map_crosslinks(self, checked_models, checked_files):
        
        # The user has selected one or multiple models and evidence
        # files, and clicked the map button: map crosslinked peptides

        # Each checked file is mapped to all checked models
        for j, evidence_file in enumerate(checked_files):
            
            # Read the file and extract the peptide pairs and search engine 
            # from it
            evidence = Evidence(evidence_file)
            input_pairs = evidence.peptide_pairs
            engine = evidence.engine

            # Display bold log message to signify which file is being
            # mapped
            self.session.logger.info(
                "<br><b>Peptide pair mapping of %s evidence file: %s</b>" 
                % (engine, evidence_file), is_html=True)
                    
            # Create a file for reference of the results to the evidence file
            # Create a code for the model with all model IDs and
            # the evidence file ID that was used
            model_ids = ",".join([
                str(model_id) for model_id in checked_models
                ])
            self.file_code = model_ids
            info_file_path = (os.path.splitext(evidence_file)[0] 
                              + "_%s.tsv" % self.file_code)
            self.info_file = InfoFile(info_file_path, engine)

            # Keep track of peptide pairs with lacking sequence
            # information
            sequence_info_lacking = 0
            
            complete_input_pairs = []

            for peptide_pair in input_pairs:
                if (peptide_pair.SequenceA == "" 
                        or peptide_pair.SequenceB == ""):
                    sequence_info_lacking += 1
                    self.info_file.add(peptide_pair.Ref, 
                                       "Sequence lacking/decoy")
                    continue
                complete_input_pairs.append(peptide_pair)
                
            # Print a message when one or multiple peptide pairs lack
            # sequence information
            if sequence_info_lacking == 1:
                print("1 Peptide pair is disregarded due to lacking/decoy "
                      "sequence")
            elif sequence_info_lacking >= 1:
                print("%s Peptide pairs are disregarded due to"
                      % str(sequence_info_lacking),
                      "lacking/decoy sequence")
            
            # Make sure that the highest score is taken in case of duplicates
            self.duplicate_scores = {}
            compare_function = self.advanced_equality_check
            peptide_pairs = list(self.deduplicate(complete_input_pairs,
                                                       compare_function))

            # Print a message stating how many peptide pairs were unique
            number_of_deduplicated = len(peptide_pairs)
            print("Unique peptide pairs: %s out of %s" 
                % (number_of_deduplicated, len(complete_input_pairs)))

            for i in range(number_of_deduplicated):
                peptide_pair = peptide_pairs[i]
                ref = peptide_pair.Ref
                if ref not in self.duplicate_scores.keys():
                    continue
                scores = self.duplicate_scores[ref]
                peptide_pair.Score = self.map_max_score(scores)
                
            self.align_peptides(peptide_pairs, checked_models)
            
            
    def align_peptides(self, peptide_pairs, checked_models):

        # Now we align all peptides to the sequences of all chains
        # open in the ChimeraX session
        #
        # The peptide sequences are compared to the sequences
        # string of all chains. The index of a crosslinked residue
        # on the sequence string is always different than its
        # "number" attribute of the corresponding Residue object in
        # ChimeraX, since numbering Residue objects commences with
        # the number 1 or higher, instead of 0. To ensure that the
        # right residue number is used in the .pb file, the first
        # residue number of each chain is stored, which is then
        # added to the index found on the sequence string.
        # Furthermore, some residues are present in the sequence
        # string, but their Residue objects are not shown in the
        # Chimerax structure. These Residue objects are of type
        # NoneType, and ChimeraX does not enable creating
        # pseudobonds between NoneType residues. To prevent adding
        # these crosslinks to the .pb file, the positions of all
        # NoneType residues is also stored.

        for model in self.session.models:
            model_id = model.id_string
            if model_id not in checked_models:
                continue
            for chain in model.chains:
                # Keep track of the number of NoneType residues at
                # the start of the sequence. These will influence
                # the first residue number         
                preceding_nonetypes = 0
                first_residue_number_found = False
                chain_sequence = chain.characters
                nonetype_positions = []
                residues = chain.residues
                for i, residue in enumerate(residues):
                    if (residue is None 
                            and not first_residue_number_found):
                        preceding_nonetypes += 1
                        nonetype_positions.append(i)
                    elif (residue is None 
                            and first_residue_number_found):
                        nonetype_positions.append(i)
                    elif (residue is not None 
                            and not first_residue_number_found):
                        first_residue_number = (
                            residue.number - preceding_nonetypes)
                        first_residue_number_found = True
                # Loop all peptide sequences over the chain to find
                # perfect alignments
                for peptide_pair in peptide_pairs:
                    peptide_sequences = peptide_pair.get_info()
                    crosslink_positions = peptide_pair.get_info("XLinkPosition")
                    letters = ["A", "B"]
                    for i, peptide_sequence in enumerate(peptide_sequences):
                        peptide_length = len(peptide_sequence)
                        for start in range(
                                len(chain_sequence) - peptide_length + 1):
                            end = start + peptide_length
                            # If the crosslinked residue is not
                            # present in the structure, the
                            # pseudobond cannot be mapped, and
                            # therefore we will disregard this
                            # alignment                                
                            crosslink_position = (
                                start + crosslink_positions[i])
                            if crosslink_position in nonetype_positions:
                                continue
                            if (chain_sequence[start:end] 
                                    == peptide_sequence):
                                alignment = Alignment(
                                    start, end, first_residue_number, 
                                    crosslink_position, model_id,
                                    chain)
                                alignments = getattr(peptide_pair, 
                                                     "Alignments" 
                                                     + letters[i])
                                alignments.append(alignment)  
                                
        self.create_pseudobonds(peptide_pairs)
            
            
    def create_pseudobonds(self, peptide_pairs):

        # Continue with creating all valid pseudobonds for all
        # peptide pairs and store them in a list. If both peptides
        # of a pair are aligned on the same chain, these alignments
        # are checked for overlap. Make separate lists for peptide
        # pairs with non-overlapping peptides and those with
        # overlapping peptides. Overlapping peptides can be
        # categorized as nonself-links and selflinks

        pbonds_unfiltered = []
        pbonds = []

        # The number of perfectly aligned peptide pairs is counted, as well
        # as the number of overlapping peptide pairs
        number_of_aligned_pairs = 0
        no_overlapping = 0

        for peptide_pair in peptide_pairs:
        # First select the peptide pairs for which both peptides
        # have alignments. Each possible pseudobond is stored as a
        # PrePseudobond object  
            ref = peptide_pair.Ref
            if not (len(peptide_pair.AlignmentsA) > 0 
                    and len(peptide_pair.AlignmentsB) > 0):
                self.info_file.add(ref, "No pseudobonds")
                continue
            number_of_aligned_pairs += 1   
            pbonds_unfiltered = [PrePseudobond(a, b, peptide_pair) 
                                 for a in peptide_pair.AlignmentsA 
                                 for b in peptide_pair.AlignmentsB]

            if len(pbonds_unfiltered) == 0:
                continue      

            # Check whether the peptide pair has pseudobonds for 
            # overlapping peptides
            peptide_pair.has_overlapping = False       
            for pb in pbonds_unfiltered:
                if pb.is_overlapping:
                    peptide_pair.has_overlapping = True
                    break
         
            # Then make a list for pseudobonds that need to be drawn, and
            # make lines for the info file
            for pb in pbonds_unfiltered:
                line = pb.line
                if (not pb.is_overlapping and not pb.is_selflink):
                    pbonds.append(pb)
                    continue
                elif (pb.is_overlapping and not pb.is_selflink):
                    self.info_file.add(ref, line, 
                                       "Overlapping (non-self)")
                elif pb.is_selflink:  
                    self.info_file.add(ref, line, 
                                       "Overlapping (self)")
                no_overlapping += 1
        
        # Print a log message stating for how many peptide pairs perfect 
        # alignments have been found
        print("Unique peptide pairs with pseudobonds: %s" 
            % number_of_aligned_pairs)
        
        info_file_path = self.info_file.path
        
        self.create_files(pbonds, info_file_path)  

        # Print a log message stating where the mapping info is stored
        print("Mapping information is stored in %s" % info_file_path)
            
            
    def create_files(self, pbonds, info_file_path):
        if len(pbonds) > 0:
            pb_file_path = os.path.splitext(info_file_path)[0] + ".pb"
            # Create a new pseudobonds model
            created_model = self.create_pseudobonds_model(pbonds,
                                                          pb_file_path)
            # Store the model and its code in the "created_models"
            # dictionary
            if created_model not in self.created_models.keys():
                self.created_models[created_model] = self.file_code
            # Show the code of this file in the pbonds menu                
            item = self.pbonds_menu.findItems(
                created_model, Qt.MatchExactly, column=0)[0]
            item.setText(1, self.file_code)
        else:
            self.info_file.create_file()
            
            
    def map_max_score(self, scores):
        
        max_score = 0
        
        for score in scores:
            if float(score) <= max_score:
                continue
            max_score = score
            
        return max_score

                   
    def create_pseudobonds_model(self, pbonds, file_path):
        
        if file_path is None:
            title = "Save pseudobonds"
            extension = "*.pb"
            file_path, _ = QFileDialog.getSaveFileName(None, 
                                                       title, "", extension)
            if file_path == "":
                return
            operation = "export"
        else:
            operation = "map"
            
        name = self.get_short_filename(file_path)
        group = self.get_pseudobonds_model(name)
        group.XMAS_made = True

        pbs_atoms_dict = self.pbs_atoms(pbonds, operation)
            
        for atoms in list(pbs_atoms_dict.keys()):
            new_pb = group.new_pseudobond(atoms[0], atoms[1])
            line = new_pb.line = self.create_pb_line(new_pb)
            peptide_pairs = new_pb.peptide_pairs = pbs_atoms_dict[atoms]
            new_pb.info_file = self.info_file
            new_pb.indices = [None] * len(peptide_pairs)
            distance = new_pb.length
            max_score = 0
            has_peptide_pairs = False
            for i, peptide_pair in enumerate(peptide_pairs):
                if peptide_pair is None:
                    continue
                has_peptide_pairs = True
                if peptide_pair.has_overlapping:
                    cat = "Overlap associated"
                else:
                    cat = "Not overlap associated"
                index = self.info_file.add(peptide_pair.Ref, line, cat, 
                                           distance)
                new_pb.indices[i] = index
                score = peptide_pair.Score
                if score <= max_score:
                    continue
                max_score = score
            if has_peptide_pairs:
                new_pb.xlinkx_score = max_score
            
        self.write_file(file_path, group, file_type=".pb")
        print("Pseudobonds are stored in %s" % file_path)
        
        if operation == "map":
            self.info_file.create_file()
        else:
            item = self.pbonds_menu.findItems(name, Qt.MatchExactly, 
                                              column=0)[0]
            item.setText(1, self.get_ids(group))

        return group.name


    def get_pseudobonds_model(self, name):
        
        group = self.pb_manager.get_group(name)
        if group.num_pseudobonds > 0:
            extension = re.search("\(\d+\)", name)
            if extension is None:
                extension = ""
                n = 1
            else:
                extension = extension.group(0)
                n = re.search("(\d+)", extension) + 1
            new_extension = "(%s).pb" % n
            name = name.replace(extension + ".pb", new_extension)
            group = self.pb_manager.get_group(name)
        self.session.models.add([group])
        group.radius = 0.5
        group.color = [255, 255, 0, 255]
        group.dashes = 8

        return group


    def create_pb_line(self, pb):

        atom1, atom2 = pb.atoms
        atom1_string = atom1.string(style="command line", omit_structure=False)
        atom2_string = atom2.string(style="command line", omit_structure=False)
        atoms_sorted = sorted([atom1_string, atom2_string])
        pb_line = atoms_sorted[0] + " " + atoms_sorted[1]

        return pb_line

    
    def pbs_atoms(self, pbs, operation="map"):

        atom_dict = {}    
        
        for i, pb in enumerate(pbs):
            if operation == "map":
                atom1, atom2 = pb.atom1, pb.atom2
            else:
                atom1, atom2 = pb.atoms
            atoms = sorted([atom1, atom2])
            atoms = tuple(atoms)
            if atoms not in atom_dict.keys():
                atom_dict[atoms] = set()
            if operation == "find shortest":
                atom_dict[atoms].add(pb.xlinkx_score)
            elif operation == "map":
                atom_dict[atoms].add(pb.peptide_pair)
            else:
                try:
                    values = pb.peptide_pairs
                    for value in values:
                        atom_dict[atoms].add(value)
                except:
                    atom_dict[atoms].add(None)                  

        return atom_dict


    def get_color(self, color_policy, value, radio_button, cmap):

        if (value > 200
            and color_policy == "score"
            and radio_button == "Gradient"):
                value = 200  

        if radio_button == "Gradient":
            rgba8_list = cmap.interpolated_rgba8([value])
            rgba8 = rgba8_list[0]
        elif radio_button == "Cut-off":
            minimum_value = self.cutoff_values["Min"][1].text()
            if minimum_value == "":
                minimum = 0
            else:
                minimum = float(minimum_value)
            maximum_value = self.cutoff_values["Max"][1].text()
            if maximum_value == "":
                maximum = float("inf")
            else:
                maximum = float(maximum_value)
            
            if (value < minimum or value > maximum):
                rgba8 = [255, 0, 0, 255]
            else:
                rgba8 = [255, 255, 0, 255]

        return rgba8


    def write_file(self, file_path, group, file_type=".pb"):

        # Write a file

        if file_type == ".pb":
            pbs = group.pseudobonds
        else:
            pbs = group

        lines = [None] * len(pbs)
        if file_type == ".pb":
            for i, pb in enumerate(pbs):
                pb = pbs[i]
                lines[i] = pb.line        
        elif file_type == "export":
            for i, pb in enumerate(pbs):
                lines[i] = pb

        lines_deduplicated = list(self.deduplicate(lines, 
                                                   self.simple_equality_check))

        created_file = open(file_path, "w")
        for line in lines_deduplicated:
            created_file.write(line + "\n")
        created_file.close()
        
        if (file_type == "export" and file_path[-3:] == ".pb"):
            text = "open \"%s\"" % file_path
            run(self.session, text, log = False)

        return len(lines_deduplicated)
        

    def deduplicate(self, lst, function):

        # Method that removes duplicate items from a list
        if function == self.simple_equality_check:
            lst.sort()
        else:
            lst.sort(key=operator.attrgetter("SequenceA", "SequenceB"))
        last = object()
        for item in lst:
            if function(item, last):
                continue
            yield item
            last = item
            
            
    def simple_equality_check(self, item, last):
        
        return item == last
    
    
    def advanced_equality_check(self, item, last):
        
        is_equal = False
        
        # The first 'last' variable is an empty object. Hence, it has no 
        # attributes
        try:
            last.SequenceA
        except:
            return False
        
        item_peptides = item.get_info()
        last_peptides = last.get_info()
        
        if item_peptides == last_peptides:
            is_equal = True
            ref = last.Ref
            if ref not in self.duplicate_scores.keys():
                self.duplicate_scores[ref] = [last.Score]
            self.duplicate_scores[ref].append(item.Score)
            self.info_file.add(item.Ref, "Duplicate of %s" % ref)
            
        return is_equal      


    def check_signal(self, item, column):

        # (De)select models depending on the item (de)selected in the
        # pbonds menu

        for model in self.session.models:
            if (model.name == item.text(column) 
                    and item.checkState(column) == Qt.Checked):
                model.selected = True
            elif (model.name == item.text(column)
                    and item.checkState(column) == Qt.Unchecked):
                model.selected = False

    def is_selection_empty(self, function):

        selection = self.get_pseudobonds()

        if selection is None:
            print("Please select pseudobonds")
        else:
            function(selection)


    def show_analyze_dialog(self, pbs):

        title = "Analyze pseudobonds in available set"        
        self.analyze_dialog = self.create_child_window(title)

        pbs_dict = self.get_pseudobonds_dictionary(pbs)

        layout = self.analyze_dialog.layout = QVBoxLayout()
        buttons_dict = {"Plot overlap": self.create_venn,
                        "Plot distances": self.create_distance_plot,
                        "Find shortest": self.find_shortest,
                        "Update distances": self.update_distances}

        names_menu = QTreeWidget()
        names_menu.setHeaderLabels(["Model name", "Chosen name"])
        edit_class = NoEditDelegate
        names_menu.setItemDelegateForColumn(0, 
                                            edit_class())
        
        # If none of the pseudobonds are linked to an info file, the
        # "Update distances" button should be disabled
        disable_update = True
        for pb in pbs:
            if not hasattr(pb, "info_file"):
                continue
            disable_update = False
            break
        
        for key in buttons_dict:
            button = QPushButton()
            button.setText(key)
            layout.addWidget(button)
            function = buttons_dict[key]
            button.clicked.connect(lambda _, f=function:
                                   self.get_names(names_menu, pbs_dict, f))
            if (key == "Update distances" and disable_update):
                button.setEnabled(False)
            
        for model in pbs_dict:
            item = QTreeWidgetItem(names_menu)
            item.setFlags(item.flags() | Qt.ItemIsEditable)
            item.setText(0, model.name)
            text = os.path.splitext(model.name)[0]
            item.setText(1, text)

        layout.insertWidget(3, QLabel("")) 
        layout.insertWidget(4, QLabel("Customize names:"))
        layout.insertWidget(5, names_menu)
        layout.insertWidget(6, QLabel("")) 

        self.analyze_dialog.ui_area.setLayout(layout)
        self.analyze_dialog.manage(None)


    def get_names(self, treewidget, pbs_dict, function):

        iterator = QTreeWidgetItemIterator(treewidget)
        names = [None]*len(pbs_dict)

        i = 0

        while iterator.value():
            item = iterator.value()
            names[i] = item.text(1)
            i += 1
            iterator += 1

        function(pbs_dict, names)
                

    def show_subset_dialog(self, pbs):

        # Through this dialog, the user can export a subset of the
        # models selected in the pbonds menu

        layout = QGridLayout()
        
        title = "Export subset of pseudobonds in the available set"
        self.subset_dialog = self.create_child_window(title)

        self.dialog_model_selector = QTreeWidget()
        self.dialog_model_selector.setHeaderLabels(["Name", "ID"])

        models = pbs.unique_structures

        for model in models:
            item = QTreeWidgetItem(self.dialog_model_selector)
            item.setText(0, model.name)
            item.setText(1, model.id_string)
            item.setCheckState(0, Qt.Checked)
            item.setFlags(item.flags() & ~Qt.ItemIsSelectable)
            item.model = model

        self.dialog_model_selector.sortItems(1, Qt.AscendingOrder)

        # Menu to select whether only intralinks, only interlinks, or
        # both need to be exported
        self.link_selector = QListWidget()
        link_types = ["Intralinks", "Chain interlinks", "Model interlinks"]
        for link_type in link_types:
            item = QListWidgetItem(self.link_selector)
            item.setText(link_type)
            item.setCheckState(Qt.Checked)
            item.setFlags(item.flags() & ~Qt.ItemIsSelectable)
        texts = []
        # If pseudobonds all belong to the same model, Model interlinks should
        # be disabled.
        if len(models) == 1:
            texts.append("Model")
        # If none of the models has multiple chains, Chain interlinks should be
        # disabled.
        multiple_chains = False
        for model in models:
            if len(model.chains) > 1:  
                multiple_chains = True 
                break
        if not multiple_chains:
            texts.append("Chain")
        for text in texts:
            item = self.link_selector.findItems("%s interlinks" % text,
                                                Qt.MatchExactly)[0]
            item.setCheckState(Qt.Unchecked)
            item.setFlags(Qt.NoItemFlags)     

        sliders = self.make_sliders(pbs, ExportSlider)
        self.distance_slider, self.score_slider = sliders
        
        label_front = QLabel("Present in at")
        self.least_or_most = QComboBox()
        self.least_or_most.addItems(["least", "most"])
        self.overlap_number = QComboBox()
        number_of_groups = len(pbs.by_group)
        overlap_list = [str(i) for i in range(1, number_of_groups + 1)]
        self.overlap_number.addItems(overlap_list)
        label_end = QLabel("out of %s pseudobond models"
                           % str(number_of_groups))
        widgets = [label_front, self.least_or_most, self.overlap_number,
                   label_end]

        overlap_layout = QHBoxLayout()
        for widget in widgets:
            overlap_layout.addWidget(widget)
        
        export_button = QPushButton("Export subset")

        self.checkbox_pb = QCheckBox(".pb")
        self.checkbox_pb.setChecked(True)
        self.checkbox_disvis = QCheckBox("DisVis")
        self.checkbox_haddock = QCheckBox("HADDOCK")
        checkboxes = {"Pb":self.checkbox_pb, "DisVis":self.checkbox_disvis,
                      "HADDOCK":self.checkbox_haddock}
        
        checkbox_layout = QHBoxLayout()
        for checkbox in checkboxes:
            current = checkboxes[checkbox]
            # Until this can be implemented:
            if current.text() == "HADDOCK":
                break
            checkbox_layout.addWidget(current)

        layout.addWidget(QLabel("Molecular models:"), 0, 0)
        layout.addWidget(self.dialog_model_selector, 1, 0)
        layout.addWidget(QLabel("Links:"), 0, 1)
        layout.addWidget(self.link_selector, 1, 1)
        layout.addWidget(QLabel(""), 2, 0)
        layout.addLayout(self.distance_slider.layout, 3, 0, 2, 2)
        layout.addLayout(self.score_slider.layout, 5, 0, 2, 2)
        layout.addWidget(QLabel(""), 7, 0)
        layout.addWidget(QLabel("Quantity:"), 8, 0)
        layout.addLayout(overlap_layout, 9, 0)
        layout.addWidget(QLabel(""), 10, 0)
        layout.addLayout(checkbox_layout, 11, 0)       
        layout.addWidget(export_button, 11, 1)

        self.subset_dialog.ui_area.setLayout(layout)
        
        export_button.clicked.connect(lambda: self.export_subset(pbs,
                                                                 checkboxes))

        self.subset_dialog.manage(None)
        self.subset_dialog.cleanup = lambda: self.display_all(pbs)


    def make_sliders(self, pbs, cls=None, dialog=None):

        maximum_distance = self.max_pseudobond_distance(pbs)
        distance_slider = cls("distance", True, 0, maximum_distance, pbs,
                              dialog)
        make_score_slider = True
        for group in pbs.by_group:
            pb = group[1][0]
            if not hasattr(pb, "xlinkx_score"):
                make_score_slider = False
                break
        if make_score_slider:
            maximum_score = self.get_maximum_score(pbs)
            score_slider = cls("score", True, 0, maximum_score, pbs, dialog)
        else:
            score_slider = cls("score", False)
        
        sliders = [distance_slider, score_slider]
        for i, slider in enumerate(sliders):
            slider.linked_slider = sliders[i - 1]

        return distance_slider, score_slider

    
    def get_maximum_score(self, pbs):

        maximum = 0

        for pb in pbs:
            if pb.xlinkx_score <= maximum:
                continue
            maximum = pb.xlinkx_score

        return int(maximum) + 1


    def max_pseudobond_distance(self, pbs):

        maximum = 0

        for pb in pbs:
            length = pb.length
            if length > maximum:
                maximum = length

        maximum = int(maximum) + 1

        return maximum


    def get_pseudobonds(self):

        all_selected_pbs = selected_pseudobonds(self.session)

        if len(all_selected_pbs) == 0:
            return None

        by_group = all_selected_pbs.by_group

        for group in by_group:
            (model, pbs) = group
            if not model.name.endswith(".pb"):
                continue
            not_selected = []
            for pb in pbs:
                if pb in all_selected_pbs:
                    continue
                not_selected.append(pb)
            for pb in not_selected:
                all_selected_pbs.remove(pb)
        
        return all_selected_pbs


    def get_pseudobonds_dictionary(self, pbs):

        pbs_dict = {}

        for group in pbs.by_group:
            (model, pbs) = group
            pbs_dict[model] = pbs

        return pbs_dict


    def export_subset(self, pseudobonds, checkboxes):

        checked = True
        
        iterator_class = QTreeWidgetItemIterator
        model_iterator = iterator_class(self.dialog_model_selector,
                                        iterator_class.Checked)
        if not model_iterator.value():
            print("Please check one or more models")
            checked = False

        links = []
        for i in range(self.link_selector.count()):
            item = self.link_selector.item(i)
            if item.checkState() == Qt.Unchecked:
                continue
            links.append(item.text())
        if len(links) == 0:
            print("Please check one or more links")
            checked = False

        if not checked:
            return

        models = []
        valid_pseudobonds = []
        
        while model_iterator.value():
            item = model_iterator.value()
            models.append(item.model)
            model_iterator += 1

        for pb in pseudobonds:
            # Ignore pseudobonds that belong to unchecked models
            in_models = True
            pb.models = set()
            a1, a2 = pb.atoms
            atoms = [a1, a2]
            for atom in atoms:
                pb.models.add(atom.structure)
            for model in pb.models:
                if model not in models:
                    in_models = False
                    break
            if not in_models:
                continue
            # Ignore pseudobonds that have an incorrect link type
            if len(pb.models) > 1:
                link_type = "Model interlinks"
            else:
                pb_chains = set()
                for atom in atoms:
                    chain = atom.residue.chain
                    pb_chains.add(chain)
                if len(pb_chains) > 1:
                    link_type = "Chain interlinks"
                else:
                    link_type = "Intralinks"
            if not link_type in links:
                continue
            # Ignore pseudobonds that have the wrong length or score
            if pb.outside_range:
                continue
                
            # Add the pseudobonds that meet all criteria        
            valid_pseudobonds.append(pb)
        
        # Check for overlap
        number_for_overlap = int(self.overlap_number.currentText())
        least_or_most = self.least_or_most.currentText()
        pb_strings = [pb.string() for pb in valid_pseudobonds]
        remove = []
        if least_or_most == "least":
            # operator.ge corresponds to ">="
            op = operator.ge
        elif least_or_most == "most":
            # operator.le corresponds to "<="
            op = operator.le
        for pb_string in pb_strings:
            if op(pb_strings.count(pb_string), number_for_overlap):
                continue
            remove.append(pb_string)
        for pb_string in remove:
            for pb in valid_pseudobonds:
                if pb.string() != pb_string:
                    continue
                valid_pseudobonds.remove(pb)
                break

        if len(valid_pseudobonds) == 0:
            print("No pseudobonds match the criteria")
            return

        if checkboxes["Pb"].isChecked():
            pb_lines = []
            for pb in valid_pseudobonds:
                if hasattr(pb, "line"):
                    pb_line = pb.line
                else:
                    pb_line = self.create_pb_line(pb)
                pb_lines.append(pb_line)
            if not checkboxes["DisVis"].isChecked():
                self.create_pseudobonds_model(valid_pseudobonds, None)
                self.subset_dialog.destroy()
            else:
                self.show_disvis_dialog(models, valid_pseudobonds, pb_lines)

        elif (not checkboxes["Pb"].isChecked()
              and checkboxes["DisVis"].isChecked()):
            self.show_disvis_dialog(models, valid_pseudobonds)


    def display_all(self, pseudobonds):

        for pb in pseudobonds:
            pb.display = True


    def save_subset(self, title, extension, lists):

        file_path, _ = QFileDialog.getSaveFileName(None, title, 
                                                   "", extension)

        if file_path == "":
            return

        file_path = os.path.splitext(file_path)[0]
        extensions = [".pb", ".txt"]

        for i, lst in enumerate(lists):
            if lst is None:
                continue
            current_path = file_path + extensions[i]
            self.write_file(current_path, lst, file_type="export")


    def show_disvis_dialog(self, models, pseudobonds, pb_lines=None):

        self.disvis_dialog = self.create_child_window("Create DisVis input "
                                                      "file")
        
        models_layout = QHBoxLayout()
        fix_layout = QVBoxLayout()
        scanning_layout = QVBoxLayout()
        distances_layout = QGridLayout()
        outer_layout = QVBoxLayout()

        chain_selection = {"Fix chain:": fix_layout,
                           "Scanning chain:": scanning_layout}
        groups = [None] * len(chain_selection)

        i = 0  
        for chain in chain_selection:
            layout = chain_selection[chain]
            layout.addWidget(QLabel(chain))
            group = QButtonGroup(outer_layout)
            for model in models:
                model_string = model.name + " (" + model.id_string + ")"
                button = QRadioButton(model_string)
                button.model = model
                group.addButton(button)
                layout.addWidget(button)
            group.buttons()[0].toggle()   
            groups[i] = group            
            i += 1 
            models_layout.addLayout(layout)
        
        distance_options = ["Minimum", "Maximum"]
        line_edits = [None] * len(distance_options)
        for i, distance in enumerate(distance_options):
            label = QLabel(distance + " distance:")
            distances_layout.addWidget(label, i, 0)
            line_edit = QLineEdit()
            line_edit.setValidator(QDoubleValidator(0.0, float("inf"), 1000))
            line_edits[i] = line_edit
            distances_layout.addWidget(line_edit, i, 1)

        def get_parameters():

            chains_keys = ["Fixed", "Scanning"]
            chains = {}
            for i, group in enumerate(groups):
                model = group.checkedButton().model
                key = chains_keys[i]
                chains[key] = model    

            distances = {}
            for i, distance in enumerate(distance_options):
                value = line_edits[i]
                distances[distance] = value

            self.create_disvis_input(chains, distances, pseudobonds, pb_lines)
            self.subset_dialog.destroy()

        ok_cancel = QDialogButtonBox(QDialogButtonBox.Ok 
                                     | QDialogButtonBox.Cancel)        
        ok_cancel.accepted.connect(get_parameters)
        ok_cancel.rejected.connect(self.disvis_dialog.destroy)  

        outer_layout.addLayout(models_layout)
        outer_layout.addLayout(distances_layout)
        outer_layout.addWidget(ok_cancel)

        self.disvis_dialog.ui_area.setLayout(outer_layout)
        self.disvis_dialog.manage(None)


    def create_disvis_input(self, chains, distances, pseudobonds,
                            pb_lines=None):

        minimum = distances["Minimum"].text()
        maximum = distances["Maximum"].text()
        lines = []

        def write_line(atom1, atom2, sort=False):
            atoms = [atom1, atom2]
            atom_strings = [None] * len(atoms)
            replace_with_space = [":", "@"]
            for i, atom in enumerate(atoms):
                string = re.search("/.+",
                                   atom.string(style="command line")).group(0)
                string = string.replace("/", "")                
                for character in replace_with_space:
                    string = string.replace(character, " ")
                atom_strings[i] = string
            if sort:
                atom_strings.sort()
            line = (atom_strings[0] + " " + atom_strings[1] + " " + minimum 
                    + " " + maximum)

            return line

        if chains["Fixed"] == chains["Scanning"]:
            chain = chains["Fixed"]
            for pb in pseudobonds:
                if len(pb.models) != 1:
                    continue
                atom1, atom2 = pb.atoms
                if (atom1 not in chain.atoms or atom2 not in chain.atoms):
                    continue
                lines.append(write_line(atom1, atom2, sort=True))                
        else:
            for pb in pseudobonds:
                if len(pb.models) == 1:
                    continue
                atom1, atom2 = pb.atoms
                atoms = [atom1, atom2]
                for atom in atoms:
                    atom.chain = object()
                    for chain in chains:
                        current_chain = chains[chain]
                        if atom in current_chain.atoms:
                            atom.chain = current_chain
                if (atoms[0].chain == chains["Fixed"] 
                    and atoms[1].chain == chains["Scanning"]):
                    lines.append(write_line(atom1, atom2))
                elif (atoms[0].chain == chains["Scanning"]
                      and atoms[1].chain == chains["Fixed"]):
                    lines.append(write_line(atom2, atom1))

        if len(lines) == 0:
            print("No pseudobonds match the criteria")
            return
        
        if pb_lines is not None:
            title = "Save pseudobonds and DisVis input"
            extension = "*.pb *.txt"
        else:
            title = "Save DisVis input"
            extension = "*.txt"

        self.save_subset(title, extension, [pb_lines, lines])

        self.disvis_dialog.destroy()
                

    def add_models(self, models):

        # The main tool window contains two treewidgets showing models;
        # the model selector for non-pseudobond structural models, and
        # the pbonds menu for pseudobond models. With this method, both
        # treewidgets of the main tool window can be filled.

        for model in models:
            self.get_structure_type(model)
            if not hasattr(model, "structure_type"):
                continue
            model_name = model.name
            # For models in the model selector, show the model's
            # "id_string" attribute as ID
            if model.structure_type == "Non-pb":
                treewidget = self.model_selector
                checkstate = Qt.Checked
                column_1_text = model.id_string
            elif model.structure_type == "Pb":
                treewidget = self.pbonds_menu
                checkstate = self.get_checkstate(model)
                # Show the model's code
                column_1_text = self.get_ids(model)
        
            model.item = QTreeWidgetItem(treewidget)
            item = model.item
            item.setText(0, model_name)
            item.setCheckState(0, checkstate)
            item.setText(1, column_1_text)
            item.setFlags(item.flags() & ~Qt.ItemIsSelectable)
            item.model = model
            
            
    def get_checkstate(self, model):
        
        if (model.get_selected() and model.get_selected(fully=True)):
            checkstate = Qt.Checked
        elif (model.get_selected() and not model.get_selected(fully=True)):
            checkstate = Qt.PartiallyChecked
        else:
            checkstate = Qt.Unchecked
            
        return checkstate        
            
            
    def get_ids(self, model):
        
        pbs = model.pseudobonds
        models = [model.id_string for model 
                           in pbs.unique_structures]
        models = sorted(models)
        ids = ",".join(models)
        
        return ids
        

    def get_structure_type(self, model):

        # Determine which treewidget the model belongs to

        # Structural, non-pseudobond models go in the model selector
        if (isinstance(model, Structure) 
                and not isinstance(model, PseudobondGroup)):
            model.structure_type = "Non-pb"
        # Models made from .pb files go in the pbonds_menu
        elif (isinstance(model, PseudobondGroup) and model.name[-3:] == ".pb"):
            model.structure_type = "Pb"          


    def trigger_handler(self):

        # Create trigger handlers for:
        # - change of selection, to change the CheckState of items in
        # the pbonds menu upon (de)selection of pseudobond models in
        # the session
        # - models being added to and removed from the session, to be
        # added and removed from the tool window

        self.triggerset = self.session.triggers

        self.change_selection_handler = self.triggerset.add_handler(
            SELECTION_CHANGED,self.selection_handler
            )
        self.add_model_handler = self.triggerset.add_handler(
            ADD_MODELS, self.model_handler
            )
        self.remove_model_handler = self.triggerset.add_handler(
            REMOVE_MODELS, self.model_handler
            )


    def selection_handler(self, trigger, trigger_data):

        # Called upon change of model selection, to adjust the
        # CheckState of models in the pbonds menu, if necessary

        for model in self.session.models:
            # Skip all models that are not in the pbonds menu
            if not model.name[-3:] == ".pb":
                continue
            item = model.item
            checkstate = self.get_checkstate(model)
            item.setCheckState(0, checkstate)


    def model_handler(self, trigger, trigger_data):

        # Called when models are added to or removed from the session.
        # The trigger_data is a list of the models that were added or
        # removed

        # If models are added, call "add_models" method
        if trigger == "add models":
            self.add_models(trigger_data)
        # If models are removed, call "remove_models" method
        else:
            self.remove_models(trigger_data)


    def remove_models(self, models):

        # Called when models are removed from the session

        for model in models:
            if not hasattr(model, "item"):
                continue
            item = model.item
            treewidget = item.treeWidget()
            root = treewidget.invisibleRootItem()
            root.removeChild(item)             


    def cleanup(self):
        
        # Called when main tool window is closed. Trigger handlers are
        # removed to prevent errors when the trigger fires

        self.triggerset.remove_handler(self.change_selection_handler)
        self.triggerset.remove_handler(self.add_model_handler)
        self.triggerset.remove_handler(self.remove_model_handler)


class Alignment:
    

    def __init__(
        self, start, end, first_residue_number, crosslink_position, model_id,
        chain):

        # Start and end position indicate the range of positions in the
        # sequence that is spanned by the peptide        
        self.start_position = start + first_residue_number
        self.end_position = end + first_residue_number
        # Position of the crosslinked residue in the sequence        
        self.crosslink_position = crosslink_position + first_residue_number
        residue = chain.residues[crosslink_position]
        atoms = residue.atoms
        for atom in atoms:
            if atom.name == "CA":
                self.atom = atom
        # String indicating on which model and chain the alignment was
        # found
        self.id_string = "#" + model_id + "/" + chain.chain_id


class PrePseudobond:

    # To avoid clashed with ChimeraX's Pseudobond class, this class is
    # named PrePseudobond    
    
    def __init__(self, alignment1, alignment2, peptide_pair):

        # The crosslink positions of the two alignments dictate
        # between which atoms the pseudobond is formed
        self.pos1 = alignment1.crosslink_position
        self.pos2 = alignment2.crosslink_position
        self.id1 = alignment1.id_string
        self.id2 = alignment2.id_string
        self.atom1 = alignment1.atom
        self.atom2 = alignment2.atom        
        self.peptide_pair = peptide_pair
        self.score = peptide_pair.Score
        # A string is created that will be one line in a .pb file
        self.line = self.create_pb_line()
        # Check for overlap between the two alignments
        self.is_overlapping = self.find_overlap(alignment1, alignment2)
        # Check whether the pseudobond is a self-link
        self.is_selflink = self.find_selflinks()

    
    def create_pb_line(self):

        pb_sorted = sorted([
            self.id1 + ":" + str(self.pos1) + "@CA",
            self.id2 + ":" + str(self.pos2) + "@CA"
            ])
        pb_line = pb_sorted[0] + " " + pb_sorted[1]
        
        return pb_line

    
    def find_overlap(self, alignment1, alignment2):

        if self.id1 != self.id2:
            is_overlapping = False
        else:
            maximum = max(alignment1.start_position, alignment2.start_position)
            minimum = min(alignment1.end_position, alignment2.end_position)
            if maximum < minimum:
                is_overlapping = True
            else:
                is_overlapping = False
        return is_overlapping

    
    def find_selflinks(self):

        is_selflink = False
        if (self.is_overlapping and self.pos1 == self.pos2):
            is_selflink = True

        return is_selflink
        

class Slider:
    
    
    def __init__(self, value_type="distance", enabled=True, minimum=None,
                 maximum=None, pbs=None):
        
        self.enabled = enabled
        self.layout = QGridLayout()
        if value_type == "distance":
            title = "Pseudobond distance (Å)"
        elif value_type == "score":
            title = "Score"
        elif value_type == "zscore":
            title = "z-score"
        label = QLabel(title + ":")
        
        self.slider = QRangeSlider(Qt.Horizontal)
        setter_min = QLineEdit()
        setter_max = QLineEdit()
        self.invert = QCheckBox("Invert")
        
        self.layout.addWidget(label, 0, 0, 1, 4)
        self.layout.addWidget(self.slider, 1, 1, 1, 2)
        self.layout.addWidget(setter_min, 1, 0)
        self.layout.addWidget(setter_max, 1, 3)
        self.layout.addWidget(self.invert, 1, 4)
        self.layout.setColumnMinimumWidth(1, 300)

        i = self.layout.count()
        self.widgets = [None]*i
        while i > 0:
            j = i - 1
            self.widgets[j] = self.layout.itemAt(j).widget()
            i -= 1 
            
        if not enabled:
            return
        
        minimum, maximum = self.scale(minimum, maximum)
        
        self.set_full_range(minimum, maximum)
        signal = self.slider.valueChanged
        self.function = None
        for pb in pbs:
            pb.outside_range = False
        within_range = lambda: self.within_range(value_type, pbs,
                                                 self.function)
        signal.connect(within_range)
        signal.connect(self.adjust_setters)
        
        self.setters = [setter_min, setter_max]
        alignments = [Qt.AlignRight, Qt.AlignLeft]
        texts = self.get_initial_texts(minimum, maximum)
        is_minimum = [True, False]
        for i, setter in enumerate(self.setters):
            setter.setAlignment(alignments[i])
            setter.setValidator(QDoubleValidator())
            setter.setText(str(texts[i]))
            setter.textEdited.connect(lambda text, m=is_minimum[i]:
                                      self.change_slider_value(text, m))
                
        self.invert.clicked.connect(within_range)
            
        # Atrribute to make sure that setters are not adjusted when slider is
        # changed due to editing of a setter   
        self.no_adjustment = False
        
        
    def get_initial_texts(self, minimum, maximum):
        
        texts = self.get_real_values([minimum, maximum])
        
        return texts
        
        
    def scale(self, minimum, maximum):
        
        self.scaling = 1
        if maximum <= 1:
            self.scaling = 1000
        minimum = minimum * self.scaling
        maximum = maximum * self.scaling + 1
        
        return minimum, maximum
    

    def set_full_range(self, minimum, maximum):
        
        # Since slider does not function when minimum != 0, create this
        # function for ZScoreSlider.
        minimum, maximum = self.get_slider_values(minimum, maximum)
        
        # After minimum and maximum values, None and None are also passed in.
        # Haven't figured out why, but for now, this is the solution:
        if minimum is None:
            return
        
        self.slider.setRange(minimum, maximum)
        self.slider.setValue((minimum, maximum))         
                
                
    def get_slider_values(self, minimum, maximum):
        
        return minimum, maximum

        
    def within_range(self, value_type, pbs, function):

        if value_type == "distance":
            distance_slider = self
            score_slider = self.linked_slider
        elif value_type == "score":
            distance_slider = self.linked_slider
            score_slider = self  
            
        if pbs is None:
            pbs = self.pbs
        
        # The operators used to determine whether the pseudobond is within or
        # outside the specified range, depend on whether the "Invert" 
        # checkbox is checked.
        if distance_slider.enabled:
            dist_op = self.get_operators(distance_slider)
            distance_range = [float(setter.text()) for setter in distance_slider.setters]
                
            for pb in pbs:
                distance = pb.length
                wrong_dist = self.check_value(dist_op, distance, distance_range)
                pb.outside_range = wrong_dist
                
        else:
            for pb in pbs:
                pb.outside_range = False
            
        if not score_slider.enabled:
            function(pbs)
            return
            
        score_op = self.get_operators(score_slider)
        score_range = [float(setter.text()) for setter in score_slider.setters]
        for pb in pbs:
            score = pb.xlinkx_score
            wrong_score = self.check_value(score_op, score, score_range)
            wrong_dist = pb.outside_range
            pb.outside_range = (wrong_dist or wrong_score) 
            
        function(pbs)
        
        
    def check_value(self, operators, value, range_lst):
        
        wrong_value = operators[0](operators[1](value, range_lst[0]),
                                   operators[2](value, range_lst[1]))
        
        return wrong_value
        
        
    def get_operators(self, slider):

        # If invert is checked, pseudobonds with a value larger than the 
        # minimum and smaller than the maximum are excluded
        if slider.invert.isChecked():
            operators = [operator.and_, operator.gt, operator.lt]
        else:
            operators = [operator.or_, operator.lt, operator.gt]
            
        return operators
        
            
    def get_real_values(self, values):
        
        if self.scaling == 1:
            dtype = int
        else:
            dtype = float
            
        return [dtype(value / self.scaling) for value in values]
            
                        
    def adjust_setters(self):
        
        # Make sure that setters are not adjusted when slider is changed due to
        # editing of a setter
        if self.no_adjustment:
            self.no_adjustment = False
            return
        
        values = self.slider.value()
        real_values = self.get_real_values(values)
        for i, setter in enumerate(self.setters):
            setter.setText(str(real_values[i]))  
    
            
    def change_slider_value(self, text, is_minimum):
        
        try:
            value = float(text)
        except:
            return
        
        self.no_adjustment = True
    
        if is_minimum:
            minimum, _ = self.get_slider_values(value, 0)
            maximum = self.slider.value()[1]
        else:
            minimum = self.slider.value()[0]
            _, maximum = self.get_slider_values(0, value)
    
        self.slider.setValue((minimum, maximum)) 


class ExportSlider(Slider):
    
    
    def __init__(self, value_type="distance", enabled=True, minimum=None,
                 maximum=None, pbs=None, *args):
        
        super().__init__(value_type, enabled, minimum, maximum, pbs)
        self.function = self.display_pseudobonds
   
        if not enabled:
            for widget in self.widgets:
                widget.setEnabled(False)
            return
        
    
    def display_pseudobonds(self, pbs):
        
        for pb in pbs:
            if pb.outside_range:
                display = False
            else:
                display = True
            pb.display = display
            


class VisualizeSlider(Slider):

    
    def __init__(self, value_type="distance", enabled=True, minimum=None,
                 maximum=None, pbs=None, dialog=None):
        
        super().__init__(value_type, enabled, minimum, maximum, pbs)
        self.value_type = value_type
        self.function = self.cutoff_style
        self.maximum = maximum
        self.pbs = pbs
        self.dialog = dialog
        self.attributes = ["color", "radius"]

        for widget in self.widgets:
            widget.setEnabled(False)

        
    def cutoff_style(self, pbs, enabled=True):
        
        dialog = self.dialog
        custom_values = dialog.custom_values
        
        for i, pb in enumerate(pbs):
            key = list(custom_values.keys())[pb.outside_range]
            settings = custom_values[key]
            for attribute in self.attributes:
                if (attribute == "color" and dialog.gradient_active 
                    and not pb.outside_range):
                    pb.color = pb.gradient_color
                    continue
                value = settings[attribute][i]
                setattr(pb, attribute, value)


class NoEditDelegate(QStyledItemDelegate):
    
    def __init__(self, parent = None):
        super().__init__(parent)

    def createEditor(self, parent, option, index):
        return None


class VisualizePseudobonds(Pseudobonds):


    def __init__(self, pbond_pointers=None):

        super().__init__(pbond_pointers)


    @property
    def dashes(self):

        return [group.dashes for group in self.groups]

	
    @dashes.setter
    def dashes(self, value):

        groups = self.groups

        if not isinstance(value, list):
            value = [value] * len(groups)

        for i, group in enumerate(groups):
            group.dashes = value[i]
            
            
class DragHandler(object):
    
    # This class isn't working yet
    
    def __init__(self, figure=None) :
        
        # Create a new drag handler and connect it to the figure's event 
        # system.
        # If the figure handler is not given, the current figure is used 
        # instead,
        
        if figure is None:
            figure = plt.gcf()
        # simple attibute to store the dragged text object
        self.dragged = None

        # Connect events and callbacks
        figure.canvas.mpl_connect("pick_event", self.on_pick_event)
        figure.canvas.mpl_connect("button_release_event", self.on_release_event)

    def on_pick_event(self, event):
        
        # Store which text object was picked and were the pick event occurs.

        print("Pick!")

        if isinstance(event.artist, Text):
            self.dragged = event.artist
            self.pick_pos = (event.mouseevent.xdata, event.mouseevent.ydata)
        return True

    def on_release_event(self, event):
        
        # Update text position and redraw

        if self.dragged is not None :
            old_pos = self.dragged.get_position()
            new_pos = (old_pos[0] + event.xdata - self.pick_pos[0],
                       old_pos[1] + event.ydata - self.pick_pos[1])
            self.dragged.set_position(new_pos)
            self.dragged = None
            plt.draw()
        return True

