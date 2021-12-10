class ZScoreSelector:


    def __init__(self):

        from PyQt5.QtWidgets import QDialog, QLabel, QLineEdit, QPushButton, QDialogButtonBox, QHBoxLayout, QVBoxLayout 

        self.main_dialog = QDialog()
        self.main_dialog.setWindowTitle("Create HADDOCK input from DisVis output")

        label = QLabel("Select a DisVis output folder:")
        line_edit = QLineEdit()
        file_button = QPushButton()
        file_button.clicked.connect(self.file_dialog(line_edit))
        ok_cancel = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        ok_cancel.accepted.connect(self.ok_clicked)
        ok_cancel.rejected.connect(self.main_dialog.close)

        file_layout = QHBoxLayout()
        file_layout.addWidget(line_edit)
        file_layout.addWidget(file_button)
        
        main_layout = QVBoxLayout()
        main_layout.addWidget(label)
        main_layout.addLayout(file_layout)
        main_layout.addWidget(ok_cancel)

        self.main_dialog.setLayout(main_layout)
        self.main_dialog.show()
        

    def file_dialog(self, line_edit):

        from PyQt5.QtWidgets import QFileDialog

        return


    def ok_clicked(self):

        return


