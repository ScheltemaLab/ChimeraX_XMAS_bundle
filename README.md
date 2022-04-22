# XMAS

The package can easily be installed from the ChimeraX toolshed
https://cxtoolshed.rbvi.ucsf.edu/apps/chimeraxxmas

Short install instructions:

- Store the 'bundle_info.xml' file and the 'src' folder in a folder of choice. The 'src' folder should contain the following files:
  - '__init__.py'
  - 'info_file.py'
  - 'mzidentml.py'
  - 'read_evidence.py'
  - 'tool.py'
  - 'z_score.py' 
  Additionally, a 'docs' and 'matplotlib_venn' folder should be present in the 'src' folder.
- In the ChimeraX application, run the command 'devel install “PATH_TO_FOLDER_OF_CHOICE”'
- The bundle can be found in the ChimeraX application under 'Tools' - 'Structure Analysis' as 'XMAS'
