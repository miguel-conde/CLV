import os
import configparser

## Example
EX_VAR = "An example global var"

# Class to store settings
class prjSettings:
    
    def __init__(self):
        pass

## FILES
prj_files = prjSettings()

prj_files.dir_root  = os.path.abspath(os.path.join(os.path.dirname(__file__), '..')) # Project root is defined by globalsettings.py location
prj_files.dir_data  = os.path.join(prj_files.dir_root, "data")
prj_files.dir_raw   = os.path.join(prj_files.dir_data, "raw")
prj_files.dir_clean = os.path.join(prj_files.dir_data, "clean")
prj_files.dir_outputs  = os.path.join(prj_files.dir_root, "outputs")

prj_files.dir_outputs_files  = os.path.join(prj_files.dir_outputs, "files")

# RAW DATA
FILE_RETAIL_DATA = os.path.join(prj_files.dir_raw, "Online Retail.xlsx")

# CLEAN DATA


# OUTPUTS

# CONFIG

prj_files.cfg_file = os.path.join(prj_files.dir_root, "config.ini")