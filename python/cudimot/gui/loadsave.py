"""
CUDIMOT GUI - Functions for loading/saving a model project
"""
import os

import yaml

import cudimot.gui.code_templates as templates

def save_project(config, generate_code=True):
    projdir = config.get("projdir", None)
    if not projdir:
        raise ValueError("Project directory not set")
    elif not os.path.exists(projdir):
        os.makedirs(projdir)
    elif not os.path.isdir(projdir):
        raise ValueError(f"{projdir} already exists and is not a directory")

    with open(os.path.join(projdir, "cudimot_config.yml"), "w") as f:
        f.write(yaml.dump(config))
    if generate_code:
        create_code(config, projdir)

def create_code(config, projdir):
    # Create modelparameters.h
    with open(os.path.join(projdir, "modelparameters.h"), "w") as f:
        f.write(templates.HEADER)
        f.write(templates.MODELPARAMETERS_H.format(**config))

    # Create modelparameters.cc
    with open(os.path.join(projdir, "modelparameters.cc"), "w") as f:
        f.write(templates.HEADER)
        f.write(templates.MODELPARAMETERS_CC.format(**config))

    # Create modelfunctions.h
    with open(os.path.join(projdir, "modelfunctions.h"), "w") as f:
        f.write(templates.HEADER)
        f.write(templates.MODELFUNCTIONS_H.format(**config))

    # Create Makefile
