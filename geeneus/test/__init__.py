# Copyright 2012 by Alex Holehouse - see LICENSE for more info
# Contact at alex.holehouse@wustl.edu

""" This is test folder for running unit tests.
"""
import sys
import os
# Add the parent directory (which holds the geeneus package)
sys.path.insert(0,os.path.abspath(__file__+"/../../../"))

import Proteome_tests
import Genome_tests
import GeneObject_tests
import GeneParser_tests
import ProteinObject_tests
import ProteinParser_tests
import Networking_tests
import UniprotAPI_tests
