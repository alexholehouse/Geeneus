#!/usr/bin/python

import sys
import os
# add the location /Geeneus/ to the Pythonpath
sys.path.insert(0,os.path.abspath(__file__+"/../../geeneus/"))

import test
import unittest

network = True
proteome = True
proteinParser = True

if network:
    suite = unittest.TestLoader().loadTestsFromTestCase(test.Networking_tests.TestNetworkingFunctions)
    unittest.TextTestRunner(verbosity=2).run(suite)

if proteome:
    suite = unittest.TestLoader().loadTestsFromTestCase(test.Proteome_tests.TestProteomeFunctions)
    unittest.TextTestRunner(verbosity=2).run(suite)

if proteinParser:
    suite = unittest.TestLoader().loadTestsFromTestCase(test.ProteinParser_tests.TestProteinParserFunctions)
    unittest.TextTestRunner(verbosity=2).run(suite)

