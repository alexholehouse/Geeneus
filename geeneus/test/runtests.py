#!/usr/bin/python

import unittest
import __init__ as test

network = False
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

