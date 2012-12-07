#!/usr/bin/python

import unittest
import __init__ as test

network = False
proteome = False
proteinParser = False
uniProtAPI = False
proteinObject = True
ALL = True



if uniProtAPI or ALL:
    suite = unittest.TestLoader().loadTestsFromTestCase(test.UniprotAPI_tests.TestUniprotAPIFunctions)
    unittest.TextTestRunner(verbosity=2).run(suite)

if network or ALL:
    suite = unittest.TestLoader().loadTestsFromTestCase(test.Networking_tests.TestNetworkingFunctions)
    unittest.TextTestRunner(verbosity=2).run(suite)

if proteome or ALL:
    suite = unittest.TestLoader().loadTestsFromTestCase(test.Proteome_tests.TestProteomeFunctions)
    unittest.TextTestRunner(verbosity=2).run(suite)

if proteinParser or ALL:
    suite = unittest.TestLoader().loadTestsFromTestCase(test.ProteinParser_tests.TestProteinParserFunctions)
    unittest.TextTestRunner(verbosity=2).run(suite)

if proteinObject or ALL:
    suite = unittest.TestLoader().loadTestsFromTestCase(test.ProteinObject_tests.TestProteinObjectFunctions)
    unittest.TextTestRunner(verbosity=2).run(suite)

