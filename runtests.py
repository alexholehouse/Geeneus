import Tests
import unittest
#import Proteome

suite = unittest.TestLoader().loadTestsFromTestCase(Tests.Proteome_tests.TestProteomeFunctions)
unittest.TextTestRunner(verbosity=2).run(suite)

suite = unittest.TestLoader().loadTestsFromTestCase(Tests.ProteinParser_tests.TestProteinParserFunctions)
unittest.TextTestRunner(verbosity=2).run(suite)

#Tests.Proteome_tests.example()
