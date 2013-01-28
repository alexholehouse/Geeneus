import unittest
from geeneus.backend import UniprotAPI


class TestUniprotAPIFunctions(unittest.TestCase):

    # Build manager object for all tests here
    def setUp(self):
        
        self.listOfIDs = ["A0AVF1", "A2A3L6", "B0I1T2"]
        
        # set listOfIDs to a big list if we want to put UniProt server through its paces
        if True:
            self.listOfIDs = ["A0AVF1", "A2A3L6", "B0I1T2", "C9JDP6", "E2RYF6", "O00116", "P01889", "P07951", "Q15714", "Q16549", "Q9Y6Z4"]


        self.UAPI = UniprotAPI.UniprotAPI()
        
    def test_rapidfire(self):

        temp = {}
        
        for ID in self.listOfIDs:
            print "ID = " + ID
            self.UAPI.getProteinObjectFromUniProt(temp, ID)
            
            
            temp[ID].get_geneID()
            temp[ID].get_protein_sequence()
            temp[ID].get_variants()
            temp[ID].get_protein_sequence_length()
            temp[ID].get_domains()
            temp[ID].get_other_accessions()
            temp[ID].get_species()
            temp[ID].get_taxonomy()
            temp[ID].get_host()
			
            temp[ID].exists()
            temp[ID].error()
            
    def test_get_name(self):

        datastore = {}
        
        self.UAPI.getProteinObjectFromUniProt(datastore, "Q9Y6Z4")
        print "\nBAD_ID should yeild a gracefull error/ignore"
        self.UAPI.getProteinObjectFromUniProt(datastore, "BAD_ID")
        self.UAPI.getProteinObjectFromUniProt(datastore, "P01889")

        self.assertEqual(datastore["Q9Y6Z4"].get_protein_name(), "Putative uncharacterized protein KIF25-AS1")
        self.assertEqual(datastore["P01889"].get_protein_name(), "HLA class I histocompatibility antigen, B-7 alpha chain")


    def test_get_protein_sequence(self):
        pass


    def test_get_variants(self):
        pass

    def test_get_geneID(self):
        pass
