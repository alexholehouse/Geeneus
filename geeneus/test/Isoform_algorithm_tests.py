import unittest
import time
from geeneus.backend import ProteinObject
from geeneus.backend import Networking
from geeneus.backend import UniprotAPI
from Bio import Entrez
import urllib2
from xml.dom.minidom import parseString

import random

class TestIsoformAlgorithm(unittest.TestCase):

    def setUp(self):

        self.nw = Networking.Networking(30)
        self.UP = UniprotAPI.UniprotAPI()
        Entrez.email = "alex.holehouse@gmail.com"
        
    def test_barrage_1(self):



        # -------------------------------------------------------------------------------
        # patho represents a set of IDs which caused failures for various reasons
        # over the system's development. They provide an initial previously pathalogical
        # dataset to test on before we hit 1000 random pulldowns.
        #
        # - P11362    really complicated, more algorithm correctness proof really 
        #             (seriously, even I'm amazed this works!)
        #
        # - P19020    closing bracket in name. Edgecase on naming conventions
        #
        # - Q13635    Has a "'" in name. Edgecase on naming conventions
        #
        # - Q14814    Also has "'" in name.
        #
        # - Q9QXS1    Straw that broke the camel's back - isoform names are full
        #             of commas. Forced a complete re-write of the name identification
        #             region, but hopefully the new code is totally robust and a *lot*
        #             easier to understand. Seriously, isoform 'PLEC-0,1C,2A,3A' - this
        #             is not an OK name!!
        #
        #
        # Other known pathalogical inputs
        # - Q9NY33    NCBI failed to get isoform '3'. Annoyingly this is because
        #             Uniprot has a '3' isoform, but NCBI doesn't. Not a lot we can
        #             do about that, really...
        # - O95467    NCBI has no splicing variant data at all, despite talking about
        #             it :-(
        #
        # - O14772    NCBI has no splicing variant data at all
        #
        # ------------------------------------------------------------------------------
        patho = ["Q8IYH5", "Q9P0K8", "Q9UK53", "P60411", "Q9NP78", "P11362", "Q9NY33", "P19020", "Q13635", "Q9QXS1", "Q14814"]
        autofail = ["Q9NY33", "O95467", "O14772"]
        pathocounter = 0

        lengthOfTest = len(patho)+10
        
        for i in xrange(0,lengthOfTest):
            
            print "On " + str(i) + " of " + str(lengthOfTest)
            print "-----------------------------------------------------------------------------------------"

            if pathocounter < len(patho):
                ID = patho[pathocounter]
                pathocounter = pathocounter+1
            else:
                ID = self.get_random_accession()
           
            print ID
            
            # Some values are just unsolvably pathalogical, so we allow the test to skip them. By unsolvable,
            # this is where a record is just missing from one of the databases - basically a human error or
            # some kind of inconsistency means we'll never get the same isoforms.
            #
            if ID in autofail:
                print "Database error (i.e. there's some reason this will *never* work - skipping"
                continue
            
            NCBIISO = self.get_NCBI_isoform(ID)
            if NCBIISO == -1:
                # If we have a problem looking up a UniProt record
                # just skip it
                continue

            print "NCBI lookup done..."

            UNIISO = self.get_Uniprot_isoform(ID)
            if UNIISO == -1:
                # If we have a problem looking up a UniProt record
                # just skip it
                continue

                        
            print "UniProt lookup done"
            
            ## useful for failure analysis
            print "NCBI dictionary"
            print NCBIISO
            print "UniProt dictionary"
            print UNIISO
            
            self.assertEqual(self.heuristic_compare(NCBIISO, UNIISO), True)
            print "\=====Isoform sequences match!===="

    def get_random_accession(self):
        # 27 000 000 IDs, so we randomly select 1 between 1 and 20 000 000
        query = "http://www.uniprot.org/uniprot/?query=reviewed:yes+AND+organism:9606&format=xml&limit=1&offset="+str(random.randint(1,20000))+"&random=yes"
        handle = urllib2.urlopen(query)
        dom = parseString(handle.read())
        return str(dom.getElementsByTagName("accession")[0].childNodes[0].nodeValue)


    def get_NCBI_isoform(self, ID):
        handle = self.nw.efetchProtein(ID)
        if not handle == -1:
            XML = Entrez.read(handle)
            PO = ProteinObject.ProteinObject(ID, XML)
            return(PO.get_isoforms())
        else:
            return -1

    def get_Uniprot_isoform(self, ID):
        datastore = {}
        PO = self.UP.getProteinObjectFromUniProt(datastore, ID)
        
        if datastore.has_key(ID):
            return(datastore[ID].get_isoforms())
        else:
            return -1
        

    def heuristic_compare(self, NCBI, UNIPROT):
        
        useRebased = False;

        if not len(NCBI) == len(UNIPROT):
            return False

        for k in UNIPROT.keys():
            try:
                seq1 = UNIPROT[k][1]
                seq2 = NCBI[k][1]
                
            except KeyError:
                print "Key mismatch"
                print "Uniprot keys are = " + str(UNIPROT.keys())
                print "NCBI key are = " + str(NCBI.keys())
                return False
                
            if not seq1 == seq2:
                return False

        return True
        
