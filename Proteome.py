# Public facing API for accessing protein information
#
# Copyright 2012 by Alex Holehouse - see LICENSE for more info
# Contact at alex.holehouse@wustl.edu

import Bio.Entrez
import Geeneus.Backend

class ProteinManager:

    def __init__(self, email, cache=True):
        self.datastore = Geeneus.Backend.ProteinParser.ProteinRequestParser(email, cache)
        if self.datastore.Error():
            self.Error = True
        else:
            self.Error = False
    
    def get_protein_name(self, ID):
        return self.datastore.get_protein_name(ID)
    
    def get_protein_sequence(self, ID):
        return self.datastore.get_sequence(ID)

    def get_raw_xml(self, ID):
        return self.datastore.get_raw_xml(ID)

    def get_variants(self, ID):
        return self.datastore.get_variants(ID)

    def get_geneID(self, ID):
        return self.datastore.get_geneID(ID)

    def get_protein_seqeuence_length(self, ID):
        print "Will return protein sequence length"

    def get_ID_type(self, ID):
        self.datastore.ID_type(ID)

    def run_translation(self, Acc):
        return self.datastore.translate_Asc2GI(Acc)

# END OF CLASS - below are general untilities

def ID_type(ProteinID):
    
    return(Geeneus.Backend.ProteinParser.ID_type(ProteinID))
 
