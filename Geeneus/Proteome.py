# Public facing API for accessing protein information
#
# Copyright 2012 by Alex Holehouse - see LICENSE for more info
# Contact at alex.holehouse@wustl.edu

import Bio.Entrez
import Geeneus.Backend

class ProteinManager:
    """ Create a manager """
    def __init__(self, email, cache=True, retry=0, timeout=20):
        self.datastore = Geeneus.Backend.ProteinParser.ProteinRequestParser(email, cache, retry, timeout)
        if self.datastore.error():
            self.error_status = True
        else:
            self.error_status = False
    
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

    def batch_get_protein_sequence(self, IDList):
        return self.datastore.batchFetch(self.datastore.get_sequence, IDList)

    def batch_get_protein_name(self, IDList):
        return self.datastore.batchFetch(self.datastore.get_name, IDList)

    def batch_get_variants(self, IDList):
        return self.datastore.batchFetch(self.datastore.get_variants, IDList)

    def purge(self):
        self.datastore.purge_data_store()

    def get_size_of_datastore():
        self.datastore.get_size()




# END OF CLASS - below are general untilities

def ID_type(ProteinID):
    return(Geeneus.Backend.ProteinParser.ID_type(ProteinID))
 
