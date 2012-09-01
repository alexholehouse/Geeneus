# Public facing API for accessing protein information
#
# Copyright 2012 by Alex Holehouse - see LICENSE for more info
# Contact at alex.holehouse@wustl.edu

import Bio.Entrez
import Geeneus.Backend

class ProteinManager:

    def __init__(self, email, cache=True):
        self.datastore = Geeneus.Backend.ProteinParser.ProteinRequestParser(email, cache)
    
    def get_protein_name(self, ID):
        return self.datastore.get_protein_name(ID)
    
    def get_protein_sequence(self, ID):
        return self.datastore.get_sequence(ID)

    def get_protein_SNPs(self, ID):
        return self.datastore.get_SNPs(ID)

    def get_GeneId(self, ID):
        return self.datastore.get_GeneId(ID)

    def get_protein_seqeuence_length(self, ID):
        print "Will return protein sequence length"

    def run_translation(self, Acc):
        return self.datastore.translate_Asc2GI(Acc)

    

    
