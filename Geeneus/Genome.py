# Public facing API for accessing gene information
#
# Copyright 2012 by Alex Holehouse - see LICENSE for more info
# Contact at alex.holehouse@wustl.edu
import Bio.Entrez
import Backend

class GeneManager:
    
    def __init__(self, email):
        self.datastore = Geeneus.Backend.GeneParser.GeneRequestParser(email)
    
    def get_gene_sequence(self, ID):
        return self.datastore.get_consensus_sequence(ID)
        
    def get_gene_coding_sequence(self, ID):
        print "This function will return a dictionary of isoform-sequence pairings"

    def get_gene_SNP(self, ID):
        print "This function will return a dictionary of SNP-number:Dictionary with SNPs for the gene ID"

    def get_gene_object(self, ID):
        print "This function will return a gene object, with various attributes associated with the gene"
