# Public facing API for accessing gene information
#
# Copyright 2012 by Alex Holehouse - see LICENSE for more info
# Contact at alex.holehouse@wustl.edu
import Bio.Entrez
import backend

class GeneManager:
    """ Manager object used to access relevant information via the below methods. The *Manager classes are the only public facing classes in Geeneus, and are the only ones which should be used to access remote information."""
    def __init__(self, email):
        """ Init function. Returns a manager object which can be quered with gene accession numbers to get relevant details"""
        self.datastore = backend.GeneParser.GeneRequestParser(email)
    
    def get_gene_sequence(self, ID):
        """Returns the DNA sequence associated with this accession number"""
        return self.datastore.get_consensus_sequence(ID)
        
    def get_gene_coding_sequence(self, ID):
        print "This function will return a dictionary of isoform-sequence pairings"

    def get_gene_SNP(self, ID):
        print "This function will return a dictionary of SNP-number:Dictionary with SNPs for the gene ID"

    def get_gene_object(self, ID):
        print "This function will return a gene object, with various attributes associated with the gene"
