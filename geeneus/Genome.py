# Public facing API for accessing gene information
#
# Copyright 2012-2015 by Alex Holehouse - see LICENSE for more info
# Contact at alex.holehouse@wustl.edu
import Bio.Entrez
import backend

class GeneManager:
    """ Manager object used to access relevant information via the below methods. The *Manager classes are the only public facing classes in Geeneus, and are the only ones which should be used to access remote information."""

    def __init__(self, email, cache=True, retry=0):
        """Returns a fully formed manager object which can be queried by the other 
        functions in this class. 

        email    Must be a valid email, as is required by the NCBI servers. For more 
                 information on NCBI usage guidelines please see 
                 [http://www.ncbi.nlm.nih.gov/books/NBK25497/].

        cache    Determines if the manager object should cache requests in memory, or 
                 the NCBI database should be queried every time a request is made. 

        retry    The number of times the networking utilities will retry on a failed 
                 connection

        timeout  The number of seconds the networking utilities wait after making a 
                 request before deciding that request has failed
        """
        
        self.datastore = backend.GeneParser.GeneRequestParser(email, cache, retry)

        if self.datastore.error():
            self.error_status = True
        else:
            self.error_status = False
    
    def get_gene_sequence(self, ID):
        """Returns the DNA sequence associated with this accession number"""
        print "+-------------------------------------------------------------------------------------+"
        print "| This function will return a gene sequence. Please wait for 0.1.4 for functionality  |"
        print "+-------------------------------------------------------------------------------------+"
        
        #return self.datastore.get_consensus_sequence(ID)
        
    def get_gene_coding_sequence(self, ID):
        print "+----------------------------------------------------------------------+"
        print "| This function will return a dictionary of isoform-sequence pairings  |"
        print "+----------------------------------------------------------------------+"

    def get_gene_SNP(self, ID):
        print "+--------------------------------------------------------------------------------------------+"
        print "| This function will return a dictionary of SNP-number:Dictionary with SNPs for the gene ID  |"
        print "+--------------------------------------------------------------------------------------------+"
    
    def get_raw_xml(self, ID):
        return self.datastore.get_raw_xml(ID)
