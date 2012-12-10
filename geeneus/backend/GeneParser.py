# Deals with cache and storing of gene data objects, as well as coordinating network access (private)
#
# Copyright 2012 by Alex Holehouse - see LICENSE for more info
# Contact at alex.holehouse@wustl.edu

from Bio import Entrez
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

import GeneObject
import Networking
import Parser as GRP

# RequestParser is a wrapper datastore layer, which facilitates non-redundant access to the
# GeneObjects. It is also where direct communication with the NCBI server occurs, but no xml
# parsing actually happens here.
class GeneRequestParser(GRP.GeneralRequestParser):
    
    def __init__(self, email, cache, retry=0, loud=True):
        try:
            GRP.GeneralRequestParser.__init__(self, email, cache, retry, loud)
            
            self.gene_datastore = {-1 : GeneObject.GeneObject(-1, [])}
            self.gene_translationMap = {-1: -1}
            self.batchableFunctions = [] # empty for now...
            
            # If we get here we're good to go!
            self.error_status = False
           
        except Exception, e:
            print eb
            print "Fatal error when creating GeneRequestParserObject"
            self.error_status = True

            
    #--------------------------------------------------------
    # PUBLIC FUNCTION
    #--------------------------------------------------------
    #
    # Check for an error in the parser object
    # Primary function for querying the gene_datastore in the RequestParser object
    # Functions which are getting information on a Gene based on GeneID should call this
    # and functionality associated with extracting information from the Gene should be
    # added to the GeneObject class.

    # Return value:
    # Will return a GeneObject, or a -1 in the case of an error
    
    def get_gene_object(self, GeneID):
        self._get_object(GeneID, self.gene_datastore, self.Networking.efetchGene, GeneObject.GeneObject)
        return self.gene_datastore[GeneID]
        
    def get_consensus_sequence(self, GeneID):
        gene = self.get_gene_object(GeneID)
        return gene.get_full_sequence()

    def get_raw_xml(self, GeneID):
        gene = self.get_gene_object(GeneID)
        return gene.get_raw_xml()

#--------------------------------------------------------
# PUBLIC FUNCTION
#--------------------------------------------------------
# Function which purges the datastore
#
    def purge_data_store(self):
        del self.gene_datastore
        self.gene_datastore = {-1 : GeneObject.GeneObject([])}

#--------------------------------------------------------
# PUBLIC FUNCTION
#--------------------------------------------------------
# Get the number of items in the datastore
#
    def get_size_of_datastore(self):

        # -1 offset because we have the -1 entry in
        return len(self.Gene_datastore)-1
   
