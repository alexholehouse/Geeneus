# Deals with cache and storing of gene data objects, as well as coordinating network access (private)
#
# Copyright 2012 by Alex Holehouse - see LICENSE for more info
# Contact at alex.holehouse@wustl.edu

from Bio import Entrez
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

import GeneObject
import Networking

# RequestParser is a wrapper datastore layer, which facilitates non-redundant access to the
# GeneObjects. It is also where direct communication with the NCBI server occurs, but no xml
# parsing actually happens here.
class GeneRequestParser:
    
    def __init__(self, email):
        
        Entrez.email = email
        self.gene_datastore = {-1 : GeneObject.GeneObject([])}
                

    # Primary function for querying the gene_datastore in the RequestParser object
    # Functions which are getting information on a Gene based on GeneID should call this
    # and functionality associated with extracting information from the Gene should be
    # added to the GeneObject class.

    # Return value:
    # Will return a GeneObject, or a -1 in the case of an error
    
    def get_gene_object(self, GeneID):
        if GeneID not in self.gene_datastore:
            
            gene_handle = Networking.efetchGene(GeneID)
            
            if gene_handle == -1:
                print "\nWARNING: Network issues, unable to download data associated with '{ID}'\n".format(ID=GeneID)
                self.gene_datastore[GeneID] = GeneObject.GeneObject(-1)
                                
            self.gene_datastore[GeneID] = GeneObject.GeneObject(Entrez.read(gene_handle))
        
        return  self.gene_datastore[GeneID]
        
    
    def get_consensus_sequence(self, GeneID):
        gene = self.get_gene_object(GeneID)
        return gene.get_full_sequence()
   
