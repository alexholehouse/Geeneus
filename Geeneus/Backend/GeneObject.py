# Provides an object-based API for a gene (private)
#
# Copyright 2012 by Alex Holehouse - see LICENSE for more info
# Contact at alex.holehouse@wustl.edu

from Bio import Entrez, Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC

# object attributes
#

class GeneObject:

#--------------------------------------------------------
#
#--------------------------------------------------------
    def __init__(self, genexml):
        
        self.__set_default_attributes(error=False)

        # check for error, and return with error flag = True
        if proteinxml == -1:        
            self._error = True
            return
        
        # check xml is a real, correct single gene
        if not self.__xml_is_OK(genexml):
            return
        
        # if we've got here we have a non-errorful, real
        # xml structure to deal with, so the object exists
        self._exists = True
        
        # Building the full sequence requires an additional network call to pull down
        # the nucleotide data. Should this fail, the __build_full_sequence() function
        # will return -1, in which case we reset the object status to error and return
        if (self.___build_full_sequence(genexml) == -1):
            self.__reset_to_error()
            return
        
        self.accession = self.__extract_accession(genexml)
        self.version = self.__extract_version(genexml)
        self.GI = self.__extract_GI(genexml)
        self.locus = self.__extract_locus(genexml)
        
#--------------------------------------------------------
#
#--------------------------------------------------------
# object gettes for attributes

    def get_full_sequence(self):
        return (self.full_sequence, self.full_sequence_c)

    def get_version(self):
        return self.version

    def get_accession(self):
        return self.accession

    def get_versioned_accession(self):
        return self.accession + "."+ self.version

    def get_gi(self):
        return self.GI
        
#--------------------------------------------------------
#
#--------------------------------------------------------
# Returns the coding sequence for the gene, as defined by
# the exons
#
    def get_coding_sequence(self, genexml):
        # stuff
        print genexml[0]["Entrezgene_locus"][0]["Gene-commentary_products"][0]["Gene-commentary_genomic-coords"][0]["Seq-loc_mix"]["Seq-loc-mix"]
        return "ok"


#--------------------------------------------------------
#
#--------------------------------------------------------
# Get gene version from XML
#
def __extract_version(self, genexml):
    return genexml[0]["Entrezgene_locus"][0]["Gene-commentary_version"]

#--------------------------------------------------------
#
#--------------------------------------------------------
# Get gene accession from XML
#
def __extract_accession(self, genexml):
    return genexml[0]["Entrezgene_locus"][0]["Gene-commentary_accession"]

#--------------------------------------------------------
#
#--------------------------------------------------------
# get gene GI from XML
#
def __extrat_GI(self, genexml):
    gene_seq_dict = genexml[0]["Entrezgene_locus"][0]["Gene-commentary_seqs"][0]["Seq-loc_int"]["Seq-interval"]
    return int(gene_seq_dict["Seq-interval_id"]["Seq-id"]["Seq-id_gi"])
#--------------------------------------------------------
#
#--------------------------------------------------------
# get gene locus from XML 
#
def __extract_locus(self, genexml):
    return genexml[0]["Entrezgene_location"][0]["Maps_display-str"]
    

def __set_default_attributes(self, error):
    self.full_sequence = ""
    self.full_sequence_c = ""
    self._exists = False
    
    # if there's an error...
    if (error):
        self._error = True
    else:
        self._error = False
    self.full_sequence_length = 0
    self.accession = ""
    self.version = ""
    self.GI = ""
    self.locus = ""

#--------------------------------------------------------
#
#--------------------------------------------------------
# Internal function which takes the genexml, and extracts the full
# coding sequence. If anything goes wrong while pulling nucleotide
# data will return -1, else returns nothing, but changes the state 
# of the objects self.full_sequence[_c] variables behind the 
# scenes #antifunctionalprogramming.
#
    def ___build_full_sequence(self, genexml):
        
        gene_seq_dict = genexml[0]["Entrezgene_locus"][0]["Gene-commentary_seqs"][0]["Seq-loc_int"]["Seq-interval"]
        
        # NB - we could get GI from self.GI, but that makes the assumption
        # that we've already run and set self.GI. Best not to do so!
        GI = int(gene_seq_dict["Seq-interval_id"]["Seq-id"]["Seq-id_gi"])
        
        # get positional and locus information (need to add 1!)
        pos_from = int(gene_seq_dict["Seq-interval_from"])+1
        pos_to = int(gene_seq_dict["Seq-interval_to"])+1

        handle = Networking.eFetchNucleotide(GI, pos_from, pos_to, 1) 
        handle_2 = Networking.eFetchNucleotide(GI, pos_from, pos_to, 2) 
        
        # if we get an error, return -1, which is dealt with in __init__ and
        # the object is set to an error state
        if (handle == -1 or handle_2 == -1):
            return -1
        
        #  given we didn't return -1 no errors, occured, so grab sequence. If t
        try:
            self.full_sequence = SeqIO.read(handle, "fasta")
            self.full_sequence_c = SeqIO.read(handle_c, "fasta")

        # TO DO - catch relevant exceptions only!`
        except:
            return -1
        

        
        
        handle.close()
        handle_c.close()

        return

#--------------------------------------------------------
#
#--------------------------------------------------------
#
#

    def __xml_is_OK(self, genexml):
        if len(proteinxml) > 1:
            print "WARNING[GeneObject.__xml_is_ok()] - GeneXML detected more than one record associated with this GI.\nThis should never happen."

        if len(genexml) == 0:
            return False
        
        
