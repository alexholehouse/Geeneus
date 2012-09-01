# Deals with cache and storing of protein data objects, as well as coordinating network access  (private)
#
# Copyright 2012 by Alex Holehouse - see LICENSE for more info
# Contact at alex.holehouse@wustl.edu

import sys
#import signal

# Biopython selective imports
from Bio import Entrez
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

import ProteinObject
import Networking

class ProteinRequestParser:

# Initialization function, set Entrez.email for calls, an ensure key -1 is set
# to a non-existant object
#
    def __init__(self, email, cache):
        Entrez.email = email
        self.cache = cache
        self.protein_datastore = {-1 : ProteinObject.ProteinObject([])}

#--------------------------------------------------------
#
#--------------------------------------------------------
#
#
    def __get_protein_object(self, ProteinID):
        if ProteinID not in self.protein_datastore or not self.cache:
            protein_handle = Networking.efetchProtein(ProteinID)
            
            # check if handle represents an error
            # PLEASE NOTE that if you use a protein ID which fails to get a protein
            # this is *not* an error (and Entrez.efetch will behave accordingly), so
            # instead in such a situtation we simply return an empty ProteinObject with
            # the exists attribute set to False, but Error is also false.
            # IF, however, protein_handle returns -1 (indicating some error) then we return
            # an empty ProteinObject with error set to True
            if protein_handle == -1:
                self.protein_datastore[ProteinID] = ProteinObject.ProteinObject(-1)

            else:
                self.protein_datastore[ProteinID] = ProteinObject.ProteinObject(Entrez.read(protein_handle))
            
        return self.protein_datastore[ProteinID]

#--------------------------------------------------------
# 
#--------------------------------------------------------
#

    def get_protein_name(self, ID):
        return (self.__get_protein_object(ID)).get_protein_name()

#--------------------------------------------------------
# 
#--------------------------------------------------------
#
    def get_sequence(self, ID):
        ProtObj = self.__get_protein_object(ID)
        return ProtObj.get_protein_sequence()
#--------------------------------------------------------
#
#--------------------------------------------------------
    def get_SNPs(self, ID):
        ProtObj = self.__get_protein_object(ID)
        return ProtObj.get_SNPs()

#--------------------------------------------------------
#
#--------------------------------------------------------
    def get_GeneId(self, ID):
        ProtObj = self.__get_protein_object(ID)
        return ProtObj.get_GeneId()


#--------------------------------------------------------
#
#--------------------------------------------------------
    def get_protein_sequence_length(self, ID):
        ProtObj = self.__get_protein_object(ID)
        return ProtObj.get_protein_sequence_length()

#--------------------------------------------------------
#
#--------------------------------------------------------
#
# Function to translate and accession value into a GI. If 
# there is more than one possible GI, then simply returs
# -1. Probably need better behaviour 

    def translate_Asc2GI(self, Accession):
        record = Networking.esearch("protein", Accession)

        IdList = record["IdList"]
        
        if len(IdList) == 1:
            return int(record["IdList"][0])
        else:
            print "There are {op} possible options".format(op=len(IdList))
            return -1
            
            
            
# +-------------------------------------------------------+
# |                    END OF CLASS                       |
# +-------------------------------------------------------+


# returns 0 for GI
# returns 1 for refseq
# returns 2 for swissprot
# returns 3 for anything else
# returns -1 for error
def ID_type(ProteinID):
    
    if (str(ProteinID).isdigit()):
        return 0
