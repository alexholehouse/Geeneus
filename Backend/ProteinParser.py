# Deals with cache and storing of protein data objects, as well as coordinating network access  (private)
#
# Copyright 2012 by Alex Holehouse - see LICENSE for more info
# Contact at alex.holehouse@wustl.edu

import sys
import re
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
            
            # if we can be sure this type of ID will not return a protein
            # because its an invald accession number
            if ID_type(ProteinID)[0] == -1:
                print "Warning - The ID {ID} is an invalid accession number, and the database will not be queried".format
                # by returning the object associated with [] we don't pollute the datastore with invalid and pointless
                # searches, we avoid queriying NCBI without a hope in hell of a hit, and we take advantage of the built in
                # bad XML behaviour without raising an error, because, technically, no error has happened, we just know
                # that the ID in question won't return protein data. It's not an error - it's just stupid.
                return ProteinObject.ProteinObject([])

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
        
        if self.protein_datastore[ProteinID].error():
            print "Warning - protein associated with {ID} could not be retrieved due to an error".format(ID=ProteinID)

        elif not self.protein_datastore[ProteinID].exists():
            print "Warning - despite searching through the database, NCBI does not appear to have a protein associated with the ID {ID}".format(ID=ProteinID)
            
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
    def get_variants(self, ID):
        ProtObj = self.__get_protein_object(ID)
        return ProtObj.get_variants()

#--------------------------------------------------------
#
#--------------------------------------------------------
    def get_geneID(self, ID):
        ProtObj = self.__get_protein_object(ID)
        return ProtObj.get_geneID()


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
            
     # returns 0 for GI
     # returns 1 for refseq
     # returns 2 for swissprot
     # returns 3 for anything else
    
    def ID_type(ProteinID):
        
        # if the ID is all digits it's a GI
        if str(ProteinID).isdigit() or re.match("gi|",ProteinID):
            return [0, "GI"]
        
        # if it begin [A|N|X|Y|Z]P then it's a refseq 
        if re.match("[ANXYZ][P]", ProteinID):
            return [1, "RefSeq"]
        
        # if it begins [O|P|Q] then it's a swissprot
        if re.match("[OPQ]", ProteinID):
            return [2, "Swissprot"]
        
        # DDBJ
        if re.match("BAA-BZZ", ProteinID) or re.match("FAA-FZZ", ProteinID) or \
                re.match("GAA-GZZ", ProteinID) or re.match("IAA-IZZ", ProteinID):
            return [3, "DDBJ"]

        # GenBank
        if re.match("AAA-AZZ", ProteinID) or re.match("AAE", ProteinID) or \
                re.match("DAA-DZZ", ProteinID) or re.match("EAA-EZZ", ProteinID) or \
                re.match("HAA-HZZ", ProteinID) or re.match("JAA-JZZ", ProteinID):  
            return [4, "GenBank"]

        
            
# +-------------------------------------------------------+
# |                    END OF CLASS                       |
# +-------------------------------------------------------+


def ID_type(ProteinID):
    # if the ID is all digits it's a GI
    if ProteinID.isdigit() or re.match("gi",ProteinID):
        return [0, "GI"]
        
    # if it begin [A|N|X|Y|Z]P then it's a refseq 
    if re.match("[ANXYZ][P]", ProteinID):
        return [1, "RefSeq"]
        
    # if it begins [O|P|Q] then it's a swissprot
    if re.match("[OPQ]", ProteinID):
        return [2, "Swissprot"]
    
    # DDBJ
    if re.match("BAA-BZZ", ProteinID) or re.match("FAA-FZZ", ProteinID) or \
            re.match("GAA-GZZ", ProteinID) or re.match("IAA-IZZ", ProteinID):
        return [3, "DDBJ"]

    # GenBank
    if re.match("AAA-AZZ", ProteinID) or re.match("AAE", ProteinID) or \
        re.match("DAA-DZZ", ProteinID) or re.match("EAA-EZZ", ProteinID) or \
        re.match("HAA-HZZ", ProteinID) or re.match("JAA-JZZ", ProteinID):  
        return [4, "GenBank"]

    if re.match("CAA-CZZ", ProteinID):
        return [5, "EMBL"]

    return [-1, "Unknown protein accession type"]
    

    
