# Deals with cache and storing of protein data objects, as well as coordinating network access  (private)
#
# Copyright 2012 by Alex Holehouse - see LICENSE for more info
# Contact at alex.holehouse@wustl.edu
import time

import sys
import re

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
        try:
            Entrez.email = email
            self.cache = cache
            self.protein_datastore = {-1 : ProteinObject.ProteinObject([])}
            self.error_status = False
            self.batchableFunctions = [self.get_sequence, self.get_protein_name, self.get_variants, self.get_geneID, self.get_protein_sequence_length]
        except: 
            print "Fatal error when creating ProteinRequestParserObject"
            self.error_status = True
            


#--------------------------------------------------------
# PUBLIC FUNCTION
#--------------------------------------------------------
#
# Check for an error in the parser object

    def error(self):
       return self.error_status


#--------------------------------------------------------
# PUBLIC FUNCTION
#--------------------------------------------------------
#
# Get the protein's name

    def get_protein_name(self, ID):
        return (self._get_protein_object(ID)).get_protein_name()

#--------------------------------------------------------
# PUBLIC FUNCTION
#--------------------------------------------------------
#
# Get the AA sequence of the protein (N-to-C)
    def get_sequence(self, ID):
        ProtObj = self._get_protein_object(ID)
        return ProtObj.get_protein_sequence()
                                       
#--------------------------------------------------------
# PUBLIC FUNCTION
#--------------------------------------------------------
#
# get a list of the single AA change varants for this
# species
    def get_variants(self, ID):
        ProtObj = self._get_protein_object(ID)
        return ProtObj.get_variants()


#--------------------------------------------------------
# PUBLIC FUNCTION
#--------------------------------------------------------
#
# get the GeneID associated with this protein
    def get_geneID(self, ID):
        ProtObj = self._get_protein_object(ID)
        return ProtObj.get_geneID()


#--------------------------------------------------------
# PUBLIC FUNCTION
#--------------------------------------------------------
#   
# Get the sequence length of the protein's AA sequence
    def get_protein_sequence_length(self, ID):
        ProtObj = self._get_protein_object(ID)
        return ProtObj.get_protein_sequence_length()


#--------------------------------------------------------
# PUBLIC FUNCTION
#--------------------------------------------------------
#
# TO DO - get linear AA sequence distnace and 3D distance
# based on PDB structure (where available)

    def get_distance_between_residues(self, ID, R1, R2):
        return 0;


#--------------------------------------------------------
# PUBLIC FUNCTION
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

#--------------------------------------------------------
# PUBLIC FUNCTION
#--------------------------------------------------------
#
# Function to get the raw XML 
    def get_raw_xml(self, ProteinID):
        if not ID_type(ProteinID)[1] == -1:
            return Entrez.read(Networking.efetchProtein(ProteinID))
        else:
            return -1


#--------------------------------------------------------
# PRIVATE FUNCTION
#--------------------------------------------------------
#
# internal function which returns a protein object either from
# the cache (if caching is on and the protein data has already 
# been downloaded) or directly from the NCBI database if caching
# is off or the protein isn't yet in the database.
#
    def _get_protein_object(self, ProteinID):
        if ProteinID not in self.protein_datastore or not self.cache:
            
            
            # if we can be sure this type of ID will not return a protein
            # because its an invald accession number
            if ID_type(ProteinID)[0] == -1:
                print "\nWarning - The ID {ID} is an invalid accession number, and the database will not be queried".format(ID=ProteinID)
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
# PRIVATE FUNCTION
#--------------------------------------------------------
# Function which takes one of the get_? functions defined
# below and an array of IDS of interest, and returns 
# tuple of ID-results while remaining INSIDE NCBI's
# fetch parameters.
#

    def batchFetch(self, function, arrayOfIDs):

        # check the function actually makes sense taking
        # a single ID element as input
        if function not in self.batchableFunctions:
            print "Warning, function cannot be run via batch"
            return

        outputList = {}

        for ID in arrayOfIDs:

            # by setting this to 0.5 we ensure that we cannot
            # go over the NCBI usage limits of more than 3 per
            # second. Best case scenario we hit 2 per second
            time.sleep(0.5)
            outputList[ID] = function(ID)

        return outputList
        

#--------------------------------------------------------
# PUBLIC FUNCTION
#--------------------------------------------------------
# Function which purges the datastore
#

    def purgeDataStore(self):
        del self.datastore
        self.protein_datastore = {-1 : ProteinObject.ProteinObject([])}
        
        
# +-------------------------------------------------------+
# |                    END OF CLASS                       |
# +-------------------------------------------------------+


#--------------------------------------------------------
# PRIVATE FUNCTION
#--------------------------------------------------------
## ProteinParser function
## General function true to the module, as does not require
## a class instance to work
def ID_type(ProteinID):

    # convert ID to all uppercase characters...
    ProteinID = ProteinID.upper()
    
    # if it begin [A|N|X|Y|Z]P_ then it's a refseq 
    if re.match("[ANXYZ][P]_", ProteinID):
        return [1, "RefSeq"]

    # if the ID is all digits it's a GI
    if ProteinID.isdigit() or re.match("GI",ProteinID):
        return [0, "GI"]
    
    # if it begins [O|P|Q] then it's a swissprot
    if re.match("[OPQ]", ProteinID) and re.match("^[A-Z0-9]+$", ProteinID):
        return [2, "Swissprot"]
    
    # Now we've removed refseq, gi and swissprot we can be a bit more discerening about 
    # what the accession should be (as now it must conform to the NCBI accesison number
    # rules). If any of the following fail then we assume the accession is malformed
    # __________________________________
    # Should be 3 letters+5digits
    # Should only contain A-Z
    # Should be less than 8 characters
    # 
    if not re.match("[A-Z][A-Z][A-Z][0-9][0-9][0-9][0-9][0-9]", ProteinID) or \
            not re.match("^[a-zA-Z0-9]+$", ProteinID) or \
            len(ProteinID) > 8:
        
        return [-1, "Unknown protein accession type"]
    
    # DDBJ
    if re.match("[B][A-Z][A-Z]", ProteinID) or re.match("F[A-Z][A-Z]", ProteinID) or \
            re.match("G[A-Z][A-Z]", ProteinID) or re.match("I[A-Z][A-Z]", ProteinID):
        return [3, "DDBJ"]

    # GenBank
    if re.match("A[A-Z][A-Z]", ProteinID) or re.match("AAE", ProteinID) or \
        re.match("D[A-Z][A-Z]", ProteinID) or re.match("E[A-Z][A-Z]", ProteinID) or \
        re.match("H[A-Z][A-Z]", ProteinID) or re.match("J[A-Z][A-Z]", ProteinID):  
        return [4, "GenBank"]

    if re.match("C[A-Z][A-Z]", ProteinID):
        return [5, "EMBL"]
    
    # if we get here something seems to have gone wrong...
    return [-1, "Unknown protein accession type"]
    

    
