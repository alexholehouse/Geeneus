# Deals with cache and storing of protein data objects, as well as coordinating network access  (private)
#
# Copyright 2012 by Alex Holehouse - see LICENSE for more info
# Contact at alex.holehouse@wustl.edu
import time
import sys
import re

# Biopython selective imports
import Bio.Entrez.Parser
from Bio import Entrez
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

import ProteinObject
import Networking
import httplib
import Parser as GRP

######################################################### 
#########################################################
class ProteinRequestParser(GRP.GeneralRequestParser):

# Initialization function, set Entrez.email for calls, an ensure key -1 is set
# to a non-existant object
#
    def __init__(self, email, cache, retry=0, loud=True):
        """"Initializes an empty requestParser object, setting the Entrez.email field and defining how many times network errors should be retried"""
        try:
            GRP.GeneralRequestParser.__init__(self, email, cache, retry, loud)
            
            self.protein_datastore = {-1 : ProteinObject.ProteinObject([])}
            self.protein_translationMap = {-1: -1}
            self.batchableFunctions = [self.get_sequence, self.get_protein_name, self.get_variants, self.get_geneID, self.get_protein_sequence_length]
            
            self.error_status = False
             
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
# TO DO - get linear AA sequence distance

    def get_distance_between_residues(self, ID, R1, R2):
        print "!!!! NOT YET IMPLEMENTED !!!!"
        return 0;


#--------------------------------------------------------
# PUBLIC FUNCTION
#--------------------------------------------------------
#
# Function to translate and accession value into a GI. If 
# there is more than one possible GI, then simply returns
# -1. Caches lookup values so we don't have to repeatedly 
# run the translation remotely

    def translate_Asc2GI(self, Accession):
        
        retry = self._build_retry_function(self.Networking.esearchProtein)
      
        if not Accession in self.protein_translationMap:
            record = -1

            while (record == -1):
                record = retry(Accession)

            # if we retry a bunch of times and it doesn't work
            if record == -2:
                print "Unable to carry out esearch for term{y}".format(y=Accession)
                return -1
                        
            # else we got some XML
            IdList = record["IdList"]

            # if there's only one GI associated with this accession number then great!
            if len(IdList) == 1:
                self.protein_translationMap[Accession] = str(record["IdList"][0])

            # if not
            else:
                print "There are {op} possible options, shown below. For PDB values, this often arises because seperate chains are treated as different proteins".format(op=len(IdList))
                for i in IdList:
                    print i
                return -1

        return self.protein_translationMap[Accession]


#--------------------------------------------------------
# PUBLIC FUNCTION
#--------------------------------------------------------
#
# Function to get the raw XML 
    def get_raw_xml(self, ProteinID):

        ProtObj = self._get_protein_object(ProteinID)
        return ProtObj.get_raw_xml()


#--------------------------------------------------------
# PUBLIC FUNCTION
#--------------------------------------------------------
# Function to get the raw XML 
# 
    def get_ID_type(self, ProteinID):
        return ID_type(ProteinID)


#--------------------------------------------------------
# PUBLIC FUNCTION
#--------------------------------------------------------
# Function which purges the datastore
#
    def purge_data_store(self):
        del self.protein_datastore
        self.protein_datastore = {-1 : ProteinObject.ProteinObject([])}

#--------------------------------------------------------
# PUBLIC FUNCTION
#--------------------------------------------------------
# Get the number of items in the datastore
#
    def get_size_of_datastore(self):
        return len(self.protein_datastore)-1

def get_protein_object(self, proteinID):

    proteinID = self._convertIfNecessary(proteinID)

    # pre-check to see if the accession matches the predefined accession format
    # rules. If it doesn't then 
    if ID_type(proteinID)[0] == -1:
        proteinID = -1
        return(self._get_object(proteinID, self.protein_datastore, self.Networking.efetchProtein, ProteinObject.ProteinObject))
    
    else:
        self._get_object(proteinID, self.protein_datastore, self.Networking.efetchProtein, ProteinObject.ProteinObject)
        return self.protein_datastore[proteinID]


#--------------------------------------------------------
# PRIVATE FUNCTION
#--------------------------------------------------------          
# Function which we can manually add cases to where converting
# from some kind of accession number to a GI would be better
# for whatever reason. Currently only true for PDB IDs
#

    def _convertIfNecessary(self, ProteinID):
        
        # if PDB
        if ID_type(ProteinID)[0] == 6:
            return self.translate_Asc2GI(ProteinID)
        
        return ProteinID

#--------------------------------------------------------
# PRIVATE FUNCTION
#--------------------------------------------------------
# Function to print a message if loud has been set to 
# true - saves littering code with if self.loud: statements
#
    def printWarning(self, message):
        if self.loud:
            print message


        
        
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
    try:
        ProteinID = ProteinID.upper()
    except AttributeError:
        return ([-1, "Unknown protein accession type"])
        

    # if it begin [A|N|X|Y|Z]P_ then it's a refseq 
    if re.match("[ANXYZ][P]_", ProteinID):
        return [1, "RefSeq"]

    # if the ID is all digits it's a GI
    if ProteinID.isdigit() or re.match("GI",ProteinID):
        return [0, "GI"]

    if re.match("[0-9][A-Z0-9][A-Z0-9][A-Z0-9]", ProteinID) or re.match("[0-9][A-Z0-9][A-Z0-9][A-Z0-9]_[A-Z0-9]", ProteinID):
        return[6, "PDB"]
    
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
        return ([-1, "Unknown protein accession type"])
    
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
    

    
