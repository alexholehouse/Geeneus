# Deals with cache and storing of protein data objects, as well as coordinating network access  (private)
#
# Copyright 2012-2015 by Alex Holehouse - see LICENSE for more info
# Contact at alex.holehouse@wustl.edu
import time
import sys
import re
import pickle

# Biopython selective imports
import Bio.Entrez.Parser
from Bio import Entrez
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from xml.dom.minidom import parseString
import ProteinObject
import Networking
import httplib
import Parser as GRP
import UniprotAPI
from DataStructures import CaseInsensitiveDict as CID

IDTYPES = {-2:"International Protein Index",\
                -1:"Unknown protein accession type", \
                0:"GI", \
                1:"RefSeq", \
                2:"UniProtKB/Swiss-Prot", \
                3:"DDBJ", \
                4:"GenBank",\
                5:"EMBL",\
                6:"PDB",\
                7:"UniProtKB/Swiss-Prot",\
                8:"Unknown protein accession type",\
                }



######################################################### 
#########################################################
class ProteinRequestParser(GRP.GeneralRequestParser):

# Initialization function, set Entrez.email for calls, an ensure key -1 is set
# to a non-existant object
#
    def __init__(self, email, cache, retry=0, loud=True, shortcut=True):
        """"Initializes an empty requestParser object, setting the Entrez.email field and defining how many times network errors should be retried"""
        try:
            GRP.GeneralRequestParser.__init__(self, email, cache, retry, loud)
        
            self.protein_datastore = {-1 : ProteinObject.ProteinObject(-1, [])}
            self.protein_translationMap = {-1: --1}
            self.batchableFunctions = [self.get_sequence, self.get_protein_name, self.get_variants, self.get_geneID, self.get_protein_sequence_length]
            self.UniprotAPI = UniprotAPI.UniprotAPI()
            
            self.shortcut = shortcut
            self.error_status = False
             
        except Exception, e: 
            print "Fatal error when creating ProteinRequestParserObject"
            self.error_status = True
            raise e



#--------------------------------------------------------
# PUBLIC FUNCTION
#--------------------------------------------------------
# return a copy of the list of keys (without -1 key)
#
    def keys(self):
        kl = []
        for i in self.protein_datastore:
            kl.append(i)

        kl.remove(-1)
        return kl

#--------------------------------------------------------
# PUBLIC FUNCTION
#--------------------------------------------------------
# return a copy of the list of keys (without -1 key)
#
    def accession_classes(self):
        return IDTYPES.values()

#--------------------------------------------------------
# PUBLIC FUNCTION
#--------------------------------------------------------
# Test if the datatstore has a key (does not trigger download
# on fail

    def has_key(self, ID):
        return self.protein_datastore.has_key(ID)

            
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
# Check for an error in the proteinObject obtained from
# the accession ID
    def protein_error(self, ID):
       ProtObj = self._get_protein_object(ID)
       return ProtObj.error()


#--------------------------------------------------------
# PUBLIC FUNCTION
#--------------------------------------------------------
#
# Check for an error in the proteinObject obtained from
# the accession ID
    def protein_exists(self, ID):
       ProtObj = self._get_protein_object(ID)
       return ProtObj.exists()

#--------------------------------------------------------
# PUBLIC FUNCTION
#--------------------------------------------------------
#
# Check for an error in the proteinObject obtained from
# the accession ID
    def get_creation_date(self, ID):
       ProtObj = self._get_protein_object(ID)
       return ProtObj.get_creation_date()


#--------------------------------------------------------
# PUBLIC FUNCTION
#--------------------------------------------------------
#
# Check for an error in the proteinObject obtained from
# the accession ID
    def get_protein_source(self, ID):
       ProtObj = self._get_protein_object(ID)
       return ProtObj.source()

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
# Get the protein sequence version
    def get_record_version(self, ID):
        return (self._get_protein_object(ID)).get_version()


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
# Get the sequence length of the protein's AA sequence
    def get_other_accessions(self, ID):
        ProtObj = self._get_protein_object(ID)
        return ProtObj.get_other_accessions()


#--------------------------------------------------------
# PUBLIC FUNCTION
#--------------------------------------------------------
#   
# Get the sequence length of the protein's AA sequence
    def get_species(self, ID):
        ProtObj = self._get_protein_object(ID)
        return ProtObj.get_species()

#--------------------------------------------------------
# PUBLIC FUNCTION
#--------------------------------------------------------
#   
# Get the sequence length of the protein's AA sequence
    def get_taxonomy(self, ID):
        ProtObj = self._get_protein_object(ID)
        return ProtObj.get_taxonomy()

#--------------------------------------------------------
# PUBLIC FUNCTION
#--------------------------------------------------------
#   
# Get the sequence length of the protein's AA sequence
    def get_host(self, ID):
        ProtObj = self._get_protein_object(ID)
        return ProtObj.get_host()

#--------------------------------------------------------
# PUBLIC FUNCTION
#--------------------------------------------------------
#   
# Get the sequence length of the protein's AA sequence
    def get_domains(self, ID):
        ProtObj = self._get_protein_object(ID)
        return ProtObj.get_domains()

#--------------------------------------------------------
# PUBLIC FUNCTION
#--------------------------------------------------------
#   
# Get the sequence length of the protein's AA sequence
    def get_gene_name(self, ID):
        ProtObj = self._get_protein_object(ID)
        return ProtObj.get_gene_name()

#--------------------------------------------------------
# PUBLIC FUNCTION
#--------------------------------------------------------
#
# Get any isoforms associated with this sequence ID

    def get_isoforms(self, ID):
        ProtObj = self._get_protein_object(ID)
        return ProtObj.get_isoforms()
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
        self.protein_datastore = {-1 : ProteinObject.ProteinObject(-1, [])}

#--------------------------------------------------------
# PUBLIC FUNCTION
#--------------------------------------------------------
# Get the number of items in the datastore
#
    def datastore_size(self):
        return len(self.protein_datastore)-1



#--------------------------------------------------------
# PUBLIC FUNCTION
#--------------------------------------------------------
# Save the datastore to file
#
    def save_datastore(self,filename='datastore'):

        with open(str(filename)+".pi",'w') as FH:
            pickle.dump(self.protein_datastore, FH)


#--------------------------------------------------------
# PUBLIC FUNCTION
#--------------------------------------------------------
# Load a datastore from file
#

    def load_datatore(self,filename):
        with open(str(filename),'r') as FH:
            self.protein_datastore = pickle.load(FH)








#--------------------------------------------------------
# PUBLIC FUNCTION
#--------------------------------------------------------
#

    def batchFetch(self, function, listOfIDs):
        outputDict = {}
        toFetch = []

        # check the function actually makes sense taking
        # a single ID element as input
        if function not in self.batchableFunctions:
            print "Warning, function cannot be run via batch"
            return

        # Firstly, we identify which, if any of these are already in the 
        # datastore, and for those which are not ignore badly formatted
        # accession numbers
        for proteinID in listOfIDs:
            if proteinID not in self.protein_datastore or not self.cache:
                proteinID = self._convertIfNecessary(proteinID)
                if not ID_type(proteinID)[0] < 0:
                    toFetch.append(proteinID)                    

        
        
        # if we're running using shortcutting we have to now split to toFetch list into a
        # shortcutable list and a non-shortcutable list, deal with each list seperatly, them
        # re-assemble them into the same order.
        if self.shortcut:
            trackerVector = []
            UniProtList = []
            NCBIList = []
            listOfXML = []
            
            # build two distinct lists and a tracking vector which describes the mapping
            # of toFetch to those two lists
            for ID in toFetch:
                if ID_type(ID)[0] == 7:
                    trackerVector.append("U")
                    UniProtList.append(ID)
                else:
                    trackerVector.append("N")
                    NCBIList.append(ID)

            # build a XMLlist using the NCBI batch lookup function
            NCBI_XML = self._get_batch_XML(NCBIList, self.Networking.efetchProtein, self.UniprotDatabaseLookup)
            
            # silently adds records to database
            self.UniprotAPI.batchFetch(UniProtList, self.protein_datastore)
            
            # now reconstruct listOfXML as a list the same order as toFetch, but with all the would-be UniProt
            # XML values equal to -1
            NCBI_counter = 0
            for i in trackerVector:
                if i == "U":
                    listOfXML.append(-1)
                if i == "N":
                    listOfXML.append(NCBI_XML[NCBI_counter])
                    NCBI_counter = NCBI_counter+1

        else:
                           
            # next we take those which are NOT in the datastore and take advantage
            # of the Biopython.Entrez' batch downloaded function. This makes a SINGLE
            # call to the server, so is a lot faster (reduces setup and teardown)
            # generates a list, each element of which is the 1:1 xml for the 
            # listOfIDs
            listOfXML = self._get_batch_XML(toFetch, self.Networking.efetchProtein, self.UniprotDatabaseLookup)
        
        
        # note we don't have to try/catch here because _get_batch_XML() guarentees
        # that each XML field is valid
        fetchCounter = 0
        
        for protein_xml in listOfXML:
            if not self.protein_datastore.has_key(toFetch[fetchCounter]) and not protein_xml == -1 :
                self.protein_datastore[toFetch[fetchCounter]] = ProteinObject.ProteinObject(toFetch[fetchCounter], [protein_xml])
            
            fetchCounter = fetchCounter+1

                        
        # finally, we build a tuple with results from everything in the input, and
        # return. Note we're ONLY doing this if that ID has already been loaded into
        # the data store. This stops the retry mechanism going through for batch, failing
        # and then just retrying again via the _get_protein_object() function which the
        # $function inevitably will call
        #
        # By checking against the datastore we also ensure we can pick up additional records
        # built by the "alternative" function, which adds to the datastore in an independent manner
        # to the for loop above 
            
        for ID in listOfIDs:
            if ID in self.protein_datastore:
                outputDict[ID] = function(ID)
                
            else:
                outputDict[ID] = function(-1)

                
        # build a case insensitive read only dictionary to output
        outFinal = CID(outputDict)

        return outFinal

#--------------------------------------------------------
# PRIVATE FUNCTION FUNCTION
#--------------------------------------------------------
#

    def _get_protein_object(self, proteinID):

        proteinID = self._convertIfNecessary(proteinID)

        # pre-check to see if the accession matches the predefined accession format
        # rules. If it doesn't then 
        if ID_type(proteinID)[0] < 0:
            if not proteinID == -1: 
                self.printWarning("\nWarning - The ID {ID} is an invalid accession number, and the database will not be queried".format(ID=proteinID))
            proteinID = -1

        else:
            # note _get_object is inherited from Parser
            self._get_object(proteinID, self.protein_datastore, self._protein_fetch_function, ProteinObject.ProteinObject, self.UniprotDatabaseLookup)
        
        return self.protein_datastore[proteinID]

#--------------------------------------------------------
# PRIVATE FUNCTION FUNCTION
#--------------------------------------------------------
# The protein fetch function allows an additional layer
# of logic when deciding how to deal with an accession value. 
#
#
    def _protein_fetch_function(self, accessionID):

               
        # if we have shortcut off then just use eFetch. If eFetch fails after a number
        # of retries we'll try UniProt if the accession is a UniProt one
        if self.shortcut == False:
            returnVal = self.Networking.efetchProtein(accessionID)
        
        # else we have shortcutting switched on, so if the accession is
        # a uniprot accession which NCBI doesn't guarentee to support we
        # automaically drop 
        else:
            if ID_type(accessionID)[0] == 7:

                
                # because we provide the getProteinObjectFromUniProt with a -1 value for 
                # the datastore, the function will return an instantiated ProteinObject
                # if possible (or -1 on failure)
                returnVal = self.UniprotAPI.getProteinObjectFromUniProt(-1, accessionID)

            else:
                returnVal = self.Networking.efetchProtein(accessionID)
               
        return returnVal
 
                


#--------------------------------------------------------
# PRIVATE FUNCTION
#--------------------------------------------------------          
# Function which we can manually add cases to where converting
# from some kind of accession number to a GI would be better
# for whatever reason. Currently only true for PDB IDs
#

    def _convertIfNecessary(self, ProteinID):
        
        # Convert to upper case for tests so we have
        # a homogenous string
        try:
            ProteinID = ProteinID.upper()
        except AttributeError:
            return ProteinID

        # if PDB translate to a GI value
        if ID_type(ProteinID)[0] == 6:
            return self.translate_Asc2GI(ProteinID)
        
        # if SwissProt/UniProt with an isoform identifier
        # chop off identifier
        if re.match("[A-Z][0-9][A-Z0-9][A-Z0-9][A-Z0-9][0-9]-[0-9]", ProteinID):
            return ProteinID[:6]
        
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


#--------------------------------------------------------
# PRIVATE FUNCTION
#--------------------------------------------------------
# Function which gives the GeneralRequestParser access to
# non NCBI based means to get a protein record (Uniprot)
#

    def UniprotDatabaseLookup(self, accessionID):
        
        IDtype = ID_type(accessionID)[0]
        
        # if shortcut is on and we're calling the alt function
        # we don't want to try UniProt *again* because it must
        # have not worked the first few times
        if self.shortcut:
            exitcodes = 2,

        # however if shortcut is off then we do want to test any
        # UniprotKB/Swissprot accession
        else:
            exitcodes = 2,7
        
        if IDtype in exitcodes:
            print "[UniProt]: Falling back and querying UniProt servers...  "
        
            # query the server through the UniprotAPI class
            self.UniprotAPI.getProteinObjectFromUniProt(self.protein_datastore, accessionID)
                    
            # NOTE this is crucial - we only return true if we can get an object 
            # from the dictionary
            try:
                self.protein_datastore[accessionID]
                print "[UniProt]: Sucess!"
                return True
            except KeyError:
                return False
        else:
            return False


    def loadfile(self, filename):
        with open(filename) as openFile:
            print "Reading in file..."
            xmlstring = openFile.read()
           

            compositeDOM = parseString(xmlstring)
           
        print "Done"
        print "Extracting protein data..."
        elementsList = compositeDOM.getElementsByTagName("entry")
        
        
        total = len(elementsList)
        counter = 0 
        
        "Parsing XML..."

        for element in elementsList:
            print "Parsing " + str(counter) +"/"+str(total)
            self.UniprotAPI.sideload_from_file(self.protein_datastore, element)
            counter = counter+1
                
        

# +-------------------------------------------------------+
# |                    END OF CLASS                       |
# +-------------------------------------------------------+


#--------------------------------------------------------
# PRIVATE FUNCTION
#--------------------------------------------------------
## ProteinParser function
## General function true to the module, as does not require
## a class instance to work
#
# List of for refrence/addition
# -2 : IPI - IPI is -2 because we can't look it up, but it 
#            is a valid type of accession. Neither NCBI nor
#            Uniprot support IPI lookups, however!
# -1 : unknown
#  0 : GI
#  1 : refseq
#  2 : Swissprot
#  3 : DDBJ
#  4 : genbank
#  5 : EMBL
#  6 : PDB
#  7 : Uniprot


def ID_type(ProteinID):
    
    # convert ID to all uppercase characters, and remove any periods 
    try:
        ProteinID = ProteinID.upper()
    except AttributeError:
        return ([-1, IDTYPES[-1]])


    # cut off . and anything after
    # "." represent versions of accession numbers, but in this context 
    # the versions a) don't change the type and b) just make all the
    # regexes way harder!

    if ProteinID.find(".") > 0:
        ProteinID = ProteinID[:ProteinID.find(".")]

        
    # if it begin [A|N|X|Y|Z]P_ then it's a refseq 
    if re.match("[ANXYZ][P]_", ProteinID):
        return [1, IDTYPES[1]]

    # if the ID is all digits it's a GI
    if ProteinID.isdigit() or re.match("GI",ProteinID):
        return [0, IDTYPES[0]]

    # is it a PDB?
    if re.match("[0-9][A-Z0-9][A-Z0-9][A-Z0-9]", ProteinID) or re.match("[0-9][A-Z0-9][A-Z0-9][A-Z0-9]_[A-Z0-9]", ProteinID):
        return[6, IDTYPES[6]]
      
    # if it begins [O|P|Q] then it is found in NCBI database
    if re.match("[OPQ]", ProteinID) and re.match("^[A-Z0-9]+$", ProteinID) and len(ProteinID) == 6:
        return [2, IDTYPES[2]]

    # else if its a uniprot/swissprot is is not necessarily found so we fall back to UniProt
    # api
    if re.match("[A-N|R-Z][0-9][A-Z0-9][A-Z0-9][A-Z0-9][0-9]", ProteinID) and len(ProteinID) == 6:
        return [7, IDTYPES[7]]

    if re.match("IPI[0-9]*", ProteinID):
        return [-2, IDTYPES[-2]]
    
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
        return [3, IDTYPES[3]]    

    # GenBank
    if re.match("A[A-Z][A-Z]", ProteinID) or re.match("AAE", ProteinID) or \
        re.match("D[A-Z][A-Z]", ProteinID) or re.match("E[A-Z][A-Z]", ProteinID) or \
        re.match("H[A-Z][A-Z]", ProteinID) or re.match("J[A-Z][A-Z]", ProteinID):  
        return [4, IDTYPES[4]]

    if re.match("C[A-Z][A-Z]", ProteinID):
        return [5, IDTYPES[5]]
    
    # if we get here something seems to have gone wrong...
    return [-1, IDTYPES[-1]]
    

    
