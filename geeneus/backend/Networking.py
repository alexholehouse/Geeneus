# Contains networking functionality  (private)
#
# Abstracts away any interaction with the eUtil tools from the user, and deals with network errors or other problems.
#
# NOTE:
# timeout code based on code from http://pguides.net/python-tutorial/python-timeout-a-function/
# -
# Copyright 2012 by Alex Holehouse - see LICENSE for more info
# Contact at alex.holehouse@wustl.edu

import datetime
import sys
import signal
import time
import urllib2
import ProteinParser

from Bio import Entrez

#--------------------------------------------------------
# Global networking timeout limits are in ProteinParser
# for fine tuning

######################################################### 
######################################################### 
# Exception class for timeouts
#
class TimeoutException(Exception):
    pass



######################################################### 
######################################################### 
# Main class to handle all the networking shinanigans
#
class Networking:    
    
    TIMEOUT=40
    def __init__(self, timeout):
        self.lastDatabaseCall = datetime.datetime.now()
        
#--------------------------------------------------------
#
#--------------------------------------------------------
# function to decorate other functions with a timeout. If the timeout is reached, causes the
# decorated function to return a -1 along with a printed warning
#    
    def timeout(timeout_time, default):
        def timeout_function(f):
            def f2(*args):

                def timeout_handler(signum, frame):
                    raise TimeoutException()

                handler = signal.signal(signal.SIGALRM, timeout_handler)
                signal.alarm(timeout_time)
                            
                try:
                    retval = f(*args)
                except TimeoutException:
                    print "\nWarning: Timeout reached after {time} seconds\n".format(time=timeout_time)
                    return default
                finally:
                    signal.signal(signal.SIGALRM, handler)
                    signal.alarm(0)
                return retval
            return f2
        return timeout_function
    

#--------------------------------------------------------
#
#--------------------------------------------------------
# function to ensure we stay with NCBI's query limit of no more than 3 per second
#
#
    def stay_within_limits(self):
        if (datetime.datetime.now() - self.lastDatabaseCall).seconds  < 1:
            if (datetime.datetime.now() - self.lastDatabaseCall).microseconds < 400000:
                time.sleep(0.5)
            
        self.lastDatabaseCall = datetime.datetime.now()

#--------------------------------------------------------
#
#--------------------------------------------------------
# efetch nucleotide sequence
#
# Function to get a live handle with nucleotide xml data. All networking
# issue should be dealt with here and abstracted totally from the user
# Decorator must decorate this function (not efetchGeneral) to avoid keyword
# conflicts
#
# A paired __internal_efNT() and efetchNucleotide() set of functions are used
# to allow decoration of one with a timeout, where the efetchNucleotide() can 
# print an error message on -1 return from EITHER the efetchGeneral function,
# or from the timeout decorator itself
#
    @timeout(TIMEOUT, -1)
    def __internal_efNT(self, GI, start, end, strand_val):
        return self.eUtilsGeneral({'function':Entrez.efetch,'db':"nucleotide", 'id':GI, 'seq_start':start, 'seq_stop':end, 'rettype':"fasta", 'strand':strand_val})
      
    def efetchNucleotide(self, GI, start, end, strand_val):
        self.stay_within_limits()
        handle = self.__internal_efNT(GI, start, end, stand_val)
    
        if (handle == -1):
            print "[NCBI]: Networking Error: Problem getting Nucleotide data for GI|{gi}".format(gi=GI)
            return -1
        else:
            return handle
#--------------------------------------------------------
#
#--------------------------------------------------------
#  efetch gene record
#
# Function to get a live handle with gene xml data. All networking
# issue should be dealt with here and abstracted totally from the user
# Decorator must decorate this function (not efetchGeneral) to avoid keyword
# conflicts
#
    @timeout(TIMEOUT, -1)
    def __internal_efG(self, GeneID):
        return self.eUtilsGeneral({'function':Entrez.efetch,'db':"gene", 'id':GeneID, 'rettype':"gene_table", 'retmode':"xml"})

    def efetchGene(self, GeneID):
        self.stay_within_limits()
        handle = self.__internal_efG(GeneID)
        if (handle == -1):
            print "[NCBI]: Network Error: Problem getting gene  data for ID: {GID}".format(GID=GeneID)
            return -1
        else:
            return handle

#--------------------------------------------------------
#
#--------------------------------------------------------
# efetch protein record
#
# Function to get a live handle with protein xml data. All networking
# issue should be dealt with here and abstracted totally from the user
# Decorator must decorate this function (not efetchGeneral) to avoid keyword
# conflicts
#
    @timeout(TIMEOUT, -1)
    def __internal_efP(self, ProteinID):
        return self.eUtilsGeneral({'function':Entrez.efetch,'db':"protein", 'id':ProteinID, 'retmode':"xml"})
       
    def efetchProtein(self, ProteinID):
        self.stay_within_limits()
        handle = self.__internal_efP(ProteinID)
        if (handle == -1):
            print "[NCBI]: Networking Error: Problem getting protein data for ID(s): {PID}".format(PID=ProteinID)
            return -1
        else:
            return handle

#--------------------------------------------------------
#
#--------------------------------------------------------
# epost protein list (not really used, but here incase we implement GI based
# asynchronous fetching in the future)
#
# Function to get post a list of IDs to the NCBI server
# for asynchronous processing. As of 23 Oct 2012 this is
# not being used, but is kept in case we add asynchronous
# epost based features in the future
#
    @timeout(TIMEOUT, -1)
    def __internal_epP(self, ProteinIDList):
        return self.eUtilsGeneral({'function':Entrez.epost,'db':"protein", 'id':",".join(ProteinIDList)})

    def epostProtein(self, ProteinIDList):
        self.stay_within_limits()
        handle = self.__internal_epP(ProteinIDList)
        if (handle == -1):
            print "[NCBI]: Networking Error: Problem ePosting ID(s): {PID}".format(PID=ProteinIDList)
            return -1
        else:
            return handle


#--------------------------------------------------------
#
#--------------------------------------------------------
# esearch
#
# Function to query the protein dtabase for the $passedTerm
# term 
#

    @timeout(TIMEOUT, -1)
    def __internal_esP(self, passedTerm):
        return self.eUtilsGeneral({'function':Entrez.esearch,'db':"protein", 'term':passedTerm})

    def esearchProtein(self, term):
        self.stay_within_limits()
        handle = self.__internal_esP(term)
        if (handle == -1):
            print "[NCBI]: Networking Error: Problem eSearching for term: {PID}".format(PID=term)
            return -1
        else:
            return handle

#--------------------------------------------------------
#
#--------------------------------------------------------
# Generic function to make some connection to the NCBI database
#
# While it would make more sense to check we're staying within the
# NCBI limits here (self.stay_within_limits()) this would mess up
# the multithreading used for the timeout decorator, so we use
# stay_within_limits() before we activate the timeout.
#
#
#
    def eUtilsGeneral(self, inputDictionary):
        function = inputDictionary["function"]
        
        del inputDictionary["function"]
        try:
            handle = function(**inputDictionary)
        except urllib2.HTTPError, err:
            print "[NCBI]: HTTP error({0}): {1}".format(err.code, err.reason)
            return -1 
        except urllib2.URLError, err:
            try:
                print "[NCBI]: URLError error({0}): {1}".format(err.code, err.reason)
            except AttributeError, err:
                print "[NCBI]: Corrupted urllib2.URLError raised"
                return -1
            return -1
        return handle



###############################################################################################
## UNIPROT NETWORKING FUNCTIONS
###############################################################################################

    @timeout(TIMEOUT, -1)
    def __internal_UniprotNR(self, queryString):
        try:
            handle = urllib2.urlopen(queryString)
        except urllib2.URLError, err:
            try:
                print "[UniProt]: URLError error({0}): {1}".format(err.code, err.reason)
            except AttributeError, err:
                print "Corrupted urllib2.URLError raised"
                return -1
            return -1
        return handle
       


    def UniProtNetworkRequest(self, accessionID):
        baseURL = 'http://www.uniprot.org/uniprot/'
        queryString = baseURL+str(accessionID)+'.xml'
        
        # probably good to set some kind of limit
        self.stay_within_limits()
        
        return self.__internal_UniprotNR(queryString)

###############################################################################################
## PFAM NETWORKING FUNCTIONS
###############################################################################################

    @timeout(TIMEOUT, -1)
    def __internal_PfamNR(self, queryString):
        
        try:
            handle = urllib2.urlopen(queryString)
        except urllib2.URLError, err:
            try:
                print "[Pfam]: URLError error({0}): {1}".format(err.code, err.reason)
            except AttributeError, err:
                print "Corrupted urllib2.URLError raised"
                return -1
            return -1
        return handle
       

    def PfamNetworkRequest(self, accessionID, limit=True):
        baseURL = "http://pfam.sanger.ac.uk/protein?output=xml&acc="
        queryString = baseURL+str(accessionID)
        
        # probably good to set some kind of limit
        self.stay_within_limits()
        
        return self.__internal_PfamNR(queryString)

            

    
