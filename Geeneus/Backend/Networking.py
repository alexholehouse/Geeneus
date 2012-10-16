# Contains networking functionality  (private)
#
# Abstracts away any interaction with the eUtil tools from the user, and deals with network errors or other problems.
#
# NOTE:
# timeout code based on code from http://pguides.net/python-tutorial/python-timeout-a-function/
# -
# Copyright 2012 by Alex Holehouse - see LICENSE for more info
# Contact at alex.holehouse@wustl.edu

import sys
import signal
import urllib2

from Bio import Entrez

#--------------------------------------------------------
# Global networking timeout limits
EFETCH_TIMEOUT = 20
ESEARCH_TIMEOUT = 20

#--------------------------------------------------------
#
#--------------------------------------------------------
# Exception class for timeouts
#
class TimeoutException(Exception):
    pass

#--------------------------------------------------------
#
#
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
@timeout(EFETCH_TIMEOUT, -1)
def __internal_efNT(GI, start, end, strand_val):
    return efetchGeneral(db="nucleotide", id=GI, 
                         seq_start=start, 
                         seq_stop=end, rettype="fasta", 
                         strand=strand_val)

def efetchNucleotide(GI, start, end, strand_val):
    handle = __internal_efNT(GI, start, end, stand_val)
    
    if (handle == -1):
        print "Networking Error: Problem getting Nucleotide data for GI|{gi}".format(gi=GI)
        return -1
    else:
        return handle

#--------------------------------------------------------
#
#--------------------------------------------------------
# Function to get a live handle with gene xml data. All networking
# issue should be dealt with here and abstracted totally from the user
# Decorator must decorate this function (not efetchGeneral) to avoid keyword
# conflicts
#
@timeout(EFETCH_TIMEOUT, -1)
def __internal_efG(GeneID):
    return efetchGeneral(db="gene", id=GeneID, rettype="gene_table", retmode="xml")

def efetchGene(GeneID):
    handle = __internal_efG(GeneID)
    if (handle == -1):
        print "Network Error: Problem getting gene  data for ID: {GID}".format(GID=GeneID)
        return -1
    else:
        return handle

#--------------------------------------------------------
#
#--------------------------------------------------------
# Function to get a live handle with protein xml data. All networking
# issue should be dealt with here and abstracted totally from the user
# Decorator must decorate this function (not efetchGeneral) to avoid keyword
# conflicts
#
@timeout(EFETCH_TIMEOUT, -1)
def __internal_efP(ProteinID):
    return efetchGeneral(db="protein", id=ProteinID,  retmode="xml")

def efetchProtein(ProteinID):
    handle = __internal_efP(ProteinID)
    if (handle == -1):
        print "Networking Error: Problem getting protein data for ID: {PID}".format(PID=ProteinID)
        return -1
    else:
        return handle

#--------------------------------------------------------
#
#--------------------------------------------------------
# Actual eFetch call to the NCBI database occurs here. Deal with network
# errors in this function, returning -1 if call fails
#
# Note the multiple tiers of try/except statements to avoid the system crashing
#
def efetchGeneral(**kwargs):
    try:
        handle = Entrez.efetch(**kwargs)
    except urllib2.HTTPError, err:
        print "HTTP error({0}): {1}".format(err.code, err.reason)
        return -1 
    except urllib2.URLError, err:
        try:
            print "URLError error({0}): {1}".format(err.code, err.reason)
        except AttributeError, err:
            print "Corrupted urllib2.URLError raised"
            return -1
        return -1
    
    return handle

#--------------------------------------------------------
#
#--------------------------------------------------------
# Function to search Entrez using an accession value
# Returns a raw protein record, or -1 if there is a problem
#
@timeout(ESEARCH_TIMEOUT, -1)
def esearch(db, Accession):
    try:
        handle = Entrez.esearch(db="protein", term=Accession)
    except urllib2.URLError, (err):
        print "Entrez.esearch error: Unable to open URL"
        return -1
    protein_record = Entrez.read(handle)

    return protein_record
