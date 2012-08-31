## contains networking functionality
## timeout code based on code from http://pguides.net/python-tutorial/python-timeout-a-function/
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
@timeout(EFETCH_TIMEOUT, -1)
def efetchNucleotide(GI, start, end, strand_val):
    return efetchGeneral(db="nucleotide", id=GI, 
                         seq_start=start, 
                         seq_stop=end, rettype="fasta", 
                         strand=strand_val)

#--------------------------------------------------------
#
#--------------------------------------------------------
# Function to get a live handle with gene xml data. All networking
# issue should be dealt with here and abstracted totally from the user
# Decorator must decorate this function (not efetchGeneral) to avoid keyword
# conflicts
#
@timeout(EFETCH_TIMEOUT, -1)
def efetchGene(GeneID):
    return efetchGeneral(db="gene", id=GeneID, rettype="gene_table", retmode="xml")

#--------------------------------------------------------
#
#--------------------------------------------------------
# Function to get a live handle with protein xml data. All networking
# issue should be dealt with here and abstracted totally from the user
# Decorator must decorate this function (not efetchGeneral) to avoid keyword
# conflicts
#
@timeout(EFETCH_TIMEOUT, -1)
def efetchProtein(ProteinID):
    return efetchGeneral(db="protein", id=ProteinID,  retmode="xml")

#--------------------------------------------------------
#
#--------------------------------------------------------
# Actual eFetch call to the NCBI database occurs here. Deal with network
# errors in this function, returning -1 if call fails
#
def efetchGeneral(**kwargs):
    try:
        handle = Entrez.efetch(**kwargs)
    except urllib2.URLError, (err):
        print "Entrez.efetch Error: Could not open URL..."
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
