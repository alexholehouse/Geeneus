# Base class which deals with caching and storing of data objects, as well as coordinating network access (private)
#
# Copyright 2012-2015 by Alex Holehouse - see LICENSE for more info
# Contact at alex.holehouse@wustl.edu

from Bio import Entrez
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

import httplib

import Networking


## This exception is imported by other classes, and is used
## when they are creating new objects and discover, to their
## horror, that the XML is semantically wrong (even if it's
## syntactically good) - e.g. if it's missing some expected
## info
class BadXMLException(Exception):
    pass


class GeneralRequestParser:
    
    def __init__(self, email, cache, retry=0, loud=True):
        try:
            Entrez.email = email
            self.Networking = Networking.Networking(30)
            self.loud = loud
            self.retry = retry
            self.cache = cache
            self.error_status = False
        
        except: 
            print "Error building generic parser object"
            self.error_status = True
    

            ###
    def error(self):
        return self.error_status
    
    
    #--------------------------------------------------------
    # Preconditions
    # ID must be a valid ID, so probably want to run through a preprocessing function before you
    # call this
    #
    # - datastore is the datastore passed from the child class, which must be an ID indexed dictionary
    #   where the value is a specific storage object type (e.g. ProteinObject, GeneObject etc)
    #
    # - fetchFunction is the Networking function which we call to
    #
    # - newObjectConstructor is a constructor for the object type we store in datastore
    #
    # - alternative is the "alternative function" which we query if NCBI services don't seem to work
    # 
    def _get_object(self, ID, datastore, fetchFunction, newObjectConstructor, alternative=False):
        
        altSuccess = False  # set flag for later
      
        if ID == -1 or ID not in datastore or not self.cache:  
           
            if ID  == -1:
                
                # by returning the object associated with [] we don't pollute the datastore with invalid and pointless
                # searches, we avoid queriying NCBI without a hope in hell of a hit, and we take advantage of the built in
                # bad XML behaviour without raising an error, because, technically, no error has happened, we just know
                # that the ID in question won't return valid XML. It's not an error - it's just stupid.
                
                return newObjectConstructor(-1, [])
                       
            xml = -1
            retry = self._build_retry_function(fetchFunction);
                     
            # retry uses a closure based approach, so this just happens a set number of times
            # - not an infinite loop if we can't get a response!
            while (xml == -1):
                xml = retry(ID)

            # if we still can't get through after retrying a number of times
            if (xml == -2):

                # If we have an alternative function to try let's try it!
                # Alternative functions have to be totally defined and initialized
                # in the calling child class. The function can, obviously, do whatever
                # you want, but the behaviour must generally;
                #
                # 1) Take just an accession ID (a datastore object is internally stored
                #    as a object variable from the child class)
                #
                # 2) Deal with everything else internally
                #
                # 3) Simply return a "True" if it was succesful, or "False" if not
                #
                # If alternative returns success it *MUST* have been able to query
                # the datastore with ID and return a relevant object.
                #
                
                if alternative:
                    altSuccess = alternative(ID)

                if not altSuccess:
                    datastore[ID] = newObjectConstructor(-1, [])
                
            # or NCBI query was succesfull
            else:
                # if the xml has been set to an object, this object MUST be
                # the already constructed *Object (e.g. ProteinObject, GeneObject etc)
                # allows additonal control in the fetch function
                if str(type(xml)) == "<type 'instance'>":
                    datastore[ID] = xml

                # else parse that XML
                else:

                    # this try catches cases where the XML, while formatted correctly,
                    # actually doesn't make sense or is lacking in an expected field
                    try:
                        datastore[ID] = newObjectConstructor(ID, xml)
                    except BadXMLException:
                        datastore[ID] = newObjectConstructor(-1, [])
                
        return datastore[ID]

    
#--------------------------------------------------------
# PRIVATE FUNCTION
#-------------------------------------------------------
# This function provides an interface with the Networking
# functionality to return a list of XML elements, each of
# which corresponds to the ID in the listofIDs
#
#
    def _get_batch_XML(self, listOfIDs, function_to_apply, alternative=False):
        
        # base case 1 (base case 2 if len == 1)
        if len(listOfIDs) == 0:
            return []

        xml = -1
        retry = self._build_retry_function(function_to_apply);

        while (xml == -1):
            xml = retry(listOfIDs)

        # If we fail implement recursive batch fetch, such that we split the list in half
        # and batch fetch each half. Do this after troubleshooting!
        if xml == -2 or not len(xml) == len(listOfIDs):

            # we only call alternative if we're down to a single accession - ie. really
            # try and get it from NCBI
            if len(listOfIDs) == 1:
                if alternative:
                    alternative(listOfIDs[0])
                # note we still return -1 so as to allow the algorithm to be maximally efficient
                return [-1]
            
            # split list in half and recursivly batch both halves
            # this serves two purposes - it lets us retry a number of times proportional
            # to the length of the list, and it provides a binary search mechanism to 
            # root out a potential rotten apple which may be causing the list to error

            xml = self._get_batch_XML(listOfIDs[:int(len(listOfIDs)/2)], function_to_apply, alternative)
            xml_2 = self._get_batch_XML(listOfIDs[int(len(listOfIDs)/2):], function_to_apply, alternative)

            xml.extend(xml_2)
            
        return xml


#--------------------------------------------------------
# PRIVATE FUNCTION
#--------------------------------------------------------          
# Builds a closure based retry function, which when called
# the first self.retry times will atempt to get the parsed
# xml for the ID in question. However, on the
# self.retry+1 time it will simply return -2 
#
#
# Return Values
#     -1 on failure
#     An initialized GeeneusObject
#     Entrez XML structure
#
    def _build_retry_function(self, function_to_apply):

        retryCounter = [0]
        numberOfRetries = self.retry

        def retry(ID):

            if retryCounter[0] < numberOfRetries+1:
                
                retryCounter[0] = retryCounter[0]+1
                
                ## if we're not on our first try
                if not retryCounter[0]-1 == 0:
                    print("Retry number " + str(retryCounter[0]) + " of " + str(numberOfRetries+1))
                
                ## return value may be a real handle or -1
                handle = function_to_apply(ID)
                
                ## if we failed return -1
                if handle == -1:
                    return -1

                # function to apply can be made to return either XML
                # or a fully formed GeeneusObject (e.g. shortcutting
                # for ProteinObject).
                try:
                    if handle.type == 'GeeneusObject':
                        return handle
                except AttributeError:
                    pass # carry on...
                
                try:
                    XML = Entrez.read(handle)
                
                # what are these errors?
                # httplib.IncompleteRead - we accidentally closed the session early, either because of a timeout on the client end (i.e. here) or because 
                #                          of some kind of server error
                # 
                # Bio.Entrez.Parser.CorruptedXMLError - Something is wrong with the XML 
                # Bio.Entrez.Parser.NotXMLError - the XML is not XML (unlikely, but worth keeping!)
                # Bio.Entrez.Parser.ValidationError - unable to validate the XML (this can be ignored, but best not to!)
                except (httplib.IncompleteRead, Entrez.Parser.CorruptedXMLError, Entrez.Parser.NotXMLError, Entrez.Parser.ValidationError), err:  
                    return -1

                return XML

            else:
                return -2;

        return retry        
