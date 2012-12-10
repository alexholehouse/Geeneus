import urllib2
from xml.dom.minidom import parseString
import ProteinObject
import Networking
import ProteinParser
import StringIO
from Bio import SeqIO
import re


class UniprotAPIException(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)    

class UniprotAPI:

    def __init__(self):   
        self.Network = Networking.Networking(40)
    
    def getProteinObjectFromUniProt(self, dataStore, accessionID):
            
        # networking and deal UniProtNetworkRequestwith XML
        dom = self.getDOMObject(accessionID)
        
        # first we check if the networking was succesful
        try:
            xml = dom.toxml()
        except AttributeError:
            print "Networking error - unable to retrieve " + accessionID + " from the UniProt server"
            return

        # If we get here then the networking portion theoretically worked
        # Next try and parse the XML for the various elements needed for the ProteinObject explcit
        # constructor
        try:
            self._isValidXML(dom)
            name = self._getProteinName(dom, accessionID)
            geneID = self._getProteinGeneID(dom, accessionID)
            geneName = self._getProteinGeneName(dom, accessionID)
            mutations = self._getProteinMutations(dom)
            sequence = self._getProteinSequence(dom, accessionID)
            creationDate = self._getCreationDate(dom, accessionID)
            other_accessions = self._getOtherAccessionValues(dom, accessionID)
            species = self._getSpecies(dom, accessionID)
            domains = self._getDomains(sequence, accessionID)
            taxonomy = self._getTaxonomy(dom, accessionID)
            isoforms = self._getIsoforms(dom, accessionID)
        except UniprotAPIException, e:
            print e
            print "Error while passing UniProt XML associated with accession + " + accessionID + ".\nSkipping..."
            
            # this means if we error at any point we don't add it to the datastore, but if all is well we do
            return 
        
        # so, if we got here were able to parse the XML ok, so jobs a good'un!
        dataStore[accessionID] = ProteinObject.ProteinObject(accessionID, xml, name, mutations, sequence, creationDate, geneID, geneName, other_accessions, species, domains, taxonomy, isoforms)

        
#--------------------------------------------------------
# PRIVATE FUNCTION
#--------------------------------------------------------
# 
    
    def getDOMObject(self, accessionID):

        handle = -1
        counter=0
        
        # crude 3-times retry. Maybe build a ProteinParser like
        # retry method, but given UniProt isn't a primary search
        # database seems like overkill at the moment
        while(handle == -1):
            handle = self.Network.UniProtNetworkRequest(accessionID)
            counter=counter+1
            if counter==3:
                return -1
        
        return parseString(handle.read())
    
#--------------------------------------------------------
# PRIVATE FUNCTION
#--------------------------------------------------------
# Internal function which gets the name of our protein from
# the DOM object
#

    def _getProteinName(self, domObject, ID):
        protein = domObject.getElementsByTagName('protein')
        self._neq1(protein, 'protein', ID)
        
        # Firstly we try and find names defined as rcommendedName->fullName
        # 
        for childNode in protein[0].childNodes:
            if childNode.nodeName == "recommendedName":
                fullName = childNode.getElementsByTagName('fullName')
            
                # check we only have a single fullName name
                self._neq1_explicit(fullName, "Multiple 'fullName' tags under the 'recommendedName' tag found in UniProt XML file associated with accession number " + ID)
                return(fullName[0].firstChild.nodeValue)        
          
        # if we get here we couldn't find a clear winner using recommendedName we try submittedName
        #
        for childNode in protein[0].childNodes:
            if childNode.nodeName == "submittedName":
                fullName = childNode.getElementsByTagName('fullName')
            
                self._neq1_explicit(fullName, "Multiple 'fullName' tags under the 'recommendedName' tag found in UniProt XML file associated with accession number " + ID)

                return(fullName[0].firstChild.nodeValue) 
        
        raise UniprotAPIException("Unable to find a protein name in the XML associated with accession " + ID)
        



#--------------------------------------------------------
# PRIVATE FUNCTION
#--------------------------------------------------------
# Internal function which gets the protien sequence
# 
# Behaviour Notes
# - We select the "protien sequence" by choosing the sequence
#   tag which has a 'length' attribute. As far as I can tell the only
#   sequence tags with 'length' attributes are those which represent the
#   protein sequence, BUT this does sort of seem like a bit of a fragile 
#   way of doing things... This is why the $found variable is included.
#   It will trip and error if we have more than one sequence tag with a 
#   length element. 
#
#
#
    def _getProteinSequence(self, domObject, ID):
        found = False
        sequencelist = domObject.getElementsByTagName('sequence')

        for sequence in sequencelist:
            if sequence.hasAttribute("length"):
                if found:
                    raise UniprotAPIException("Found multiple 'sequence' tags which have a 'length' attribute. Means we are unable to select a unique sequence. This is a flaw in the XML parser - please subit a bug report!" + ID)
                    
                retval = str(sequence.firstChild.toxml()).translate(None, '\n')
                found = True

        return retval
                
#--------------------------------------------------------
# PRIVATE FUNCTION
#--------------------------------------------------------
# Internal function which gets the a list of dictionaries of protien mutations
# congruent with the ProteinObject variants list format
#

    def _getProteinMutations(self, domObject):
        mutations = []

        featureList = domObject.getElementsByTagName('feature')

        for feature in featureList:
          
            if feature.attributes['type'].nodeValue == 'sequence variant':
                try:
                    tempDict = {'Mutant' : str(feature.getElementsByTagName('variation')[0].firstChild.toxml())}
                except IndexError:
                    # if we find mutations that lack a varation field, skip over it
                    continue
                        
                # for now we're only looking for single mutants, so if we find something else the reset,
                # discard and continue
                if len(tempDict['Mutant']) > 1:
                    del(tempDict)
                    continue
                
                tempDict['Original'] = str(feature.getElementsByTagName('original')[0].firstChild.toxml())

                # build as complete a description as possible, but some of these fields may be missing...
                try:
                    desc = str(feature.attributes['description'].nodeValue)    
                except KeyError:
                    desc = ""

                try:
                    idDesc = str(feature.attributes['id'].nodeValue)
                    idDesc = " id="+idDesc
                except KeyError:
                    idDesc = ""

                tempDict['Notes'] = desc + idDesc
                try:
                    tempDict['Location'] = str(feature.getElementsByTagName('location')[0].getElementsByTagName('position')[0].attributes["position"].nodeValue)
                except IndexError:
                    # Sometimes mutations have a begin and end tag instead of an internal position tag. If this is the case then the 
                    # position tag getElementByTagNames request comes up as an empty list triggering an IndexError.
                    # Given we're only looking at single mutations here we would ignore a region like this anyway. We just continue over and 
                    #
                    # Note that in a future version we might re-write this whole method to support many types of mutations, but at the moment we're
                    # just going for single mutations
                    continue
                    
                tempDict['Type'] = "Single"

                mutations.append(tempDict)
                

        return mutations
                
                
#--------------------------------------------------------
# PRIVATE FUNCTION
#--------------------------------------------------------
# Internal function which gets the name of the gene associated with this protein
#
# Behaviour notes
# - Stop searching for genes after the first primary gene is found (i.e. ignore later primary genes. 
#   May want to add warning if multiple primary genes are found, although this shouldn't happen.
#
# - If no primary genes are found we automatically select the first gene element ($fallBack) as the
#   gene of interest
#
# - This may not be known/recorded
##
    def _getProteinGeneName(self, domObject, ID):

        fallBack = ""
    
        geneList = domObject.getElementsByTagName('gene')
        
        # if there is no "gene" tag then set Gene name to an empty string
        try:
            self._neq1(geneList, 'gene', ID)
        except UniprotAPIException,e :
            print e
            return fallBack
        
        
        for subElement in geneList[0].getElementsByTagName('name'):
            
            # if fallBack is unset, set it (means fallBack is always set to the first real gene)
            if fallBack == "":
                fallBack = subElement.firstChild.toxml()
                
            # if we find the primary gene
            if subElement.attributes['type'].nodeValue == 'primary':
                return (subElement.firstChild.toxml())

        # if we get here no primary gene was found, so default to the first gene. Make sure we actually found 
        # a name somewhere in all that XML though!
        self._isEmpty(fallBack, "Unable to find a name tag associated with the UniProt XML for accession " + ID)
                
        return(fallBack)


#--------------------------------------------------------
# PRIVATE FUNCTION
#--------------------------------------------------------
# Internal function which gets the ID of the gene associated
# with this accession value
#
# Behaviour notes
# - Stop searching for genes after the first GeneID is found 
#   (there should only be one anyway)
#
# - If no IDs are found we just return an empty string
#
# - This may not be known/recorded
##

    def _getProteinGeneID(self, domObject, ID):
        elementList = domObject.getElementsByTagName("entry")[0].childNodes
        for element in elementList:
            if element.nodeType == element.ELEMENT_NODE and element.getAttribute("type") == "GeneID":
                return element.getAttribute("id") 
        return ""

    


#--------------------------------------------------------
# PRIVATE FUNCTION
#--------------------------------------------------------
# Internal function which gets the date the record was created
#

    def _getCreationDate(self, domObject, ID):
        retval = ""
        entry = domObject.getElementsByTagName('entry')
        self._neq1(entry, "entry", ID)
        
        retval = entry[0].attributes['created'].nodeValue
        
        # raise an exception if we didn't find a creation date
        self._isEmpty(retval, "Unable to find a creation date in the UniProt XML string associated with accession number " + ID)
        
        return retval

#--------------------------------------------------------
# PRIVATE FUNCTION
#--------------------------------------------------------
# Internal function which gets a list of alternative accession 
# values associated with this record
#

    def _getOtherAccessionValues(self, domObject, ID):
        
        accList = []

        elementList = domObject.getElementsByTagName("entry")[0].childNodes
        for element in elementList:
            
            # ________________________________________________________________
            # look for root XML based "assession" values
            if element.nodeName == 'accession':
                tempID = element.childNodes[0].nodeValue
                accList.append((tempID, ProteinParser.ID_type(tempID)[1]))

            # ________________________________________________________________
            # also look for a 'name' XML value, add as Swissprot as we do in
            # the equivalent function in ProteinObject.py 
            elif element.nodeName == 'name':
                accList.append((element.childNodes[0].nodeValue, "Swissprot"))

                
            # ________________________________________________________________
            # Look for EMBL protein sequences
            elif element.nodeType == element.ELEMENT_NODE and element.getAttribute("type") == "EMBL":
                
                # EMBL entires need to be further searched for their protien sequence ID values
                for subElement in element.childNodes:
                    if subElement.nodeType == subElement.ELEMENT_NODE and subElement.getAttribute('type') == "protein sequence ID":
                        tempID = subElement.getAttribute('value')
                        accList.append((tempID, ProteinParser.ID_type(tempID)[1]))

                        
            # ________________________________________________________________
            # Finally also look for refseq and IPI are a little accession valyes
            elif element.nodeType == element.ELEMENT_NODE and element.getAttribute("type") in ["IPI", "RefSeq"]:
                tempID = element.getAttribute('id')
                accList.append((tempID, ProteinParser.ID_type(tempID)[1]))

        return list(set(accList))

#--------------------------------------------------------
# PRIVATE FUNCTION
#--------------------------------------------------------
# Internal function which gets the species from which
# this protein is derived
#

    def _getSpecies(self, domObject, ID):
        elementList = domObject.getElementsByTagName("entry")[0].childNodes
        
        for element in elementList:
            if element.nodeName == "organism":
                for subElement in element.childNodes:
                    if subElement.nodeType == subElement.ELEMENT_NODE and subElement.getAttribute("type") == "scientific":
                        return (subElement.childNodes[0].nodeValue)

#--------------------------------------------------------
# PRIVATE FUNCTION
#--------------------------------------------------------
# Internal function which gets a list of dictionaries which
# hold the pfam domains. These domains are built by querying
# Pfam database directly because domain information is not
# stored on 
#

    def _getDomains(self, sequence, ID):

        #===========================================================
        # internal/local function - avoid
        # class-space polution
        def getPfamDOM(ID):
            handle = self.Network.PfamNetworkRequest(ID, False)
        
            try:
                domObject = parseString(handle.read())
            except Exception, e:
                print e 
                raise UniprotAPIException("Error when converting Pfam XML to DOM object")
        
            return domObject
        #===========================================================
        
        domainList = []
        tempDict = {}
        
        try:
            PfamDomObject = getPfamDOM(ID)
            self._isValidXML(PfamDomObject)
        except UniprotAPIException, e:
            print e
            raise UniprotAPIException("Error when accessing Pfam database with accession " + ID)

        elementList = PfamDomObject.getElementsByTagName("entry")[0].childNodes

        for element in elementList:
            # Ensure the Pfam sequence and our Uniprot sequence match, or the domain locations
            # will miss align!
            if element.nodeName == 'sequence':

                if not len(element.childNodes[0].nodeValue) == len(sequence):
                    print "Length difference between sequences"
                    print "Pfam sequence (" + str(len(element.childNodes[0].nodeValue)) +") = " + str(element.childNodes[0].nodeValue).lower()
                    print "Uniprot sequence ("  + str(len(sequence)) +") = " + str(sequence).lower()
                    
                    raise UniprotAPIException("Pfam and Uniprot sequences fail to match for accession " + ID)

            elif element.nodeName == 'matches':
                for match in element.childNodes:
                    if match.nodeName == 'match':
                        for location in match.childNodes:
                            if location.nodeName == 'location':
                                tempDict["start"] = str(location.getAttribute("start"))
                                tempDict["stop"] = str(location.getAttribute("end"))
                                tempDict["type"] = "pfam"
                                tempDict["accession"] = ID
                                tempDict["label"] = str(match.getAttribute("id"))
                                domainList.append(tempDict)
                                tempDict = {}
        return domainList

    def _getTaxonomy(self, domObject, ID):
        elementList = domObject.getElementsByTagName("entry")[0].childNodes

        taxon = []

        for element in elementList:
            if element.nodeName == "organism":
                for subElement in element.childNodes:
                    if subElement.nodeName == "lineage":
                        for taxonVal in subElement.childNodes:
                            if taxonVal.nodeName == "taxon":
                                taxon.append(taxonVal.childNodes[0].nodeValue)

        return taxon


    def _getIsoforms(self, domObject, ID):
        isoforms = {}
        isoformIDs = []
        isoformDomList = domObject.getElementsByTagName("isoform")
        for isoformDom in isoformDomList:
            for subElement in isoformDom.childNodes:
                if subElement.nodeName == "id":
                    isoformIDs.append(subElement.childNodes[0].nodeValue)

        if len(isoformIDs) > 0:
            isoformRaw = self.Network.UniProtBatchIsoformNetworkRequest(isoformIDs)

            if isoformRaw == -1:
                raise UniprotAPIException("Unable to carry out isoform sequence lookup through UniProt for accession " + ID)
            
            isoformList = list(SeqIO.parse(StringIO.StringIO(isoformRaw), 'fasta'))
            
            for record in isoformList:
                
                # first try and fine a number after the word isoform in the description
                # Originally I pulled the number after the dash as the isoform number (e.g. Q9NP78-5). However, it turns
                # out this number doesn't necessarily point to it's corresponding isoform, case in point, Q9NP78-5 actually
                # refers to isoform 4.
                #
                # If this extraction method fails we do go and pull the accession reference - it's right, but potentially
                # misleading!

                try:
                    boundary = str(record.description).find("Isoform")+8
                    
                    if boundary == 7:
                        # If we've actually pulled down a non-isoform sequence then
                        # let's just skip over it. Note boundary would == 7 because a failure
                        # to match = -1 and then we add 8!
                        continue
                    endBoundary = record.description[boundary:].find("of") 
                    
                    # it seems that literally every isoform is defined by Isoform <name> of ....
                    # so we just extract everything in the <name> portion to use as the name key
                    # we may have to add additional rules here...
                    # Weirdly, using .find gives a much cleaner and more robust solution
                    # than a regex...
                    #
                    isoformName = str(record.description[boundary:endBoundary+boundary-1])
                    
                    
                # this is a bit bad (catching every exception) but means if any part of the above
                # process goes wrong we default to exracting from the ID    
                except Exception, e:
                    print e 
                    print record.description
                    print "Error when parsing isoform name using description (above), falling back to identifiers..."
                    isoformName = (record.id[record.id.find("|")+1:record.id.rfind("|")])

                
                isoforms[isoformName] = record.seq.tostring().lower()

            return isoforms
                
        else:
            return {}
            
            
#--------------------------------------------------------
# PRIVATE FUNCTION
#--------------------------------------------------------
# ERROR HANDLING FUNCTIONS
    
    # does the  $nodelist have more than 1 $element tag?    
    def _neq1(self, nodelist, element, ID ):
        if nodelist.length > 1:
            raise UniprotAPIException("Multiple '" + element + "' tags found in UniProt XML string associated with accession number " + ID ) 
        if nodelist.length < 1:
            raise UniprotAPIException("No '" + element + "' tags found in UniProt XML string associated with accession number " + ID ) 

    # does the $nodelist have more than 1 $element tag - show explicit message
    def _neq1_explicit(self, nodelist, message):
        if nodelist.length != 1:
            raise UniprotAPIException(message) 

    # is a string empty - show explicit message
    def _isEmpty(self, string, message):
        if string == "":
            raise UniprotAPIException(message)

    def _isValidXML(self, domObject):
        if not len(domObject.getElementsByTagName("entry")) == 1:
            raise UniprotAPIException("Invalid XML")






        

        


            
        
            
            

        



        


