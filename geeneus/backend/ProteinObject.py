# Provides an object-based API for a protein (private)
# Part of the Geeneus package
#
# Copyright 2012 by Alex Holehouse, with significant contributions
# from Matt Matlock - see LICENSE for more info
# Contact at alex.holehouse@wustl.edu

from Bio import Entrez, Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC

import string
import re


import ProteinParser
from Utilities import show_warning, show_error, show_status

# ==========================================================================================
# Object attributes
# 
# ------------------------------------------------------------------------------------------
# self.exists                 Does the object associated with this ID exist in the database
# self.error                  ALWAYS false in a ProteinObject (true in ProteinErrorObject)
# self.raw_XML                Unprocessed XML string
# -------------------------------------------------------------------------------------------
# self.name                   Protein name
# self.database               String defining which database the record was obtained from ("NCBI" 
#                             or "UniProt")
# self.version                Record version (if unversioned defaults to 1)
# self.sequence               Protein amino acid sequence
# self.sequence_create_date   Date the sequence was entered into the protein database
# self.sequence_length        Number of amino acids residues
# self.other_accessions       List of other accesison values
# self.geneID                 NCBI GeneID for the protein, should you want to lookup Gene information
# self.gene_name              Gene name (may be an empty string)
# self.species                Species of origin
# self.taxonomy               Ordered taxonomy string list
# self.host                   Name of viral host organism ('N/A' if not applicable)
# self.domain_list            List of pfam defined domains
# self.isoforms               List of dictionaries, each dictionary corresponding to an isoform.
#                             Dictionaries are keyed by the isoform ID (e.g. Q12345-4) and each 
#                             value is a 2 element list. Element 1 is the isoform name, and 
#                             element 2 is the sequence.
# self.protein_variants       List of dictionaries, each dictionary corresponding to a unique
#                             variant. Each dictionary has location, mutation and notes





######################################################### 
######################################################### 
# Exception class for ProteinObject errors
#
class ProteinObjectException(BaseException):
    pass


######################################################### 
######################################################### 
# Main class for dealing with Protein objects. Can construct
# them from an explicit initializer, or can build them
# by parsing NCBI derived XML
#
class ProteinObject:


#--------------------------------------------------------
# PUBLIC GETTER FUNCTIONS
#--------------------------------------------------------
# Getter functions should be used rather than direct member access - try and
# maintain some encapsulation.
#

    def exists(self):
        return self._exists

    def error(self):
        return self._error

    def get_raw_xml(self):
        return self.raw_XML
    
    def get_protein_name(self):
        return self.name

    def source(self):
        return self.database
    
    def get_version(self):
        return self.version

    def get_protein_sequence(self):
        return self.sequence

    def get_creation_date(self):
        return self.sequence_create_date

    def get_protein_sequence_length(self):
        return self.sequence_length
    
    def get_other_accessions(self):
        return self.other_accessions

    def get_geneID(self):
        return self.geneID

    def get_gene_name(self):
        return self.gene_name

    def get_species(self):
        return self.species

    def get_taxonomy(self):
        return self.taxonomy

    def get_host(self):
        return self.host

    def get_domains(self):
        return self.domains

    def get_isoforms(self):
        return self.isoforms

    def get_variants(self):
        return self.protein_variants

#--------------------------------------------------------
# PUBLIC FUNCTION
#--------------------------------------------------------
# Object initializer
# Multiple dispatch because no function overloading is allowed
# in Python. We need to be able to create ProteinObjects from
# both NCBI XML, but also from pre-extracted values from
# other sources (e.g. UniProt XML)
#
    def __init__(self, *args):

        # metavalue
        self.type = "GeeneusObject"

        # if we have a single xml string
        if len(args) == 2:
            self.__init_1(*args)

        # if we're passing in pre-parsed data
        else:
            self.__init_2(*args)


#--------------------------------------------------------
# PRIVATE FUNCTION
#--------------------------------------------------------
# Initializes all the objects attributes to default values before
# populating with proteinxml based data. If no xml data is available
# default values are not overwritten, so rather than exceptions being
# raised on request for non-existant information default values are, and
# exists remains set to False
#
    def __init_1(self, accession, proteinxml):
       
        # set the default values (these are kept for empty/
        # error calls
        
        
        self.accession = accession
        self.sequence = None
        self._exists = False
        self._error = False
        self.sequence_create_date= None
        self.protein_variants = None
        self.geneID = None
        self.sequence_length = None
        self.name = None
        self.other_accessions = None
        self.species = None
        self.taxonomy = None
        self.host = None
        self.domains = None
        self.gene_name = None
        self.isoforms = None
        self.raw_XML = None
        self.database = None
        self.version = None
        
        if proteinxml == -1:        
            self._error = True
            return

        if not self._xml_is_OK(proteinxml):

            # for intrest, lets set the raw XML value if possible
            # this is useful because sometimes the XML is valid, but there
            # is some aspect which makes it invalid for our purposes (e.g.
            # its for mRNA!)
            try:
                self.raw_XML = proteinxml[0]
            except IndexError, e:
                return
            return

        self._exists = True

        # Now we set the rest of the values using the parsed XML
        # Note this is kept in a KeyError try/except block 
        # to catch any malformed XML edge cases.
        #

        try:
            self.database = 'NCBI'
            self.raw_XML = proteinxml[0]
            self.name = self._extract_protein_name(proteinxml[0]["GBSeq_definition"])
            self.version = self._extract_version(proteinxml[0]) 
            self.sequence = proteinxml[0]["GBSeq_sequence"].lower()
            self.sequence_create_date = proteinxml[0]["GBSeq_create-date"]
            self.sequence_length = len(self.sequence)
            self.geneID = self._extract_geneID(proteinxml[0]["GBSeq_source-db"], proteinxml[0]["GBSeq_feature-table"])
            self.gene_name = self._extract_gene_name(proteinxml[0]["GBSeq_feature-table"])
            
            self.other_accessions = self._extract_other_accessions(proteinxml[0])
            self.species = self._extract_species(proteinxml[0])
            self.taxonomy = self._extract_taxonomy_string(proteinxml[0]['GBSeq_taxonomy'])
            self.host = self._extract_host(proteinxml[0]["GBSeq_feature-table"])


            self.domains = self._extract_domain_list(proteinxml[0]["GBSeq_feature-table"])
            
            self.isoforms = self._extract_isoforms(proteinxml[0], accession, self.sequence)

            self.protein_variants = self._extract_variant_features(proteinxml[0]["GBSeq_feature-table"])        
        except KeyError, e:
            print "ERROR when building ProteinObject using accession " + accession + " (NCBI XML)"
            print e
            raise e
            

#--------------------------------------------------------
# PRIVATE FUNCTION
#--------------------------------------------------------
# Autocreate initializer
# If you already have the relevant data you can build a ProteinObject directly.
# Useful for non NCBI based record construction (e.g. UniProt)
#
    def __init_2(self, accession, version, xml, name, mutations, sequence, creation_date, geneID, gene_name, other_accessions, species, domains, taxonomy, isoforms, database, host):
        
        self.accession = accession
        self._exists = True
        self._error = False
        self.version = version
        self.sequence = sequence
        self.sequence_create_date = creation_date
        self.protein_variants = mutations
        self.geneID = geneID
        self.name = name
        self.raw_XML = xml
        self.sequence_length = len(self.sequence)        
        self.other_accessions = other_accessions
        self.species = species
        self.taxonomy = taxonomy
        self.domains = domains
        self.gene_name = gene_name
        self.isoforms = isoforms
        self.database = database # presumably always "UniProt"...
        self.host = host 


#--------------------------------------------------------
# PRIVATE FUNCTION
#--------------------------------------------------------
# Function to check that the xml we've downloaded is good, and represents a 
# viable protein xml structure. Additional tests may be added here as we find more
# edge cases! Returns FALSE if there's a problem, TRUE otherwise
#

    def _xml_is_OK(self, proteinxml):
        if len(proteinxml) > 1:
            show_warning(" [ProteinObject._xml_is_ok()] - ProteinXML detected more than one record associated with this GI.\nThis should never happen.")
            return False
        
        # Nothing in XML - so return an empty-initiailized object with exists = 0
        if len(proteinxml) == 0:
            return False

        # Check that we're really dealing with protein (despite specifying db="protein"
        # on the efetch call, when a GI is used other databases seem to be searched too)
        try:
            if not (proteinxml[0]["GBSeq_moltype"] == "AA"):
                return False
        except TypeError:
            return False

        return True


#--------------------------------------------------------   
#
#--------------------------------------------------------
# Function to get variant features. Easy to extended should extra
# variant data be needed, but for the function creates a list of n
# dictionaries (where n = number of variants in protein xml data)
# and each dictonary contains variant location, mutation and notes.
#
    def _extract_variant_features(self, featurelist):

        #===========================================================
        # Function to actually build the mutation. At the moment only
        # supports single mutations although the plan is to expand that
        # out in the future
        #
        def buildMutationEntry(feature, acc, sequence):
            
            # get location of variant feature
            location = -1
            defString = -1
            mutation= {}
            
            location = feature["GBFeature_location"]

            
            
            # get mutation defString
            for feature_subsection in feature["GBFeature_quals"]:
                if not feature_subsection.has_key("GBQualifier_value"):
                    continue
                
                # somewhat niave approach - take first one we find...
                if feature_subsection["GBQualifier_value"].find("->") > -1 or feature_subsection["GBQualifier_value"].find("Missing") > -1 :
                    defString = feature_subsection["GBQualifier_value"]
                    break
                
                
            # if we found a location but failed to build a new defString raise an exception
            # as this suggests a failing in the way the defString is parsed
            if not location == -1 and defString == -1:
                raise ProteinObjectException("Flaw in how we parse the Variant defining string - no missing or -> found but location was defined")
                
            
            # first we parse the location information
            
            try:
                location = int(location)
            except ValueError:
                try:
                    location = int(location.split("..")[0])
                except Exception, e:
                    print "Fundemental flaw while parsing location information for " + acc
                    raise e
                
                
            # now we parse the actual variant 
            if defString.find("Missing") > -1:
                mutation["variant"] = ""
                mutation["original"] = sequence[location-1].upper()
                mutation["mutant"] = "-"
                mutation["type"] = "Deletion"
                mutation["notes"] = defString[7:] # So defString is "Missing (..." so we just get (...
                mutation["location"] = location
                return mutation
                
            # note that at the moment it seems like we can either have "Missing" residues or replacements
            # ( -> ) 
            else:
                mutation["location"] = location
                mutation["variant"] = defString[:defString.find("(")-1]
                mutation["original"] = defString.split(" -> ")[0]
                defString = defString.split(" -> ")[1]

                # the mutation is the string of A-Z or spaces after the " -> " symbol. This is 
                # because we can end this string with a variety of characters
                mutation["mutant"] = re.search("[A-Z ]*", defString).group()

                # we define the notes as the rest of the defString
                mutation["notes"] = defString[len(mutation["mutant"]):]
                
                # finally, strip spaces from the mutation and original strings
                mutation["original"] = mutation["original"].replace(" ","")
                mutation["mutant"] = mutation["mutant"].replace(" ","")
                
                ## Now we've built the mutations we define the type based on the
                ## change being made
                
                # insertion
                if len(mutation['mutant']) > len(mutation['original']):
                    if mutation['original'] in mutation['mutant']:
                        mutation["type"] = "Insertion"
                    else:
                        mutation["type"] = "Insertion & substitution"
                    return mutation

                # deletion/subsitution
                if len(mutation['mutant']) < len(mutation['original']):
                    mutation["type"] = "Deletion & substitution"
                    
                        
                    return mutation
                
                # substitution
                if len(mutation['mutant']) == len(mutation['original']):
                    
                    if len(mutation['mutant']) == 1:
                        mutation["type"] = "Substitution (single)"
                        return mutation
                    
                    if len(mutation['mutant']) == 2:
                        mutation["type"] = "Substitution (double)"
                        return mutation
                    
                    else:
                        mutation["type"] = "Substitution (" + str(len(mutation['mutant'])) + ")"
                        return mutation
                    
                # if we get here then need to raise and exception
                raise ProteinObjectException("Couldn't identify the mutation type for accession " + acc)
            
        #===========================================================

        variant_list = []
        
        for feature in featurelist:            
            if not feature.has_key("GBFeature_quals"):
                continue            
            for feature_subsection in feature["GBFeature_quals"]:               
                if not feature_subsection.has_key("GBQualifier_value"):
                    continue
                if feature_subsection["GBQualifier_value"] == "Variant":
                    mutDictTemp = buildMutationEntry(feature, self.accession, self.sequence)
                    
                    # if mutDictTemp = False we didn't find an appropriate
                    # mutation match. Else add it.
                    if mutDictTemp:
                        variant_list.append(mutDictTemp)
        
        if len(variant_list) == 0:
            return []
        else:
            return variant_list
        
        
#--------------------------------------------------------
# PRIVATE FUNCTION
#--------------------------------------------------------
# Function to extract the geneID from the protein data for 
# use in getting gene information from the Genome class if
# needed
#
    def _extract_geneID(self, GBSeq_source_db, featuretable):
        source = str(GBSeq_source_db) + " " + str(featuretable)
    
        # See if there's a geneID in the DB source data, and if not
        # return an empty string. NB some proteins won't have a GeneID
        geneID_location = source.find("GeneID:")
        if (geneID_location == -1):
            return ""
        
        # Given we found a "GeneID:" label in the text, we cut out the
        # value subsequent to the tag before the next comma, and return
        # that as the GeneID
        geneID_location = geneID_location+7
        geneID_end_location = re.search("\W", source[geneID_location:]).regs[0][0] + geneID_location

        geneID = source[geneID_location:geneID_end_location]
        
        return geneID


    def _extract_protein_name(self, nameString):
        
        start = nameString.find("=")
        if not start < 0:
            return nameString[start+1:]
        else:
            return nameString
        

#--------------------------------------------------------
# PRIVATE FUNCTION
#--------------------------------------------------------
# Extract the gene name from the source XML
#
    def _extract_gene_name(self, feature_table):

        gene_name = ""
        
        for f in feature_table:

            # look for explicit gene fields - if we find one return
            # automatically
            if  f["GBFeature_key"]  ==  'gene':
                gene_name  =  self._get_qualifier (["Gene","gene"],  f['GBFeature_quals'])
                
                # avoids accidentally returning None if nothing is found
                if gene_name:
                    return gene_name

            # look for coding sequence names - alternative to gene
            # name but may be overridden if an explicit gene name also exists
            if  f [ 'GBFeature_key' ]  ==  'CDS' :
                gene_name = self._get_qualifier( 'gene', f['GBFeature_quals'])
                
        # if gene_name is set to None then return an empty string instead
        if not gene_name:
            gene_name = ""
        
        return gene_name
                

#--------------------------------------------------------
# PRIVATE FUNCTION
#--------------------------------------------------------
# Build a list of dictionaries containing domains
#

    def _extract_domain_list(self, feature_table):

        domainList = []

        for f in feature_table:
            tempDomainDictionary = {}
            
            # check first
            if not f.has_key('GBFeature_quals') or not f.has_key('GBFeature_key'):
                continue

            try:
                note_val = self._get_qualifier('note', f['GBFeature_quals'])
            except KeyError:
                continue
            if f['GBFeature_key'] == 'Region' and note_val != None and note_val.find('pfam') > 0:
                
                tempDomainDictionary["type"] = "pfam"
                tempDomainDictionary["label"] = self._get_qualifier('region_name', f['GBFeature_quals'])
                tempDomainDictionary["accession"] = f['GBFeature_intervals'][0]['GBInterval_accession']
                tempDomainDictionary["start"] = int(f['GBFeature_intervals'][0]['GBInterval_from'])
                tempDomainDictionary["stop"] = int(f['GBFeature_intervals'][0]['GBInterval_to'])
                
                # add the domainDictionary element
                domainList.append(tempDomainDictionary)

        return domainList


    def _extract_version(self, xml):

        # if no version is provided assume this is version 1 
        # (for older records the initial version was not set to 1)
        try:
            version = xml['GBSeq_accession-version'].split(".")[1]
        except (KeyError, IndexError):
            version = "1"

        return version

        

#--------------------------------------------------------
#
#--------------------------------------------------------
# Internal function to get the isoform sequences associated
# with a specific protein. To minimize the number of networking
# calls, it does this by parsing the splicing variant annotation
# information and reconstructing a new protein sequence based on
# that information.
#
# 
#
#
    def _extract_isoforms(self, xml, ID, sequence):
        

        #===========================================================
        # Local exception, because isoforms are serius buznuz
        #
        class IsoformException(BaseException):
            pass
        
        #===========================================================
        #
        # Function which takes the "note" string from a splicing variant
        # record and determines the isoform names to which that splicing
        # variant applies.
        # 
        # The splicing variant "note" structure is typically something like
        #
        #    "<isoform change>' (in isoform 1 [and isoform 2 and isoform 3....]"
        #
        # Where 1, 2 and 3 represent the relevant isoform names (i.e. this
        # could be) 
        #
        #    "1..34 missing (in isoform small and isoform short and isoform pointless"
        #
        # This method pulls out those names, returning a list of them. 
        # 
           
        def getRelevantIsoforms(defString):
            isoforms = []   

            # locations is now a list of everywhere 'isoform' is found
            # in the defString
            locations = [m.start() for m in re.finditer('isoform', defString)]
            
            for i in xrange(0,len(locations)):
                
                # if we're at the last isoform in the list (i.e. if we have 1)
                # then 0 = 1-1
                if i == len(locations)-1 :
                    end = defString.rfind(")")
                    name = defString[locations[i]+8:end]

                # we're inside the isoform
                else:
                    
                    end1 = defString[locations[i]+8:].find('and isoform')
                    end2 = defString[locations[i]+8:].find('isoform')
                    
                    # if the next end is also the "and isoform" end 
                    if end1+4 == end2:
                        end = end1 - 1 #-1 to cut off space
                    else:
                        end = end2 - 2 #-2 to cut off comma and space
                        
                    name = defString[locations[i]+8:locations[i]+8+end]
                    
                isoforms.append(name)
            return isoforms
            
        #
        # Function which takes the 'note' string (as defined above) along with the 
        # location of that splicing variant and encodes it into a list with the 
        # following structure. 
        # 
        # [TYPE OF SPLICE VARIANT, START, END, ORIGINAL VALUE, NEW VALUE]
        # Note that
        # - type is one of "missing" or "replacement"
        # - if "replacement", original and new values represent the switcharoo
        # - if "missing" then original and new values are empty, we just cut out the region 
        #
        # The location parameter parameter can be a single point, or a range, as defined
        # in the XML

        def getSpliceEvent(defString, location):
                        
            defString = defString.lower()

            # grab location information
            try:
                locationList = [int(location["GBInterval_from"])-1, int(location["GBInterval_to"])]
            except KeyError:
                # some regions are defined as  a single point
                locationList = [int(location["GBInterval_point"])-1, int(location["GBInterval_point"])] 
            
            # if this defString defines a missing event   
            if defString.find("missing") > -1:
                return["missing", locationList[0], locationList[1], "", ""]
            
            # if the defstring defines a replacement event
            if defString.find("->") > -1:

                # to get the string of AA either side of -> we have to do some
                # quick regexs/reformatting
                defString = defString.replace(" ","")
                
                old = str(defString.split("->")[0])
                new = re.search("[a-zA-Z]*", defString.split("->")[1]).group()

                return ["replacement", locationList[0], locationList[1], old, new]

            ########################################
            # if we get here this XML is, technically, maformed. HOWEVER, we can try and salvage by adding 
            # 'exception' handling

            show_warning("The isform line \n\n[" + defString + "] \nis malformed\n\nTrying to parse regardless...\n")

            # Exception 1 use of ">" instead of "->"
            if defString.find(">") > -1:

                

                # to get the string of AA either side of -> we have to do some
                # quick regexs/reformatting
                defString = defString.replace(" ","")
                
                old = str(defString.split(">")[0])
                new = re.search("[a-zA-Z]*", defString.split(">")[1]).group()

                return ["replacement", locationList[0], locationList[1], old, new]

            
            raise ProteinObjectException("Isoform 'note string' had neither 'missing' nor '->' in it, and\neven specific known exception handling couldn't help!\nParse error!") 
        
        #
        # Function to ensure the semantics of the defined splicing
        # variants make sense.
        # It's worth noting that as of writing (and maybe 100 000 accessions tested) no 
        # set has ever not met these constraints, which gives us confidence that the
        # algorithm approach is valid. It does, however, rely on the quality of the NCBI
        # annotations, so considering this a seperate constraints test method seemed appropriate
        # rather than risking silent errors.
        #
        # The general idea is that no two splicing regions should overlap.
        # We only pass the ID and isoform parameters to create more expressive
        # error messages.
        #
        # The constraintsList is progressivly built. For each isoform in turn, we 
        # cycle through every splicing event, and see if that isoform was involved in
        # the splicing event. If it was, we check that the region defined by that splicing
        # event does not overlap with any of the regions already defined in the
        # constraintsList, and if it does not, we then add this region to the constraintsList
        # This means that for each isoform we can check that none of the regions overlap
        #
        
        
        
        def checkAndBuildConstraints(constraintsList, eventDetails, ID, isoform):
            
            returnConst = constraintsList

            for constraint in constraintsList:
                # if lower value of our new constraint is smaller than
                # the upper value of another constraint these regions must overlap
               
                if eventDetails[1] < constraint[2] and eventDetails[1] > constraints[1]:
                    raise IsoformException("Isoform definition overlap (type 1) on accession "+ID+" between isoform " + isoform +" and "+constraint[0])
                
                # if the upper bound on our new constrain is larger than
                # the lower bound on an existing constaint those must also 
                # overlap
                if eventDetails[2] > constraint[1] and eventDetails[2] < constraint[2]:
                    raise IsoformException("Isoform definition overlap (type 2) on accession "+ID+" between isoform " + isoform +" and "+constraint[0])
                    
            returnConst.append([isoform, eventDetails[1], eventDetails[2]])
            
            return returnConst

        
        # BuildIsoIDDictionary is the first thing done, and involves taking the
        # the raw xml and building a dictionary which maps the isoform ID (e.g. Q12345-2)
        # to the isoform's name 
        #
        def buildIsoIDDictionary(xml):
            
            isoID = {}
            try:
                comments = xml["GBSeq_comment"]
            except KeyError:
                return {}
            start = comments.find("[ALTERNATIVE PRODUCTS]")
            if start == -1:
                return {}
            
            stop = comments[start+23:].find("; [")
            
            # safety check for edge cases
            if stop == -1:
                stop = comments[start+23:].find(";[")
                if stop == -1:
                    stop = comments[start+23:].find("[")
                    if stop == -1:
                        stop = len(comments)-(start+23)
            
            # the apString is the string which contains all the 
            # altenative product isoformID<->name mappings
            apString = comments[start+23:start+23+stop]
            
            # first find the first "Name=" delimiter
            startOfName = apString.find("Name=")+5
            
            # we now just cycle through the apString, redefining it on every cycle
            # as the original string minus the region with the preceding name-ID mapping
            while startOfName > 4:
                              
                # get name
                endOfName = apString[startOfName:].find(";")+startOfName
                name = apString[startOfName:endOfName]
                
                # update the apString
                apString = apString[endOfName:]
                
                # get isoform ID
                startOfID = apString.find("IsoId=")+6
                endOfID = apString[startOfID:].find(";")+startOfID
                ID = apString[startOfID:endOfID]
                
                # add into dictionary
                isoID[name] = ID
                
                # update the apString again!
                apString = apString[endOfID:]
                
                # finally, reset the startOfName again and repeat until
                # we can't find any more names
                startOfName = apString.find("Name=")+5
            return isoID
 
        
        #===========================================================
        
        # MAIN FUNCTION BEGINS HERE
                
        # pull out the feature table
        ft = xml["GBSeq_feature-table"]
        
        # build the IsoID - name conversion dictionary
        nametoIsoID = buildIsoIDDictionary(xml)

        isoformList = []
        splicingEvents = []

        for i in ft:
            if i["GBFeature_key"] == "Region":
                if self._get_qualifier("region_name", i["GBFeature_quals"]) == "Splicing variant":

                    
                    defString = self._get_qualifier("note", i["GBFeature_quals"])

                    try:
                        splicingEvents.append((getRelevantIsoforms(defString), getSpliceEvent(defString, i["GBFeature_intervals"][0])))
                    except IsoformException, e:
                        print e
                        raise IsoformException("Error while getting isoform data for accession " + ID) 

        # now we've built a list of tuples of the form ([isoform numbers], [description]) we have to do each isoform in sequence
        # Our isoform formation assumes that each isoform description talks about refrence in relation to reference sequence (isoform 1)
        # e.g. if isform 2 has G->A at 299 and has 100->200 missing you would not remove 100-200 and the change the *newly indexed* 299 
        # which may not even exist, but instead would change 299 in the original, or 199 in the new one.
        # 
        # This has semantic difficulties. Notably, we have to keep track of the offsets generated by our changes so as to change
        # the correct locations in the future.

        # First identify the number of different isoforms
        for event in splicingEvents:
            isoformList.extend(event[0])

        isoformList = list(set(isoformList))
        isoformList.sort()
        isoformSequenceList = {}

        for isoform in isoformList:
            
            # offsetVector provides a mapping between the basic sequence indices
            # and the progressive splice variant indices
            offsetVector = [0]*len(sequence)
            seqLen = len(offsetVector)

            # constraints is used to ensure we have non-overlapping splicing regions
            constraints = []

            # initialize the isoform sequence to the primary sequence at first 
            isoformSequenceList[isoform] = sequence
                   
            # first loop checks that the splicing instructions actually make sense
            # eg [A->G at 235] and [missing 100-300] doesnt make sense, so would raise an exeption
            for event in splicingEvents:
                # if our isoform is affected by the event
                if isoform in event[0]:
                    constraints = checkAndBuildConstraints(constraints, event[1],ID, isoform)
                              
            # having checked our constraints are not broken, we can construct the isoform sequence
            for event in splicingEvents:
                if isoform in event[0]:

                    if event[1][0] == "missing":
                        
                        # note we -1 because the position indicies
                        # start on 1, not 0

                        start = event[1][1]
                        stop = event[1][2]
                        deltaOffset = (start-stop)
                        
                        # we enclose this in a try/except block because if we're off on
                        # or indexing offsetVector[start] might be a None
                        try:
                            sect1 = isoformSequenceList[isoform][:offsetVector[start]+start]
                            sect2 = isoformSequenceList[isoform][offsetVector[stop-1]+stop:]
                                                    
                        except Exception, e:
                            print e
                            raise IsoformException("ERROR when cutting out splice variant regions [missing]")

                        isoformSequenceList[isoform] = sect1+sect2

                        # set values in removed region to None, such that should we ever mess up
                        # we avoid any silent creeping errors
                        for i in xrange(start, stop):
                            offsetVector[i] = None

                        # for all values greater than stop we now update our offsetVector
                        # to include the correct offset 
                        for i in xrange(stop, seqLen):
                            offsetVector[i] = offsetVector[i] + deltaOffset
            
                    if event[1][0] == "replacement":
                        start = event[1][1]
                        stop = event[1][2] 

                        # i.e. if we replace a 5-AA region with a 1 AA region
                        # then we get 1-5 = -4 deltaOffset
                        deltaOffset = len(event[1][4]) - len(event[1][3])
                        
                        try:
                            sect1 = isoformSequenceList[isoform][:offsetVector[start]+start]
                            sect2 = isoformSequenceList[isoform][offsetVector[stop-1]+stop:]
                        except Exception, e:
                            print e
                            raise IsoformException("ERROR when cutting out splice variant regions [replacement]")   
                        
                        isoformSequenceList[isoform] = sect1+event[1][4]+sect2
                                                
                        for i in xrange(start, stop):
                            offsetVector[i] = None
            
                        for i in xrange(stop, seqLen):
                            offsetVector[i] = offsetVector[i] + deltaOffset
                            

        isoformReturnVal = {}
        salvage=0
        for isoformName in isoformSequenceList:
            
            try:
                isoformReturnVal[nametoIsoID[isoformName]] = [isoformName, isoformSequenceList[isoformName].lower()]
            except KeyError:                                
                missing_data=[isoformName, isoformSequenceList[isoformName].lower()]                       
                salvage=salvage+1

                
        # salvage works in the case where a single isoform was skipped and a single isoform name remains, meaning
        # chance are they should actually go together
        if salvage == 1:
            isoformSequenceList_keys = isoformSequenceList.keys()
            nametoIsoID_keys = nametoIsoID.keys()

            # remove 1
            try:
                nametoIsoID_keys.remove("1")
            except Exception:
                pass

            if len(set(isoformSequenceList_keys + nametoIsoID_keys))-len(nametoIsoID_keys) == 1:
                
                for i in nametoIsoID_keys:
                    if i not in isoformSequenceList_keys:
                        try:
                            isoformReturnVal[nametoIsoID[i]] = missing_data
                            show_warning("Malformed isoform data lead to an inconsistency, but the gremlins think they fixed it... [ID=" + self.accession + "]")
                        except KeyError:
                            show_warning("Malformed isoform data has lead to an inconsistency - skipping that isoform [ID=" + self.accession + "]")
                            # OH COME ON!!

        if salvage > 1:
            show_warning("Malformed isoform data has lead to an inconsistency - skipping several isoforms [ID=" + self.accession + "]")

            
                
            
                
    

        return isoformReturnVal
         
    
                    
      
                
    
#--------------------------------------------------------
#
#--------------------------------------------------------
# Returns a nice string defining the taxonomy of the protein species
#

    def _extract_taxonomy_string(self, taxonomy):
        return [ t.strip().lower() for t in taxonomy.split(';') ]

#--------------------------------------------------------
#
#--------------------------------------------------------
# Returns the species
#

    def _extract_species(self, xml):
        return xml['GBSeq_organism']


#--------------------------------------------------------
#
#--------------------------------------------------------
# Returns the proteins host species in the format <scientific> (<common>)
# if it's a viral protein,  else "N/A"
#

    def _extract_host(self, featureTable):
        for i in featureTable:
            for subElement in i["GBFeature_quals"]:
                if subElement['GBQualifier_name'] == "host":
                    return subElement["GBQualifier_value"]
                
        # if we don't find no host
        return "N/A"


#--------------------------------------------------------
#
#--------------------------------------------------------
# Returns the unique set of other protein accessions associated
# with this protein
#
    
    def _extract_other_accessions(self, xml):

        #===========================================================
        # internal/local function - avoid
        # class-space polution
        def parse_sequence_id(seqid):
            result = []
            
            if seqid.find('sp') == 0:
                _, sp_acc, sp_locus = seqid.split("|")
                #i = sp_acc.rfind(".")
                #if i > 0:
                #    sp_acc = sp_acc[0:i]
        
                seqid = sp_acc
                result.append(('Swissprot-Locus', sp_locus))

            else:
                # clip 'ref|...|' or dbj|...| wrapper off those 
                # which have it
                if seqid.find('ref') == 0 or seqid.find('dbj') == 0:
                    seqid = seqid[4:-1]
                    
                # clip gi| from front
                elif seqid.find('gi|') == 0 or seqid.find('GI|') == 0 :
                    seqid = seqid[3:]

                else:
                    splitSeqID = seqid.split("|")
                    if len(splitSeqID) == 3:
                        seqid = splitSeqID[1]
                        
                    
            # now we filter such that only accessions of known types
            # are stored [0.1.7 updated]

            seq_type = ProteinParser.ID_type(seqid)
            if seq_type[0] > -1:
                result.append((seq_type[1], seqid))
        
            return result
        #===========================================================

        # build source from which we parse
        acc_ids_to_parse = set(xml['GBSeq_other-seqids'])
        acc_ids_to_parse.add(self.accession)
        acc_ids_to_parse.add(xml['GBSeq_primary-accession'] )
        
        prot_accessions = []
        for seqid in acc_ids_to_parse:
            prot_accessions.extend(parse_sequence_id(seqid))

        return list(set(prot_accessions))


    def _get_qualifier(self, listOfNames, feature_quals):
        for q in feature_quals:

            # adds robustness for badly formatted XML
            try:
                if q['GBQualifier_name'] in listOfNames:
                    val = q['GBQualifier_value']
                    return val
            except KeyError:
                continue

    
            
