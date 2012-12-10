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

# Object attributes
# 
#
# self.sequence - protein sequence
# self.exists - does the object associated with this ID exist in the database
# self.error - ALWAYS false in a ProteinObject (true in ProteinErrorObject)
# self.sequence_create_date - date the sequence was entered into the protein database
# self.protein_variants - list of dictionaries, each dictionary corresponding to a unique
#                     variant. Each dictionary has location, mutation and notes
# self.sequence_length
# self.geneID - NCBI GeneID for the protein, should you want to lookup Gene information
# self.species - Species of origin
# self.taxonomy - Taxonomy string list
# self.domain_list - list of pfam defined domains
# self.gene_name - gene name
# self.other_accessions - list of other accesison values

class ProteinObjectException(BaseException):
    pass

class ProteinObject:

#--------------------------------------------------------
# Getter functions should be used rather than direct member access - try and
# maintain some encapsulation.
#
    def get_geneID(self):
        return self.geneID

    def get_protein_sequence(self):
        return self.sequence

    def get_variants(self):
        return self.protein_variants

    def get_protein_sequence_length(self):
        return self.sequence_length

    def get_protein_name(self):
        return self.name

    def get_raw_xml(self):
        return self.raw_XML

    def get_other_accessions(self):
        return self.other_accessions

    def get_species(self):
        return self.species

    def get_taxonomy(self):
        return self.taxonomy

    def get_domains(self):
        return self.domains

    def get_gene_name(self):
        return self.gene_name

    def get_isoforms(self):
        return self.isoforms
        
    def exists(self):
        return self._exists

    def error(self):
        return self._error


#-------------------------------------------------------
#=======================================================



#--------------------------------------------------------
# Object initializer
# Multiple dispatch because no function overloading is allowed
# in Python. We need to be able to create ProteinObjects from
# both NCBI XML, but also from pre-extracted values from
# other sources (e.g. UniProt XML)
#
    def __init__(self, *args):

        # if we have a single xml string
        if len(args) == 2:
            self.__init_1(*args)

        # if we're passing in pre-parsed data
        else:
            self.__init_2(*args)
        
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

        #proteinxml = proteinxml[0]
        
        self.accession = accession
        self.sequence = ""
        self._exists = False
        self._error = False
        self.sequence_create_date= "01-JAN-1900"
        self.protein_variants = []
        self.geneID = ""
        self.sequence_length = 0
        self.name = ""
        self.other_accessions = []
        self.species = ""
        self.taxonomy = ""
        self.domains = ""
        self.gene_name = ""
        self.isoforms = {}

        
        if proteinxml == -1:        
            self._error = True
            return

        if not self._xml_is_OK(proteinxml):
            return

        self._exists = True

        # Now we set the rest of the values using the parsed XML
        self.raw_XML = proteinxml[0]
        self.sequence = proteinxml[0]["GBSeq_sequence"]
        self.sequence_length = len(self.sequence)
        self.sequence_create_date = proteinxml[0]["GBSeq_create-date"]
        self.protein_variants = self._extract_variant_features(proteinxml[0]["GBSeq_feature-table"])        
        self.geneID = self._extract_geneID(proteinxml[0]["GBSeq_source-db"], proteinxml[0]["GBSeq_feature-table"])
        self.name = self._extract_protein_name(proteinxml[0]["GBSeq_definition"])
        self.other_accessions = self._extract_other_accessions(proteinxml[0])
        self.species = self._extract_species(proteinxml[0])
        self.taxonomy = self._extract_taxonomy_string(proteinxml[0]['GBSeq_taxonomy'])
        self.domains = self._extract_domain_list(proteinxml[0]["GBSeq_feature-table"])
        self.gene_name = self._extract_gene_name(proteinxml[0]["GBSeq_feature-table"])
        self.isoforms = self._extract_isoforms(proteinxml[0]["GBSeq_feature-table"], accession, self.sequence)


#--------------------------------------------------------
# Autocreate initializer
# If you already have the relevant data you can build a ProteinObject directly. Useful for non NCBI based record
# access (e.g. UniProt)
#
    def __init_2(self, accession, xml, name, mutations, sequence, creation_date, geneID, gene_name, other_accessions, species, domains, taxonomy, isoforms):
        
        self.accession = accession
        self._exists = True
        self._error = False
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

#--------------------------------------------------------
#
#--------------------------------------------------------
# Function to check that the xml we've downloaded is good, and represents a 
# viable protein xml structure. Additional tests may be added here as we find more
# edge cases! Returns FALSE if there's a problem, TRUE otherwise

    def _xml_is_OK(self, proteinxml):
        if len(proteinxml) > 1:
            print "WARNING [ProteinObject._xml_is_ok()] - ProteinXML detected more than one record associated with this GI.\nThis should never happen."
            return False
        
        # Nothing in XML - so return an empty-initiailized object with exists = 0
        if len(proteinxml) == 0:
            return False

        # Check that we're really dealing with protein (despite specifying db="protein"
        # on the efetch call, when a GI is used other databases seem to be searched too...
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

        variant_list = []
                
        for feature in featurelist:            

            if not feature.has_key("GBFeature_quals"):
                continue
            
            for feature_subsection in feature["GBFeature_quals"]:
                
                if not feature_subsection.has_key("GBQualifier_value"):
                    continue
            
                # look for single variant
                featurematch_single = re.match("[QWERTYIPASDFGHKLCVNM] -> [QWERTYIPASDFGHKLCVNM]",feature_subsection["GBQualifier_value"])
                
                # look for double variant
                featurematch_double = re.match("[QWERTYIPASDFGHKLCVNM][QWERTYIPASDFGHKLCVNM] -> [QWERTYIPASDFGHKLCVNM][QWERTYIPASDFGHKLCVNM]",feature_subsection["GBQualifier_value"])
                
                if featurematch_single or featurematch_double:
                    break

            # deal with singles
            if featurematch_single:
                temp_dic = {"Variant" : featurematch_single.string[:6]}
                temp_dic["Original"] = featurematch_single.string[:1]
                temp_dic["Mutant"] = featurematch_single.string[5:6]
                temp_dic["Type"] = "Single"
                temp_dic["Notes"] = featurematch_single.string[7:]
                temp_dic["Location"] = feature["GBFeature_location"]
                variant_list.append(temp_dic)
                del(temp_dic)
                
            # deal with doubles
            if featurematch_double:
                temp_dic = {"Variant" : featurematch_double.string[:8]}
                temp_dic["Original"] = featurematch_double.string[:2]
                temp_dic["Mutant"] = featurematch_double.string[6:8]
                temp_dic["Type"] = "Double"
                temp_dic["Notes"] = featurematch_double.string[9:]
                temp_dic["Location"] = feature["GBFeature_location"]
                variant_list.append(temp_dic)
                del(temp_dic)

        if len(variant_list) == 0:
            return []
        else:
            return variant_list
        
        
#--------------------------------------------------------
#
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
#
#--------------------------------------------------------
# Extract the gene name from the source XML
#
    def _extract_gene_name(self, feature_table):
        gene_name = ""

        for f in feature_table:

            # look for explicit gene fields - if we find one return
            # automatically
            if  f [ 'GBFeature_key' ]  ==  'gene' :
                gene_name  =  self._get_qualifier ( 'Gene' ,  f [ 'GBFeature_quals' ])
                return gene_name

            # look for coding sequence names - alternative to gene
            # name but may be overridden if an explicit gene name also exists
            if  f [ 'GBFeature_key' ]  ==  'CDS' :
                gene_name = self._get_qualifier( 'gene', f['GBFeature_quals'])

        return gene_name
                

#--------------------------------------------------------
#
#--------------------------------------------------------
# Build a list of dictionaries containing domains
#

    def _extract_domain_list(self, feature_table):

        domainList = []

        for f in feature_table:
            tempDomainDictionary = {}
            note_val = self._get_qualifier('note', f['GBFeature_quals'])
            if f['GBFeature_key'] == 'Region' and note_val != None and note_val.find('pfam') > 0:
                
                tempDomainDictionary["type"] = "pfam"
                tempDomainDictionary["label"] = self._get_qualifier('region_name', f['GBFeature_quals'])
                tempDomainDictionary["accession"] = f['GBFeature_intervals'][0]['GBInterval_accession']
                tempDomainDictionary["start"] = int(f['GBFeature_intervals'][0]['GBInterval_from'])
                tempDomainDictionary["stop"] = int(f['GBFeature_intervals'][0]['GBInterval_to'])
                
                # add the domainDictionary element
                domainList.append(tempDomainDictionary)

        return domainList

#--------------------------------------------------------
#
#--------------------------------------------------------
# Returns a nice string defining the taxonomy of the protein species
#
    def _extract_isoforms(self, ft, ID, sequence):
        

        #===========================================================
        # Local exception, because isoforms are serius buznuz
        #
        class IsoformException(BaseException):
            pass
        
        #===========================================================
        #
        # Function which takes the "note" string from a splicing variant
        # record and determins the isoform numbers to which that splicing
        # variant applies.
        # 
        # The splicing variant "note" structure is typically something like
        #    '<isoform change>' (in isoform 1 [and 2 and 3....]'
        # Where 1, 2 and 3 represent the relevant isoform numbers
        # This method pulls out those numbers, returning a list of them
        # 
        def getRelevantIsoforms(defString):
            isoforms = []   
            locations = [m.start() for m in re.finditer('isoform', defString)]
            
            for i in xrange(0,len(locations)):

                # REGEX breakdown
                # For clairty we can break this REGEX down
                # ([\w| |\-|\+]*)   match any number of alphanumeric characters or 
                #                   white space or - or +
                # ( and|,|\)))      until you get to " and" or "," or ")" in which 

                options = []
                
                # 5 isoform 6 isoform  --> 5
                try:
                    isoSearch = re.search("(^[\w| |\.|\-|\+]*)(?= isoform)", defString[locations[i]+8:]).group()
                    options.append(isoSearch)
                except AttributeError:
                    pass
               
                # 2 and 4 .... -> 2
                try:
                    andSearch = re.search("(^[\w| |\.|\-|\+]*)(?= and)", defString[locations[i]+8:]).group()
                    options.append(andSearch)
                except AttributeError:
                    pass
                    
                # 2, 3 and ... -> 2
                try:
                    commaSearch = re.search("(^[\w| |\.|\-|\+]*)(?=,)", defString[locations[i]+8:]).group()
                    options.append(commaSearch)
                except AttributeError:
                    pass

                # 2) var=DB.... -> 2
                try:
                    parenSearch = re.search("(^[\w| |\.|\-|\+]*)(?=\))", defString[locations[i]+8:]).group()
                    options.append(parenSearch)
                except AttributeError:
                    pass
                
                try:
                    name = options[0]
                except IndexError:
                    raise IsoformException("ERROR while trying to identify isoform names - no names conformed to rule scheme")
                
                # finally now we have a list of possible names, we find which of these is the shortest and select it 
                for ID in options:
                    if len(ID) < len(name):
                        name = ID

                isoforms.append(name)

            return isoforms
        
        #
        # Function which takes the 'note' string (as defined above) along with the 
        # domain bounds and encodes it into a list with the following structure
        # 
        # [TYPE OF SPLICE VARIANT, START, END, ORIGINAL VALUE, NEW VALUE]
        # Note that
        # - type is one of "missing" or "replacement"
        # - if "replacement", original and new values represent the switcharoo
        # - if "missing" then original and new values are empty, we just cut out the region 

        def getSpliceEvent(defString, location):
            
            defString = defString.lower()

            # grab location information
            try:
                locationList = [int(location["GBInterval_from"])-1, int(location["GBInterval_to"])]
            except KeyError:
                # some regions are defined as  a single point
                locationList = [int(location["GBInterval_point"])-1, int(location["GBInterval_point"])] 
            
            # missing event   
            if defString.find("missing") > -1:
                return["missing", locationList[0], locationList[1], "", ""]
            
            # replacement event
            if defString.find("->") > -1:

                # to get the string of AA either side of -> we have to do some
                # quick regexs/reformatting
                defString = defString.replace(" ","")
                
                old = str(defString.split("->")[0])
                new = re.search("[a-zA-Z]*", defString.split("->")[1]).group()

                return ["replacement", locationList[0], locationList[1], old, new]
            
            raise ProteinObjectException("Note string had neither missing nor -> in it - parse error!") 
        
        # Possible overlao
        #  1 2 3 4 5 6
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

        
        #===========================================================

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
    
        return isoformSequenceList 
         
    
                    
      
                
    
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
                i = sp_acc.rfind(".")
                if i > 0:
                    sp_acc = sp_acc[0:i]
        
                result.append(('Swissprot', sp_acc))
                result.append(('Swissprot', sp_locus))

            else:
                # clip 'ref|...|' or dbj|...| wrapper off those 
                # which have it
                if seqid.find('ref') == 0 or seqid.find('dbj') == 0:
                    seqid = seqid[4:-1]

                seq_type = ProteinParser.ID_type(seqid)[1]
                result.append((seq_type, seqid))
        
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


    def _get_qualifier(self, name, feature_quals):
        qualifier = [ q for q in feature_quals if q['GBQualifier_name'] == name ]
        if len(qualifier) == 1:
            return qualifier[0]['GBQualifier_value']



    
            
