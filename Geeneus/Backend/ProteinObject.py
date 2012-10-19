# Provides an object-based API for a protein (private)
#
# Copyright 2012 by Alex Holehouse - see LICENSE for more info
# Contact at alex.holehouse@wustl.edu


from Bio import Entrez, Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC

import string
import re

# object attributes
# self.sequence - protein sequence
# self.exists - does the object associated with this ID exist in the database
# self.error - ALWAYS false in a ProteinObject (true in ProteinErrorObject)
# self.sequence_create_date - date the sequence was entered into the protein database
# self.protein_variants - list of dictionaries, each dictionary corresponding to a unique
#                     variant. Each dictionary has location, mutation and notes
# self.sequence_length
# self.geneID - NCBI GeneID for the protein, should you want to lookup Gene information

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
        
    def exists(self):
        return self._exists

    def error(self):
        return self._error


#-------------------------------------------------------
#=======================================================
#
#--------------------------------------------------------
# Object initializer. 
#
# Initializes all the objects attributes to default values before
# populating with proteinxml based data. If no xml data is available
# default values are not overwritten, so rather than exceptions being
# raised on request for non-existant information default values are, and
# exists remains set to False
#
    def __init__(self, proteinxml):

        # set the default values (these are kept for empty/
        # error calls
        
        self.sequence = ""
        self._exists = False
        self._error = False
        self.sequence_create_date= "01-JAN-1900"
        self.protein_variants = []
        self.geneID = 0
        self.sequence_length = 0
        self.name = ""
        
        if proteinxml == -1:        
            self._error = True
            return

        if not self._xml_is_OK(proteinxml):
            return
        
        self._exists = True

        # Now we set the rest of the values using the parsed XML
        self.raw_XML = proteinxml
        self.sequence = proteinxml[0]["GBSeq_sequence"]
        self.sequence_length = len(self.sequence)
        self.sequence_create_date = proteinxml[0]["GBSeq_create-date"]
        self.protein_variants = self._extract_variant_features(proteinxml[0]["GBSeq_feature-table"])        
        self.geneID = self._extract_geneID(proteinxml[0]["GBSeq_source-db"])
        #self.name = self._extract_names[proteinxml[0]["GBSeq_definition"]]

#--------------------------------------------------------
#
#--------------------------------------------------------
# Function to check that the xml we've downloaded is good, and represents a 
# viable protein xml structure. Additional tests may be added here as we find more
# edge cases! Returns FALSE if there's a problem, TRUE otherwise

    def _xml_is_OK(self, proteinxml):
        if len(proteinxml) > 1:
            print "WARNING [ProteinObject._xml_is_ok()] - ProteinXML detected more than one record associated with this GI.\nThis should never happen."
        
        # Nothing in XML - so return an empty-initiailized object with exists = 0
        if len(proteinxml) == 0:
            return False
        
        # Check that we're really dealing with protein (despite specifying db="protein"
        # on the efetch call, when a GI is used other databases seem to be searched too...
        if not (proteinxml[0]["GBSeq_moltype"] == "AA"):
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
    def _extract_geneID(self, GBSeq_source_db):
        source = str(GBSeq_source_db)
    
        # See if there's a geneID in the DB source data, and if not
        # return -1
        geneID_location = string.find(source, "GeneID:")
        if (geneID_location == -1):
            return "No GeneId Found"
        
        # Given we found a "GeneID:" label in the text, we cut out the
        # value subsequent to the tag before the next comma, and return
        # that as the GeneID
        geneID_location = geneID_location+7
        geneID_end_location = string.find(source[geneID_location:], ",") + geneID_location
        geneID = source[geneID_location:geneID_end_location]
        
        return(int(geneID))

