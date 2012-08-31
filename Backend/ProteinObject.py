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
# self.protein_snps - list of dictionaries, each dictionary corresponding to a unique
#                     SNP. Each dictionary has location, mutation and notes
# self.sequence_length
# self.GeneID - NCBI GeneID for the protein, should you want to lookup Gene information

class ProteinObject:

#--------------------------------------------------------
# Getter functions should be used rather than direct member access - try and
# maintain some encapsulation.
#
    def get_GeneId(self):
        return self.GeneId

    def get_protein_sequence(self):
        return self.sequence

    def get_SNPs(self):
        return self.protein_snps

    def get_protein_sequence_length(self):
        return self.sequence_length

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
        self.protein_snps = []
        self.GeneId = 0
        self.sequence_length = 0
        
        if proteinxml == -1:        
            self._error = True
            return

        if not self.__xml_is_OK(proteinxml):
            return
                
        self._exists = True

        # Now we set the rest of the values using the parsed XML
        self.sequence = proteinxml[0]["GBSeq_sequence"]
        self.sequence_length = len(self.sequence)
        self.sequence_create_date = proteinxml[0]["GBSeq_create-date"]
        self.protein_snps = self.__extract_snp_features(proteinxml[0]["GBSeq_feature-table"])        
        self.GeneId = self.__extract_GeneId(proteinxml[0]["GBSeq_source-db"])

#--------------------------------------------------------
#
#--------------------------------------------------------
# Function to check that the xml we've downloaded is good, and represents a 
# viable protein xml structure. Additional tests may be added here as we find more
# edge cases! Returns FALSE if there's a problem, TRUE otherwise

    def __xml_is_OK(self, proteinxml):
        if len(proteinxml) > 1:
            print "WARNIN [ProteinObject.__xml_is_ok()] - ProteinXML detected more than one record associated with this GI.\nThis should never happen."
        
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
#
#--------------------------------------------------------
# Function to get SNP features. Easy to extended should extra
# SNP data be needed, but for the function creates a list of n
# dictionaries (where n = number of SNPs in protein xml data)
# and each dictonary contains SNP location, mutation and notes.
#
    def __extract_snp_features(self, featurelist):

        SNP_list = []
                
        for feature in featurelist:
            for feature_subsection in feature["GBFeature_quals"]:
                featurematch = re.match("[QWERTYIPASDFGHKLCVNM] -> [QWERTYIPASDFGHKLCVNM]*",feature_subsection["GBQualifier_value"])
                if featurematch:
                    variant = True
                    break
            if featurematch:
                temp_dic = {"SNP" : featurematch.string[:6]}
                temp_dic["Notes"] = featurematch.string[7:]
                temp_dic["Location"] = int(feature["GBFeature_location"])
                SNP_list.append(temp_dic)
                del(temp_dic)

        if len(SNP_list) == 0:
            return [0]
        else:
            return SNP_list
        
#--------------------------------------------------------
#
#--------------------------------------------------------
# Function to extract the GeneId from the protein data for 
# use in getting gene information from the Genome class if
# needed
#
    def __extract_GeneId(self, GBSeq_source_db):
        source = str(GBSeq_source_db)
    
        # See if there's a GeneId in the DB source data, and if not
        # return -1
        GeneId_location = string.find(source, "GeneID:")
        if (GeneId_location == -1):
            return "No GeneId Found"
        
        # Given we found a "GeneID:" label in the text, we cut out the
        # value subsequent to the tag before the next comma, and return
        # that as the GeneID
        GeneId_location = GeneId_location+7
        GeneId_end_location = string.find(source[GeneId_location:], ",") + GeneId_location
        GeneId = source[GeneId_location:GeneId_end_location]
        
        return(int(GeneId))
