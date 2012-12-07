# Public facing API for accessing protein information
#
# Copyright 2012 by Alex Holehouse - see LICENSE for more info
# Contact at alex.holehouse@wustl.edu

import Bio.Entrez
import geeneus.backend

class ProteinManager:
    
    def __init__(self, email, cache=True, retry=0):
        """Returns a fully formed manager object which can be queried by the other 
        functions in this class. 

        email    Must be a valid email, as is required by the NCBI servers. For more 
                 information on NCBI usage guidelines please see 
                 [http://www.ncbi.nlm.nih.gov/books/NBK25497/].

        cache    Determines if the manager object should cache requests in memory, or 
                 the NCBI database should be queried every time a request is made. 

        retry    The number of times the networking utilities will retry on a failed 
                 connection

        timeout  The number of seconds the networking utilities wait after making a 
                 request before deciding that request has failed
        """

        self.datastore = geeneus.backend.ProteinParser.ProteinRequestParser(email, cache, retry)
        if self.datastore.error():
            self.error_status = True
        else:
            self.error_status = False


    def has_key(self, ID):
        """ Check if the datastore has a protein ID loaded or not. Does
            not preload the ID in question
        """
        return self.datastore.has_key(ID)

    def keys(self):
        """ Return a list of the keys in the datastore (note this hides 
        the always present '-1' key
        """
        return self.datastore.keys()
    
    def get_protein_name(self, ID):
        """ Returns the name of the protein """
        return self.datastore.get_protein_name(ID)
    
    def get_protein_sequence(self, ID):
        """ Return the protein's primary amino acid sequence as a string """
        return self.datastore.get_sequence(ID)

    def get_raw_xml(self, ID):
        """ Return the raw XML associated with this accession value """
        return self.datastore.get_raw_xml(ID)

    def get_variants(self, ID):
        """ Return a list of dictionaries, where each dictionary is a mutation 
            dictionary with six keys - 

            Location   Position in the primary sequence
            Original   Original amino acid 
            Mutant     Mutant amino acid
            Type       Type of mutation (double or single)
            Variant    Single term summary showing "Original -> Mutant"
            Notes      Annotation notes from download 
         """
            
        return self.datastore.get_variants(ID)

    def get_geneID(self, ID):
        """ Proteins are also associated with specific genes. This returns 
            the gene ID associated with this accession number
        """
        return self.datastore.get_geneID(ID)

    def get_gene_name(self, ID):
        """ Get the name of the gene associated with this protien (may be empty if 
            no true gene is associated)
        """
        return self.datastore.get_gene_name(ID)

    def get_taxonomy(self, ID):
        return self.datastore.get_taxonomy(ID)

    def get_domains(self, ID):
        return self.datastore.get_domains(ID)

    def get_species(self, ID):
        return self.datastore.get_species(ID)

    def get_other_accessions(self, ID):
        return self.datastore.get_other_accessions(ID)
    
    def get_protein_sequence_length(self, ID):
        """ Returns an integer equal to the length of the protein's primary sequence """
        return len(self.datastore.get_sequence(ID))

    def get_ID_type(self, ID):
        """ Returns a two position list, where list[0] is an exit code 
            and list[1] is the name of the type of accession number 
        """
        return self.datastore.get_ID_type(ID)

    def run_translation(self, Acc):
        """ Translates an alphanumeric accession number to a GI number """
        return self.datastore.translate_Asc2GI(Acc)

    def batch_get_protein_sequence(self, IDList):
        """ Allows for primary sequence retrieval en-masse using a list of IDs
            as input. Returns a dictionary, indexed by accession number, of
            the various sequences. Makes a single request to the NCBI server """
        return self.datastore.batchFetch(self.datastore.get_sequence, IDList)

    def batch_get_protein_name(self, IDList):
        """ Get the names of a list of IDs """
        return self.datastore.batchFetch(self.datastore.get_protein_name, IDList)

    def batch_get_variants(self, IDList):
        """ Allows variants retrieval en-mass using a list of IDs as input.
            Returns a dictionary of lists, where each dictionary entry is 
            indexed by a protein ID from the list and each list is made up of
            a dictionaries which contian variant information """
        return self.datastore.batchFetch(self.datastore.get_variants, IDList)

    def purge(self):
        """ Wipe the managers memory store. Only relevant if you're worried
            about memory or running a script as a daemon which is constantly
            pulling down data for a one time use """
        self.datastore.purge_data_store()

    def get_size_of_datastore(self):
        """ Get the number of items in the internal datastore """
        return self.datastore.get_size_of_datastore()

    
    
