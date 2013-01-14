# geeneus
### Remote protein record access made simple

**NOTE: This github page hosts the development version, but not the distribution version. To get and install the current stable version of `geeneus` use `pip`**

       sudo pip install geeneus

**The command above will install the current stable release. **

### Introduction
**geeneus** is designed as a simple to use, robust and reliable Python API to obtain biological record information (currently only protein data is supported). The tool primarily uses NCBI's protein database, but falls back on the EBI UniProt database in a totally seamless fashion if needed.

The motivation comes from the simple fact that when I began working with NCBI's databases I wanted an interface which allowed me to do:

       interface.get_protein_sequence(accession_number)

and would just return the sequence associated with that accession number. This didn't exist, so I decided to create it.

The primary focus of **geeneus** from day 1 has been ease of use. NCBI allows access to their records through an Entrez based RESTful API called [eUtils](http://www.ncbi.nlm.nih.gov/books/NBK25500/). However, this can be complicated to set up, and to people who are less used to networking or programming can pose a major barrier to access. [Biopython](http://biopython.org/) goes some way [to help with this](http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc98), but still requires the user to parse XML, deal with networking, read handles, and lots of other things which are just more work.

Considering this, my goals were to;

* Parse the returned XML - the user should *never* have to deal with XML unless they want to (note the raw XML associated with a record is still easily available!)

* Deal with all the networking errors - the user should be able to turn their internet connection off while running and the system won't crash (it will probably stop getting results though)

* Abstract any of the complexity to create a uniform and easy to use, API which returns built-in types. This means the various functions only return strings, lists, tuples and dictionaries.

### Installation
The simplest way to is to just use `pip` (note there's no need to download anything here, `pip` downloads the source from the package index database and then installs it)

       pip install geeneus # (may need to be sudo)

More installation options can be found [here](http://pypi.python.org/pypi/Geeneus/)

### Usage

For now only protein record access usage is provided, as the gene record functionality is still in development. However, the backend provides an extendable framework for abstracting networking, caching records and error handling, such that the addition of additional data types is reasonably straight forward (i.e. the only additional code needed is record type specific XML parsing)

**geeneus** as a package contains a module for each type of record you might want to parse (protein, gene etc). To set up for protein records, you do;

       from geeneus import Proteome
       manager = Proteome.ProteinManager("your.emailaddress@email.com")

And from there, the NCBI/UniProt protein data is at your fingertips. Say we want the sequence for the *protein sprouty homolog 4 isoform 2*. This protein has the accession number "`NP_001120968`", so we simply do

       manager.get_protein_sequence("NP_001120968")

and we're greeted with the sequence (formatting for easy reading here)

       'meppipqsap ltpnsvmvqp lldsrmshsr lqhpltilpi dqvktshve ndyidnpsla 
        lttgpkrtrg gapelaptpa rcdqdvthhw isfsgrpssv ssssstssd qrlldhmapp          
        pvadqaspra vriqpkvvhc qpldlkgpav ppeldkhfll ceacgkckc kecasprtlp 
        scwvcnqecl csaqtlvnyg tcmclvqgif yhctneddeg scadhpcsc srsnccarws
        fmgalsvvlp cllcylpatg cvklaqrgyd rlrrpgcrck htnsvicka asgdaktsrp 
        dkpf'

For the full range of functions available see below, or try `help(manager)` or `help(geneeus.Proteome)`. 

## Functions
Below there is a brief reference list for the available functions. The '*datastore*' being mentioned here refers to the internal storage structure that holds the data as it's parsed and stored.

### List of access modules
    Proteome # for accessing protein data
    Genome   # for accessing gene information (currently non-functional)

### List of Proteome functions

#### Protein Manager initializer parameters (default in parentheses)
    email                 # required for NCBI database access

    cache (True)          # True or False: If set to True then the
                            datastore cache's records, while if False
                            a new record is downloaded on every                             
                            request. Unless you have a very specific
                            reason it is highly recommended to
                            keep this as True

    retry (0)             # Number of retries the NCBI networking
                            protocols use if access fails

    uniprotShortut (True) # If true we check UniProt first and
                            exclusivly for some record types (see
                            the section "Shortcutting" below)

#### Introspective functions (queries relating to the datastore)
    has_key(ID)         # check if an ID is currently cached in 
                          the datastore (does not trigger downloading
                          on False)

    keys()              # returns a list of all the protein IDs in the              
                          datastore

    datastore_size()    # get the number of records in the datastore

    purge()             # delete all data from the datastore

    error(ID)           # See if the record associated with this ID
                          caused an error (discussed below in more
                          detail)

    exists(ID)          # See if the protein record associated with 
                          this ID exists (discussed below in more  
                          detail)

    run_translation(ID) # function to translate an accession to a GI 
                          (uses a call to the NCBI lookup server)

#### Protein record functions (queries relating to the record)
    get_record_source(ID)        # determine which database this 
                                   record was downloaded from (will 
                                   return either 'NCBI' or 'UniProt')
                            
    get_record_creation_date(ID) # get the date the record was created
                                   in the database (not the datastore) 

    get_record_version(ID)       # get the record version (discussed     
                                   in more detail below)

    get_ID_type(ID)              # returns a 2 position list where  
                                   [0] is an exit code and [1] is 
                                   the name of the ID type (if it 
                                   meets the criterion for any 
                                   specific type

#### Protein information functions (queries relating to the protein)
    get_raw_xml(ID)                 # get an XML string for the record

    get_protein_name(ID)            # get the protein name

    get_protein_sequence(ID)        # get the amino acid sequence

    get_protein_sequence_length(ID) # get number of amino acids in 
                                      sequence

    get_geneID(ID)                  # if the protein is associated 
                                      with a specific gene ID then 
                                      this returns that ID

    get_gene_name(ID)               # if the protein is associated 
                                      with a specific gene then 
                                      this returns that gene's name

    get_variants(ID)                # get a list of variant 
                                      dictionaries (covered in detail 
                                      in a later section)

    get_isoforms(ID)                # get a list of isoforms indexed
                                      by their isoform reference 
                                      (covered in detail in a later
                                      section)

    get_domains(ID)                 # get a list of the PFAM 
                                      identified domains (covered in
                                      detail later)

    get_other_accessions(ID)        # get a list of other accessions
                                      which also point at this record

    get_taxonomy(ID)                # get the ordered taxonomy list of
                                      the species this protein  
                                      originates from

    get_species(ID)                 # get the species of origin

#### Batch protein information functions
    batch_get_protein_sequence(List) # returns a list which  
                                       corresponds with the input list
                                       of protein sequences

    batch_get_protein_name(List)     # returns a list which  
                                       corresponds with the input list
                                       of protein names

    batch_get_variants(List)         # returns a list which  
                                       corresponds with the input list
                                       of protein variants lists


#### Genome functions
    Currently untested, so best to ignore for the time being...
    
## Complex return types
Domains, variants, isoforms and other_accessions queries return complex structures (i.e. not a string or an integer). To understand what is being returned we briefly summarize them here. We also have a quick discussion on version numbers, UniProt shortcutting, and on `error` and `exists` statuses.

#### Domains
A domain query returns a list of domain dictionaries, where each dictionary has the following key value pairs;

    start        # domain start location
    stop         # domain stop location
    type         # domain type (currently always `pfam`, but left in for 
                   potential upgrades to different domain recognition 
                   approaches in the future)
    label        # label information
    accession    # accession associated with the domain

#### Variants
*UPDATE IN 0.1.6: Previously he key names were capitalized (e.g. 'Location', 'Original' etc). To add consistency with other complex types keys are now all lower case*

A variant query returns a list of variant dictionaries, where each dictionary has the following key value pairs;

    location   # start location of the variant being reported
    original   # original amino acid(s) at location
    mutant     # new amino acid(s) at location
    variant    # a convenient X -> Y string for easy visualization of
                 what the variant is
    type       # Describes the type of variant 
    notes      # any annotation provided in the reference database as
                 well as mutation references where possible

For `type` there are a number of possible values;

    "Deletion"                 # represent regions or amino acids
                                 which are missing (e.g. AGDDT -> -)
    "Deletion & substitution"  # represents the situation where the
                                 mutant is shorter and but other amino
                                 acids are added (e.g. AGH -> AK)
    "Insertion"                # represent situations where the
                                 mutated version is longer than the
                                 original but the original is still  
                                 present (e.g A -> AKL)
    "Insertion & substitution" # represent situations where the
                                 mutated version is longer than the 
                                 original and we lose the original 
                                 (e.g. A -> GKL)
    "Substitution (single)"    # Single amino acid switch (e.g. A->G)
    "Substitution (double)"    # Double amino acid switch 
                                 (e.g. AK -> GL)
    "Substitution (<X>)"       # Greater than 2 substitution where <X>
                                 is the length of the exchange 
                                 (e.g. if <X> = 4 then AKHI -> CDYW)



One potentially confusing is the use of DNA-typical vocabulary (insertion, deletion and subsitution) when talking about variant changes. It was decided that these provide appropriate and easy to understand terms, even if they typically refer to changes in DNA, not amino acid sequence. It's important to remember that a single substitution of an amino acid does not necessarily correspond to a SNP, and these descriptions refer *exclusively* to changes to the amino acid sequence, not the underlying DNA sequence.

#### Isoforms
Isoform queries return dictionary, where each key value is keyed by the isoform ID (which has the form "accession-X" e.g. Q12345-3). The value is a 2 position list, where position 0 is the isoform name/reference number and position 1 is the full isoform sequence. 

For example, a protein with one isoform (Q12345-2) might return the following dictionary

       {'Q1234-2':['2','GTAGHKLPKKLRSDF']}

Isoforms present some difficulties. Isoforms are defined in NCBI records in two different ways.

One set are based on feature annotation. A canonical record will have a number of features which describe the change in the various isoforms relative to that canonical sequence. Those isoforms don't have their own records - they exist as annotations in a canonical record and their sequences are represented as Q12345-X. **geeneus** can reconstruct these isoform amino acid sequences based on the annotations, which is how it builds the NCBI isoform lists (UniProt isoforms are each individually defined and linked to the reference record, so no reconstruction is required). One additional complication is that each isoform has both an isoform ID (Q12345-X) and an isoform name (isoform Y). While often X and Y are the same this is not always the case, sometimes Y is a word name (such a small, short, beta etc). More unhelpfully, sometimes Y is a number, but a *different* number to the isoform ID X value. 

The other set exist as stand alone records, often with "isoform" in their name. These records unhelpfully includes no references to other isoforms, or to the canonical sequence/record from which they relate to. Similarly, canonical sequences have no reference to their "isoform" records. As such there is no way to include these in the list of isoforms (although they do have their own accession values as a result of existing as their own records). It is therefore up to the user to identify such isoforms. Theoretically, you could query the entire list of cross-referenced accessions and compare names to try and identify isoforms, but this would be a hugely expensive (in terms of network traffic) and there's no guarantee any/all isoforms would be included in the cross reference section.

**NB:** If you want to request a specific isoform sequence using the dash notation, doing 

       manager.get_protein_sequence("Q12345-5") 

will return **the canonical sequence**, not the isoform sequence.

Instead, you need to do

       manager.get_isoforms("Q12345-5")["Q12345-5"][1]
       # or
       manager.get_isoforms("Q12345")["Q12345-5"][1]

This will do the following

* Get the record associated with `Q12345`, then get the isoforms
* Pull out the isoform data associated with `Q12345-5`
* As discussed above, the isoform data is a tuple of (`name`, `sequence`) so position 1 gives you the sequences

It would probably be wise to wrap this in a `try/except` block in case the isoform ID is missing, such as

       try:
          manager.get_isoforms("Q12345-5")["Q12345-5"][1]
       except KeyError:
          print "Isoform ID not found"



#### Other accession
Other accession queries returns a list of tuples which define, (`type of accession`, `accession`).  The type of accession will be one of the following;

* "UniProtKB/Swiss-Prot"
* "RefSeq"
* "GI"
* "PDB"
* "International Protein Index"
* "DDBJ"
* "GenBank"
* "EMBL"
* "Unknown accession type"

The way **geeneus** classifies these accession types is based on a set of hardcoded regular expressions. If you expect one accession type to be classified as something it's not, this may be an error in either how these values are parsed or the defining regular expressions, so please submit a bug report!

## Record versions
Accession versioning is done by appending a period followed by a version number, e.g. Q12345.1 or Q12345.2 would be two different versions of the same record. Versioning represents updates made to records, typically as new information becomes available.

One potential source of problems is that NCBI records obtained directly through the eUtils interface (as opposed to through the website) do not contain any information on related versions. This means that **geeneus** is unable to give this information either.

Querying a non-versioned accession (e.g. Q12345 or NP_1234567) will give the most up-to-date record associated with that accession, while querying a versioned value (Q12345.3 or NP_1234567.5) will give that specific version. However, there is no way to know if any specific versioned record is the most up-to-date record, or access previous records. This is not *necessarily* a problem, it's just worth being aware that if you query with explicit version numbers this may not give the most up to date version. 

Note that GI numbers are unique for each different version, so deal with versioning in a different manner. The `version` number returned here refers to the non-GI accession version, where available. If no explicit version is available then we assume the version is 1.

## Shortcutting
Shortcutting allows **geeneus** to bypass the NCBI servers entirely for UniProtKB/Swiss-Prot accession values which NCBI doesn't guarantee to support (e.g. all of them except those beginning with O, P or Q) and go directly to the UniProt servers, both in batch mode and individual mode.

This leads to a *massive* increase in performance when dealing with large number of accessions which meet this criterion as we avoid countless database misses, retry-rounds, and the eventual fallback to UniProt and just go straight there from the word go. 

## Error and Exists status
`error` and `exists` represent two flags to help deal with problems.	If a networking request fails, this sets a records `error` flag to true, and any of the other methods return `None`. In this case, `exist` is also set to False, as the record does not exist in the datastore. 

Sometimes we may obtain a record correctly and without error, but find that it does not represent a valid protein record. If this happens we have not triggered an error (so error=False), but the record does not exists, in which case `exist` is also set to False. In these instances all methods *except* the `get_raw_xml()' method return `None`, but `get_raw_xml()` will return the xml for further inspection.

*UPDATE IN 0.1.6: Previously if an error=True or exists=False the methods would return empty strings, lists of dictionaries. This can be dangerously ambiguous, so the behavior was changed to return None (equivalent to NA, or unknown)*

## Design Decisions

### UniProt fallback
**NB** *The discussion below is more relevant when* `uniprotShortcut` *is set to* `False`.

NCBI guarantees support for UniProt/Swissprot IDs which begin with O, P or Q. However, it also offers some support for other types of UniProt IDs (i.e. those which begin with other letters). There is sometimes a discrepancy between those which are available through the website and through the eUtils interface, where the eUtils lookup fails on an accession that should succeed. To deal with this, **geeneus** can use UniProt IDs and fall back on the UniProt servers. Dealing with UniProt calls instead of NCBI calls is more expensive in terms of network traffic because it requires an addition batch request to the PFAM servers to define the domains. The UniProt records only hold references to PFAM domains, but lack the necessary details to build informative domain data structures. However, like NCBI, UniProt servers do offer a batch request mode, which is utilized in `batch_get`  methods.

If `uniprotShortcut` is set to `True` then we default to the UniProt servers instead of NCBI for those accessions not beginning with O, P or Q. However, **geeneus** still provides UniProt fallback for those O/P/Q records, should they fail on NCBI lookup.

The user is totally oblivious to this behavior - **geeneus** provides an entirely uniform access to the information regardless of its source. To check which database a record has come from you can do

       manager.get_record_source(ID)

Which returns the remote database type from which the record associated with accession ID was obtained. Note that if ID is not already in the datastore this will trigger the record to be fetched.

### Caching 
Local and temporary caching is an important feature of **geeneus** which helps make it an ideal tool for interactive data exploration. The manager object builds up a local data structure, and by default caches the records it fetches from remote database. The upshot of this from the user's perspective is that if I run

    manager.get_protein_sequence("NP_001120968")
    manager.get_protein_sequence("NP_001120968")

The second call doesn't query the database, but just reads off the cached value. This caching behavior can be turned off on by setting `cache=False` when initializing the ProteinManager object.


### Batch queries 
A key design decision was how to deal with batch queries.

The eUtils recommended approach for making large (100+ IDs) queries is to initially ePost a list of those queries. The ePost operation sends this list to the ENTREZ servers, returning a `WebEnv` value and a `QueryKey` value. These two can then be used with an `eFetch` to go to the sever and get the result of the list submitted previously. The difficulty is that this list *must* be made up of UIDs (unique identifiers) which for proteins means GI numbers. If you only have an accession value (as is common) the only way to get this GI number is to query the database, and this `eSearch` operation can *not* be done in batch. Essentially, this is a chicken and egg problem - to get the GIs we need to do an `ePost` based batch query we have to run serial database queries.

This means that using `ePost`/`eFetch` would be great if you had a list of UIDs, but in practice accession numbers are a lot more common and useful, and the mapping of *n* accession values to UIDs would require *n* calls the server anyway. 

To get around this, I use concatenated `eFetch` calls for batch queries, whereby a single call is submitted with a list of IDs. This is a fast and stable way to get around this problem, and, so far, as shown no issues with lists up to 100 IDs long. The potential issue is that the HTTP GET request being made here literally gets longer as we add accession values, so this represents a top limit in terms of networking protocols. However, I implemented a recursive cascading retry mechanism which halves the list and retries each half, so should a list be too long it should only result in two calls instead of one.

### Robustness
A primary goal with **geeneus** was to create an API which is robust to input. By this, it should be able to handle case insensitivity, convert accession values where necessary, and correctly recognize valid accession values while rejecting irrelevant ones to minimize server burden. For accession filtering, we use regular expressions to ensure the only accession values which we query could be real values (based on NCBI's [accession rules](http://www.ncbi.nlm.nih.gov/Sequin/acc.html). While PDBs don't fall into this category, we allow translation between PDBID and GI, although often a chain identifier is required as the NCBI protein database treats separate chains as separate records.

All the networking is dealt with in a highly modular fashion, and network failure tolerance is a priority.


## Background and licence
This code was developed by [Alex Holehouse](http://holehouse.org) at [Washington University in Saint Louis](http://www.wustl.edu/) as part of the [Naegle lab](http://naegle.wustl.edu/people/lab_members.html). It is licensed under the the GNU General Public License (GPL-2.0). For more information see LICENCE.

### Acknowledgement
A truly massive hat tip to Matt Matlock for countless suggestions, discussions and constant feedback regarding geeneus in a production setting. Matt has significantly shaped this code for the better, and it would be nothing like it is today without his input (the need for isoform support, domain support, the idea for shortcutting etc etc).
