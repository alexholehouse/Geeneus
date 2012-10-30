# Geeneus
#### NCBI database access made simple

**NOTE: This github page hosts the development version, but not the distribution version. To get and install Geeneus, go to [the PyPI site](http://pypi.python.org/pypi/Geeneus/) and get the installable source from there.**

### Introduction
Geeneus is designed as a simple to use, robust and reliable API to NCBI's various databases. Right now it provides in depth access for protein records, with limited gene record functionality. However, over time I plan to add more database types and more in depth functionality.

The motivation comes from the simple fact that when I began working with NCBI's databases I wanted an interface which allowed me to do:

    interface.get_protein_sequence(accession_number)

and would just return the sequence associated with that accession number. This didn't exist, so I decided to create it.

The primary focus of Geeneus from day 1 has been ease of use. NCBI allows access to their records through an ENTREZ based RESTful API called [eUtils](http://www.ncbi.nlm.nih.gov/books/NBK25500/). However, this can be complicated to set up, and to people who are less used to networking or programming can pose a major barrier to access. [Biopython](http://biopython.org/) goes some way [to help with this](http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc98), but still requires the user to parse XML, deal with networking, read handles, and lots of other things which are just more work.

Considering this, my goals were to;

* Parse the returned XML - the user should *never* have to deal with XML unless they want to

* Deal with all the networking errors - you should be able to turn your internet connection off while running and the system doesn't crash (it will probably stop getting results though)

* Abstract any of the complexity to create a uniform, easy to use, Python built-in types based API. This means the various functions only return strings, lists and dictionaries.

### Installation
The simplest way to is to just use `pip` (note there's no need to download anything here, `pip` downloads the source from the package index database and then installs it)

    pip install Geeneus # (may need to be sudo)

More installation options can be found [here](http://pypi.python.org/pypi/Geeneus/)

### Usage

For now only protein record access usage is provided, as the gene record functionality is still in development.

Geeneus as a package contains a module for each type of record you might want to parse (protein, gene etc). To set up for protein records, you do;

    from geeneus import Proteome
    manager = Proteome.ProteinManager("your.emailaddress@email.com")

And from there, the NCBI protein data is at your fingertips. Say we want the sequence for the *protein sprouty homolog 4 isoform 2*. This protein has the accession number NP_001120968, so we simply do

    manager.get_protein_sequence("NP_001120968")

and we're greeted with the sequence

    'meppipqsapltpnsvmvqplldsrmshsrlqhpltilpidqvktshvendyidnpslalttgpkrtrggapelaptparcdqdvthhwisfsgrpssvs
     sssstssdqrlldhmapppvadqaspravriqpkvvhcqpldlkgpavppeldkhfllceacgkckckecasprtlpscwvcnqeclcsaqtlvnygtc
     mclvqgifyhctneddegscadhpcscsrsnccarwsfmgalsvvlpcllcylpatgcvklaqrgydrlrrpgcrckhtnsvickaasgdakt
     srpdkpf'

For the full range of functions available, try `help(manager)` or help(geneeus.Proteome) 

#### List of access modules
    Proteome # for accessing protein data
    Genome # for accessing gene information

#### List of Proteome functions
    get_protein_name(ProteinID)
    get_protein_sequence(ProteinID)
    get_raw_xml(ProteinID)
    get_variants(ProteinID)
    get_geneID(ProteinID)
    get_protein_sequence_length(ProteinID)
    get_ID_type(ProteinID)
    run_translation(accession number)
    batch_get_protein_sequence([List of IDs])
    batch_get_protein_name([List of IDs])
    batch_get_variants([List of IDs])
    purge()
    get_size_of_datastore()

#### Genome functions
    Currently untested, so best to ignore for the time being...

### Design Decisions

#### Caching 
There were a number of design decisions which were made during the projects development, and no doubt will continue to be made. The manager object builds up a local data structure, and by default caches requests it makes to the database. The upshot of this from the user's perspective is that if I run

    manager.get_protein_sequence("NP_001120968")
    manager.get_protein_sequence("NP_001120968")

The second call doesn't query the database, but just reads off the cached value. This caching behavior can be turned off on by setting `cache=False` when initializing the ProteinManager object.

#### Batch queries 
A key design decision was how to deal with batch queries.

The eUtils recommended approach for making large (100+ IDs) queries is to initially ePost a list of those queries. The ePost operation sends this list to the ENTREZ servers, returning a `WebEnv` value and a `QueryKey` value. These two can then be used with an `eFetch` to go to the sever and get the result of the list submitted previously. The difficulty is that this list *must* be made up of UIDs (unique identifiers) which for proteins means GI numbers. If you only have an accession value (as is common) the only way to get this GI number is to query the database, and this `eSearch` operation can *not* be done in batch. Essentially, this is a chicken and egg problem - to get the GIs we need to do an `ePost` based batch query we have to run serial database queries.

This means that using `ePost`/`eFetch` would be great if you had a list of UIDs, but in practice accession numbers are a lot more common and useful, and the mapping of *n* accession values to UIDs would require *n* calls the server anyway. 

To get around this, I use concatenated `eFetch` calls for batch queries, whereby a single call is submitted with a list of IDs. This is a fast and stable way to get around this problem, and, so far, as shown no issues with lists up to 100 IDs long. The potential issue is that the HTTP GET request being made here literally gets longer as we add accession values, so this represents a top limit in terms of networking protocols. However, I implemented a recursive cascading retry mechanism which halves the list and retries each half, so should a list be too long it should only result in two calls instead of one.

#### Robustness
A primary goal with Geeneus was to create an API which is robust to input. By this, it should be able to handle case insensitivity, convert accession values where necessary, and correctly recognize valid accession values while rejecting irrelevant ones to minimize server burden. For accession filtering, we use regular expressions to ensure the only accession values which we query could be real values (based on NCBI's [accession rules](http://www.ncbi.nlm.nih.gov/Sequin/acc.html). While PDBs don't fall into this category, we allow translation between PDBID and GI, although often a chain identifier is required as the NCBI protein database treats separate chains as separate records.

All the networking is dealt with in a highly modular fashion, and network failure tolerance is a priority.

### Background and licence
This code was developed by [Alex Holehouse](http://holehouse.org) at [Washington University in Saint Louis](http://www.wustl.edu/) as part of the [Naegle lab](http://naegle.wustl.edu/people/lab_members.html). It is licensed under the the GNU General Public License (GPL-2.0). For more information see LICENCE.
