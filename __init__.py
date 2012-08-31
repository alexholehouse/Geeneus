# Copyright 2012 by Alex Holehouse.  All rights reserved. 
# Contact at alex.holehouse@wustl.edu
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

""" Acts as a easy-to-use wrapper around the Entrez package, for an intuitive
and simple way to get sequence information from NCBI. Abstracts NCBI's
eUtils from the user, providing a uniform interface. The absolute goal here is 
not to provide a customizable or complex system, but to provide a limited but
very easy to use small set of tools for getting the necessary data.

Variables:
email      Set the Entrez email parameter (default is none)
tool       Set the Entrez tool parameter (default is biopython)

Internal Modules - call using Geeneus.Genome or Geeneus.Proteome
genome     Functions for dealing with genome data (e.g. getting sequences,
           coding regions)
proteome    Functions for dealing with protein data (e.g. protein sequence)

"""

import Proteome
import Genome

tool = "biopython"
email = None

