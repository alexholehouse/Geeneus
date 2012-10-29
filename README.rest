Geeneus
=========

Overview
----------
Geeneus is a very simple API for accessing data from the NCBI databases.
It currently supports primarily protein records (i.e. queries to the "protein" database) but has some early stage gene record functionality too. We plan to progressively add functionality

The general idea is that once a manager object has been created, all the complexity associated with parsing and networking should be totally abstracted from the user. To achieve this goal, the manager object has a set of functions which can be called and will return relevant data in Python basic data types only. For example, sequences are strings, sets are in dictionaries and lists of dictionaries, etc.

Installation
-------------
By far, the easiest way is to use pip to directly install. Running::

    sudo pip install Geeneus


Will install Geeneus with the biopython dependency.

Geeneus can also be installed from source using either `easy_install`, `pip` or just the`setup.py` script. See the github page for details on this.

------------

Usage
----------

Once installed, general usage is as follows::

    #!/usr/bin/env python

    from geeneus import Proteome
    manager = Proteome.ProteinManager("your.emailaddress@email.com")
    manager.get_protein_name("accession number")

NCBI requires that requests are accompanied by an email address. Geeneus is built such that it is impossible to exceed the eUtils guideline queries-per-minute rate if you use a single `manager` object. However, there are other guideliness which may be broken. As a user, you are recommended to review `NCBI's usage guidelines
<http://www.ncbi.nlm.nih.gov/books/NBK25497/>`_ and ensure you act accordingly.

For more information regarding the possible functions try::

    help(geeneus)
    help(geeneus.Proteome)

------------

More information
----------------

`For information on this project, including underlying design principles just click here 
<http://rednaxela.github.com/Geeneus>`_.

Requirements
-------------
Geeneus requires `biopython <http://biopython.org/DIST/docs/install/Installation.html>`_. Initially we're assuming biopython version 1.6, although earlier versions haven't been tested. To put it another way, we've tested on 1.6 and all's well. Earlier versions may also work, but you're on your own.


About
----------
This code was developed by `Alex Holehouse
<http://holehouse.org/>`_ at `Washington University in Saint Louis
<http://www.wustl.edu/>`_ as part of the `Naegle lab
<http://naegle.wustl.edu/>`_. It is licensed under the the GNU General Public License (GPL-2.0). For more information see LICENCE.
