Geeneus
=========

Overview
----------
Geeneus is a very simple API for accessing data from the NCBI databases.
It current has extensive support for protein information and limited support for gene information.

The general idea is that once a manager object has been created, all the complexity associated with parsing and networking should be totally abstracted from the user. To achieve this goal, the manager object has a set of functions which can be called and will return relevant data in Python basic data types only. For example, sequences are strings, sets are in dictionaries and lists of dictionaries, etc.

-------------

Installation
-------------

By far the easiest way to install Geeneus is to use pip to directly install. Running

    sudo pip install Geeneus

Will install Geeneus with the biopython dependency.

Geeneus can also be installed from source using either `easy_install`, `pip` or just the`setup.py` script.

#### PIP and easy_install from source
Both of these will also install the Biopython dependency if needed.

    sudo pip install Geeneus-0.1.0.tar.gz
    # or you can use
    sudo easy_install Geeneus-0.1.0.tar.gz


#### setup.py from source
Note that sdist setup tools do not support the `install_requires` option in the `setup.py` file, so if you install directly using `setup.py` this will generate a warning message like

    /usr/lib/python2.7/distutils/dist.py:267: UserWarning: Unknown distribution option: 'install_requires' 
    warnings.warn(msg)

If you already have biopython installed this shouldn't be a problem. If you do not please install biopython first!

    tar xvf Geeneus-0.1.0.tar.gz 
    cd Geeneus-0.1.0
    sudo python setup.py install

------------


Usage
------------

Once install, general usage is as follows;

    #!/usr/bin/env python

    from geeneus import Proteome
    manager = Proteome.ProteinManager("your.emailaddress@email.com")
    manager.get_protein_name("accession number")

For more information regarding the possible functions try

    help(geeneus)
    help(geeneus.Proteome)

------------

More information
----------
[For information on this project, including underlying design principles just click here](http://rednaxela.github.com/Geeneus)


Requirements
----------
geeneus [requires biopython](http://biopython.org/DIST/docs/install/Installation.html). Initially we're assuming 1.6, although earlier versions haven't been tested. To put it another way, we've tested on 1.6 and all's well. Earlier versions may also work, but you're on your own.


About
----------
This code was developed by [Alex Holehouse](http://holehouse.org) at [Washington University in Saint Louis](http://www.wustl.edu/) as part of the [Naegle lab](http://naegle.wustl.edu/people/lab_members.html). It is licenced under the the GNU General Public License (GPL-2.0). For more information see LICENCE.

