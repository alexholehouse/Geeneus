# Copyright 2012 by Alex Holehouse - see LICENSE for more info
# Contact at alex.holehouse@wustl.edu

""" This is the Ontology Tools module which allows for access to disease ontology data. This acts as a stand alone section of Geeneus, so users should have
the option of including it or not. Uses sqlalchamy and automatically downloads a DO database for local access, so if either of these things are an issue then
be aware
"""
import DOParser

try:
    import sqlalchamy
except ImportError:
    print "\n########"
    print "ERROR: The OntologyTools module requires that you have sqlalchamy installed in location on your Python path."
    print "As a result, OntologyTools won't work. This will not effect any other part of Geeneus"
    print "########\n"






