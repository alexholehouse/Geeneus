#!/bin/bash

echo "Delete existing deployment\n(note yes will delete EVERYTHING, not just\nGenome.py"
rm -ri deployment/Genome.py

rm -r deployment/

mkdir deployment/
mkdir deployment/Backend/

cp Geeneus/Genome.py deployment/
cp Geeneus/Proteome.py deployment/
cp Geeneus/__init__.py deployment/
cp Geeneus/LICENSE deployment/
cp Geeneus/Backend/*.py deployment/Backend/
