#!/bin/bash
 ./build_deployment.sh

cd ../deployment/
python setup.py sdist
cd dist/
unpack Geeneus-0.1.1.tar.gz
cd Geeneus-0.1.1/
sudo python setup.py install