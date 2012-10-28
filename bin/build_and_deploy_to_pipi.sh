#!/bin/bash
 ./build_deployment.sh

cd ../deployment/
python setup.py sdist
python setup.py sdist upload