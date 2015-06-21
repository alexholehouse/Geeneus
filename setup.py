from distutils.core import setup

setup(
    name='Geeneus',
    version='0.1.9',
    author='Alex Holehouse',
    author_email='alex.holehouse@gmail.com',
    packages=['geeneus', 'geeneus.test', 'geeneus.backend'],
    scripts=[],
    url='http://alexholehouse.github.com/Geeneus/',
    license='LICENSE.txt',
    description='Simple API for NCBI database access',
    long_description=open('README.txt').read(),
    install_requires=[
        "Biopython >= 1.60",
        "requests >= 0.14.2"   
        ],
    )
