from setuptools import setup, find_packages
import os, re

def get_version():
    '''setup.py must never import the package it builds: pip runs the build in an isolated
        environment that does not contain the dependencies, so importing CanSNPerTree died with
        "ModuleNotFoundError: No module named ete3" and neither pip install . nor pip install -e .
        could get past reading the metadata
    '''
    init = os.path.join(os.path.dirname(os.path.abspath(__file__)), "CanSNPer2", "__init__.py")
    with open(init) as version_file:
        found = re.search(r"__version__\s*=\s*['\"]([^'\"]+)['\"]", version_file.read())
    if not found:
        raise RuntimeError("No __version__ was found in CanSNPer2/__init__.py")
    return found.group(1)

setup(
    name="CanSNPer2",
    version=get_version(),
    url="https://git-int.foi.se/bioinfo/cansnper2",
    description="CanSNPer2: A toolkit for SNP-typing bacterial genomes.",

    # Author details
    author='David Sundell',
    author_email='david.sundell@foi.se',

    license=' GNU GENERAL PUBLIC LICENSE version 3',

    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        "Operating System :: OS Independent"
    ],
    python_requires='>=3.6',

    keywords="Bioinformatics SNP-typing sequence-data",

    install_requires=['ete3','PyQt5'],
    packages=find_packages(exclude=['contrib', 'docs', 'test*']),
    #py_modules=['CanSNPer2.modules.ParseXMFA',"CanSNPer2.modules.DatabaseConnection","CanSNPer2.modules.CanSNPer2","CanSNPer2.modules.NewickTree"],
    entry_points={'console_scripts': [
                    'CanSNPer2=CanSNPer2.CanSNPerTree:main',
                    'CanSNPer2-database=CanSNPer2.SNPDatabase:main',
                    'CanSNPer2-download=CanSNPer2.DownloadGenomes:main',
					'CanSNPer2-test=CanSNPer2.selftest:main'
    ]})
