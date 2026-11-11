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

    # setuptools >=77 wants the license as an SPDX expression and deprecates the
    # "License ::" classifier; the file in the package folder stays the license text
    license='GPL-3.0-or-later',

    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
        'Programming Language :: Python :: 3.14',
        "Operating System :: OS Independent"
    ],
    # the ete3 3.1.3 import of cgi (removed in python 3.13) capped CanSNPer2 below 3.13 for
    # years (issue #184); the package no longer imports ete3/ete4 at all, it writes a
    # self-contained html tree viewer instead, and runs on 3.14
    python_requires='>=3.11',

    keywords="Bioinformatics SNP-typing sequence-data",

    # no third party runtime imports any more: ete3, PyQt5, six and numpy were all tree
    # drawing dependencies (issues #180/#182: "pip install -e CanSNPer2" had to pull ete3 and
    # its imports in, and environments whose ete3 lacked them died at import time). With the
    # html viewer in modules/NewickTree.py the whole class of failures is gone; pytest is
    # development only and stays in requirements.txt instead
    install_requires=[],
    packages=find_packages(exclude=['contrib', 'docs', 'test*']),
    #py_modules=['CanSNPer2.modules.ParseXMFA',"CanSNPer2.modules.DatabaseConnection","CanSNPer2.modules.CanSNPer2","CanSNPer2.modules.NewickTree"],
    entry_points={'console_scripts': [
                    'CanSNPer2=CanSNPer2.CanSNPerTree:main',
                    'CanSNPer2-database=CanSNPer2.SNPDatabase:main',
                    'CanSNPer2-download=CanSNPer2.DownloadGenomes:main',
					'CanSNPer2-test=CanSNPer2.selftest:main'
    ]})
