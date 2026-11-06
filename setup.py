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
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        "Operating System :: OS Independent"
    ],
    # ete3 3.1.3 imports cgi (ete3/webplugin/webapp.py), which was removed in Python 3.13,
    # so importing ete3 - and therefore CanSNPer2 - is impossible on 3.13 and newer; the
    # upper bound makes pip refuse such an environment with a clear Requires-Python error
    # instead of failing at import time (issue #184). The pinned numpy==1.24.4 in
    # requirements.txt only has wheels up to 3.11, so a built conda package should target
    # 3.8-3.11 and needs numpy>=1.26 before 3.12 can be supported.
    python_requires='>=3.8,<3.13',

    keywords="Bioinformatics SNP-typing sequence-data",

    # ete3 is the only third party module the package imports directly, but ete3 itself
    # imports six and numpy at import time, and "pip install -e CanSNPer2" does not read
    # requirements.txt. Environments whose ete3 was installed without its dependencies
    # therefore died with "No module named 'six'" and then "No module named 'numpy'" on
    # the first CanSNPer2 call (issue #180, #182)
    install_requires=['ete3','PyQt5','six','numpy'],
    packages=find_packages(exclude=['contrib', 'docs', 'test*']),
    #py_modules=['CanSNPer2.modules.ParseXMFA',"CanSNPer2.modules.DatabaseConnection","CanSNPer2.modules.CanSNPer2","CanSNPer2.modules.NewickTree"],
    entry_points={'console_scripts': [
                    'CanSNPer2=CanSNPer2.CanSNPerTree:main',
                    'CanSNPer2-database=CanSNPer2.SNPDatabase:main',
                    'CanSNPer2-download=CanSNPer2.DownloadGenomes:main',
					'CanSNPer2-test=CanSNPer2.selftest:main'
    ]})
