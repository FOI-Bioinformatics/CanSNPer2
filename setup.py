from setuptools import setup, find_packages
from CanSNPer2.CanSNPerTree import __version__

setup(
    name="CanSNPer2",
    version=__version__,
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

    install_requires=['ete3','flextaxd'],
    packages=find_packages(exclude=['contrib', 'docs', 'test*']),
    #py_modules=['CanSNPer2.modules.ParseXMFA',"CanSNPer2.modules.DatabaseConnection","CanSNPer2.modules.CanSNPer2","CanSNPer2.modules.NewickTree"],
    entry_points={'console_scripts': [
                    'CanSNPer2=CanSNPer2.CanSNPerTree:main',
                    'CanSNPer2-database=CanSNPer2.SNPDatabase:main',
                    'CanSNPer2-download=CanSNPer2.DownloadGenomes:main',
					'CanSNPer2-test=CanSNPer2.selftest:main'
    ]})
