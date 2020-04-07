from setuptools import setup, find_packages
from CanSNPer2.CanSNPerTree import __version__

setup(
    name="CanSNPer2",
    version=__version__,
    url="https://git-int.foi.se/bioinfo/cansnper2",
    description="CanSNPer2: A toolkit for SNP-typing bacterial genomes.",
    license="GPLv3'",

    # Author details
    author='David Sundell',
    author_email='david.sundell@foi.se',

    keywords="Bioinformatics SNP-typing sequence-data",
    classifiers=[
        'Development Status :: 5 - Beta',
        'License :: OSI Approved :: GPL',
        'Programming Language :: Python :: 3'
        ],
    install_requires=['ete3','flextaxd'],
    packages=find_packages(exclude=['contrib', 'docs', 'test*']),
    #py_modules=['CanSNPer2.modules.ParseXMFA',"CanSNPer2.modules.DatabaseConnection","CanSNPer2.modules.CanSNPer2","CanSNPer2.modules.NewickTree"],
    entry_points={'console_scripts': [  'CanSNPer2=CanSNPer2.CanSNPerTree:main',
                                        'CanSNPer2-database=CanSNPer2.SNPDatabase:main',
                                        'CanSNPer2-download=CanSNPer2.DownloadGenomes:main'
                                        ]})
