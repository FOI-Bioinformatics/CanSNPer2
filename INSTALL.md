## Installing CanSNPer2
CanSNPer2 is written in Python and needs Python 3.11 or later (3.11 to 3.14 are
listed in setup.py). It uses the Python standard library only, so there are no
Python packages to install besides CanSNPer2 itself. An aligner program is
needed as well, see the dependencies section below.

### Bioconda
With an activated Bioconda channel (see [Set up channels](https://bioconda.github.io/index.html#set-up-channels)), 
install with:
`conda install -c bioconda -c conda-forge cansnper2`
and update with:
`conda update cansnper2`

This installs the single command `cansnper2` together with its default aligner
progressiveMauve; the optional aligners minimap2 and mummer4 can be installed
alongside it and selected with `--aligner`.

Remark that Bioconda only supports MacOS and Linux systems. 

### Python Package Index
CanSNPer2 is a pure-stdlib package, so pip from a checkout or a release archive is
all it needs (the bioconda package is the primary route):
`python -m pip install .        # or pipx install .`

The aligner programs need to be installed separately. See the dependencies section below.

### Latest source code
To install the latest source code, start by cloning this repo, change into the CanSNPer2 directory and
install it with pip. CanSNPer2 has no Python dependencies.

```
git clone https://github.com/FOI-Bioinformatics/CanSNPer2.git
cd CanSNPer2
pip install .
```

The aligner programs need to be installed separately. See the dependencies section below.

## Testing the installation
When the dependencies are all installed CanSNPer2 can be run as it is from the 
shell.

To begin with, try:
```
cansnper2 --help
```

If there are no errors and you see the help text printed, CanSNPer2 is working 
correctly and you can go on to do your analysis.

From within the CanSNPer2 directory, without any install step, the same works as
`python -m CanSNPer2.cli --help`.

More on how to run CanSNPer2 in its various modes can be found in the README 
file that came with the distribution.

CanSNPer2 has been tested on Linux and MacOS (OSX), but feel free to tinker with setting up 
a working version on your favorite OS. Mauve is available on other platforms.

## Dependencies
CanSNPer2 needs Python 3.11 or later and at least one aligner program. It has
no Python package dependencies (ETE3, PyQt5 and NumPy are not used; the tree
view is written as a self-contained html file).

CanSNPer2 aligns each query against the reference genomes with one of these
programs, chosen with `--aligner` (default progressiveMauve):

[progressiveMauve](http://darlinglab.org/mauve/mauve.html)  
The progressiveMauve binary must be in the PATH, or be given with
`--aligner-path` (the program itself or the folder that holds it).

[minimap2](https://github.com/lh3/minimap2)  
Used with `--aligner minimap2`.

[MUMmer4](https://github.com/mummer4/mummer)  
Used with `--aligner mummer4`; it needs both `nucmer` and `delta-filter`.

`cansnper2 selftest` reports which of these programs are found in PATH.
