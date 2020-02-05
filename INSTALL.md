##Installing CanSNPer2
CanSNPer2 is written in Python and requires Python installed as well as several 
dependencies listed below. The software is developed and tested with Python 
version 3.7 but earlier (3.5, 3.6) may work just as well. Due to other 
dependencies other versions of Python are not supported.

### Bioconda
With an activated Bioconda channel (see [Set up channels](https://bioconda.github.io/index.html#set-up-channels)), 
install with:
`conda install cansnper2`
and update with:
`conda update cansnper2`

This will install CanSNPer2 with all Python and other dependencies (except progressivemauve).

Remark that Bioconda only supports MacOS and Linux systems. 

### Python Package Index
An alternative way to install CanSNPer2 and Python dependencies is using PyPI with: 
`pip install cansnper2`

Non-Python dependencies needs to be installed manuelly. See dependency secion below.

### Latest source code
To install the latest source code, start by cloning this repo, change into the CanSNPer2 directory and
run the setup file. This will install CanSNPer2 and the Python dependecies.

```
git clone https://github.com//CanSNPer2.git
cd CanSNPer2
python setup.py
```

Non-Python dependencies needs to be installed manuelly. See dependency secion below.

## Testing the installation
When the dependencies are all installed CanSNPer2 can be run as it is from the 
shell.

To begin with, try:
```
CanSNPer2 --help
```

If there are no errors and you see the help text printed, CanSNPer2 is working 
correctly and you can go on to do your analysis.

If all dependencies are installed you can also run CanSNPer2 without running the
`setup.py`. This is done by running `python CanSNPer2 -h` from within the CanSNPer2
directory.

More on how to run CanSNPer2 in its various modes can be found in the README 
file that came with the distribution.

CanSNPer2 has been tested on Linux and MacOS (OSX), but feel free to tinker with setting up 
a working version on your favorite OS. Mauve is available on other platforms.

##Dependencies
Software that must be installed before running CanSNPer2:

[Python (3.7.X)](http://www.python.org/getit/)  
Most Linux distributions come with Python installed. However we recommend using 
[Conda](http://conda.pydata.org/docs/) to setup a flexible Python environment. 

[ETE3](http://ete.cgenomics.org/)  
ETE3 has a number of additional dependencies, listed in their install 
notes. Most notably, there are several dependencies that are not needed 
for CanSNPer2, but they may raise warnings as ETE3 is loaded. Note that
Qt < 5 is needed to run ETE3.

[NumPy](http://www.numpy.org/)  
Simple install instructions are available for this package.

[progressiveMauve](http://darlinglab.org/mauve/mauve.html)  
The progressiveMauve binary must be in the PATH or specifically set in 
the CanSNPer2.conf file.

