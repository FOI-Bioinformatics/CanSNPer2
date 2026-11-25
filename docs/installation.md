# Installation

This page describes how to install CanSNPer2 and the external programs it runs, and how to check
the installation. It is written for anyone who sets up CanSNPer2 on a workstation, a server or a
compute cluster.

CanSNPer2 is a Python package without third-party Python dependencies. It needs Python 3.11 or
later and, depending on the task, the external programs in the table below.

| Program | Bioconda package | Needed for |
| --- | --- | --- |
| minimap2 | `minimap2` | typing assemblies with the default aligner (`--aligner minimap2`), and typing sequencing reads |
| samtools (1.10 or later) | `samtools` | typing sequencing reads |
| progressiveMauve | `progressivemauve` | typing assemblies with `--aligner progressiveMauve` |
| MUMmer4 (`nucmer` and `delta-filter`) | `mummer4` | typing assemblies with `--aligner mummer4` |

At least one aligner is needed to type assemblies. Building databases, fetching schemes and
downloading reference genomes need no external program. How the aligners differ is described in
[method.md](method.md) and [validation.md](validation.md).

## Install with Bioconda

The commands below create a conda environment with CanSNPer2, minimap2 and samtools, which
covers typing of assemblies with the default aligner and typing of sequencing reads:

```sh
mamba create -n cansnper2 -c conda-forge -c bioconda cansnper2 minimap2 samtools
mamba activate cansnper2
```

`conda` can be used in place of `mamba` with the same arguments. The alternative aligners can be
added to the environment:

```sh
mamba install -n cansnper2 -c conda-forge -c bioconda mummer4
mamba install -n cansnper2 -c conda-forge -c bioconda progressivemauve
```

Check the installation with `cansnper2 --version`, which should print 2.6.0 or later, and with
`cansnper2 selftest` (see [Check the installation](#check-the-installation)).

Notes on the Bioconda package:

- The recipe is maintained in the
  [bioconda-recipes](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/cansnper2)
  repository, not in this repository. The external programs are named explicitly in the
  commands above, so that they are installed whatever the recipe lists.
- progressiveMauve is available from Bioconda for Linux (x86-64) only. On macOS, type
  assemblies with minimap2 (the default) or MUMmer4 (`--aligner mummer4`).

## Install from source

For development, or in environments without conda, CanSNPer2 can be installed with pip into
any environment with Python 3.11 or later (versions 3.11 to 3.14 are tested):

```sh
git clone https://github.com/FOI-Bioinformatics/CanSNPer2.git
cd CanSNPer2
python -m pip install .
```

pip installs the `cansnper2` command. The external programs are installed separately, for
example from Bioconda into the same conda environment:

```sh
mamba install -c conda-forge -c bioconda minimap2 samtools mummer4
```

Each program must be found in `PATH`. For assemblies, the aligner program can instead be given
with `cansnper2 type --aligner-path`, as the program itself or the folder that holds it.

To work on the code itself, see [developer.md](developer.md).

## Platforms

CanSNPer2 is tested on Linux and run on macOS. Native Windows is not supported: the downloads of
reference genomes and parts of `cansnper2 type` run in worker processes started with the `fork`
method, which Windows does not provide, and read typing runs the mapping commands
through `bash`.

## Check the installation

`cansnper2 selftest` checks that the Python modules of CanSNPer2 can be imported and reports
which external programs are found in `PATH`. In an environment with minimap2 and samtools, and
without progressiveMauve and MUMmer4, it prints (the path depends on the installation):

```text
minimap2: /home/user/miniforge3/envs/cansnper2/bin/minimap2 (--aligner minimap2, default)
nucmer: not found (--aligner mummer4)
delta-filter: not found (needed by nucmer)
progressiveMauve: not found (--aligner progressiveMauve)
usable aligners: minimap2
read typing: usable (minimap2 and samtools)
```

The lines mean:

| Line | Meaning |
| --- | --- |
| `<program>: <path> (--aligner <name>)` | the aligner program was found at that path and is selected with `--aligner <name>`; `default` marks the aligner used when `--aligner` is not given |
| `<program>: not found (...)` | the program is not in `PATH` |
| `delta-filter: ... (needed by nucmer)` | MUMmer4 is usable only when both `nucmer` and `delta-filter` are found |
| `usable aligners: ...` | the aligners that `cansnper2 type` can run for assemblies |
| `read typing: usable (minimap2 and samtools)` | sequencing reads can be typed |
| `read typing: not usable, <program> not found` | read typing is not possible until the named programs are installed |

The exit status is 0 when at least one aligner is usable. It is 1, with an `ERROR:` line, when
no aligner is found or when CanSNPer2 itself cannot be imported. A missing read typing program
does not change the exit status. `selftest` does not check the versions of the programs.

`cansnper2 --version` prints the installed version, and `cansnper2 --help` lists the commands.

## Upgrade an installation

With Bioconda, update the package in its environment:

```sh
mamba update -n cansnper2 -c conda-forge -c bioconda cansnper2
```

From source, update the clone (`git pull`) and run `python -m pip install .` again.

Version 2.6.0 replaced the earlier commands by the single command `cansnper2`, and databases
built by versions before 2.4.0 must be converted before they can be used. Both are described in
[upgrading.md](upgrading.md).
