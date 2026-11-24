__version__ = "2.6.0"

COMMANDS_HELP = """
available CanSNPer2 commands (installed together):
  cansnper2            one command with a verb for every job: cansnper2 --help
    type               type query assemblies against a CanSNP database
    selftest           check the installation: imports and aligner programs
    db                 create, import, export, info, convert, modify, remove
    scheme             list schemes of CanSNPer2-data, get builds their database
    ref                get the reference genomes a database needs from NCBI
  CanSNPer2, CanSNPer2-database, CanSNPer2-download and CanSNPer2-test do the
  same jobs with their old arguments and keep working; cansnper2 is the
  interface to use from now on.
Run any command with --help for its arguments; examples in README.md.
"""

