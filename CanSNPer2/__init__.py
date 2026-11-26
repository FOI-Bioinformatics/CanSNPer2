__version__ = "2.6.0"

COMMANDS_HELP = """
CanSNPer2 since 2.6.0 is one command with a verb for every job:
  cansnper2 type      type assemblies or reads against a scheme database
  cansnper2 selftest  check the installation: imports and aligner programs
  cansnper2 db        create, import, export, info, convert, modify, remove
  cansnper2 scheme    list schemes of CanSNPer2-data, get builds their database
  cansnper2 ref       get the reference genomes a database needs from NCBI
The commands CanSNPer2, CanSNPer2-database, CanSNPer2-download and
CanSNPer2-test are no longer installed in 2.6.0; the verbs above do the
same jobs, with hyphenated spellings of the options.
Run cansnper2 <command> --help for its arguments; examples in README.md.
"""
