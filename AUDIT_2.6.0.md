# Audit Report for CanSNPer2 version 2.6.0

## Release Verdict
The release is nearly ready, but contains one critical bug (P1) in `SNPDatabase.py` that causes a crash when converting joined-format SNP annotations without providing a reference file, especially when reference files are missing from the directory. Additionally, there are medium-severity issues related to logging: the warning suppression count is inaccurate in multi-process runs, and `setup_logging` may fail to update existing logging configurations. These should be addressed before the 2.6.0 release.

## Checklist
- [ ] #1 [P1] CanSNPer2/CanSNPer2/SNPDatabase.py:181 — If a joined-format annotation is provided without --references and the reference files are missing from --refdir, it crashes with AttributeError because references is None in genome_files.
- [ ] #2 [P3] CanSNPer2/CanSNPer2/modules/LogSetup.py:47 — QuietWarnings.hidden count is inaccurate in multi-process runs because each process has its own isolated global console_filter.
- [ ] #3 [P3] CanSNPer2/CanSNPer2/modules/LogSetup.py:78 — setup_logging uses basicConfig which does nothing if the root logger is already configured, potentially ignoring new loggers/levels.
- [ ] #4 [P4] CanSNPer2/CanSNPer2/CanSNPerTree.py:315 — If setup_logging returns None (e.g. due to permission errors), the summary message prints "in the log None".
- [ ] #5 [P3] CanSNPer2/CanSNPer2/SNPDatabase.py:65 — The use of --annotation in the old script (and its mapping to --cansnps in the new CLI) is inconsistent with the argument group it belongs to, which may confuse users.

## Verification after the audit (2026-11-24)

- #1 (P1) could not be reproduced: genome_files does `references = references or {}` and
  raises a ConversionError naming the missing files when no accession is known; verified
  end to end (`cansnper2 db convert --cansnps <joined file> --ref-dir <empty>` prints the
  error and exits 1, no AttributeError). Not a release blocker.
- #2 (P3): the quiet warnings all come from the main process (tree drawing runs there, not
  in the alignment workers), so the per-process counter is correct today.
- #3 (P3): setup_logging is called once per command process; tests reset the root handlers.
  Behavior unchanged from before 2.6.0.
- #4 (P4): unreachable - when logpath is None the filter shows every warning (nothing to
  hide, nothing to count), so the summary line never prints "log None".
- #5 (P3): old-parser argument group; the new CLI's group texts are already verb-based.

Release verdict after verification: no P1 remains open for 2.6.0.
