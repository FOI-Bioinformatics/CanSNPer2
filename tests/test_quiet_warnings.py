"""Warnings the terminal holds back, and the summary line that replaces them.

One warning per sample about SNP positions without alignment data made the
results overview look like an error screen: the warning belongs in the log,
the terminal gets one line with the count, and --verbose shows the warnings
themselves again."""
import logging

from CanSNPer2.modules import LogSetup


def reset_logging():
	for handler in list(logging.getLogger().handlers):
		logging.getLogger().removeHandler(handler)
		handler.close()


def test_quiet_warnings_stay_in_the_log(tmp_path, capsys):
	reset_logging()
	logpath = LogSetup.setup_logging(str(tmp_path), "quiet-test")
	logging.getLogger("quiet").warning("hidden from the window", extra={"quiet_console": True})
	logging.getLogger("quiet").warning("loud to the window")
	read = capsys.readouterr()
	captured = read.out + read.err
	assert "hidden from the window" not in captured
	assert "loud to the window" in captured
	assert LogSetup.quiet_warning_count() == 1
	with open(logpath) as f:
		logged = f.read()
	assert "hidden from the window" in logged and "loud to the window" in logged
	reset_logging()


def test_verbose_runs_show_the_quiet_warnings(tmp_path, capsys):
	reset_logging()
	LogSetup.setup_logging(str(tmp_path), "verbose-test", level=logging.INFO)
	logging.getLogger("quiet").warning("shown when verbose", extra={"quiet_console": True})
	read = capsys.readouterr()
	captured = read.out + read.err
	assert "shown when verbose" in captured
	assert LogSetup.quiet_warning_count() == 0
	reset_logging()


def test_quiet_warnings_are_not_hidden_without_a_log_file(tmp_path, capsys):
	'''a log file that cannot be written is reported, and then every warning has to
		reach the terminal because there is no log left to point the user at'''
	reset_logging()
	blocker = tmp_path / "logs"
	blocker.write_text("a file cannot be a log folder")
	logpath = LogSetup.setup_logging(str(blocker), "no-log-test")
	assert logpath is None
	logging.getLogger("quiet").warning("no log to hide in", extra={"quiet_console": True})
	read = capsys.readouterr()
	captured = read.out + read.err
	assert "no log to hide in" in captured
	assert "cannot be used" in captured
	assert LogSetup.quiet_warning_count() == 0
	reset_logging()
