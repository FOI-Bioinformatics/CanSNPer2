import logging
import os
import sys
from datetime import datetime

'''
One place for the log file naming of the commands (issue #217). Every verb passes its own
tag (cansnper2-type, cansnper2-db-info, cansnper2-scheme-get, ...), so two verbs never append
to each other's log file, and a name that is already taken steps through -HHMM and then
.2, .3 ... so two runs inside the same minute get separate files.
'''


def log_tag(prog):
	'''the log file tag of a command: its name with the verbs, cansnper2 scheme get gives
		cansnper2-scheme-get'''
	return "-".join(prog.split())


def get_logpath(logdir, tag, now=None):
	'''First log file name for this command that is not already taken'''
	now = now or datetime.now()
	logdir = logdir.rstrip("/")+"/"
	if not os.path.exists(logdir):
		os.makedirs(logdir, exist_ok=True)
	base = "{logdir}{tag}-{date}".format(logdir=logdir, tag=tag, date=now.strftime("%b-%d-%Y"))
	logpath = base+".log"
	if not os.path.exists(logpath):
		return logpath
	## no colon in the time: it is not allowed in file names on Windows
	stem = "{base}-{hour:02d}{minute:02d}".format(base=base, hour=now.hour, minute=now.minute)
	logpath = stem+".log"
	part = 1
	while os.path.exists(logpath):
		part += 1
		logpath = "{stem}.{part}.log".format(stem=stem, part=part)
	return logpath


class QuietWarnings(logging.Filter):
	'''A warning can ask not to be repeated on the console with
		extra={"quiet_console": True}: it always goes into the log file, the filter
		holds it back from the terminal and counts it, so the command can end with
		one line saying that there were warnings and where to read them. Run with
		--verbose or --debug and the warnings show up as they always did; when
		there is no log file, nothing may be hidden from the user.'''
	def __init__(self, show_quiet):
		super().__init__()
		self.show_quiet = show_quiet
		self.hidden = 0

	def filter(self, record):
		if getattr(record, "quiet_console", False) and not self.show_quiet:
			self.hidden += 1
			return False
		return True


console_filter = None


def quiet_warning_count():
	'''the warnings the console held back since the last setup_logging'''
	return console_filter.hidden if console_filter is not None else 0


def setup_logging(logdir, tag, level=logging.WARNING, now=None):
	'''Log to the console and to the log file of this command, and return that path.
		When the log folder cannot be made or written (a read-only working directory), one
		warning is written to stderr, the command logs to stderr only and None is returned.'''
	global console_filter
	console = logging.StreamHandler()
	handlers = []
	logpath = None
	try:
		logpath = get_logpath(logdir, tag, now=now)
		handlers.append(logging.FileHandler(logpath))
	except OSError as e:
		logpath = None
		print("Warning: the log folder {logdir} cannot be used ({error}), logging to stderr only".format(
			logdir=logdir, error=e.strerror or e), file=sys.stderr)
	console_filter = QuietWarnings(show_quiet=(level <= logging.INFO or logpath is None))
	console.addFilter(console_filter)
	handlers.append(console)
	logging.basicConfig(
			level=level,
			format="%(asctime)s %(module)s [%(levelname)-5.5s]  %(message)s",
			handlers=handlers)
	return logpath
