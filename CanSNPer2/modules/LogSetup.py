import logging
import os
from datetime import datetime

'''
One place for the log file naming that the three commands used to copy paste (issue #217).
Every command passes its own tag, so CanSNPer2-database and CanSNPer2-download can no
longer append to each other's log file, and a name that is already taken steps through
-HH:MM and then .2, .3 ... so two runs inside the same minute get separate files.
'''


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
	stem = "{base}-{hour:02d}:{minute:02d}".format(base=base, hour=now.hour, minute=now.minute)
	logpath = stem+".log"
	part = 1
	while os.path.exists(logpath):
		part += 1
		logpath = "{stem}.{part}.log".format(stem=stem, part=part)
	return logpath


def setup_logging(logdir, tag, level=logging.WARNING, now=None):
	'''Log to the console and to the log file of this command, and return that path'''
	logpath = get_logpath(logdir, tag, now=now)
	logging.basicConfig(
			level=level,
			format="%(asctime)s %(module)s [%(levelname)-5.5s]  %(message)s",
		    handlers=[
		        logging.FileHandler(logpath),
		        logging.StreamHandler()
		    ])
	return logpath
