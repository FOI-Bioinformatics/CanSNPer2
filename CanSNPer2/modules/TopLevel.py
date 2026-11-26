"""
The handler of unexpected errors at the top level of every command.

A known problem with an input (a missing file, a bad scheme or database) ends in one message
and exit status 1 where it is found. Any other exception used to end the command with a
Python traceback. run_main() turns it into one line on stderr,

    CanSNPer2 stopped: <type>: <message>; run with --debug for details

and exit status 2. With --debug the traceback is printed as well. KeyboardInterrupt and
SystemExit are not exceptions of this kind and pass unchanged. The traceback is also
written to the log file of the command, when it has one.
"""

import logging
import sys
import traceback

## the exit status of a command stopped by an unexpected error; known errors exit with 1
UNEXPECTED_EXIT = 2


def log_to_files(message):
    """write message with the current traceback to the log files only: the console already
    has the message from run_main"""
    record = logging.LogRecord("CanSNPer2", logging.ERROR, __file__, 0, message, None, sys.exc_info())
    for handler in logging.getLogger().handlers:
        if isinstance(handler, logging.FileHandler):
            try:
                handler.handle(record)
            except Exception:
                pass


def run_main(main, argv=None):
    """Run main() and return what it returns. An unexpected exception becomes one line on
    stderr and exit status UNEXPECTED_EXIT; argv (default sys.argv[1:]) is searched for --debug"""
    try:
        return main()
    except Exception as e:
        argv = sys.argv[1:] if argv is None else argv
        message = "CanSNPer2 stopped: {type}: {message}".format(type=type(e).__name__, message=e)
        log_to_files(message)
        if "--debug" in argv:
            traceback.print_exc()
            print(message, file=sys.stderr)
        else:
            print(message + "; run with --debug for details", file=sys.stderr)
        sys.exit(UNEXPECTED_EXIT)
