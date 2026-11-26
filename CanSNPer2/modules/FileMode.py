"""The permission bits of files that are written through tempfile.mkstemp.

mkstemp creates its file with mode 0600, so an output that is written to a temporary
file and then renamed into place could only be read by its owner. A converted database
or a reference folder is often shared within a group, so such files get the mode that
open() gives a new file: 0666 masked by the umask of the process."""

import os


def default_file_mode():
    """the mode open() gives a new file: 0666 masked by the current umask. The umask can only
    be read by setting it, so it is set and restored again at once; the value set meanwhile
    is the most restrictive one, so a file made by another thread in that moment is not
    opened up"""
    umask = os.umask(0o077)
    os.umask(umask)
    return 0o666 & ~umask


def set_default_mode(path):
    """give path the mode of a file made by open() instead of the 0600 of mkstemp"""
    os.chmod(path, default_file_mode())
