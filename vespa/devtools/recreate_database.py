#!/usr/bin/env python

"""A semi-official script that deletes the existing database and rebuilds it
from scripts."""

# Python modules

import os
import sys
import logging
import optparse
import subprocess
import tempfile

# 3rd party modules

# Our modules
import vespa.common.create_database as create_database
import vespa.common.util.logging_ as util_logging
import vespa.common.util.import_ as util_import
import vespa.common.util.config as util_config
import vespa.common.util.misc as util_misc
import vespa.common.util.db as util_db

INSTRUCTIONS = \
"""This semi-official script (re)creates Vespa's database, setting it back to
"factory defaults". If there's already an existing Vespa database, it is
destroyed in the process without any chance of recovery.

There are two possible arguments; both are optional.

If you specify -r or --restore, this script will attempt to restore all of
the custom objects in your database. If the database or export format has
changed, this might fail.

If you specify a filename, this script will write the database to that
filename instead of the default Vespa database file.

"""

msg = \
"""This tool is semi-deprecated until further notice as it was created for 
use in Vespa versions < 0.5.0 and did not include functionality for dealing
with either PulseProjecs (version < 0.6.0) or PulseDesigns (version > 0.6.0)

Returning without executing ... Sorry!
"""

print(msg)
return

# bjs - 2021-03-07, temp fix to keep people from using this

parser = optparse.OptionParser(usage=INSTRUCTIONS)
parser.add_option("-r", "--restore", action="store_true", dest="restore")

(options, args) = parser.parse_args()

restore = options.restore

if not args:
    # User didn't supply a filename. This is normal. Instead I use the standard
    # database path & name.
    data_dir = util_misc.get_data_dir()
    # The database name is in vespa.ini
    config = util_config.VespaConfig()
    target_filename = os.path.join(data_dir, config["database"]["filename"])
elif len(args) == 1:
    # User supplied a filename so I write the database there.
    target_filename = args[0]
else:
    # len(args) > 1
    print(INSTRUCTIONS)

    print("Wrong number of arguments.")

    sys.exit(-1)

# restore_files tracks the temp files to which a given type has been exported
restore_files = {}

if options.restore:
    restorable_objects = ("metabolites", "pulse_sequences", "experiments")

    # export all object types to temp files
    for object_type in restorable_objects:
        print("Attempting to export %s..." % object_type)
        fd, filename = tempfile.mkstemp()

        # Save this filename so I can import from it (and clean it up) later
        restore_files[object_type] = filename

        command = ["python", "export.py", object_type]
        subprocess.call(command, stdout=fd)


logger = logging.getLogger(util_logging.Log.DATABASE)

# Set up a handler that just prints to the screen.
stream_handler = logging.StreamHandler()

logger.addHandler(stream_handler)

create_database.create_database(target_filename, logger)

if options.restore:
    db = util_db.Database(target_filename)

    for object_type in restore_files:
        print("Attempting to restore %s..." % object_type)

        import_filename = restore_files[object_type]

        # This code is kind of ugly. Isn't there a better way to do this
        # than hardcoding a bunch of names and types?
        if object_type == "metabolites":
            importer = util_import.MetaboliteImporter(import_filename, db)
        elif object_type == "experiments":
            importer = util_import.ExperimentImporter(import_filename, db)
        if object_type == "pulse_sequences":
            importer = util_import.PulseSequenceImporter(import_filename, db)

        importer.go()

        # Clean up
        os.remove(import_filename)

    print("\nDone restoring!")


