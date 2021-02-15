#!/usr/bin/env python

"""A semi-official script that exports objects from the Vespa database into
Vespa's standard XML export format.
"""

# Python modules

import os
import sys

# 3rd party modules

# Our modules
import vespa.common.util.db as util_db
import vespa.common.util.misc as util_misc
import vespa.common.util.config as util_config
import vespa.common.util.export as util_export


COMMENT = """This is the master file for populating the database."""

INSTRUCTIONS = \
"""This script exports objects from the Vespa database into Vespa's standard 
XML export format. Output is to stdout; if you want to export to a file,
redirect the output as in the example below.

Typical invocation:
   python export.py  [object] > target_filename.xml

Where [object] is one of experiments, metabolites, or pulse_sequences.

When exporting experiments, their simulations are always included.

The comment included in the export file is hardcoded in this script. It is:
-----------------------------------------------------------------------------
%s
-----------------------------------------------------------------------------
""" % COMMENT


def export_experiments(db, filename):
    # There's no function for fetching all experiments, so I have to fetch
    # experiment previews and then each individual experiment.
    experiments = db.fetch_experiment_previews()
    experiments = [db.fetch_experiment(experiment.id) for experiment \
                                                        in experiments]
    util_export.export(filename, experiments, db, COMMENT)


def export_metabolites(db, filename):
    metabolites = db.fetch_metabolites()
    util_export.export(filename, metabolites, db, COMMENT)


def export_pulse_sequences(db, filename):
    pulse_sequences = db.fetch_pulse_sequences()
    util_export.export(filename, pulse_sequences, db, COMMENT)

# dispatcher maps each potential command line param to the appropriate 
# function
dispatcher = { "experiments"        : export_experiments, 
               "metabolites"        : export_metabolites,
               "pulse_sequences"    : export_pulse_sequences,
             }
                          

if len(sys.argv) == 2:
    function = dispatcher.get(sys.argv[1])
    
    if not function:
        print(INSTRUCTIONS, file=sys.stderr)
        sys.exit(-1)
else:
    # Wrong # of arguments
    print(INSTRUCTIONS, file=sys.stderr)
    sys.exit(-1)


# Get the database path & name
data_dir = util_misc.get_data_dir()
# The database name is in vespa.ini
config = util_config.VespaConfig()
db_path = os.path.join(data_dir, config["database"]["filename"])

db = util_db.Database(db_path, True)

function(db, sys.stdout)
