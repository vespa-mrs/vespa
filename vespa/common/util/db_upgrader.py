# Python modules

import shutil
import logging
import os
import xml.etree.cElementTree as ElementTree
import sqlite3 as sqlite

# 3rd party modules
import numpy

# Our modules
import vespa.common.util.db as db_module
from vespa.common.util.db import _Fetch
import vespa.common.constants as constants
import vespa.common.util.misc as util_misc
import vespa.common.util.time_ as util_time
import vespa.common.util.import_ as util_import
import vespa.common.util.logging_ as util_logging
import vespa.common.rfp_machine_specs as rfp_machine_specs

# deprecated import vespa.common.rfp_machine_settings as rfp_machine_settings
# deprecated from vespa.common.rfp_pulse_design import _convert_project_to_design


"""
This very specialized module is an extension of our standard database module
that adds a few functions for upgrading the database. It's meant to be used by
init code which will upgrade the database to the most current version so
that mainstream code never has to worry about dealing with old versions.

To add a new upgrade method, scroll down and examine the template (in
a big comment) and/or copy & paste from an existing method.
"""

RESOURCES_PATH = os.path.join(util_misc.get_vespa_install_directory(), 
                              "common", "resources", "db_upgrader")


# ---------------------------- Philosophy ----------------------------
# There are two ways would can do database upgrades: stepwise and in jumps. 
# The difference is that stepwise upgrades only increase the version one step 
# at a time, e.g. 1==>2, 2==>3, 3==>4, etc. whereas jumps upgrade the database 
# from any version to current in one fell swoop, e.g. 1==>4.

# We only do stepwise upgrades to keep our lives simple. The advantage is that
# we write individual functions for each upgrade step and they never have to
# be revisited. The disadvantage is that it's less efficient than jumping when
# Vespa finds a database that needs to be upgraded more than one step. For
# instance, if Vespa finds a database that's at version 1 while the current
# version is 4, three serial upgrades will happen (1==>2, 2==>3, and 3==>4).

# By contrast, a jump upgrade would handle the scenario above by upgrading the
# database from 1==>4 in one shot. The downside is that someone needs to write
# the code for that. That means that if the current database is version N,
# then N-1 functions need to exist. For example if N == 4, we'd need functions
# to upgrade 1==>4, 2==>4 and 3==>4 to cover all possible scenarios.

# Furthermore, every time a new database version is released, *each function
# would need to be touched*. For instance, if version 5 was the new version,
# we'd add a function to upgrade 4==>5 AND change these functions:
# 1==>4, 2==>4 and 3==>4
# to these:
# 1==>5, 2==>5 and 3==>5

# To summarize, stepwise upgrades can be inefficient in certain scenarios but
# those scenarios are uncommon. The benefit of using stepwise upgrades is that
# each new database version requires writing just one new function.

# Jump upgrades never fall into the inefficiency trap that stepwise upgrades
# sometimes do. However, when we introduce new database version N then we need
# to write a new function and edit (and debug and test) N-1 existing
# functions.
# --------------------------------------------------------------------------

        
class DatabaseUpgrader(db_module.Database):
    """A database class expressly for upgrading the database version if
    needed. Instantiate and call upgrade() to force an upgrade to the most
    recent version. If no upgrade is needed, calling upgrade() is harmless.
    """
    def __init__(self, filename):
        db_module.Database.__init__(self, filename)


    def fetch_version(self):
        """Returns the database version number as an int."""
        sql = """SELECT
                    database_version
                 FROM
                    vespa
              """
        version = self._execute(sql, None, _Fetch.SINGLETON)

        return int(version)


    def update_version(self, version):
        """Writes the version param (which must be an int) to the database."""
        sql = """UPDATE
                    vespa
                 SET
                    database_version = ?
              """
        self._execute(sql, int(version), _Fetch.NONE)


    def upgrade(self):
        """Upgrades the database to the most recent version. If no upgrade is
        needed, calling this is harmless.
        """
        database_version = self.fetch_version()
        
        upgrade_happened = False

        while database_version < constants.DATABASE_VERSION:
            upgrade_happened = True
            # Upgrade from database_version to database_version + 1
            
            logger = logging.getLogger(util_logging.Log.DATABASE)
            # Set up a handler for this log
            logger.addHandler(util_logging.DatabaseFileHandler())
            
            # Write some metadata to the log
            logger.info("Vespa Database Upgrade Log")
            logger.info("-" * 60)
            logger.info("Starting...")
            
            logger.info("Current database version is %d" % database_version)

            # Copy the current database as a CYA. The filename to which it's 
            # copied contains the db version and a timestamp.
            target = self._filename + (".v%d" % database_version)
            target += util_time.now(".%Y%m%d.%H%M%S")
            
            logger.info("Backing up current database to %s" % target)
            shutil.copyfile(self._filename, target)

            # Call the appropriate upgrade method for this version
            method = getattr(self, "_upgrade_%d" % database_version)
            method(logger)
            
            logger.info("Done!")

            database_version = self.fetch_version()

        if upgrade_happened:
            # Clean up. Since we have enabled SQLite's auto-vacuuming, this
            # isn't strictly necessary but it might help performance esp. 
            # after an upgrade which may have involved lots of changes to 
            # the database.
            self._execute("VACUUM")


    # ---------------------------------------------------------------------
    # Upgrade methods follow. Each is named "_upgrade_N" where N is the
    # version of the current database, i.e. the version you're upgrading FROM.
    # e.g to upgrade from 1==>2, call self._upgrade_1(logger).

    # Template for upgrading N==>N+1, using N==2 as an example --
    
    # def _upgrade_2(self, logger):
    #     # Upgrades the database from version 2 to 3
    #     # Specific changes --
    #     # - Changed all foos to bars
    #     # - Added table whizzo_quality_assortment
    #     # - Removed index from pepper.pot
    #     # - etc.
    #
    #     self.begin_transaction()
    #
    #     # Do important stuff here
    #     #logger.info("Doing important stuff...")
    #     pass
    #
    #     # Bump the version number
    #     self.update_version(3)
    #
    #     self.commit()

    def _upgrade_1(self, logger):
        # Upgrades the database from version 1 to 2
        # Version 1 was used for Vespa <= 0.1.6 which was Simulation only
        # Version 2 was used for Vespa >= 0.2.0 which was Simulation + RFPulse

        msg = 'Vespa only updates from Database v10, you have v'+self.fetch_version()+', please do new install.'
        raise RuntimeError(msg)

        # # The change from 1==>2 is just adding RFPulse-specific tables.
        # self.begin_transaction()
        #
        # # Do important stuff here
        # logger.info("Adding RFPulse-specific tables...")
        #
        # sql = open(os.path.join(RESOURCES_PATH, "1_to_2.sql")).read()
        #
        # self.executescript(sql)
        #
        # logger.info("Adding default machine settings templates...")
        #
        # filename = os.path.join(RESOURCES_PATH, "1_to_2_machine_settings_templates.xml")
        #
        # self._insert_machine_settings_templates(filename)
        #
        # # PS This code used to work but now causes a problem. It will only be
        # # invoked for very old databases; see here:
        # # http://scion.duhs.duke.edu/vespa/project/ticket/41
        # # We'll never use this code again and ordinarily I prefer to delete
        # # dead code, but in this case I'll leave it here as sample code.
        # # logger.info("Importing pulse projects...")
        # # filename = os.path.join(RESOURCES_PATH, "1_to_2_pulse_projects.xml")
        # # importer = util_import.PulseProjectImporter(filename, self)
        # # importer.go(False)
        #
        # # Bump the version number
        # self.update_version(2)
        #
        # self.commit()


    def _upgrade_2(self, logger):
        # Upgrades the database from version 2 to 3
        # Database version 2 was used in Vespa 0.2.0 and 0.2.1.
        # Database version 3 was first used in Vespa 0.2.2.

        msg = 'Vespa only updates from Database v10, you have v'+self.fetch_version()+', please do new install.'
        raise RuntimeError(msg)

        # # The change from 2==>3 moves dims info out of the simulations table,
        # # moves spectrum data into the simulations table, discards
        # # spectrum_lines, changes experiment_metabolites from a view to a
        # # real table. Also, spectrum data is now stored in BLOBs which
        # # makes for a much more compact database.
        #
        # self.begin_transaction()
        #
        # # Important stuff starts here...
        #
        # # We're going to replace the simulations table with a table of the
        # # same name but a different structure. For now, we rename the
        # # "old" (V2) table.
        # logger.info("Renaming simulations table...")
        # sql = """ALTER TABLE
        #             simulations
        #          RENAME TO
        #             simulations_v2
        #       """
        # self._execute(sql)
        #
        # # We're going to create a table called experiment_metabolites; the
        # # view will get in the way.
        # logger.info("Dropping experiment_metabolites view...")
        # sql = """DROP VIEW
        #             experiment_metabolites
        #       """
        # self._execute(sql)
        #
        # # OK, now we're done with the prep work. Let's create some new stuff.
        # logger.info("Creating new tables...")
        # sql = """CREATE TABLE experiment_dims (
        #             id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
        #             experiment_id TEXT NOT NULL REFERENCES experiments (id),
        #             dim1 REAL NOT NULL,
        #             dim2 REAL NOT NULL,
        #             dim3 REAL NOT NULL,
        #             UNIQUE (experiment_id, dim1, dim2, dim3)
        #             )
        #         """
        # self._execute(sql)
        #
        # sql = """CREATE TABLE simulations (
        #             metabolite_id INTEGER NOT NULL REFERENCES metabolites (id),
        #             dims_id INTEGER NOT NULL REFERENCES experiment_dims (id),
        #             started TIMESTAMP DEFAULT NULL,
        #             completed TIMESTAMP DEFAULT NULL,
        #             ppms BLOB DEFAULT NULL,
        #             areas BLOB DEFAULT NULL,
        #             phases BLOB DEFAULT NULL,
        #             PRIMARY KEY (metabolite_id, dims_id)
        #             )
        #         """
        # self._execute(sql)
        #
        # sql = """CREATE TABLE experiment_metabolites (
        #             experiment_id TEXT NOT NULL,
        #             metabolite_id TEXT NOT NULL,
        #             PRIMARY KEY (experiment_id, metabolite_id)
        #             )
        #         """
        # self._execute(sql)
        #
        # # Move data from the V2 to the V3 tables
        # logger.info("Breaking up simulations...")
        # # The V2 simulations table has been renamed to simulations_v2. The
        # # V3 simulations table is just called "simulations". For each row
        # # in V2, I create a corresponding entry in the tables "simulations"
        # # and "experiment_dims".
        # sql = """SELECT
        #             *
        #          FROM
        #             simulations_v2
        #       """
        # rows = self._execute(sql)
        #
        # # Track which metabs are associated with which experiments so I can
        # # populate the experiment_metabolites table after rebuilding all the
        # # simulations.
        # experiment_metabolite_map = { }
        # # In this map I track the ids of the rows in the experiment_dims
        # # table. The combination of (experiment, dim1, dim2, dim3) comprises
        # # a primary key (i.e. is unique), and so I use that as the key to
        # # this dict. Each key points to id of the corresponding row in
        # # experiment_dims.
        # experiment_dims_map = { }
        #
        # for simulation in rows:
        #     # Update the metabs map
        #     experiment_id = simulation["experiment_id"]
        #     metabolite_id = simulation["metabolite_id"]
        #     if experiment_id not in experiment_metabolite_map:
        #         # For each experiment we track the list of associated metabs.
        #         # Each metab only need be listed once, so we use a set.
        #         experiment_metabolite_map[experiment_id] = set()
        #
        #     experiment_metabolite_map[experiment_id].add(metabolite_id)
        #
        #     # Create & populate an entry in the dims map, if necessary.
        #     key = (experiment_id, simulation["dim1"], simulation["dim2"],
        #                           simulation["dim3"])
        #     if key not in experiment_dims_map:
        #         # Create the experiment_dims entry
        #         sql = """INSERT INTO
        #                      experiment_dims
        #                         (experiment_id, dim1, dim2, dim3)
        #                  VALUES
        #                      (?, ?, ?, ?)
        #               """
        #         dims_id = self._execute(sql, key, _Fetch.LAST_ID)
        #         experiment_dims_map[key] = dims_id
        #     else:
        #         dims_id = experiment_dims_map[key]
        #
        #     # fetch the spectrum from the V2 table.
        #     ppms, areas, phases = self._fetch_spectrum(simulation["id"])
        #
        #     ppms = db_module._numpy_array_to_blob(ppms)
        #     areas = db_module._numpy_array_to_blob(areas)
        #     phases = db_module._numpy_array_to_blob(phases)
        #
        #     # Create the simulation itself, including the spectra
        #     sql = """INSERT INTO
        #                  simulations
        #                      (metabolite_id, dims_id,
        #                       started, completed,
        #                       ppms, areas, phases)
        #              VALUES
        #                  (?, ?, ?, ?, ?, ?, ?)
        #           """
        #     params = (metabolite_id, dims_id,
        #               simulation["started"], simulation["completed"],
        #               ppms, areas, phases
        #              )
        #     self._execute(sql, params, _Fetch.NONE)
        #
        # # We're done (permanently) with the old simulations table and
        # # with spectrum_lines.
        # self._execute("DROP TABLE simulations_v2")
        # self._execute("DROP TABLE spectrum_lines")
        #
        # # Now populate experiment_metabolites
        # logger.info("Populating experiment_metabolites...")
        # sql = """INSERT INTO
        #              experiment_metabolites
        #                  (experiment_id, metabolite_id)
        #          VALUES
        #              (?, ?)
        #       """
        # for experiment_id, metabolite_ids in experiment_metabolite_map.items():
        #     for metabolite_id in metabolite_ids:
        #         self._execute(sql, (experiment_id, metabolite_id))
        #
        #
        # # The indices on the V2 simulations table got dropped along with the
        # # table. Now that the inserting is done, I recreate indices.
        # sql = """CREATE INDEX
        #             simulations_metabolite_id_index
        #          ON
        #             simulations (metabolite_id)
        #       """
        # self._execute(sql)
        #
        # sql = """CREATE INDEX
        #             experiment_dims_experiment_id_index
        #          ON
        #             experiment_dims (experiment_id)
        #       """
        # self._execute(sql)
        #
        # sql = """CREATE INDEX
        #             simulations_dims_id_index
        #          ON
        #             simulations (dims_id)
        #       """
        # self._execute(sql)
        #
        # # Bump the version number
        # self.update_version(3)
        #
        # self.commit()
        
        
    def _upgrade_3(self, logger):
        # Upgrades the database from version 3 to 4
        # Database version 3 was used in Vespa 0.2.2 - 0.2.4.
        # Database version 4 was first used in Vespa 0.2.5.

        msg = 'Vespa only updates from Database v10, you have v'+self.fetch_version()+', please do new install.'
        raise RuntimeError(msg)

        # # Specific changes --
        # # - Added the table gaussian_pulse_parameters
        #
        # self.begin_transaction()
        #
        # logger.info("Adding table gaussian_pulse_parameters...")
        # sql = """CREATE TABLE gaussian_pulse_parameters (
        #             id INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
        #             tip_angle REAL NOT NULL  DEFAULT 0.0,
        #             time_steps INTEGER NOT NULL  DEFAULT 0,
        #             duration REAL NOT NULL  DEFAULT 0.0,
        #             bandwidth REAL NOT NULL  DEFAULT 0.0,
        #             filter_type TEXT NOT NULL  DEFAULT '''',
        #             filter_application REAL NOT NULL  DEFAULT 0.0
        #             )
        #         """
        # self._execute(sql)
        #
        # # Bump the version number
        # self.update_version(4)
        #
        # self.commit()


    def _upgrade_4(self, logger):
        # Upgrades the database from version 4 to 5
        # Database version 4 was used in Vespa 0.2.5.
        # Database version 5 was first used in Vespa 0.2.6.

        msg = 'Vespa only updates from Database v10, you have v'+self.fetch_version()+', please do new install.'
        raise RuntimeError(msg)

        # # Specific changes --
        # # - Added the table randomized_pulse_parameters
        # # - Added the table import_pulse_parameters
        #
        # self.begin_transaction()
        #
        # logger.info("Adding table randomized_pulse_parameters...")
        # sql = """CREATE TABLE randomized_pulse_parameters (
        #             id INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
        #             time_steps INTEGER NOT NULL DEFAULT 0,
        #             duration REAL NOT NULL DEFAULT 0.0
        #          )
        #       """
        # self._execute(sql)
        #
        # logger.info("Adding table import_pulse_parameters...")
        # sql = """CREATE TABLE import_pulse_parameters (
        #             id INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
        #             file_path TEXT NOT NULL  DEFAULT '''',
        #             comment TEXT DEFAULT NULL,
        #             dwell_time REAL NOT NULL  DEFAULT 0.0,
        #             max_intensity REAL NOT NULL  DEFAULT 0.0,
        #             file_format TEXT NOT NULL  DEFAULT ''''
        #             )
        #       """
        # self._execute(sql)
        #
        # # Bump the version number
        # self.update_version(5)
        #
        # self.commit()


    def _upgrade_5(self, logger):
        # Upgrades the database from version 5 to 6
        # Database version 5 was used in Vespa 0.2.6.
        # Database version 6 was first used in Vespa 0.3.0.

        msg = 'Vespa only updates from Database v10, you have v'+self.fetch_version()+', please do new install.'
        raise RuntimeError(msg)

        # # Specific changes --
        # # - Added the columns scale_factor and use_max_intensity to
        # #   the table import_pulse_parameters
        # # - Added the table pulse_sequence_pulse_projects
        #
        # logger.info("Adding column scale_factor...")
        # sql = """ALTER TABLE import_pulse_parameters
        #             ADD COLUMN
        #                 scale_factor REAL NOT NULL  DEFAULT 0.0
        #       """
        # self._execute(sql)
        #
        # logger.info("Adding column use_max_intensity...")
        # sql = """ALTER TABLE import_pulse_parameters
        #             ADD COLUMN
        #                 use_max_intensity BOOLEAN NOT NULL  DEFAULT '0'
        #       """
        # self._execute(sql)
        #
        # # Using max intensity used to be the default (it was the only
        # # option, actually) but now it's not. For backwards compatibility,
        # # we ensure that existing projects have this flag set to true.
        # logger.info("Updating existing projects...")
        # sql = """UPDATE import_pulse_parameters
        #             SET use_max_intensity = '0'
        #       """
        # self._execute(sql)
        #
        # logger.info("Adding table pulse_sequence_pulse_projects...")
        # sql = """CREATE TABLE pulse_sequence_pulse_projects (
        #             pulse_sequence_id TEXT NOT NULL  REFERENCES pulse_sequences (id),
        #             pulse_project_id TEXT NOT NULL  REFERENCES pulse_projects (id),
        #             progression INTEGER NOT NULL
        #         );
        #       """
        # self._execute(sql)
        #
        #
        # # Bump the version number
        # self.update_version(6)
        #
        # self.commit()


    def _upgrade_6(self, logger):
        # Upgrades the database from version 6 to 7
        # Database version 6 was used in Vespa 0.3.0.
        # Database version 7 was first used in Vespa 0.3.1

        msg = 'Vespa only updates from Database v10, you have v'+self.fetch_version()+', please do new install.'
        raise RuntimeError(msg)

        # # The change from 6 ==> 7 was all for support of optimal control in
        # # RFPulse.
        # # Specific changes --
        # # - Added the tables ocn_states, ocn_parameters, deltab1_points,
        # #   and ocn_residual_errors
        # # - Added the column ocn_state_id to the table results
        #
        # logger.info("Adding table ocn_states...")
        # sql = """CREATE TABLE ocn_states (
        #             id INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
        #             multiplier REAL NOT NULL ,
        #             met_max_iterations BOOLEAN NOT NULL ,
        #             met_residual_error BOOLEAN NOT NULL ,
        #             met_differential_error BOOLEAN NOT NULL ,
        #             met_max_time BOOLEAN NOT NULL ,
        #             met_increasing_error BOOLEAN NOT NULL ,
        #             run_time REAL NOT NULL ,
        #             iterations INTEGER NOT NULL ,
        #             decreases INTEGER NOT NULL
        #           );
        #       """
        # self._execute(sql)
        #
        # logger.info("Adding table deltab1_points...")
        # sql = """CREATE TABLE deltab1_points (
        #             ocn_state_id INTEGER NOT NULL  REFERENCES ocn_states (id),
        #             progression INTEGER NOT NULL ,
        #             real_amplitude REAL NOT NULL ,
        #             imaginary_amplitude REAL NOT NULL
        #           );
        #       """
        # self._execute(sql)
        #
        # logger.info("Adding table ocn_parameters...")
        # sql = """CREATE TABLE ocn_parameters (
        #             id INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
        #             pulse_type TEXT NOT NULL  DEFAULT '''',
        #             phase_type TEXT NOT NULL  DEFAULT '''',
        #             gradient_refocusing_value REAL NOT NULL  DEFAULT 0.0,
        #             tip_angle REAL NOT NULL  DEFAULT 0.0,
        #             bandwidth REAL NOT NULL  DEFAULT 0.0,
        #             step_size_multiplier REAL NOT NULL ,
        #             step_size_modification TEXT NOT NULL  DEFAULT '''',
        #             excite_band_points INTEGER NOT NULL  DEFAULT 0,
        #             b1_immunity_range REAL NOT NULL  DEFAULT 0.0,
        #             steps INTEGER NOT NULL  DEFAULT 0,
        #             b1_maximum REAL NOT NULL ,
        #             limit_sar BOOLEAN NOT NULL  DEFAULT '0',
        #             sar_factor REAL NOT NULL  DEFAULT 0.0,
        #             error_increase_tolerance REAL NOT NULL  DEFAULT 0.0,
        #             max_iteration_check BOOLEAN NOT NULL  DEFAULT '0',
        #             max_iterations INTEGER NOT NULL  DEFAULT 0,
        #             residual_error_check BOOLEAN NOT NULL  DEFAULT '0',
        #             residual_error_tolerance REAL NOT NULL  DEFAULT 0.0,
        #             differential_error_check BOOLEAN NOT NULL  DEFAULT '0',
        #             differential_error_tolerance REAL NOT NULL  DEFAULT 0.0,
        #             halt_if_error_increasing BOOLEAN NOT NULL  DEFAULT '0',
        #             halt_on_max_time BOOLEAN NOT NULL  DEFAULT '0',
        #             max_time REAL NOT NULL  DEFAULT 0.0,
        #             enforce_symmetry BOOLEAN NOT NULL  DEFAULT '0'
        #         );
        #
        #       """
        # self._execute(sql)
        #
        # logger.info("Adding table ocn_residual_errors...")
        # sql = """CREATE TABLE ocn_residual_errors (
        #             ocn_state_id INTEGER NOT NULL  REFERENCES ocn_states (id),
        #             value REAL NOT NULL ,
        #             progression INTEGER NOT NULL
        #           );
        #       """
        # self._execute(sql)
        #
        # logger.info("Adding column ocn_state_id...")
        # sql = """ALTER TABLE results
        #             ADD COLUMN
        #                 ocn_state_id INTEGER NOT NULL  DEFAULT 0 REFERENCES ocn_states (id)
        #       """
        # self._execute(sql)
        #
        # # Bump the version number
        # self.update_version(7)
        #
        # self.commit()


    def _upgrade_7(self, logger):
        # Upgrades the database from version 7 to 8
        # Database version 7 was used from Vespa 0.3.1 to 0.8.4.
        # Database version 8 was first used in Vespa 0.8.5  

        msg = 'Vespa only updates from Database v10, you have v'+self.fetch_version()+', please do new install.'
        raise RuntimeError(msg)

        # # The change from 7 ==> 8 was where we added the Pulse application
        # # in alpha release to test if it could replace the RFPulse app.
        # # Specific changes --
        # # The change from 7==>8 is just adding Pulse-specific tables and
        # # a few default machine specs, Create transform kernels, and
        # # example Pulse Designs.
        # self.begin_transaction()
        #
        # # Do important stuff here
        # logger.info("Adding Pulse-specific tables...")
        #
        # sql = open(os.path.join(RESOURCES_PATH, "7_to_8.sql")).read()
        #
        # self.executescript(sql)
        #
        # logger.info("Adding default machine specs templates...")
        # filename = os.path.join(RESOURCES_PATH, "7_to_8_machine_specs_templates.xml")
        # self._insert_machine_specs_templates(filename)
        #
        # # Note. BJS, as of _upgrade_9, the pulse_design table has an additional
        # #  attribute in it for gyromagnetic_nuclei. Here in _upgrade_7, we have
        # #  not added that yet. That's not an issue for these transform kernels
        # #  except that they do not take advantage of the gyromagnetic nuclei
        # #  value in their algorithms. The new kernels added in _upgrade_9 do
        # #  make use of that value. So, since these kernels are deprecated, I
        # #  choose not to import them here since the new ones will be added in
        # #  the _upgrade_9 step.
        # #
        # # logger.info("Importing pulse transform kernels ...")
        # # filename = os.path.join(RESOURCES_PATH, "7_to_8_transform_kernels.xml")
        # # importer = util_import.TransformKernelImporter(filename, self)
        # # importer.go(False)
        #
        # # Note. BJS, as of _upgrade_9, the pulse_design table has an additional
        # #  attribute in it for gyromagnetic_nuclei. Here in _upgrade_7, we have
        # #  not added that yet. So the 7_to_9_pulse_designs in the xml file will
        # #  not load in properly, causing an exception.  If the user is going
        # #  all the way to _upgrade_9, they will get a set of pulse designs at
        # #  that point, so I'm removing this step here.
        # #
        # # logger.info("Importing example pulse designs ...")
        # # filename = os.path.join(RESOURCES_PATH, "7_to_8_pulse_designs.xml")
        # # importer = util_import.PulseDesignImporter(filename, self)
        # # importer.go(False)
        #
        # # Bump the version number
        # self.update_version(8)
        #
        # self.commit()


    def _upgrade_8(self, logger):
        # Upgrades the database from version 8 to 9
        # Database version 8 was first used in Vespa 0.8.5
        # Database version 9 was first used in Vespa 0.8.6
        # The change from 8 ==> 9 was where we fully implement the Pulse
        # application and use it to replace the RFPulse app.
        #

        msg = 'Vespa only updates from Database v10, you have v'+self.fetch_version()+', please do new install.'
        raise RuntimeError(msg)

        # # Specific changes --
        # # The change from 8==>9 does not change the footprint of the
        # # database, it just transfers all existing RFPulse results into
        # # Pulse results. The Pulse versions are just Imports of the final
        # # RF waveform in the RFPulse result.
        #
        # #
        # # 1. Import new Pulse transform kernel we need to convert RFPulse results
        # #
        #
        # self.begin_transaction()
        #
        # logger.info("Importing pulse transform kernels ...")
        # filename = os.path.join(RESOURCES_PATH, "8_to_9_transform_kernels.xml")
        # importer = util_import.TransformKernelImporter(filename, self)
        # importer.go(False)
        #
        # self.commit()
        #
        # #
        # # 2. Transfer all RFPulse project results to Pulse design results
        # #
        #
        # self.begin_transaction()
        #
        # logger.info("Converting all RFPulse pulse projects into Pulse pulse design objects ...")
        #
        # projects = self.fetch_pulse_project_previews()
        # for project in projects:
        #     # check if project already in Pulse table, if not convert and insert
        #     if self.count_pulse_designs(project.id):
        #         # Don't bother to import this; it's already in the database.
        #         logger.info("Ignoring pulse design %s, name %s, already in database" % (project.id, project.name))
        #     else:
        #         # fetch the next project and get a representation of its final waveform
        #         rfpulse_project = self.fetch_pulse_project(project.id)
        #         # embedded conversion into the pulse_design module
        #         pulse_design = _convert_project_to_design(self, rfpulse_project)
        #
        #         # as of db version 10, we have a new column in the pulse_design
        #         # table, so we can not call the insert_pulse_design() method.
        #         # We have created a backwards compatible method instead for
        #         # this one use.
        #         #
        #         # self.insert_pulse_design(pulse_design)
        #         self.insert_rfpulse_convert_pulse_design(pulse_design)
        #         logger.info(pulse_design.comment)
        #
        #
        #
        # # 3. Transfer Simulation referrer data between Pulse/RFPulse referrer tables
        # #
        # # Do this in three steps
        # # a. get all data from pulse_sequence_pulse_projects table
        # # b. remove all entries in pulse_sequence_pulse_designs table with a
        # #    pulse_design_id == pulse_project_id just in case there is some
        # #    older version with different pulses referred to
        # # c. insert all entries from pulse_sequence_pulse_projects table into
        # #    the pulse_sequence_pulse_designs table
        #
        # sql = """SELECT
        #             pulse_sequence_id, pulse_project_id, progression
        #          FROM
        #             pulse_sequence_pulse_projects
        #          WHERE
        #             1=1
        #       """
        # all_sequences = [(row[0], row[1], row[2]) for row in self._execute(sql)]
        #
        # all_ids = [item[0] for item in all_sequences]
        # sql = """DELETE FROM
        #             pulse_sequence_pulse_designs
        #          WHERE
        #             pulse_sequence_id = ?"""
        # for id in all_ids:
        #     self._execute(sql, id, _Fetch.NONE)
        #
        # sql = """INSERT INTO
        #             pulse_sequence_pulse_designs
        #                 (pulse_sequence_id, pulse_design_id, progression)
        #          VALUES
        #             (?, ?, ?)
        #       """
        # for row in all_sequences:
        #     self._execute(sql, (row[0], row[1], row[2]))
        #
        #
        # # Bump the version number
        # self.update_version(9)
        #
        # self.commit()


    def _upgrade_9(self, logger):
        # Upgrades the database from version 9 to 10
        #
        # - Database version 9 was first used in Vespa 0.8.6  
        # - Database version 10 was firt used in Vespa 0.9.4
        #
        # The change from 9 ==> 10 was where we made Pulse gamma values
        #  user selectable. Which required an additional entry in the
        #  master parameters table.
        #

        msg = 'Vespa only updates from Database v10, you have v'+self.fetch_version()+', please do new install.'
        raise RuntimeError(msg)

        # # Specific changes --
        # # We had to change all the kernels to make use of the system gamma value
        # # Examples files default to 1H, but we need a new Example Pulse Design
        # #  to demo the user selection
        # #
        #
        # #
        # # 1. Import new entry into master parameters
        # #
        #
        # self.begin_transaction()
        #
        # logger.info("Adding column gyromagnetic_nuclei ...")
        # sql = """ALTER TABLE pulse_designs
        #             ADD COLUMN
        #                 gyromagnetic_nuclei TEXT NOT NULL DEFAULT '1H'
        #       """
        # self._execute(sql)
        #
        # self.commit()
        #
        # #
        # # 2. Import new Pulse kernels and example files
        # #
        #
        # self.begin_transaction()
        #
        # logger.info("Importing pulse transform kernels ...")
        # filename = os.path.join(RESOURCES_PATH, "9_to_10_transform_kernels.xml")
        # importer = util_import.TransformKernelImporter(filename, self)
        # importer.go(False)
        #
        # logger.info("Importing example pulse designs ...")
        # filename = os.path.join(RESOURCES_PATH, "9_to_10_pulse_designs.xml")
        # importer = util_import.PulseDesignImporter(filename, self)
        # importer.go(False)
        #
        # # Bump the version number
        # self.update_version(10)
        #
        # self.commit()

  
    def _upgrade_10(self, logger):
        # Upgrades the database from version 10 to 11
        #
        # - Database version 10 was first used in Vespa 0.9.4
        # - Database version 11 was first used in Vespa 0.9.5
        #
        # The change from 10 ==> 11 was to fix a bug in the base template for
        # importing pulse waveforms from a file. Dwell times were being truncated
        # to nearest usec. Changed import type from int to float. This imports a 
        # new kernel with an updated version number and leaves the old buggy one
        # in place.
        #
        # Specific changes --
        #
        
        #
        # 1. Import new Pulse kernel
        #
         
        self.begin_transaction()
        
        logger.info("Importing pulse transform kernels ...")
        filename = os.path.join(RESOURCES_PATH, "10_to_11_transform_kernels.xml")
        importer = util_import.TransformKernelImporter(filename, self)
        importer.go(False)
        
        # Bump the version number
        self.update_version(11)
        
        self.commit()        


    def _upgrade_11(self, logger):
        # Upgrades the database from version 11 to 12
        #
        # - Database version 11 was first used in Vespa 0.9.4
        # - Database version 12 was first used in Vespa 0.9.5  
        #
        # The change from 11 ==> 12 was where we made Pulse transform kernels
        #  able to be deprecated. Which required an additional entry in the
        #  transform_kernels table. 
        #
        # The change from 11 ==> 12 was to fix a bug in a pulse kernel that 
        # caused an error in the latest version of Python/Numpy. Previous 
        # versions of numpy threw a warning for this bug, but newer version
        # throws an exception. So I had to fix an exisiting transform_kernel
        # that was 'public' and figure some way to swap this out. I decided 
        # to create a 'deprecated' field in the transform_kernels table that
        # would allow the dialog for browsing kernels to choose whether to
        # display the 'deprecated' kernels or not. This leaves the older or
        # buggier kernel in the database, but keeps the GUI from being 
        # messy.
        
        #
        # 1. Import new entry into transform_kernels table
        #
        
        self.begin_transaction()
        
        logger.info("Adding column deprecated to transform_kernels table ...")
        sql = """ALTER TABLE transform_kernels
                    ADD COLUMN
                        deprecated BOOLEAN NOT NULL  DEFAULT '0'
              """
        self._execute(sql)

        self.commit()

        #
        # 2. Mark out-of-date tranform kernels as deprecated = True
        #
        
        # set the deprecate field to True on a number of transform_kernels. Some had bugs
        #   but most were just pre-multinuc. This does not delete them, just allows the 
        #   gui to know if they should be displayed or not.
        
        logger.info("Setting deprecate flag on pulse transform kernels ...")
        
        uuid_list = ["9f7f9bae-8d38-49ca-8df0-ea71d9deea2c",        # Matpulse SLR multinuc v1 buggy
                     "563cf919-d3f0-41ae-8651-0914e0a5d9bf",        # Dinesh GOIA multinuc v1 buggy
                     "3d58c135-f9ff-44a7-a8b0-2b73ce9f7e2f",        # Dinesh GOIA 1H
                     "1278ab4f-bcd9-48fd-aff8-6d97e3062181",        # Bassi 1H
                     "c1960e53-0ef7-422c-8600-f457c033b107",        # simple sinc 1H
                     "814e1d0e-7068-4583-9110-72f0ccf4914e",        # Import from file no Gradient 1H
                     "1930f3f6-322a-48d3-811d-279ba0dce514",        # Matpulse Gaussian 1H
                     "d6b9efe0-1d85-43f4-871c-3b06aa637ec6",        # Matpulse Hyp Sec 1H
                     "c004c294-1a57-472c-96b3-aa65e9fb9383",        # Matpulse Interp-Rescale 1H
                     "ffc9c755-ad1f-4d37-af19-cb5993bb525d",        # Matpulse SLR 1H
                     ]
        
        for id_ in uuid_list:
            self.begin_transaction()

            sql = """UPDATE
                        transform_kernels
                     SET
                        deprecated = 1
                     WHERE
                        id = ?
                     """
            self._execute(sql, id_, _Fetch.NONE)

            self.commit()
            
        #
        # 3. Import updated Pulse kernels and updated example files
        #
        
        self.begin_transaction()
        
        logger.info("Importing pulse transform kernels ...")
        filename = os.path.join(RESOURCES_PATH, "11_to_12_transform_kernels.xml")
        importer = util_import.TransformKernelImporter(filename, self)
        importer.go(False)

        logger.info("Importing example pulse designs ...")
        filename = os.path.join(RESOURCES_PATH, "11_to_12_pulse_designs.xml")
        importer = util_import.PulseDesignImporter(filename, self)
        importer.go(False)

        # Bump the version number
        self.update_version(12)

        self.commit()        
        

    ##################     Helper  functions    ######################


    # def _insert_machine_settings_templates(self, filename):
    #     """
    #     This is a helper function for _upgrade_1().
    #
    #     Machine settings templates don't participate in import/export. I wrote
    #     this special function to populate the machine_settings_templates table
    #     from an XML file that I wrote by hand although it conforms to VIFF.
    #
    #     """
    #     # Read the XML file
    #     tree = ElementTree.ElementTree(file=filename)
    #
    #     root = tree.getroot()
    #
    #     machine_settings_elements = root.findall("machine_settings")
    #     ids = [ ]
    #     for e in machine_settings_elements:
    #         machine_settings = rfp_machine_settings.MachineSettingsTemplate(e)
    #         ids.append(self.insert_machine_settings(machine_settings))
    #
    #     # There must always be a default template.
    #     self.mark_machine_settings_template_as_default(ids[0])


    def _insert_machine_specs_templates(self, filename):
        """
        This is a helper function for _upgrade_7().
    
        Machine specs templates don't participate in import/export. I wrote
        this special function to populate the machine_specs_templates table
        from an XML file that I wrote by hand although it conforms to VIFF.
        
        """
        # Read the XML file
        tree = ElementTree.ElementTree(file=filename)
    
        root = tree.getroot()
    
        machine_specs_elements = root.findall("machine_specs")
        ids = [ ]
        for e in machine_specs_elements:
            machine_specs = rfp_machine_specs.MachineSpecsTemplate(e)
            ids.append(self.insert_machine_specs(machine_specs))
       
        # There must always be a default template. 
        self.mark_machine_specs_template_as_default(ids[0])


    def _fetch_spectrum(self, simulation_id):
        """This is a helper function for _upgrade_2()."""
        # Fetch the spectrum lines for this simulation
        sql = """SELECT
                    ppm, area, phase
                 FROM
                    spectrum_lines
                 WHERE
                    simulation_id = ?
                 ORDER BY
                    ppm
              """
        rows = self._execute(sql, simulation_id)
        
        if rows:
            # Here we use the nifty Python idiom zip(*rows) for matrix 
            # transposition as demonstrated in the example below:
            # >>> rows = [ (1, 2, 3), (1, 2, 3), (1, 2, 3), (1, 2, 3), ]
            # >>> print zip(*rows)
            # [(1, 1, 1, 1), (2, 2, 2, 2), (3, 3, 3, 3)]
        
            # In this case it will turn N rows of (ppm, area, phase) into
            # a 3-tuple of (ppms, areas, phases) in which each list has N 
            # items. That's perfect except that I want numpy arrays instead 
            # of lists.
            return [numpy.array(item, float) for item in zip(*rows)]
        else:
            # This is the unusual condition where a simulation has no 
            # spectrum lines associated. This can happen in two ways. Either 
            # it produced none (rare but possible) or it's an experiment that
            # was exported with no results and then re-imported.
            return [numpy.array([], float), 
                    numpy.array([], float), 
                    numpy.array([], float), ]
                    

    
    
    
