# Python modules

import re
import os
import tempfile
import shutil
import xml.etree.cElementTree as ElementTree
import decimal

# 3rd party modules

# Our modules
import vespa.common.constants as constants
import vespa.common.util.db as util_db
import vespa.common.util.import_ as util_import
import vespa.common.util.misc as util_misc
import vespa.common.rfp_machine_settings as rfp_machine_settings
import vespa.common.rfp_machine_specs as rfp_machine_specs

# _ISOTOPE_REGEX matches isotope strings that are already in our preferred 
# format. (See normalize_isotope_name() in util_misc.)
# This is used by _sort_isotope_names().
_ISOTOPE_REGEX = re.compile("""
    ^              # String starts here
    (\d{1,2})      # Match & capture 1 or 2 digits
    ([A-Z]{1,2})   # Match & capture 1 or 2 letters
    $              # String must end here
                           """, re.VERBOSE)

def create_database(target_filename, logger):
    """A utility function that deletes the existing database and creates a
    new one based on the XML export files that ship with the app.
    
    You can think of function as "reset to factory default".
    """
    # Write some metadata to the log
    logger.info("Vespa Database Creation Log")
    logger.info("-" * 60)
    logger.info("Starting...")
    
    # Here's a slight hack for ya...this code relies on/assumes that the
    # the resources directory is immediately below where this file lives.
    resources_path = os.path.join(os.path.split(__file__)[0], "resources")
    logger.info("Resources path is %s" % resources_path)

    logger.info("Database version is %d" % constants.DATABASE_VERSION)
    
    # mkstemp() returns an open file descriptor. I don't want to write to 
    # the file via the file descriptor so I close it immediately.
    fd, temp_filename = tempfile.mkstemp()
    os.close(fd)
    
    logger.info("Creating new database %s..." % temp_filename)
    db = util_db.Database(temp_filename)

    # First things first -- set auto vaccum on. SQLite insists that we do this
    # before any tables are created.
    db.executescript("PRAGMA auto_vacuum = 1;")

    # These are the SQL scripts that (re)create the database. 
    logger.info("Executing create_tables.sql...")
    filename = os.path.join(resources_path, "create_tables.sql")
    db.execute_script_from_file(filename)
    
    logger.info("Executing create_views.sql...")
    filename = os.path.join(resources_path, "create_views.sql")
    db.execute_script_from_file(filename)

    logger.info("Populating database version...")
    # We write the version into the database so that future versions of 
    # Vespa will know how to upgrade the database when the format changes.
    # The code below does a slight end run around our goal of keeping all SQL 
    # confined to the db.py modules. Shhhh, don't tell anyone!
    sql = "INSERT INTO vespa (database_version) VALUES (%d);" % \
                                                    constants.DATABASE_VERSION
    db.executescript(sql)

    logger.info("Creating isotopes...")
    filename = os.path.join(resources_path, "isotopes.xml")
    _insert_isotopes(db, filename)
    
    logger.info("Creating B0 bins...")
    filename = os.path.join(resources_path, "b0_bins.xml")
    _insert_b0_bins(db, filename)
    
    # The order in which we import metabs, pulse seqs and experiments doesn't
    # matter in theory. However, in practice this is a bit of defensive 
    # programming. That's because experiments.xml contains definitions of
    # metabs & pulse seqs. If a metab and/or pulse seq changes it's easy 
    # enough to remember to update the appropriate XML file. However, it's 
    # equally easy to forget to update experiments.xml which may also contain
    # the metab or pulse seq. Thus the definitions in experiments.xml are
    # sometimes out of date. Incidences of this have been rare, but it has
    # happened.
    # As long as we import the metabs and pulse seqs first, their definitions
    # in experiments.xml (which are more likely to be wrong) will be ignored 
    # which buys us extra time to discover/notice that they're out of sync.
    
    # Ignore part of above.  Need to import Pulse bits first in case I need
    # to convert pulse project to pulse design
    
    logger.info("Importing transform kernels...")
    filename = os.path.join(resources_path, "transform_kernels.xml")
    importer = util_import.TransformKernelImporter(filename, db)
    importer.go(False)

    logger.info("Importing machine specs...")
    filename = os.path.join(resources_path, "machine_specs_templates.xml")
    _insert_machine_specs_templates(db, filename)

    logger.info("Importing pulse designs...")
    filename = os.path.join(resources_path, "pulse_designs.xml")
    importer = util_import.PulseDesignImporter(filename, db)
    importer.go(False)
    
    logger.info("Importing pulse sequences...")
    filename = os.path.join(resources_path, "pulse_sequences.xml")
    importer = util_import.PulseSequenceImporter(filename, db)
    importer.go(False)
    
    logger.info("Importing metabolites...")
    filename = os.path.join(resources_path, "metabolites.xml")
    importer = util_import.MetaboliteImporter(filename, db)
    importer.go(False)
    
    logger.info("Importing experiments...")
    filename = os.path.join(resources_path, "experiments.xml")
    importer = util_import.ExperimentImporter(filename, db)
    importer.go(False)

#    logger.info("Importing machine settings...")
#    filename = os.path.join(resources_path, "machine_settings_templates.xml")
#    _insert_machine_settings_templates(db, filename)
#
#    logger.info("Importing pulse projects...")
#    filename = os.path.join(resources_path, "pulse_projects.xml")
#    importer = util_import.PulseProjectImporter(filename, db)
#    importer.go(False)

    # I create indices last -- inserts go faster without indices to update.
    logger.info("Executing create_indices.sql...")
    filename = os.path.join(resources_path, "create_indices.sql")
    db.execute_script_from_file(filename)

    # Close the db and ensure I've dropped all references to it before I copy 
    # it and remove the tempfile.
    db.close()

    del db

    if os.path.exists(target_filename):
        logger.info("Deleting target (%s)..." % target_filename)
        os.remove(target_filename)

    logger.info("Copying %s to %s..." % (temp_filename, target_filename))
    shutil.copyfile(temp_filename, target_filename)

    # Clean up the temp file
    logger.info("Removing %s..." % temp_filename)
    os.remove(temp_filename)

    logger.info("Done!")


def _insert_b0_bins(db, filename):
    """This is a helper function for create_database().
    
    B0 bins are a bit different than everything else. They're one of the 
    few database tables that aren't part of export/import. I wrote this 
    special function to populate the table from a special XML file that looks
    like a Vespa export but is actually just handcoded.
    """
    # Read the XML file
    tree = ElementTree.ElementTree(file=filename)
    
    root = tree.getroot()
    
    # The XML file contains a bin width. Our database wants a lower and 
    # upper bound for the bin, so we divide the width by two. To avoid
    # floating point noise, we use a Decimal object.
    width = decimal.Decimal(root.findtext("bin_width")) / 2
    
    bins = root.find("bins").findall("bin")
    
    bins = [decimal.Decimal(element.text) for element in bins]
    
    # bins is now a list of Decimals like [45.0, 64.0, 124.0, ...]
    # They represent the bin centers.
    
    # Build a script to insert them
    sql = """INSERT INTO
                b0_bins (left, center, right)
             values
                (%f, %f, %f);
          """
    lines = [sql % (center - width, center, center + width) for center in bins]
    
    db.executescript("\n".join(lines))
    

def _insert_isotopes(db, filename):
    """This is a helper function for create_database().
    
    Isotopes are a bit different than everything else. They're one of the 
    few database tables that aren't part of export/import. I wrote this 
    special function to populate the table from a special XML file that looks
    like a Vespa export but is actually just handcoded.

    In a pinch, should it prove that we've overlooked an isotope that someone
    wants to use, s/he could simply add it to the XML file and recreate the
    database.
    """
    # Read the XML file
    tree = ElementTree.ElementTree(file=filename)
    
    root = tree.getroot()
    
    isotopes = root.findall("isotope")
    
    isotopes = [element.text for element in isotopes]
    
    # isotopes is now a list of strings like ['1H', '2H', ... ]
    
    # Ensure there are no duplicates.
    isotopes = set(isotopes)
    
    # Ensure they're properly sorted. 
    isotopes = _sort_isotope_names(isotopes)
    
    # Build a script to insert them
    sql = """INSERT INTO
                isotopes (name, display_order)
             values
                ('%s', %d);
          """
    lines = [sql % (isotope, i) for i, isotope in enumerate(isotopes)]
    
    db.executescript("\n".join(lines))
    

#def _insert_machine_settings_templates(db, filename):
#    """This is a helper function for create_database().
#    
#    Machine settings templates don't participate in import/export. I wrote
#    this special function to populate the machine_settings_templates table
#    from an XML file that I wrote by hand although it follows the Vespa
#    interchange format.
#    """
#    # Read the XML file
#    tree = ElementTree.ElementTree(file=filename)
#    
#    root = tree.getroot()
#    
#    machine_settings_elements = root.findall("machine_settings")
#    ids = [ ]
#    for e in machine_settings_elements:
#        machine_settings = rfp_machine_settings.MachineSettingsTemplate(e)
#        ids.append(db.insert_machine_settings(machine_settings))
#       
#    # There must always be a default template. 
#    db.mark_machine_settings_template_as_default(ids[0])


def _insert_machine_specs_templates(db, filename):
    """This is a helper function for create_database().
    
    Machine specs templates don't participate in import/export. I wrote
    this special function to populate the machine_specs_templates table
    from an XML file that I wrote by hand although it follows the Vespa
    interchange format.
    """
    # Read the XML file
    tree = ElementTree.ElementTree(file=filename)
    
    root = tree.getroot()
    
    machine_specs_elements = root.findall("machine_specs")
    ids = [ ]
    for e in machine_specs_elements:
        machine_specs = rfp_machine_specs.MachineSpecsTemplate(e)
        ids.append(db.insert_machine_specs(machine_specs))
       
    # There must always be a default template. 
    db.mark_machine_specs_template_as_default(ids[0])
    

def _sort_isotope_names(isotopes):
    """Given an iterable (e.g. list) of valid isotope names, returns them 
    in sorted order. 

    This function is necessary because a naive alphabetic sort returns the 
    incorrect order. e.g. "10B" sorts before "2H".
    """
    isotopes = [util_misc.normalize_isotope_name(name) for name in isotopes]

    # Python sorts them perfectly as long as it knows that the number
    # is actually a number, so I turn the list of strings into a list of
    # tuples of (int, string).
    for i, isotope in enumerate(isotopes):
        # My regex matching is simplified by the facts that (a) this function
        # requires valid names, so I know that the regex will find a match,
        # and (b) I normalized them above so I know that the regex will 
        # return match groups of (number, letter, None, None).
        isotope = _ISOTOPE_REGEX.match(isotope).groups()[:2]
        isotopes[i] = (int(isotope[0]), isotope[1])

    # Now sorting is trivial
    isotopes.sort()

    return ["%d%s" % isotope for isotope in isotopes]

