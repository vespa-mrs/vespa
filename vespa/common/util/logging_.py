# Python modules
import logging
import os

# 3rd party modules

# Our modules
import vespa.common.util.time_ as util_time
import vespa.common.util.misc as util_misc

# Set up the app-wide logger
logger = logging.getLogger()
logger.setLevel(logging.INFO)


class Log(object):
    # These constants are arbitrary and may change.
    # However, bool(NONE) is always guaranteed to be False.
    NONE = None
    IMPORT = "import"
    DATABASE = "database"
    
# Below are handler classes for different logs. Handlers can do a bunch of
# things; these handlers mostly use the default behavior of 
# logging.FileHandler. They mainly exist to control the name of the log 
# file.


class NullHandler(logging.FileHandler):
    """Logs to /dev/null. Handy for when you need a handler that won't 
    write anything."""
    def __init__(self):
        self._filename = os.devnull        
        logging.FileHandler.__init__(self, self._filename)

    def __get_filename(self): return self._filename
    filename = property(__get_filename)



class ImportFileHandler(logging.FileHandler):
    def __init__(self):
        # Create a unique filename for this import.
        self._filename = "vespa_import."
        self._filename += util_time.filename_timestamp()
        self._filename += ".txt"
        
        self._filename = os.path.join(util_misc.get_data_dir(), "logs", 
                                      self._filename)
        
        logging.FileHandler.__init__(self, self._filename)

        formatter = logging.Formatter('%(asctime)s %(levelname)-8s: %(message)s')
        self.setFormatter(formatter)

    def __get_filename(self): return self._filename
    filename = property(__get_filename)
    

class DatabaseFileHandler(logging.FileHandler):
    def __init__(self):
        # Create a filename for the database log.
        self._filename = "vespa_database."
        self._filename += util_time.filename_timestamp()
        self._filename += ".txt"
        
        self._filename = os.path.join(util_misc.get_data_dir(), "logs", 
                                      self._filename)
                
        logging.FileHandler.__init__(self, self._filename)

        formatter = logging.Formatter('%(asctime)s %(levelname)-8s: %(message)s')
        self.setFormatter(formatter)

    def __get_filename(self): return self._filename
    filename = property(__get_filename)
    
