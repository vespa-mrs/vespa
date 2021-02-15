# Python modules

# 3rd party modules

# Our modules
import vespa.common.util.db as common_util_db



class Database(common_util_db.Database):
    """
    A database class that looks and works just like common.util.db.Database.
    This exists to provide methods that are specific to Analysis. In
    addition, this caches metabolites by default. (If you don't know what
    that means, don't worry about it.)

    """
    def __init__(self, filename, cache_metabolites=True):
        common_util_db.Database.__init__(self, filename, cache_metabolites)
        


