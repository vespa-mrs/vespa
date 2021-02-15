# Python modules

import sqlite3 as sqlite
import decimal
import re
import xdrlib

# 3rd party modules
import numpy

# Our modules
import vespa.common.util.misc as util_misc
import vespa.common.util.time_ as util_time
import vespa.common.constants as constants
import vespa.common.mrs_experiment as mrs_experiment
import vespa.common.mrs_j_coupling as mrs_j_coupling
import vespa.common.mrs_metabolite as mrs_metabolite
import vespa.common.mrs_pulse_sequence as mrs_pulse_sequence
import vespa.common.mrs_simulation as mrs_simulation
import vespa.common.mrs_spin as mrs_spin
import vespa.common.experiment_preview as experiment_preview
import vespa.common.transform_kernel_preview as transform_kernel_preview
import vespa.common.rfp_pulse_design as rfp_pulse_design
import vespa.common.rfp_pulse_design_preview as rfp_pulse_design_preview
import vespa.common.rfp_transform as rfp_transform
import vespa.common.rfp_rf_result as rfp_rf_result
import vespa.common.rfp_transform_kernel as rfp_transform_kernel
import vespa.common.rfp_machine_specs as rfp_machine_specs
import vespa.common.rfp_master_parameters as rfp_master_parameters
#import vespa.common.rfp_pulse_project as rfp_pulse_project
#import vespa.common.pulse_project_preview as pulse_project_preview
#import vespa.common.rfp_transformation_factory as rfp_transformation_factory
#import vespa.common.rfp_result as rfp_result
#import vespa.common.rfp_create_slr as rfp_create_slr
#import vespa.common.rfp_create_hs as rfp_create_hs
#import vespa.common.rfp_create_gaussian as rfp_create_gaussian
#import vespa.common.rfp_create_randomized as rfp_create_randomized
#import vespa.common.rfp_create_import as rfp_create_import
#import vespa.common.rfp_interpolate_rescale as rfp_interpolate_rescale
#import vespa.common.rfp_root_reflection as rfp_root_reflection
#import vespa.common.rfp_optimal_control_nonselective as rfp_ocn
#import vespa.common.rfp_ocn_state as rfp_ocn_state
#import vespa.common.rfp_machine_settings as rfp_machine_settings
import vespa.common.mrs_prior as mrs_prior
import vespa.common.mrs_prior_metabolite as mrs_prior_metabolite
import vespa.common.util.fileio as util_fileio

#######################     Table of Contents     ########################
#
# This long file consists mostly of one class, the Database class which 
# is a sqlite.Connection object with lots of Vespa-specific methods added
# for manipulating the database.
#
# The Database class is organized into sections for "private" (i.e.
# internal use only) methods and public methods. 

# The public methods are further divided into methods for controlling
# transactions (just a few of these), common methods (not primarily used by
# just one Vespa app) and one category each for Simulation and RFPulse
# methods. Technically, there's nothing that restricts any Vespa app from
# using a Simulation or RFPulse method. But the methods in those sections are
# used primarily by the app after which they're named, and we tend to think of
# e.g. fetch_experiment() as a Simulation method even though it's also used 
# by Analysis.

# The private methods are similarly categorized. There's Simulation and
# RFPulse private methods, and the _execute method which is used by just
# about all the other methods.

# In addition to the Database class, this file contains a few module-level
# items that are not for public use. It also contains the _BetterRow class
# which is again, not for public use. The database code uses it heavily.



class _Fetch(object):
    # These constants are arbitrary and may change in the future.
    # However, bool(NONE) is always guaranteed to be False and 
    # bool(ANYTHING_ELSE) is always guaranteed to be True.
    NONE      = 0
    LAST_ID   = 1
    SINGLETON = 2
    ONE       = 3
    COLUMN    = 4
    ALL       = 5
    CURSOR    = 6


def _numpy_array_to_blob(data):
    """Given a numpy array, returns a BLOB ready for insertion into a 
    SQLite table. The BLOB contains the array represented in XDR 
    "Variable-Length Array" format:    
    http://tools.ietf.org/html/rfc4506#section-4.13
    
    The BLOB doesn't preserve shape or data type information; it's up
    to the caller to do that. The reasoning is that all BLOBs in a given 
    column will have the same shape and data type, and it wastes space to
    repeat the same information for each row.
    
    Note that this format is deliberately simple. It's also language neutral.
    It's convenient for reading into numpy, but it's also readable from C,
    Fortran, etc. in case someone wants to read our database without Python.
    """    
    data_type = util_fileio.get_data_type(data)
    
    # Flatten the array into a 1d list. ravel() is conceptually the same as 
    # flatten() but sometimes faster.
    data = data.ravel().tolist()

    p = xdrlib.Packer()

    if constants.DataTypes.is_complex(data_type):
        # XDR doesn't understand complex numbers so I expand these into 
        # (real, imaginary) pairs.
        data = util_fileio.expand_complexes(data)

    # Decide which packing function to use
    if data_type in (constants.DataTypes.BOOL, ):
        function = p.pack_bool
    if data_type in (constants.DataTypes.BYTE, constants.DataTypes.INT32):
        function = p.pack_int
    if data_type in (constants.DataTypes.INT64, ):
        function = p.pack_hyper
    elif data_type in (constants.DataTypes.FLOAT32, constants.DataTypes.COMPLEX64):
        function = p.pack_float
    elif data_type in (constants.DataTypes.FLOAT64, constants.DataTypes.COMPLEX128):
        function = p.pack_double

    p.pack_array(data, function)
    
    # I tried using zlib.compress() on this data, but the arrays I was 
    # working with were short (1-25 elements) and the compression rarely 
    # accomplished anything except to waste a little time.
    
    return sqlite.Binary(p.get_buffer())
    
 
def _blob_to_numpy_array(data, data_type, shape=None):
    """Given a SQLite BLOB written by _numpy_array_to_blob(), returns the
    numpy array contained in that BLOB. If the caller provides a shape tuple,
    the returned array is shaped appropriately.
    """
    # FYI, Python's sqlite3 module translates BLOB columns into Python's 
    # rarely-used "buffer" type.
    
    # This code is similar to util_fileio.decode_xdr(), but it's simpler.
    # This code reads XDR variable-length arrays in which the length of the
    # array is encoded in the XDR. decode_xdr() reads fixed-length XDR arrays
    # in which one must know the length of the array in advance.
    data_type = constants.DataTypes.any_type_to_internal(data_type)
    
    p = xdrlib.Unpacker(data)

    if data_type in (constants.DataTypes.FLOAT32, constants.DataTypes.COMPLEX64):
        unpack_function = p.unpack_float
    elif data_type in (constants.DataTypes.FLOAT64, constants.DataTypes.COMPLEX128):
        unpack_function = p.unpack_double
    elif data_type in (constants.DataTypes.BOOL, ):
        unpack_function = p.unpack_bool
    elif data_type in (constants.DataTypes.BYTE, constants.DataTypes.INT32,):
        unpack_function = p.unpack_int
    elif data_type in (constants.DataTypes.INT64,):
        unpack_function = p.unpack_hyper
    else:
        raise ValueError("Unknown data type '%s'" % data_type)    
    
    data = p.unpack_array(unpack_function)
    
    p.done()
    
    if constants.DataTypes.is_complex(data_type):
        # Knit the (real, imaginary) pairs back together into complex numbers.
        data = util_fileio.collapse_complexes(data)
    
    data_type = constants.DataTypes.any_type_to_numpy(data_type)
    data = numpy.fromiter(data, data_type)
    
    if shape:
        data.shape(shape)
    
    return data
    

class Database(sqlite.Connection):
    """A specialized class for connecting to our SQLite database.
    Most of the functions you'll need are of the form fetch_xxx() and
    return one or more objects (Metabolite, Experiment, etc.)

    Any function that takes an id as a parameter assumes the object exists
    and will fail gracelessly if passed a non-existent id.

    All of the fetch_xxx functions that return objects (e.g. fetch_experiment,
    fetch_metabolites) return fully-fledged objects. In other words, if you
    fetch an metabolite, the returned object will be populated with its
    associated spins and J couplings as well.
    """

    def __init__(self, filename, cache_metabolites=False):
        sqlite.Connection.__init__(self, filename,
                                   detect_types=sqlite.PARSE_DECLTYPES)
        self._filename = filename
        self.cache_metabolites = cache_metabolites

        # This puts us in autocommit mode.
        self.isolation_level = None
        
        # Register the adapters & converters for my custom objects.
        # Note that Python has built-in adapters & coverters for columns
        # with type "date" and "timestamp".
        # I was originally concerned with case-sensitivity here with regard
        # to column types and converter names. The Python doc is silent
        # on the topic, but in the code it's clear that converter names are
        # forced to upper case as are column type names (prior to comparison).
        # e.g. --
        # >>> sqlite3.register_converter("bunnies", lambda: None)
        # >>> sqlite3.converters.keys()
        # ['DATE', 'TIMESTAMP', 'BUNNIES']
        sqlite.register_adapter(bool, _adapt_boolean)
        sqlite.register_converter("BOOLEAN", _convert_boolean)

        # We use our own row class rather than the default provided by the
        # sqlite module.
        self.row_factory = _BetterRow

        # I explicitly turn off foreign key support. In SQLite prior to
        # 3.6.19, this is a no-op. In versions >= 3.6.19, this is probably
        # off by default but we want to make sure. I think our code is
        # FK-friendly (deletes are done in the correct order, etc.) but I'm
        # not sure and that's not well-tested, so I'd prefer to have FKs
        # turned off for everyone.
        self._execute("PRAGMA foreign_keys = OFF")

        # The metab cache always exists. When self.cache_metabolites is
        # False, the cache is always empty.
        self._metabolite_cache = { }

        if self.cache_metabolites:
            sql = """SELECT
                        count(*)
                     FROM
                        sqlite_master
                     WHERE
                        type = 'table' AND
                        tbl_name = 'metabolites'
                  """
            if self._execute(sql, None, _Fetch.SINGLETON):
                # The metabolites table exists.
                self._update_metabolite_cache()
            #else:
                # The metabs table doesn't exist. That's not an error because
                # this module can be used on a blank database.


    @property
    def filename(self): 
        return self._filename


    #####################     Private methods     #####################

    def _execute(self, sql, parameters=(), fetch=_Fetch.ALL):
        """Executes SQL and provides some help with the results.
        Parameters are --
        sql: the parameterized SQL to execute
        params: the params (if any) referenced by the SQL. Single params
        can be passed as-is, multiple params should be in a tuple or list.
        fetch: one of the values defined in the Fetch class. Here's what each
        one means.
            _Fetch.ALL is for queries that return zero-to-many rows with
                multiple columns (N x M).
                e.g. select * from employees

            _Fetch.ONE means that the function will return a single item. That
                item is either a row, or None if the query returns no results.
                It's appropriate for queries that return one row with
                multiple columns (1 x M).
                e.g. select * from employees where unique_id = 42

            _Fetch.COLUMN is a special case of FETCH_ALL. It's appropriate for
                queries that return many rows with one column (N x 1).
                e.g.  select name from employees

                In this case sqlite returns a list of single-item tuples
                which is more awkward than it needs to be. _Fetch.COLUMN
                reduces that two-dimensional list to a single dimension.

            _Fetch.SINGLETON means that the function will return a single
                value (1 x 1)
                e.g. select count(*) from employees

            _Fetch.LAST_ID is appropriate for INSERT statements where the
                id of the newly-inserted row is of interest. It returns
                the result of cursor.lastrowid.

            _Fetch.NONE is appropriate for functions that don't return any
                rows. Inserts, updates and deletes are good examples. (But
                see also LAST_ID for INSERT statements.)

            _Fetch.CURSOR returns a cursor object that you can use to iterate
                through results.
        """

        # If the caller passed a simple object, turn it into an iterable.
        if isinstance(parameters, str) or    \
           isinstance(parameters, int) or           \
           isinstance(parameters, int) or          \
           isinstance(parameters, float) or         \
           isinstance(parameters, complex):
            parameters = (parameters, )
        elif parameters is None:
            parameters = ()
        elif not util_misc.is_iterable(parameters):
            raise TypeError("'parameters' must be None, a string, number, or iterable (e.g. list or tuple)")
        #else:
            # parameters already what it needs to be
        cursor = self.cursor()

        cursor.execute(sql, parameters)

        # Retrieve column names so that returned rows can implement .keys()
        # It's perhaps a little unintuitive to generate that information here,
        # but it's very efficient since each row refers to the same list of
        # column names.
        column_names = ( )
        if cursor.description:
            column_names = [column_info[0] for column_info in cursor.description]

        result = None

        if fetch == _Fetch.NONE:
            # nothing to do
            pass
        elif fetch == _Fetch.LAST_ID:
            # return the id of the newly-inserted row.
            result = cursor.lastrowid
        elif fetch == _Fetch.SINGLETON:
            # return one object (a float, int, string, None, etc.)
            result = cursor.fetchone()
            if result:
                result = result[0]
        elif fetch == _Fetch.ONE:
            # return one row
            result = cursor.fetchone()
            if result:
                result._column_names = column_names
        elif fetch == _Fetch.ALL:
            # return many rows
            result = cursor.fetchall()
            for row in result:
                row._column_names = column_names
        elif fetch == _Fetch.COLUMN:
            # Reduce the list of lists to a simple list of values.
            result = cursor.fetchall()

            result = [row[0] for row in result]
        elif fetch == _Fetch.CURSOR:
            result = cursor
        else:
            raise ValueError

        return result



    #################   Public  Transaction  Methods    #################

    def begin_transaction(self):
        self.isolation_level = "DEFERRED"


    def commit(self):
        sqlite.Connection.commit(self)
        # This puts us in autocommit mode.
        self.isolation_level = None


    def rollback(self):
        self.rollback()
        # This puts us in autocommit mode.
        self.isolation_level = None


    ###########    Miscellaneous/Common  Public  Methods    #############
    ##################   (in alphabetic order, no less!)   #################

    # These are called from various parts of Vespa. They don't "belong" to
    # any one particular app.

    def execute_script_from_file(self, filename):
        """Executes the contents of a SQL script contained in a file. Note
        that this stubbornly refuses to respect transactions regardless of
        whether they are embedded in the script or invoked through the
        isolation_level attribute. It's a bad idea, therefore, to use this
        method for large changes.
        """
        self.executescript(open(filename, "r").read())


    def find_unique_name(self, an_object, detail):
        """Given an object and a detail string, returns a name guaranteed 
        to be unique (i.e. otherwise unused) in the database. The object must
        be one of the types in the COUNTING_FUNCTION_MAP defined in this
        method. Right now, those are --
            - rfp_pulse_design.PulseDesign
            - rfp_transform_kernel.TransformKernel
            - rfp_machine_specs.MachineSpecsTemplate 
            - mrs_metabolite.Metabolite
            - mrs_pulse_sequence.PulseSequence
            - mrs_experiment.Experiment
        
        The detail string is typically "clone" or "import" and is added to the
        name (if necessary) to make it unique.
        
        If object's name is already unique, this function returns it 
        unchanged.
        
        This only guarantees that the name is unique at the time the
        function completes. In a context where multiple threads or 
        processes call this function, the "unique" name might not be unique
        anymore by the time the caller saves it to the database. This isn't a
        practical concern in the current (Mar 2011) incarnation of Vespa.
        """
        # COUNTING_FUNCTION_MAP connects a type to the method that knows how
        # to count that type in the database.
        COUNTING_FUNCTION_MAP = { 
            rfp_pulse_design.PulseDesign : self.count_pulse_designs,
            rfp_transform_kernel.TransformKernel : self.count_transform_kernels,
            rfp_machine_specs.MachineSpecsTemplate : self.count_machine_specs_templates,
            mrs_metabolite.Metabolite : self.count_metabolites,
            mrs_pulse_sequence.PulseSequence : self.count_pulse_sequences,
            mrs_experiment.Experiment : self.count_experiments,
        }
        
        name = an_object.name
        date = util_time.now(util_time.ISO_DATE_FORMAT)

        # Finding a unique name is a generic algorithm but the function that
        # decides whether or not the name is unique is specific to the object
        # passed. Here I find the function that counts the names for the
        # type in question.
        counting_function = None
        for key in COUNTING_FUNCTION_MAP:
            if isinstance(an_object, key):
                counting_function = COUNTING_FUNCTION_MAP[key]
                break
                
        if not counting_function:
            # Developer error -- no counting function for this type is 
            # defined in the COUNTING_FUNCTION_MAP
            raise TypeError("Object must be a known type")
                
        while counting_function(name=name):
            if date in name:
                # The date is already in the name. Append a number too.
                # First, we see if the name already ends in _N, where N
                # is one or more digits (e.g. "_1" or "_99").
                match = re.search(r"(_\d+)$", name)
                if match:
                    digits = match.group(0)
                    # Remove the _N from the end of the name and prepare
                    # to append _(N + 1)
                    name = name[:-len(digits)]
                    digits = digits.replace("_", "")
                    i = int(digits) + 1
                else:
                    i = 1
                
                name = "%s_%d" % (name, i)
            else:
                name = "%s_%s_%s" % (name, detail, date)

        # Now the name is unique.
        
        return name


    def mark_public(self, mrs_objects):
        """Given a single instance or list of objects, sets the "is_public" 
        flag for each object in the list. The object must be one of the types
        in the TABLE_NAME_MAP defined in this method. Right now, those are --
            - rfp_pulse_design.PulseDesign
            - rfp_transform_kernel.TransformKernel
            - mrs_metabolite.Metabolite
            - mrs_pulse_sequence.PulseSequence
            - mrs_experiment.Experiment
            
        In addition, all relevant subobjects will be marked public, too.
        For instance, if mrs_objects contains an experiment, and the 
        experiment refers to a pulse sequence, and that pulse sequence 
        refers to a pulse project, then the experiment, pulse sequence, and
        pulse project will all be marked public. 
            
        If mrs_objects is a list, the objects don't all have to be of the 
        same type. The list can contain, for example, both experiments and
        metabs.

        Making something public is a permanent, one-way operation. There is
        no mark_private() or privatize().
        """
        # TABLE_NAME_MAP associates a type to the name of the table that must
        # be altered in order to mark it public.
        TABLE_NAME_MAP = { 
            rfp_pulse_design.PulseDesign        : "pulse_designs",
            rfp_transform_kernel.TransformKernel: "transform_kernels",
            mrs_metabolite.Metabolite           : "metabolites",
            mrs_pulse_sequence.PulseSequence    : "pulse_sequences",
            mrs_experiment.Experiment           : "experiments",
        }

        if not util_misc.is_iterable(mrs_objects):
            # The caller passed a single object; make it a list.
            mrs_objects = [mrs_objects]

        sql = """UPDATE
                    %s
                 SET
                    is_public = 1
                 WHERE id = ?"""
                 
        # For metabs, I keep track of which are affected so I can update
        # the cache all in one fell swoop.
        affected_metabolite_ids = [ ]

        for mrs_object in mrs_objects:
            table_name = ""
            # Find the table name associated with this type of object
            for key in TABLE_NAME_MAP:
                if isinstance(mrs_object, key):
                    table_name = TABLE_NAME_MAP[key]
                    break
                    
            if not table_name:
                # Developer error -- no table name for this type is 
                # defined in the TABLE_NAME_MAP
                raise TypeError("Vespa Internal Error - object must be a known type")

            # Mark the object as public
            self._execute(sql % table_name, mrs_object.id, _Fetch.NONE)

            # Mark sub-objects (referenced objects) as public, as appropriate
            if isinstance(mrs_object, mrs_experiment.Experiment):
                self.mark_public(mrs_object.metabolites)
                if mrs_object.pulse_sequence:
                    self.mark_public(mrs_object.pulse_sequence)  
            elif isinstance(mrs_object, mrs_pulse_sequence.PulseSequence):
                self.mark_public(mrs_object.pulse_projects)
            elif isinstance(mrs_object, mrs_metabolite.Metabolite):
                affected_metabolite_ids.append(mrs_object.id)


        if affected_metabolite_ids and self.cache_metabolites:
            # Ensure that the cache is up-to-date
            self._update_metabolite_cache(affected_metabolite_ids)


    #################   Analysis  Public  Methods    #################
    ##################   (in alphabetic order, no less!)   #################

    def fetch_prior(self, experiment_id, dim_id):
        """Given an experiment id and an id from the experiment_dims table, 
        constructs and returns an MrsPrior object based on that.
        """
        prior = mrs_prior.Prior()

        preview = self.fetch_experiment_preview(experiment_id)
    
        prior.source = 'experiment'
        prior.source_id = experiment_id
        prior.comment = preview.comment
        prior.nucleus = preview.isotope
    
        # Fetch the simulations associated with this dimension.
        # There will be one for each metab in the experiment. 
        simulations = self._fetch_simulations_by_dim(dim_id)

        for simulation in simulations:
            prior_metabolite = mrs_prior_metabolite.PriorMetabolite()
            prior_metabolite.dims = [simulation.metabolite.name] + simulation.dims
            prior_metabolite.ppms = simulation.ppms
            prior_metabolite.areas = simulation.areas
            prior_metabolite.phases = simulation.phases
            prior.metabolites[prior_metabolite.name] = prior_metabolite
            
        return prior

    #################   Analysis  Private  Methods    #################
    ##################   (in alphabetic order, no less!)   #################

    def _fetch_simulations_by_dim(self, dims_id):
        """Given an id from the experiment_dims table, returns a sorted 
        list of the simulations (Simulation objects) associated with 
        that dim. There will be one for each metabolite in the simulation.
        """
        sql = """SELECT
                   simulations.metabolite_id AS metabolite_id,
                   simulations.ppms AS ppms,
                   simulations.areas AS areas,
                   simulations.phases AS phases,
                   metabolites.name AS metabolite_name,
                   experiment_dims.dim1 AS dim1,
                   experiment_dims.dim2 AS dim2,
                   experiment_dims.dim3 AS dim3
                 FROM
                    simulations, experiment_dims, metabolites
                 WHERE
                    experiment_dims.id = ?                       AND
                    simulations.dims_id = experiment_dims.id     AND
                    metabolites.id = simulations.metabolite_id
                 ORDER BY
                    dim3, dim2, dim1, metabolite_name
              """
        rows = self._execute(sql, dims_id)

        # Cache the metabs referenced here in a dict
        metabolites = { }
        for row in rows:
            id_ = row["metabolite_id"]
            if id_ not in metabolites:
                metabolites[id_] = self.fetch_metabolite(id_)        

        simulations = [ ]
        for row in rows:
            # If we modify this row a little, we can pass it directly to 
            # the Simulation object's ctor.
            row = dict(row)
            row["metabolite"] = metabolites[row["metabolite_id"]]
            ndims = list(range(1, constants.RESULTS_SPACE_DIMENSIONS))
            row["dims"] = [row["dim%d" % i] for i in ndims]
            
            if row["ppms"] or row["areas"] or row["phases"]:
                row["ppms"] = _blob_to_numpy_array(row["ppms"], 
                                                   constants.DataTypes.FLOAT64)
                row["areas"] = _blob_to_numpy_array(row["areas"], 
                                                    constants.DataTypes.FLOAT64)
                row["phases"] = _blob_to_numpy_array(row["phases"], 
                                                    constants.DataTypes.FLOAT64)

            simulations.append(mrs_simulation.Simulation(row))
            
        return simulations




    #################   Simulation  Public  Methods        #################
    ##################   (in alphabetic order, no less!)   #################

    # I call these "Simulation" methods because they're mostly used by 
    # Simulation. However, they can be used by any Vespa app.

    def count_experiments(self, id_=None, name=None):
        """Given an id or a name, returns the number of experiments matching
        that criterion. Exactly one of id_ or name must be supplied.
        """
        if id_:
            column_name = "id"
            param = id_
        else:
            column_name = "name"
            param = name

        sql = """SELECT
                    count(*)
                 FROM
                    experiments
                 WHERE
                    %s = ?
              """ % column_name
        return self._execute(sql, param, fetch=_Fetch.SINGLETON)


    def count_metabolites(self, id_=None, name=None):
        """Given an id or a name, returns the number of metabolites matching
        that criterion. Exactly one of id_ or name must be supplied.
        """
        if id_:
            column_name = "id"
            param = id_
        else:
            column_name = "name"
            param = name

        sql = """SELECT
                    count(*)
                 FROM
                    metabolites
                 WHERE
                    %s = ?
              """ % column_name
        return self._execute(sql, param, fetch=_Fetch.SINGLETON)


    def count_pulse_sequences(self, id_=None, name=None):
        """Given an id or a name, returns the number of pulse sequences
        matching that criterion. Exactly one of id_ or name must be supplied.
        """
        if id_:
            column_name = "id"
            param = id_
        else:
            column_name = "name"
            param = name

        sql = """SELECT
                    count(*)
                 FROM
                    pulse_sequences
                 WHERE
                    %s = ?
              """ % column_name
        return self._execute(sql, param, _Fetch.SINGLETON)


    def fetch_b0_bins(self, only_in_use=False):
        """Returns a list of b0 values in the database. Values that fall into
        standard bins (as defined in the table b0_bins) are returned as
        the bin value (e.g. 64.1 is returned as 64.0). Other values are
        returned unaltered.

        If only_in_use is True, then only values from the experiments table
        are returned. When only_in_use is False, those values are combined
        with common bin values (from b0_bins).
        """
        if only_in_use:
            union = ""
        else:
            # add the full list of bins
            union = """ UNION
                            SELECT
                                center AS binned_b0
                            FROM
                                b0_bins
                    """

        # Subtle point -- when the UNION clause is present, the DISTINCT
        # clause applies to the UNION-ed query results rather than just the
        # result of the first query. Very convenient as that's exactly
        # what I want it to do.
        sql = """SELECT DISTINCT
                    IFNULL(center, b0) AS binned_b0
                 FROM
                    experiments
                 LEFT OUTER JOIN
                    b0_bins ON
                        (left  <= b0) AND
                        (right >= b0)
                 %s
                 ORDER BY
                    binned_b0;
              """ % union
        return self._execute(sql, None, _Fetch.COLUMN)


    def fetch_experiment(self, experiment_id):
        """Given an experiment id, returns the associated experiment object."""
        sql = """SELECT
                    *
                 FROM
                    experiments
                 WHERE
                    id = ?
              """
        row = self._execute(sql, experiment_id, fetch=_Fetch.ONE)

        experiment = mrs_experiment.Experiment(row)

        # Get the associated pulse sequence, if necessary.
        if row["pulse_sequence_id"]:
            experiment.pulse_sequence = \
                        self.fetch_pulse_sequence(row["pulse_sequence_id"])

            # Fetch the experiment's params. An experiment won't have params
            # if there's no associated pulse sequence.
            sql = """SELECT
                        value
                     FROM
                        experiment_user_static_parameters
                     WHERE
                        experiment_id = ?
                     ORDER BY
                        display_order
                  """
            experiment.user_static_parameters = \
                            self._execute(sql, experiment.id, _Fetch.COLUMN)

        # Get the metabs
        metabolites = self.fetch_metabolites_by_experiment(experiment.id)

        experiment.metabolites = metabolites

        # Create a dict of this experiment's metabs keyed by metab id. (I'll
        # need it in a minute.)
        metabolites = dict( [(metabolite.id, metabolite) for metabolite 
                                                         in metabolites] )

        # Get the dims
        sql = """SELECT
                    dim1, dim2, dim3
                 FROM
                    experiment_dims 
                 WHERE
                    experiment_id = ?
                 ORDER BY 
                    dim3, dim2, dim1
              """
        rows = self._execute(sql, experiment_id)

        # Here we use the nifty Python idiom zip(*rows) for matrix 
        # transposition as demonstrated in the example below:
        # >>> rows = [  (1, 2, 3), (1, 2, 3), (1, 2, 3), (1, 2, 3), ]
        # >>> print zip(*rows)
        # [(1, 1, 1, 1), (2, 2, 2, 2), (3, 3, 3, 3)]
        dims = list(zip(*rows))
        # Filter out duplicates & sort what's left
        experiment.dims = [sorted(list(set(dim))) for dim in dims]

        # Fetch simulations for this experiment
        sql = """SELECT
                   simulations.started AS started,
                   simulations.completed AS completed,
                   simulations.metabolite_id AS metabolite_id,
                   simulations.ppms AS ppms,
                   simulations.areas AS areas,
                   simulations.phases AS phases,
                   metabolites.name AS metabolite_name,
                   experiment_dims.dim1 AS dim1,
                   experiment_dims.dim2 AS dim2,
                   experiment_dims.dim3 AS dim3
                 FROM
                    simulations, experiment_dims, metabolites
                 WHERE
                    experiment_dims.experiment_id = ?            AND
                    simulations.dims_id = experiment_dims.id     AND
                    metabolites.id = simulations.metabolite_id
                 ORDER BY
                    dim3, dim2, dim1, metabolite_name
              """
        rows = self._execute(sql, experiment.id)
        # results look like this --
        # id    dim1    dim2    dim3    metabolite_id
        # ----  -----   -----   -----   ------------------------------------
        # 169   0.01    1.0     0.00    a1b9f07b-3665-4ce8-ba4a-5f454baf9681
        # 170   0.01    1.2     0.00    f03efaa3-d686-41db-8df2-8a2e450f596b
        # 171   0.02    1.0     0.00    a1b9f07b-3665-4ce8-ba4a-5f454baf9681
        # 172   0.02    1.2     0.00    f03efaa3-d686-41db-8df2-8a2e450f596b
        # 173   0.03    1.0     0.00    a1b9f07b-3665-4ce8-ba4a-5f454baf9681
        # 174   0.03    1.2     0.00    f03efaa3-d686-41db-8df2-8a2e450f596b
        # 175   0.04    1.0     0.00    a1b9f07b-3665-4ce8-ba4a-5f454baf9681
        # 176   0.04    1.2     0.00    f03efaa3-d686-41db-8df2-8a2e450f596b
        # 176   0.05    1.0     0.00    a1b9f07b-3665-4ce8-ba4a-5f454baf9681
        # 176   0.05    1.2     0.00    f03efaa3-d686-41db-8df2-8a2e450f596b

        experiment.simulations = [ ]
        for row in rows:
            # If we modify this row a little, we can pass it directly to 
            # the Simulation object's ctor.
            row = dict(row)
            row["metabolite"] = metabolites[row["metabolite_id"]]
            ndims = list(range(1, constants.RESULTS_SPACE_DIMENSIONS))
            row["dims"] = [row["dim%d" % i] for i in ndims]
            
            if row["ppms"] or row["areas"] or row["phases"]:
                row["ppms"] = _blob_to_numpy_array(row["ppms"], 
                                                   constants.DataTypes.FLOAT64)
                row["areas"] = _blob_to_numpy_array(row["areas"], 
                                                    constants.DataTypes.FLOAT64)
                row["phases"] = _blob_to_numpy_array(row["phases"], 
                                                    constants.DataTypes.FLOAT64)

            experiment.simulations.append(mrs_simulation.Simulation(row))

        return experiment


    def fetch_experiment_dims(self, experiment_id):
        """Given an experiment id, returns all of the dims (pulse seq loop 
        values) associated with that experiment.
        
        The dims are returned in a dict that maps the database ids to the
        a tuple like (1.0, .33, 42.0).
        
        The dim tuple always has a length of RESULTS_SPACE_DIMENSIONS - 1
        (currently and possibly forever 3). For pulse seqs that don't use
        every loop, non-existent loop values are represented as 0.0.
        """
        sql = """SELECT
                    id, dim1, dim2, dim3
                 FROM
                    experiment_dims 
                 WHERE
                    experiment_id = ?
                 ORDER BY 
                    dim3, dim2, dim1
              """
        rows = self._execute(sql, experiment_id)
        
        # # The construct (row[1], row[2], row[3]) could be replaced by 
        # # row[1:] if our _BetterRow objects supported slicing.
        return dict( [(row["id"], (row[1], row[2], row[3])) for row in rows])


    def fetch_experiment_preview(self, experiment_id):
        """Given an experiment id, returns the associated ExperimentPreview
        object.
        """

        sql = """SELECT DISTINCT
                    experiments.id AS id,
                    experiments.name AS name,
                    experiments.is_public AS is_public,
                    experiments.comment AS comment,
                    experiments.b0 AS b0,
                    experiments.isotope AS isotope,
                    experiments.pulse_sequence_id AS pulse_sequence_id,
                    pulse_sequences.name AS pulse_sequence_name
                 FROM
                    experiments, pulse_sequences
                 WHERE
                    experiments.pulse_sequence_id = pulse_sequences.id AND
                    experiments.id = ?
              """

        row = self._execute(sql, experiment_id, _Fetch.ONE)
        
        row = dict(row)
        row["metabolites"] = self.fetch_metabolites_by_experiment(experiment_id)
        
        return experiment_preview.ExperimentPreview(row)
        


    def fetch_experiment_previews(self, b0=None, isotope=None):
        """Returns a list of experiment previews.

           The results are optionally filtered by the params b0 and isotope.
           If b0 is present, results are returned that are +/-1 int(b0).
        """
        # The query below is written in its simplest form which doesn't
        # filter on b0 or isotope. If either of those params are present,
        # the relevant SQL & params are added.
        where_clause = ""
        params = [ ]
        if b0:
            # Check to see if this b0 falls into a standard bin.
            sql = """SELECT
                        left, right
                     FROM
                        b0_bins
                     WHERE
                        ? >= left AND
                        ? <= right"""
            row = self._execute(sql, (b0, b0), _Fetch.ONE)
            if row:
                left, right = row
            else:
                # This doesn't fall into any standard bins. I create an
                # ad-hoc bin.
                drift = decimal.Decimal(str(constants.AD_HOC_B0_BIN_DRIFT))
                left = (decimal.Decimal(str(b0)) - drift).to_eng_string()
                right = (decimal.Decimal(str(b0)) + drift).to_eng_string()

            where_clause += " AND (b0 >= ? AND b0 <= ?) "
            params.append(left)
            params.append(right)

        if isotope:
            where_clause += """ AND experiments.isotope = ? """
            params.append(isotope)

        sql = """SELECT DISTINCT
                    experiments.id AS id,
                    experiments.name AS name,
                    experiments.is_public AS is_public,
                    experiments.comment AS comment,
                    experiments.b0 AS b0,
                    experiments.isotope AS isotope,
                    experiments.pulse_sequence_id AS pulse_sequence_id,
                    metabolites.id AS metabolite_id,
                    pulse_sequences.name AS pulse_sequence_name
                 FROM
                    experiments, experiment_metabolites, metabolites,
                    pulse_sequences
                 WHERE
                    experiment_metabolites.experiment_id = experiments.id AND
                    experiment_metabolites.metabolite_id = metabolites.id AND
                    experiments.pulse_sequence_id = pulse_sequences.id
                    %s
                 ORDER BY
                    experiments.name, experiments.id, metabolites.name
              """ % where_clause

        rows = self._execute(sql, params)

        # I get back N rows per experiment where N is the # of unique
        # metabolites. The experiment data repeats. I order the results by
        # both id and name so that if names aren't unique, the results will
        # still be grouped correctly.

        # Sample output below with some columns removed and id truncated.

        # id        name                       metabolite_name
        # --------  -------------------------  -----------------
        # 14f32a89  Test PRESS Ideal no loops  choline-truncated
        # 14f32a89  Test PRESS Ideal no loops  aspartate
        # 14f32a89  Test PRESS Ideal no loops  lactate
        # 14f32a89  Test PRESS Ideal no loops  glutamate
        # 14f32a89  Test PRESS Ideal no loops  glutamine
        # 14f32a89  Test PRESS Ideal no loops  creatine
        # 14f32a89  Test PRESS Ideal no loops  myo-inositol
        # 14f32a89  Test PRESS Ideal no loops  n-acetylaspartate
        # 32a69fd8  Example OnePulse Data      choline-truncated
        # 32a69fd8  Example OnePulse Data      aspartate
        # 32a69fd8  Example OnePulse Data      lactate
        # 32a69fd8  Example OnePulse Data      glutamate
        # 32a69fd8  Example OnePulse Data      glutamine
        # 32a69fd8  Example OnePulse Data      creatine
        # 32a69fd8  Example OnePulse Data      myo-inositol
        # 32a69fd8  Example OnePulse Data      n-acetylaspartate
        # 32a69fd8  Example OnePulse Data      gaba
        # 32a69fd8  Example OnePulse Data      taurine
        previews = [ ]
        # I create an empty preview to make the first iteration of the loop
        # work properly.
        current_preview = experiment_preview.ExperimentPreview()
        for row in rows:
            if current_preview.id != row["id"]:
                # This is an experiment I haven't seen before.
                current_preview = experiment_preview.ExperimentPreview(row)
                previews.append(current_preview)

            metabolite = self.fetch_metabolite(row["metabolite_id"])
            current_preview.metabolites.append(metabolite)

        return previews


    def check_experiment_id(self, experiment_id):
        """Given an experiment id, returns True if an experiment with that
        id exists in the database. Returns False is not.
        
        """
        sql = """SELECT
                    id, dim1
                 FROM
                    experiment_dims 
                 WHERE
                    experiment_id = ?
              """
        rows = self._execute(sql, experiment_id)
    
        if not rows: 
            return False
        else:
            return True
        


    def fetch_isotopes(self, only_in_use=False):
        """Returns a sorted list of isotopes with no duplicates. The param
        only_in_use controls whether the returned list contains only isotopes
        that are in use or all known isotopes.
        """
        if only_in_use:
            sql = """SELECT DISTINCT
                        isotope
                     FROM
                        metabolite_spins, isotopes
                     WHERE
                        isotope = isotopes.name
                     ORDER BY
                        isotopes.display_order
                  """
        else:
            sql = """SELECT
                        name
                     FROM
                        isotopes
                     ORDER BY
                        display_order
                  """

        return self._execute(sql, fetch=_Fetch.COLUMN)


    def fetch_metabolite(self, metabolite_id, include_experiment_names=False):
        """Returns the metabolite with the given id.

        If include_experiment_names is True, the metabolite will be given
        an experiment_names object which is a list of the names of the
        experiments that reference the metabolite.
        """
        if self.cache_metabolites:
            # I fetch the metabolite from the cache.
            metabolite = self._metabolite_cache.get(metabolite_id, None)
        else:
            # We're not using the cache.
            metabolite = self._fetch_metabolite(metabolite_id)

        if metabolite and include_experiment_names:
            sql = """SELECT
                        experiments.name
                     FROM
                        experiment_metabolites, experiments
                     WHERE
                        experiment_metabolites.metabolite_id = ? AND
                        experiment_metabolites.experiment_id = experiments.id
                     ORDER BY
                        experiments.name, experiments.id
                  """
            metabolite.experiment_names = self._execute(sql, metabolite.id,
                                                        _Fetch.COLUMN)

        return metabolite


    def fetch_metabolites(self, isotope=None, include_deactivated=True,
                          include_experiment_names=False):
        """Returns a list of 0 or more metabolite objects sorted by name.
        They can optionally be filtered by isotope and whether or not they've
        been taken out of service.

        include_experiment_names has the same effect as for
        fetch_metabolite().
        """
        # I select the metabs from the database in order to let the database
        # do filtering and sorting for me. However, I select only the ids
        # and then I use them to get the full metab objects out of my cache.
        params = [ ]
        where_clause = ""
        from_clause = ""
        if not include_deactivated:
            where_clause += " AND deactivated IS NULL "

        if isotope:
            where_clause += """
                    AND metabolite_spins.metabolite_id = metabolites.id AND
                    metabolite_spins.isotope = ?
                            """
            from_clause += ", metabolite_spins "
            params.append(isotope)

        # The metab ordering here needs to be kept in sync with that of the
        # mrs_metabolite.Metabolite.__lt__(). They're currently both
        # case-insensitive.
        sql = """SELECT DISTINCT
                    metabolites.id AS id
                 FROM
                    metabolites %s
                 WHERE
                    1 = 1
                    %s
                 ORDER BY
                    lower(name)
               """ % (from_clause, where_clause)

        ids = self._execute(sql, params, _Fetch.COLUMN)

        return [self.fetch_metabolite(id_, include_experiment_names) for id_ in ids]


    def fetch_metabolites_by_experiment(self, experiment_id):
        sql = """SELECT
                    id
                 FROM
                    metabolites 
                 INNER JOIN
                    experiment_metabolites
                 ON
                    metabolites.id = experiment_metabolites.metabolite_id AND
                    experiment_metabolites.experiment_id = ?
                 ORDER BY 
                    name
              """
        ids = self._execute(sql, experiment_id, _Fetch.COLUMN)

        return [self.fetch_metabolite(id_) for id_ in ids]
        

    def fetch_metabolites_by_name(self, name, include_deactivated=True):
        """Returns a sorted list of the metabolites with a given name."""
        where = ""
        if not include_deactivated:
            where = " AND deactivated IS NULL "

        sql = """SELECT
                    id
                 FROM
                    metabolites
                 WHERE
                    name = ?
                    %s
                 ORDER BY
                    metabolites.created
              """ % where

        ids = self._execute(sql, name, _Fetch.COLUMN)

        return [self.fetch_metabolite(id_) for id_ in ids]



    # The functions fetch_pulse_sequence() and fetch_pulse_sequences() are
    # wrappers that do very little code-wise, but they keep the database API
    # for pulse sequences consistent with the APIs for other objects (i.e.
    # one function that accepts an id as a param and returns one object,
    # another function that accepts no params and returns a list of all
    # such objects).
    def fetch_pulse_sequence(self, pulse_sequence_id,
                             include_experiment_names=False):
        """Returns the pulse sequence with the given id.

        If include_experiment_names is True, the returned pulse sequence
        gains an additional attribute called experiment_names. It's a list
        (possibly empty) of the names of the experiments that use the
        pulse sequence.
        """
        return self._fetch_pulse_sequences(id_=pulse_sequence_id,
                       include_experiment_names=include_experiment_names)[0]


    def fetch_pulse_sequences(self, include_experiment_names=False):
        """Returns a list of all pulse sequences.

        If include_experiment_names is True, the returned pulse sequences
        gain an additional attribute called experiment_names. It's a list
        (possibly empty) of the names of the experiments that use the
        pulse sequence.
        """

        return self._fetch_pulse_sequences(include_experiment_names=include_experiment_names)


    def fetch_pulse_sequences_by_name(self, name,
                                      include_experiment_names=False):
        """Returns a list of all pulse sequences.

        If include_experiment_names is True, the returned pulse sequences
        gain an additional attribute called experiment_names. It's a list
        (possibly empty) of the names of the experiments that use the
        pulse sequence.
        """
        return self._fetch_pulse_sequences(name=name,
                            include_experiment_names=include_experiment_names)


    def insert_experiment(self, experiment, use_transaction=True):
        """Given a populated experiment object, creates database rows
        representing the new experiment and its associated simulations
        and spectra.

        If the id is not populated, an id will be assigned.
        
        If use_transaction=True, the database will wrap the inserts in a 
        transaction. There's no reason to set this to False unless you're
        calling this from code which already has a transaction open. 

        It is assumed that the referenced pulse sequence (if any) and
        metabolites exist, that all fields in the experiment are valid and
        that database constraints (e.g. unique name) will not be violated.

        When this function completes, the id fields of the simulations
        are populated.
        """
        if use_transaction:
            self.begin_transaction()

        if not experiment.id:
            experiment.id = util_misc.uuid()

        if not experiment.created:
            experiment.created = util_time.now()

        pulse_sequence_id = experiment.pulse_sequence.id \
                                        if experiment.pulse_sequence else None

        sql = """INSERT INTO
                    experiments
                        (id, name, is_public, created, investigator, comment,
                         b0, isotope,
                         peak_search_ppm_low, peak_search_ppm_high,
                         blend_tolerance_ppm, blend_tolerance_phase,
                         pulse_sequence_id)
                 VALUES
                    (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
              """
        params = (experiment.id, experiment.name, experiment.is_public,
                  experiment.created,
                  experiment.investigator, experiment.comment,
                  experiment.b0, experiment.isotope,
                  experiment.peak_search_ppm_low,
                  experiment.peak_search_ppm_high,
                  experiment.blend_tolerance_ppm,
                  experiment.blend_tolerance_phase,
                  pulse_sequence_id
                 )
        self._execute(sql, params, _Fetch.NONE)

        # Add parameters
        self._insert_experiment_user_static_parameters(experiment.id, 
                                           experiment.user_static_parameters)

        # Add metab references
        sql = """INSERT INTO
                    experiment_metabolites
                        (experiment_id, metabolite_id)
                 VALUES
                    (?, ?)
              """
        for metabolite in experiment.metabolites:
            self._execute(sql, (experiment.id, metabolite.id), _Fetch.NONE)

        # Add dimensions
        sql = """INSERT INTO
                    experiment_dims
                        (experiment_id, dim1, dim2, dim3)
                 VALUES
                    (?, ?, ?, ?)
              """
        for dims in experiment.all_indices:
            self._execute(sql, [experiment.id] + list(dims))

        self._insert_simulations(experiment.id, experiment.simulations)

        if use_transaction:
            self.commit()

        return experiment.id


    def insert_metabolite(self, metabolite):
        """Given a populated metabolite object, creates database rows
        representing the new metabolite and its associated spins and
        J Couplings.

        If the id is not populated, an id will be assigned.

        It is assumed that all fields in the metabolite are valid and
        that database constraints (e.g. unique metabolite name) will not
        be violated.

        When this function completes, the id fields of the metabolite,
        spins & j_couplings are populated.
        """
        self.begin_transaction()

        if not metabolite.id:
            metabolite.id = util_misc.uuid()

        if not metabolite.created:
            metabolite.created = util_time.now()

        sql = """INSERT INTO
                    metabolites (id, name, is_public, created, creator,
                                 comment, deactivated)
                 VALUES
                    (?, ?, ?, ?, ?, ?, ?)
              """
        params = (metabolite.id, metabolite.name, metabolite.is_public,
                  metabolite.created,
                  metabolite.creator, metabolite.comment,
                  metabolite.deactivated)

        self._execute(sql, params, _Fetch.NONE)

        # Add spins
        sql = """INSERT INTO
                    metabolite_spins (metabolite_id, isotope, chemical_shift,
                                      display_order)
                 VALUES
                    (?, ?, ?, ?)
              """
        for i, spin in enumerate(metabolite.spins):
            params = (metabolite.id, spin.isotope, spin.chemical_shift, i)

            spin.id = self._execute(sql, params, _Fetch.LAST_ID)

        # Add J couplings
        sql = """INSERT INTO
                    j_couplings (value, spin1_id, spin2_id)
                 VALUES
                    (?, ?, ?)
              """
        for j_coupling in metabolite.j_couplings:
            params = (j_coupling.value, j_coupling.spin1.id, j_coupling.spin2.id)

            j_coupling.id = self._execute(sql, params, _Fetch.LAST_ID)

        self.commit()

        if self.cache_metabolites:
            self._update_metabolite_cache()

        return metabolite.id


    def insert_pulse_sequence(self, pulse_sequence):
        """Saves the pulse_sequence described by the param. The
        pulse_sequence must not exist in the database.

        If the id is not populated, an id will be assigned.

        It is assumed that all fields in the pulse sequence are valid and
        that database constraints (e.g. unique name) will not be violated.

        When this function completes, the id fields of the parameters
        are populated.
        """
        self.begin_transaction()

        if not pulse_sequence.id:
            pulse_sequence.id = util_misc.uuid()

        if not pulse_sequence.created:
            pulse_sequence.created = util_time.now()

        # Insert the pulse_sequence
        sql = """INSERT INTO
                    pulse_sequences
                        (id, name, is_public, created, creator, comment,
                         sequence_code, binning_code)
                 VALUES
                    (?, ?, ?, ?, ?, ?, ?, ?)
              """
        params = (pulse_sequence.id, pulse_sequence.name,
                  pulse_sequence.is_public,
                  pulse_sequence.created, pulse_sequence.creator,
                  pulse_sequence.comment,
                  pulse_sequence.sequence_code, pulse_sequence.binning_code)
        self._execute(sql, params, _Fetch.NONE)

        # Insert the loops
        sql = """INSERT INTO
                    pulse_sequence_loops
                        (pulse_sequence_id, label, display_order)
                 VALUES
                    (?, ?, ?)
              """
        for i, label in enumerate(pulse_sequence.loop_labels):
            if label:
                params = (pulse_sequence.id, label, i)
                self._execute(sql, params, _Fetch.NONE)

        # Insert the params
        sql = """INSERT INTO
                    pulse_sequence_user_static_parameters
                        (pulse_sequence_id, type, name, default_value,
                         display_order)
                 VALUES
                    (?, ?, ?, ?, ?)
              """
        for i, parameter in enumerate(pulse_sequence.user_static_parameters):
            params = (pulse_sequence.id, parameter.type, parameter.name,
                      parameter.default, i)

            parameter.id = self._execute(sql, params, _Fetch.LAST_ID)

        # Insert the pulse_design references
        sql = """INSERT INTO
                    pulse_sequence_pulse_designs
                        (pulse_sequence_id, pulse_design_id, progression)
                 VALUES
                    (?, ?, ?)
              """
        for i, pulse_design in enumerate(pulse_sequence.pulse_projects):
            self._execute(sql, (pulse_sequence.id, pulse_design.id, i))
            
        self.commit()


    #################   Simulation  Private  Methods    #################
    ##################   (in alphabetic order, no less!)   #################

    def _fetch_metabolite(self, metabolite_id, metabolite=None):
        # An internal-use only function for fetching the metabolite from
        # the database. If the metabolite param is not None, this function
        # doesn't create a new metab to hold the values fetched from the
        # database but uses the metab in the parameter.

        # I construct a metabolite from three queries. First I get the basic
        # metab info, then the spins, and then the J coupling data.
        sql = """SELECT
                    *
                 FROM
                    metabolites
                 WHERE
                    id = ?
              """
        row = self._execute(sql, metabolite_id, _Fetch.ONE)
        
        assert(bool(row))
        
        # If the caller didn't pass a metab to populate, I create a new one.
        if metabolite:
            metabolite.inflate(row)
        else:
            metabolite = mrs_metabolite.Metabolite(row)

        # Fetch spins
        sql = """SELECT
                    *
                 FROM
                    metabolite_spins
                 WHERE
                    metabolite_id = ?
                 ORDER BY
                    display_order
              """
        rows = self._execute(sql, metabolite_id)
        
        metabolite.spins = [mrs_spin.Spin(row) for row in rows]

        # OK, now I have the basic metab info + spins.

        # I build a dict of my spin objects, keyed by id. This will come in
        # handy twice. First, it gives me an easy way to provide the ids
        # I'm looking for to the query below. Second, when I get my results
        # back and I want to associate each JCoupling object with the
        # right spins, this dict gives me an easy way to look them up.
        spins = dict([ (spin.id, spin) for spin in metabolite.spins ])

        # Here I build a list of the param placeholders that I feed to
        # the SQL below.
        placeholders = _build_placeholders(len(spins))

        sql = """SELECT
                    j_couplings.id AS id, value, spin1_id, spin2_id
                 FROM
                    j_couplings, metabolite_spins MS1, metabolite_spins MS2
                 WHERE
                    (j_couplings.spin1_id IN (%s) OR
                     j_couplings.spin2_id IN (%s)
                    )                                   AND
                    MS1.id = j_couplings.spin1_id       AND
                    MS2.id = j_couplings.spin2_id
                 ORDER BY
                    MS1.display_order, MS2.display_order
              """ % (placeholders, placeholders)
        rows = self._execute(sql, (list(spins.keys()) + list(spins.keys())) )

        metabolite.j_couplings = [ ]
        for row in rows:
            j_coupling = mrs_j_coupling.JCoupling(row)
            # Here's where I make use of the dictionary I built above.
            j_coupling.spin1 = spins[row["spin1_id"]]
            j_coupling.spin2 = spins[row["spin2_id"]]

            metabolite.j_couplings.append(j_coupling)

        return metabolite


    def _fetch_pulse_sequences(self, id_=None, name=None,
                               include_experiment_names=False):
        """Returns a list of all pulse sequences, or a single-item list
        if the id_ param is populated.

        If include_experiment_names is True, the returned pulse sequences
        gain an additional attribute called experiment_names. It's a list
        (possibly empty) of the names of the experiments that use the
        pulse sequence.
        """
        # I fetch the pulse sequence[s] in three queries -- one for the
        # pulse sequence, one for its params, and one for its loop values. 
        # If include_experiment_names is True, I run a fourth query to 
        # fetch them.
        where = ""
        params = None
        if id_:
            where = " AND pulse_sequences.id = ? "
            params = id_
        if name:
            where = " AND pulse_sequences.name = ? "
            params = name

        # The SQL below contains a dummy WHERE condition to ensure that the
        # query is always syntactically valid even if I don't add anything 
        # to it.
        sql = """SELECT
                    pulse_sequences.*
                 FROM
                    pulse_sequences
                 WHERE
                    1 = 1  %s
                 ORDER BY
                    name
                 """ % where
        rows = self._execute(sql, params)
        
        pulse_sequences = \
                    [mrs_pulse_sequence.PulseSequence(row) for row in rows]
        
        for pulse_sequence in pulse_sequences:
            # Fetch pulse sequence params. (There might not be any.)
            sql = """SELECT
                        *
                     FROM
                        pulse_sequence_user_static_parameters
                     WHERE
                        pulse_sequence_id = ?
                     ORDER BY
                        display_order
                     """
            rows = self._execute(sql, pulse_sequence.id)

            # Some special attention is required here. Unlike most of our
            # tables/objects, there's not a 1:1 map between the column and 
            # attribute names for pulse sequence params. The table 
            # pulse_sequence_user_static_parameters contains a column called 
            # default_value while the UserStaticParameter class has an 
            # attribute called default. We can't use the latter as a column 
            # name because it's a reserved word in SQL, and neither can we 
            # use it as a column name alias. 
            # So here I rename the column in the returned row. Remember 
            # that the rows returned by self._execute() are _BetterRow objects
            # and need to be turned into proper Python dicts before I 
            # manipulate them.
            rows = [dict(row) for row in rows]
            
            for row in rows:
                row["default"] = row["default_value"]
                del row["default_value"]

            pulse_sequence.user_static_parameters = \
                [mrs_pulse_sequence.UserStaticParameter(row) for row in rows]
            
            # Grab the loops. (There might not be any.)
            sql = """SELECT
                        label
                     FROM
                        pulse_sequence_loops
                     WHERE
                        pulse_sequence_id = ?
                     ORDER BY
                        display_order
                    """
            pulse_sequence.loop_labels = self._execute(sql, pulse_sequence.id,
                                                       _Fetch.COLUMN) 
            
            # Fetch pulse designs. (There might not be any.)
            sql = """SELECT
                        pulse_design_id
                     FROM
                        pulse_sequence_pulse_designs
                     WHERE
                        pulse_sequence_id = ?
                     ORDER BY
                        progression
                  """
            ids = self._execute(sql, pulse_sequence.id, _Fetch.COLUMN)
         
            pulse_designs = [self.fetch_pulse_design(id_) for id_ in ids]
            
            pulse_sequence.pulse_projects = pulse_designs

            if include_experiment_names:
                sql = """SELECT
                            name
                         FROM
                            experiments
                         WHERE
                            pulse_sequence_id = ?
                         ORDER BY
                            name"""
                pulse_sequence.experiment_names = \
                    self._execute(sql, pulse_sequence.id, _Fetch.COLUMN)

        return pulse_sequences


    def _insert_experiment_user_static_parameters(self, experiment_id, 
                                                  parameters):
        """Given an experiment id and an iterable of parameter strings,
        attaches those parameters to the experiment.
        
        This function does not use a transaction. The caller should.
        """
        # Ensure that the experiment id exists.
        assert(bool(self.count_experiments(experiment_id)))

        sql = """INSERT INTO
                    experiment_user_static_parameters 
                        (experiment_id, value, display_order)
                 VALUES
                    (?, ?, ?)
              """
        for i, parameter in enumerate(parameters):
            self._execute(sql, (experiment_id, parameter, i), _Fetch.NONE)


    def _insert_simulations(self, experiment_id, simulations):
        """Given an experiment id and an iterable of mrs_simulation objects,
        attaches those simulations to the experiment.
        
        This function does not use a transaction. The caller should.
        """
        # Most of our code assumes that the caller "gets it right" and 
        # doesn't attempt nonsensical actions such as adding simulations
        # to an experiment id that isn't in the database. However, this
        # function gets special treatment in the form of a check to ensure 
        # that the experiment id exists.
        assert(bool(self.count_experiments(experiment_id)))
        
        # There will be one simulation associated with each dimension 
        # for each metab (e.g. 7 metabs = 7 sims for each dimension). Rather
        # than looking up the dims_id repeatedly, I cache them all here.
        sql = """SELECT
                    id, dim1, dim2, dim3
                 FROM
                    experiment_dims
                 WHERE
                    experiment_id = ?
              """
        rows = self._execute(sql, experiment_id)
        
        # Reorganize rows as a list of 2-tuples of ( (dims), dims_id )
        rows = [((row["dim1"], row["dim2"], row["dim3"]), row["id"]) for
                                                                row in rows]
                                                                
        # Now we can make a dict keyed by the dims tuple which points to the
        # id associated with the dim.
        dims_map = dict(rows)
        
        sql = """INSERT INTO
                    simulations 
                        (metabolite_id, dims_id, started, completed, 
                         ppms, areas, phases)
                 VALUES
                    (?, ?, ?, ?, ?, ?, ?)
              """
        for simulation in simulations:
            dims_id = dims_map[tuple(simulation.dims)]
            
            params = [simulation.metabolite.id, dims_id,
                      simulation.started, simulation.completed]
            
            if simulation.ppms.size or simulation.areas.size or \
               simulation.phases.size:
                params.append(_numpy_array_to_blob(simulation.ppms))
                params.append(_numpy_array_to_blob(simulation.areas))
                params.append(_numpy_array_to_blob(simulation.phases))
            else:
                params += [None, None, None]
            
            self._execute(sql, params, _Fetch.NONE)


    def _update_metabolite_cache(self, update_ids=[ ]):
        # Updates the metabolite cache. Update_ids can be a list or a
        # single metabolite id. If any update ids are supplied, the data
        # for those metabs is re-fetched from the database. This is how
        # metab edits are handled.
        # If no ids are supplied, all metabs in the database but not in the
        # cache are added to the cache. This is how the cache is made current
        # on initialization, imports and when new metabs are added via the
        # GUI.
        # Deleted metabs are removed from the cache manually.

        if isinstance(update_ids, str):
            update_ids = [update_ids]

        if not update_ids:
            # The caller doesn't know what to update, which means new metabs
            # were added. This also applies during init when all metabs
            # are "new" to the cache.
            # In this case I examine the database to see which ids are not
            # in my cache.
            # In case the cache is empty, I create one fake id in the list
            # that's guaranteed to not match anything.
            cache_ids = list(self._metabolite_cache.keys()) or ["NULL"]
            placeholders = _build_placeholders(len(cache_ids))
            sql = """SELECT
                        id
                     FROM
                        metabolites
                     WHERE
                        id NOT IN (%s)
                  """ % placeholders
            update_ids = self._execute(sql, cache_ids, _Fetch.COLUMN)

        for id_ in update_ids:
            # Subtle but important -- if a metab is already represented in
            # the cache, I want to update that object rather than replace
            # it. This ensures that all references to that metab refer to
            # the updated (correct) version.
            metabolite = self._metabolite_cache.get(id_, None)
            self._metabolite_cache[id_] = self._fetch_metabolite(id_, metabolite)



   #  #################    RFPulse  Public  Methods          #################
   #  ##################   (in alphabetic order, no less!)   #################
   #
   #  # I call these "RFPulse" methods because they're mostly used by
   #  # RFPulse. However, they can be used by any Vespa app.
   #
   #  def count_machine_settings_templates(self, id_=None, name=None):
   #      """Given an id or a name, returns the number of machine settings
   #      templates  matching that criterion. Exactly one of id_ or name must be
   #      supplied.
   #      """
   #      if id_:
   #          column_name = "id"
   #          parameter = id_
   #      else:
   #          column_name = "name"
   #          parameter = name
   #
   #      sql = """SELECT
   #                  count(*)
   #               FROM
   #                  machine_settings
   #               WHERE
   #                  %s = ?
   #            """ % column_name
   #      return self._execute(sql, parameter, fetch=_Fetch.SINGLETON)
   #
   #
   #  def count_pulse_projects(self, id_=None, name=None):
   #      """Given an id or a name, returns the number of pulse projects
   #      matching that criterion. Exactly one of id_ or name must be supplied.
   #      """
   #      if id_:
   #          column_name = "id"
   #          parameter = id_
   #      else:
   #          column_name = "name"
   #          parameter = name
   #
   #      sql = """SELECT
   #                  count(*)
   #               FROM
   #                  pulse_projects
   #               WHERE
   #                  %s = ?
   #            """ % column_name
   #      return self._execute(sql, parameter, fetch=_Fetch.SINGLETON)
   #
   #
   #
   #  def delete_machine_settings(self, ids):
   #      """Given a single machine settings id or a list of them, deletes the
   #      machine settings. It doesn't matter if any/all of the machine settings
   #      are templates.
   #      """
   #      # Make ids a list if it isn't one already
   #      if not util_misc.is_iterable(ids, False):
   #          ids = (ids, )
   #
   #      sql = """DELETE FROM
   #                  machine_settings
   #               WHERE
   #                  id IN (%s)
   #             """ % _build_placeholders(len(ids))
   #
   #      self._execute(sql, ids, _Fetch.NONE)
   #
   #
   #  def exists_pulse_project_in_database(self, ppid):
   #      # Lets first make sure it's there!
   #      sql = """SELECT
   #                  *
   #               FROM
   #                  pulse_projects
   #               WHERE
   #                  id = ?
   #            """
   #      return bool(self._execute(sql, ppid, fetch=_Fetch.ONE))
   #
   #
   #  def exists_transform_kernel_in_database(self, ppid):
   #      # Lets first make sure it's there!
   #      sql = """SELECT
   #                  *
   #               FROM
   #                  transform_kernels
   #               WHERE
   #                  id = ?
   #            """
   #      return bool(self._execute(sql, ppid, fetch=_Fetch.ONE))
   #
   #
   # def fetch_default_machine_settings_template(self):
   #     """Returns the default machine settings template. (There's always one)."""
   #     sql = """SELECT
   #                 *
   #              FROM
   #                 machine_settings
   #              WHERE
   #                 is_template = 1 AND
   #                 is_default  = 1
   #           """
   #
   #     row = self._execute(sql, None, _Fetch.ONE)
   #
   #     return rfp_machine_settings.MachineSettingsTemplate(row)
   #
   #
   # def fetch_machine_settings(self, id_):
   #     """Given the id of a machine settings object, returns that object. It
   #     doesn't matter if it's a template or not.
   #     """
   #
   #     sql = """SELECT
   #                 *
   #              FROM
   #                 machine_settings
   #              WHERE
   #                 id = ?
   #           """
   #
   #     row = self._execute(sql, id_, _Fetch.ONE)
   #
   #     if row["is_template"]:
   #         return rfp_machine_settings.MachineSettingsTemplate(row)
   #     else:
   #         return rfp_machine_settings.MachineSettings(row)
   #
   #
   # def fetch_machine_settings_templates(self):
   #     """Returns a list of machine settings templates ordered by name"""
   #     sql = """SELECT
   #                 *
   #              FROM
   #                 machine_settings
   #              WHERE
   #                 is_template = 1
   #              ORDER BY
   #                 name
   #           """
   #
   #     rows = self._execute(sql)
   #
   #     return [rfp_machine_settings.MachineSettingsTemplate(row) for row in rows]
   #
   #
   # def fetch_pulse_project(self, pulse_project_id):
   #     """Given a pulse project id, returns the associated pulse_project object."""
   #     sql = """SELECT
   #                 *
   #              FROM
   #                 pulse_projects
   #              WHERE
   #                 id = ?
   #           """
   #     row = self._execute(sql, pulse_project_id, fetch=_Fetch.ONE)
   #
   #     # Will inflate row into a pulse_project.
   #     pulse_project = rfp_pulse_project.PulseProject(row)
   #
   #     if row["master_parameters_id"]:
   #         master_parameters_id = row["master_parameters_id"]
   #         pulse_project.master_parameters = self._fetch_master_parameters(master_parameters_id)
   #
   #     if row["machine_settings_id"]:
   #         machine_settings_id = row["machine_settings_id"]
   #         pulse_project.machine_settings = self.fetch_machine_settings(machine_settings_id)
   #
   #     # Fetch transformations for this pulse_project
   #     sql =   """SELECT
   #                id, progression, transformation_type, parameters_id, result_id
   #              FROM
   #                 transformations
   #              WHERE
   #                 pulse_project_id = ?
   #              ORDER BY
   #                 progression
   #           """
   #     trows = self._execute(sql, pulse_project_id, fetch=_Fetch.ALL)
   #
   #     transformations = [ ]
   #
   #     calc_resolution = pulse_project.master_parameters.calc_resolution
   #
   #     for row in trows:
   #         trans_type = row["transformation_type"]
   #         trans_type = constants.TransformationType.get_type_for_value(trans_type, 'db')
   #         parameters = self._fetch_parameters(trans_type, row["parameters_id"])
   #
   #         transformation = \
   #             rfp_transformation_factory.create_transformation(trans_type, parameters)
   #         result_id = row["result_id"]
   #         if result_id:
   #             transformation.result = self._fetch_result(result_id)
   #         transformations.append(transformation)
   #
   #     pulse_project.transformations = transformations
   #
   #     pulse_project.referrers = self.fetch_pulse_project_referrers(pulse_project_id)
   #
   #     return pulse_project
   #
   #
   # def fetch_pulse_project_preview(self, id_):
   #     """Returns a preview of the pulse project matching the given UUID"""
   #     previews = self._fetch_pulse_project_previews(id_)
   #
   #     return (previews[0] if previews else None)
   #
   #
   # def fetch_pulse_project_previews(self):
   #     """Returns a list of pulse project previews."""
   #     return self._fetch_pulse_project_previews()
   #
   #
   # def fetch_pulse_project_referrers(self, id_):
   #     """Given a pulse project id, returns a (possibly empty) list of the
   #     ids of pulse sequences that use this pulse project.
   #     """
   #     sql = """SELECT
   #                 id, name
   #              FROM
   #                 pulse_sequences
   #              WHERE
   #                 id IN
   #                     (SELECT
   #                         pulse_sequence_id
   #                      FROM
   #                         pulse_sequence_pulse_projects
   #                      WHERE
   #                         pulse_project_id = ?
   #                     )
   #           """
   #     return [(row[0], row[1]) for row in self._execute(sql, id_)]
   #
   #
   #  def insert_pulse_project(self, pulse_project, use_transactions=True,
   #                           results_to_skip=[ ]):
   #      """Given a populated pulse_project object, creates database rows
   #      representing the new pulse_project and its associated transformations
   #      (and results - a.k.a. waveforms, profiles and gradients).
   #
   #      results_to_skip is a list of result objects. Any result present in
   #      this list will be skipped; i.e. not written to the database. This
   #      allows RFPulse to ignore out-of-sync results without destroying them.
   #      """
   #      if use_transactions:
   #          self.begin_transaction()
   #
   #      # These next two if's should not be required
   #      if not pulse_project.id:
   #          pulse_project.id = util_misc.uuid()
   #      if not pulse_project.created:
   #          pulse_project.created = util_time.now()
   #
   #      master_parameters_id = self._insert_master_parameters(pulse_project.master_parameters)
   #      machine_settings_id = self.insert_machine_settings(pulse_project.machine_settings)
   #
   #      sql = '''INSERT INTO
   #                  pulse_projects
   #                      (id, is_public, name, creator, created,
   #                      comment, machine_settings_id, master_parameters_id)
   #               VALUES
   #                  (?, ?, ?, ?, ?, ?, ?, ?)
   #            '''
   #
   #      sql_params = (pulse_project.id, pulse_project.is_public ,
   #                    pulse_project.name, pulse_project.creator,
   #                    pulse_project.created, pulse_project.comment,
   #                    machine_settings_id, master_parameters_id)
   #
   #      self._execute(sql, sql_params, _Fetch.NONE)
   #
   #      # Add Transformations
   #      self._insert_transformations(pulse_project.id,
   #                                   pulse_project.transformations,
   #                                   results_to_skip)
   #
   #      if use_transactions:
   #          self.commit()
   #
   #      return pulse_project.id
   #
   #
   #  def insert_machine_settings(self, machine_settings):
   #      # Templates are a little different from the settings associated with
   #      # a pulse project.
   #      is_template = isinstance(machine_settings,
   #                               rfp_machine_settings.MachineSettingsTemplate)
   #
   #      sql = '''INSERT INTO
   #                  machine_settings
   #                      (name, is_template, is_default,
   #                       machine_type,
   #                       field_strength, max_b1_field,
   #                       zero_padding, min_dwell_time,
   #                       dwell_time_increment,
   #                       gradient_raster_time,
   #                       gradient_slew_rate,
   #                       gradient_maximum)
   #               VALUES
   #                  (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
   #            '''
   #
   #      name = machine_settings.name if is_template else None
   #      is_default = machine_settings.is_default if is_template else False
   #
   #      machine_type = machine_settings.machine_type
   #      # If this is one of our predefined machine types, use its associated
   #      # constant. Otherwise, this is freeform text entered by the user so
   #      # we leave it as-is.
   #      if machine_type in constants.MachineType.ALL:
   #          machine_type = machine_type["db"]
   #
   #      sql_params = (name, is_template, is_default, machine_type,
   #                    machine_settings.field_strength,
   #                    machine_settings.max_b1_field,
   #                    machine_settings.zero_padding,
   #                    machine_settings.min_dwell_time,
   #                    machine_settings.dwell_time_increment,
   #                    machine_settings.gradient_raster_time,
   #                    machine_settings.gradient_slew_rate,
   #                    machine_settings.gradient_maximum)
   #
   #      return self._execute(sql, sql_params, _Fetch.LAST_ID)
   #
   #
   #  def insert_last_transformation(self, pulse_project, use_transactions=True):
   #
   #      progression = len(pulse_project.transformations)
   #      transformation = pulse_project.transformations[progression-1]
   #
   #      if use_transactions:
   #          self.begin_transaction()
   #
   #      parameters_id = self._insert_parameters(transformation.type,
   #                                              transformation.parameters)
   #      result_id = self._insert_result(transformation.result)
   #
   #      sql = '''INSERT INTO
   #                  transformations
   #                      (pulse_project_id,
   #                       progression,
   #                       transformation_type,
   #                       parameters_id,
   #                       result_id)
   #               VALUES
   #                  (?, ?, ?, ?, ?)
   #            '''
   #
   #      sql_params = (pulse_project.id, progression,
   #                    transformation.type['db'],
   #                    parameters_id, result_id)
   #
   #      transformation_id = self._execute(sql, sql_params, _Fetch.LAST_ID)
   #
   #      if use_transactions:
   #          self.commit()
   #
   #      return transformation_id
   #
   #
   #  def mark_machine_settings_template_as_default(self, id_):
   #      self.begin_transaction()
   #
   #      # Ensure no template is marked as default
   #      sql = """UPDATE
   #                  machine_settings
   #               SET
   #                  is_default = 0
   #               """
   #      self._execute(sql, None, _Fetch.NONE)
   #
   #      # Now mark one of them as default.
   #      sql = """UPDATE
   #                  machine_settings
   #               SET
   #                  is_default = 1
   #               WHERE
   #                  id = ?
   #               """
   #      self._execute(sql, id_, _Fetch.NONE)
   #
   #      self.commit()
        

    # #################   RFPulse  Private  Methods          #################
    # ##################   (in alphabetic order, no less!)   #################
    #
    # def _fetch_gradient(self, gradient_id):
    #
    #     sql = """SELECT
    #                 *
    #              FROM
    #                 gradients
    #              WHERE
    #                 id = ?
    #           """
    #
    #     row = self._execute(sql, gradient_id, _Fetch.ONE)
    #
    #     gradient = rfp_result.Gradient(row)
    #
    #     # Retrieve the individual data points, and recreate the gradient waveform(s).
    #     sql = '''SELECT
    #                 time_point, gradient_value, f2_value
    #              FROM
    #                 gradient_waveforms
    #              WHERE
    #                 gradient_id = ?
    #              ORDER BY
    #                 time_point
    #           '''
    #
    #     rows = self._execute(sql, gradient_id, fetch=_Fetch.ALL)
    #
    #     # Here we use the nifty Python idiom zip(*rows) for matrix
    #     # transposition as demonstrated in the example below:
    #     # >>> rows = [  (1, 2, 3), (1, 2, 3), (1, 2, 3), (1, 2, 3), ]
    #     # >>> print zip(*rows)
    #     # [(1, 1, 1, 1), (2, 2, 2, 2), (3, 3, 3, 3)]
    #     time_values, g2_values, f2_values = list(zip(*rows))
    #
    #     if rows:
    #         gradient.set_gradient_waveform(g2_values, time_values)
    #         if any(f2_values):
    #             gradient.set_f2(f2_values)
    #
    #     return gradient
    #
    #
    # def _fetch_gaussian(self, gaussian_id):
    #     """
    #     Given the gaussian parameters ID, return the associated object.
    #     """
    #
    #     sql = """SELECT
    #                 *
    #              FROM
    #                 gaussian_pulse_parameters
    #              WHERE
    #                 id = ?
    #           """
    #     row = self._execute(sql, gaussian_id, fetch=_Fetch.ONE)
    #
    #     # We turn the row into a true dict in order to modify it a little.
    #     row = dict(row)
    #
    #     row["time_points"] = row["time_steps"]
    #
    #     row['filter_type'] = \
    #         constants.FilterType.get_type_for_value(row['filter_type'], 'db')
    #
    #     return rfp_create_gaussian.GaussianPulseParameters(row)
    #
    #
    # def _fetch_hs(self, hs_id):
    #     """
    #     Given the hyperbolic secant parameters ID, return the associated object.
    #     """
    #
    #     sql = """SELECT
    #                 *
    #              FROM
    #                 hs_pulse_parameters
    #              WHERE
    #                 id = ?
    #           """
    #     row = self._execute(sql, hs_id, fetch=_Fetch.ONE)
    #
    #     # We turn the row into a true dict in order to modify it a little.
    #     row = dict(row)
    #
    #     row["time_points"] = row["time_steps"]
    #
    #     row['filter_type'] = \
    #         constants.FilterType.get_type_for_value(row['filter_type'], 'db')
    #
    #     return rfp_create_hs.HyperbolicSecantPulseParameters(row)
    #
    #
    # def _fetch_import(self, id_):
    #     """Given the id to the params of a create by import pulse, returns
    #     the associated object."""
    #     sql = """SELECT
    #                 *
    #              FROM
    #                 import_pulse_parameters
    #              WHERE
    #                 id = ?
    #           """
    #
    #     row = self._execute(sql, id_, fetch=_Fetch.ONE)
    #
    #     # We turn the row into a true dict in order to modify it a little.
    #     row = dict(row)
    #
    #     row['file_format'] = \
    #         constants.ImportFileFormat.get_type_for_value(row['file_format'], 'db')
    #
    #     return rfp_create_import.CreateImportPulseParameters(row)
    #
    #
    # def _fetch_interpolate_rescale(self, interpolate_rescale_id):
    #     """
    #     Given the interpolate rescale parameters ID, return the associated
    #     interpolation, filtering, rescaling object.
    #     """
    #
    #     sql = """SELECT
    #                 *
    #              FROM
    #                 interpolate_rescale_parameters
    #              WHERE
    #                 id = ?
    #           """
    #     row = self._execute(sql, interpolate_rescale_id, fetch=_Fetch.ONE)
    #
    #     return rfp_interpolate_rescale.InterpolateRescaleParameters(row)
    #
    #
    # def _fetch_master_parameters(self, id_):
    #     """Given a master_parameters_id, returns the associated
    #     master_parameters object."""
    #     sql = """SELECT
    #                 *
    #              FROM
    #                 master_parameters
    #              WHERE
    #                 id = ?
    #           """
    #     row = self._execute(sql, id_, fetch=_Fetch.ONE)
    #
    #     # We turn the row into a true dict in order to modify it a little.
    #     row = dict(row)
    #
    #     # Convert this to a standard constant.
    #     row["pulse_bandwidth_type"] = \
    #         constants.PulseConvention.get_type_for_value(row["pulse_bandwidth_type"], "db")
    #
    #     # Will inflate row into a MasterParameters.
    #     return rfp_master_parameters.MasterParameters(row)
    #
    #
    # def _fetch_ocn(self, pulse_id):
    #
    #     """Given an optimal control (non-selective) pulse_id,
    #        return the associated OCNParameters object."""
    #
    #     sql = """SELECT
    #                 *
    #              FROM
    #                 ocn_parameters
    #              WHERE
    #                 id = ?
    #           """
    #     row = self._execute(sql, pulse_id, fetch=_Fetch.ONE)
    #
    #     # We turn the row into a true dict in order to modify it a little.
    #     row = dict(row)
    #
    #     # These get converted into proper constants
    #     row['pulse_type'] = \
    #         constants.UsageType.get_type_for_value(row['pulse_type'], 'db')
    #
    #
    #     # This is a bit of a hack as we have a "type" that spans two types.
    #     phase_type = row['phase_type']
    #     row['phase_type'] = \
    #         constants.PhaseType.get_type_for_value(phase_type, 'db')
    #     if not row['phase_type']:
    #         constants.NonCoalescedPhaseSubtype.get_type_for_value(phase_type, 'db')
    #
    #     row["step_size_modification"] = \
    #         constants.StepSizeModification.get_type_for_value(row['step_size_modification'], 'db')
    #
    #     return rfp_ocn.OCNParameters(row)
    #
    #
    # def _fetch_ocn_state(self, ocn_state_id):
    #     # First retrieve all the individual data points for deltab1.
    #     sql = '''SELECT
    #                 real_amplitude, imaginary_amplitude
    #              FROM
    #                 deltab1_points
    #              WHERE
    #                 ocn_state_id = ?
    #              ORDER BY
    #                 progression
    #           '''
    #
    #     rows = self._execute(sql, ocn_state_id)
    #
    #     deltab1 = [complex(real, imaginary) for real, imaginary in rows]
    #     deltab1 = numpy.array(deltab1, dtype=numpy.complex)
    #
    #     # Fetch residual errors
    #     sql = '''SELECT
    #                 value
    #              FROM
    #                 ocn_residual_errors
    #              WHERE
    #                 ocn_state_id = ?
    #              ORDER BY
    #                 progression
    #           '''
    #
    #     residual_errors = self._execute(sql, ocn_state_id, _Fetch.COLUMN)
    #
    #     # Now get all of the scalar values
    #     sql = """SELECT
    #                 *
    #              FROM
    #                 ocn_states
    #              WHERE
    #                 id = ?
    #           """
    #     row = self._execute(sql, ocn_state_id, fetch=_Fetch.ONE)
    #
    #     # We turn the row into a true dict in order to modify it a little.
    #     row = dict(row)
    #
    #     row["deltab1"] = deltab1
    #     row["residual_errors"] = residual_errors
    #
    #     return rfp_ocn_state.OCNState(row)
    #
    #
    # def _fetch_parameters(self, trans_type, p_id):
    #
    #     if trans_type == constants.TransformationType.CREATE_SLR:
    #         return self._fetch_slr(p_id)
    #
    #     elif trans_type == constants.TransformationType.CREATE_GAUSSIAN:
    #         return self._fetch_gaussian(p_id)
    #
    #     elif trans_type == constants.TransformationType.CREATE_RANDOMIZED:
    #         return self._fetch_randomized(p_id)
    #
    #     elif trans_type == constants.TransformationType.CREATE_IMPORT:
    #         return self._fetch_import(p_id)
    #
    #     elif trans_type == constants.TransformationType.CREATE_HYPERBOLIC_SECANT:
    #         return self._fetch_hs(p_id)
    #
    #     elif trans_type == constants.TransformationType.INTERPOLATE_RESCALE:
    #         return self._fetch_interpolate_rescale(p_id)
    #
    #     elif trans_type == constants.TransformationType.OCN:
    #         return self._fetch_ocn(p_id)
    #
    #     elif trans_type == constants.TransformationType.ROOT_REFLECTION:
    #         return self._fetch_root_reflect(p_id)
    #
    #     else:
    #         return None
    #
    #
    # def _fetch_pulse_project_previews(self, id_=None):
    #     """Returns a list of pulse project previews. If the id_ parameter is
    #     populated, only projects with that UUID are returned. Of course there
    #     will be at most one matching project, but it's still returned in a
    #     list.
    #     """
    #     parameters = [ ]
    #     where_clause = ""
    #     if id_:
    #         where_clause = " AND id = ? "
    #         parameters = (id_, )
    #
    #     sql = """SELECT
    #                 id, name, creator, created, is_public, comment
    #              FROM
    #                 pulse_projects
    #              WHERE
    #                 1 = 1  %s
    #              ORDER BY
    #                 name
    #           """ % where_clause
    #
    #     rows = self._execute(sql, parameters)
    #
    #     previews = [pulse_project_preview.PulseProjectPreview(row) for row
    #                                                                in rows]
    #
    #     for preview in previews:
    #         preview.referrers = self.fetch_pulse_project_referrers(preview.id)
    #
    #     return previews
    #
    #
    #
    # def _fetch_randomized(self, id_):
    #     """Given the id to the params of a randomized create pulse, returns
    #     the associated object."""
    #     sql = """SELECT
    #                 *
    #              FROM
    #                 randomized_pulse_parameters
    #              WHERE
    #                 id = ?
    #           """
    #     row = self._execute(sql, id_, fetch=_Fetch.ONE)
    #
    #     # We turn the row into a true dict in order to modify it a little.
    #     row = dict(row)
    #
    #     row["time_points"] = row["time_steps"]
    #
    #     return rfp_create_randomized.RandomizedPulseParameters(row)
    #
    #
    # def _fetch_result(self, result_id):
    #
    #     sql = '''SELECT
    #                 *
    #              FROM
    #                 results
    #              Where
    #                 id = ?
    #           '''
    #
    #     row = self._execute(sql, result_id, fetch=_Fetch.ONE)
    #
    #     created = row["created"]
    #
    #     if row["gradient_id"]:
    #         gradient = self._fetch_gradient(row["gradient_id"])
    #     else:
    #         gradient = None
    #
    #     if row["ocn_state_id"]:
    #         ocn_state = self._fetch_ocn_state(row["ocn_state_id"])
    #     else:
    #         ocn_state = None
    #
    #     # Retrieve all the individual data points, and recreate the waveform.
    #     sql = '''SELECT
    #                 time_point, real_amplitude, imaginary_amplitude
    #              FROM
    #                 rf_waveforms
    #              WHERE
    #                 result_id = ?
    #              ORDER BY
    #                 time_point
    #           '''
    #
    #     rows = self._execute(sql, result_id)
    #
    #     time_array = []
    #     complex_array = []
    #     for row in rows:
    #         time_array.append(row["time_point"])
    #         complex_array.append(complex(row["real_amplitude"],
    #                                      row["imaginary_amplitude"]))
    #
    #     rf = rfp_result.Waveform({"waveform"        : complex_array,
    #                               "waveform_x_axis" : time_array})
    #
    #     return rfp_result.Result({"rf"        : rf,
    #                               "gradient"  : gradient,
    #                               "ocn_state" : ocn_state,
    #                               "created"   : created})
    #
    #
    # def _fetch_root_reflect(self, rr_id):
    #
    #     """
    #     Given the root_reflect_parameters ID, return the associated object.
    #     """
    #
    #     sql = """SELECT
    #                 *
    #              FROM
    #                 root_reflect_parameters
    #              WHERE
    #                 id = ?
    #           """
    #     row = self._execute(sql, rr_id, fetch=_Fetch.ONE)
    #
    #     # We turn the row into a true dict in order to modify it a little.
    #     d = dict(row)
    #
    #     d['anorm'] = complex(d['anorm_real'], d['anorm_imaginary'])
    #     d['bnorm'] = complex(d['bnorm_real'], d['bnorm_imaginary'])
    #
    #
    #     # A roots.
    #
    #     sql = """SELECT
    #                 *
    #              FROM
    #                 a_roots
    #              WHERE
    #                 root_reflect_id = ?
    #           """
    #     arows = self._execute(sql, rr_id, fetch=_Fetch.ALL)
    #
    #     aroots = []
    #     aroots_flipped = []
    #     for arow in arows:
    #         aroots_flipped.append(arow["was_flipped"])
    #         aroots.append(complex(arow["aroot_real"], arow["aroot_imaginary"]))
    #
    #     d['aroots'] = aroots
    #     d['aroots_flipped'] = aroots_flipped
    #
    #
    #     # B roots
    #
    #     sql = """SELECT
    #                 *
    #              FROM
    #                 b_roots
    #              WHERE
    #                 root_reflect_id = ?
    #           """
    #     brows = self._execute(sql, rr_id, fetch=_Fetch.ALL)
    #
    #     broots = []
    #     broots_flipped = []
    #     for brow in brows:
    #         broots_flipped.append(brow["was_flipped"])
    #         broots.append(complex(brow["broot_real"], brow["broot_imaginary"]))
    #
    #     d['broots'] = broots
    #     d['broots_flipped'] = broots_flipped
    #
    #     return rfp_root_reflection.RootReflectionParameters(d)
    #
    #
    # def _fetch_slr(self, pulse_id):
    #
    #     """Given an slr pulse_id,
    #        return the associated SlrPulseParameters object."""
    #
    #     sql = """SELECT
    #                 *
    #              FROM
    #                 slr_pulse_parameters
    #              WHERE
    #                 id = ?
    #           """
    #     row = self._execute(sql, pulse_id, fetch=_Fetch.ONE)
    #
    #     # We turn the row into a true dict in order to modify it a little.
    #     row = dict(row)
    #
    #     # This value gets munged a little.
    #     row["time_points"] = row["time_steps"]
    #
    #     # These get converted into proper constants
    #     row['nc_phase_subtype'] = \
    #         constants.NonCoalescedPhaseSubtype.get_type_for_value(row['nc_phase_subtype'], 'db')
    #
    #     row["use_remez"] = \
    #             (row['slr_filter_type'] == constants.SLRFilterType.REMEZ['db'])
    #
    #     return rfp_create_slr.SLRPulseParameters(row)
    #
    #
    # def _insert_gradient(self, gradient):
    #     # Currently this code is only called from within a transaction, so
    #     # it doesn't need to open a transaction itself.
    #     if not gradient:
    #         return 0
    #
    #     sql = '''INSERT INTO
    #                 gradients
    #                     (linear_gradient_value,
    #                      refocused_gradient,
    #                      frequency_offset)
    #              VALUES
    #                 (?, ?, ?)
    #           '''
    #
    #     sql_params = (gradient.linear_gradient,
    #                   gradient.refocusing_gradient,
    #                   gradient.frequency_offset)
    #
    #     gradient_id = self._execute(sql, sql_params, _Fetch.LAST_ID)
    #
    #     f2_len = len(gradient.f2_frequency_waveform)
    #     ii = 0
    #
    #     sql = '''INSERT INTO
    #                gradient_waveforms
    #                   (gradient_id,
    #                    time_point,
    #                    gradient_value,
    #                    f2_value)
    #              VALUES
    #                   (?, ?, ?, ?)
    #          '''
    #     for tms in gradient.time_points:
    #         f2 = 0.0
    #         if f2_len > 0:
    #            f2 = gradient.f2_frequency_waveform[ii]
    #
    #         gt = gradient.time_points[ii]
    #         gw = gradient.gradient_waveform[ii]
    #
    #         sql_params = (gradient_id, gt, gw, f2)
    #         self._execute(sql, sql_params, _Fetch.NONE)
    #
    #         ii += 1
    #
    #     return gradient_id
    #
    #
    # def _insert_gaussian(self, parameters):
    #
    #     sql = '''INSERT INTO
    #                 gaussian_pulse_parameters
    #                     (tip_angle,
    #                      time_steps,
    #                      duration,
    #                      bandwidth,
    #                      filter_type,
    #                      filter_application)
    #              VALUES
    #                 (?, ?, ?, ?, ?, ?)
    #           '''
    #
    #     sql_params = (parameters.tip_angle,
    #                   parameters.time_points,
    #                   parameters.duration,
    #                   parameters.bandwidth,
    #                   parameters.filter_type['db'],
    #                   parameters.filter_application)
    #
    #     return self._execute(sql, sql_params, _Fetch.LAST_ID)
    #
    #
    # def _insert_hyperbolic_secant(self, parameters):
    #
    #     sql = '''INSERT INTO
    #                 hs_pulse_parameters
    #                     (total_rotation,
    #                      time_steps,
    #                      dwell_time,
    #                      was_bandwidth_specified,
    #                      quality_cycles,
    #                      power_n,
    #                      sharpness_mu,
    #                      filter_type,
    #                      filter_application)
    #              VALUES
    #                 (?, ?, ?, ?, ?, ?, ?, ?, ?)
    #           '''
    #
    #     sql_params = (parameters.total_rotation,
    #                   parameters.time_points,
    #                   parameters.dwell_time,
    #                   parameters.was_bandwidth_specified,
    #                   parameters.quality_cycles,
    #                   parameters.power_n,
    #                   parameters.sharpness_mu,
    #                   parameters.filter_type['db'],
    #                   parameters.filter_application)
    #
    #     return self._execute(sql, sql_params, _Fetch.LAST_ID)
    #
    #
    # def _insert_import(self, parameters):
    #
    #     sql = '''INSERT INTO
    #                 import_pulse_parameters
    #                     (file_path,
    #                      comment,
    #                      file_format,
    #                      dwell_time,
    #                      use_max_intensity,
    #                      max_intensity,
    #                      scale_factor
    #                     )
    #              VALUES
    #                 (?, ?, ?, ?, ?, ?, ?)
    #           '''
    #
    #     params = (parameters.file_path,
    #               parameters.comment,
    #               parameters.file_format['db'],
    #               parameters.dwell_time,
    #               parameters.use_max_intensity,
    #               parameters.max_intensity,
    #               parameters.scale_factor
    #               )
    #
    #     return self._execute(sql, params, _Fetch.LAST_ID)
    #
    #
    # def _insert_interpolate_rescale(self, parameters):
    #
    #     sql = '''INSERT INTO
    #                 interpolate_rescale_parameters
    #                     (do_interpolate,
    #                      interpolation_factor,
    #                      new_dwell_time,
    #                      do_rescaling,
    #                      angle)
    #              VALUES
    #                 (?, ?, ?, ?, ?)
    #           '''
    #
    #     sql_params = (parameters.do_interpolate,
    #                   parameters.interpolation_factor,
    #                   parameters.new_dwell_time,
    #                   parameters.do_rescaling,
    #                   parameters.angle)
    #
    #     return self._execute(sql, sql_params, _Fetch.LAST_ID)
    #
    #
    # def _insert_master_parameters(self, master_parameters):
    #
    #     sql = '''INSERT INTO
    #                 master_parameters
    #                     (calc_resolution, pulse_bandwidth_type)
    #              VALUES
    #                 (?, ?)
    #           '''
    #
    #     pulse_bandwidth_type = master_parameters.pulse_bandwidth_convention['db']
    #
    #     sql_params = (master_parameters.calc_resolution, pulse_bandwidth_type)
    #
    #     return self._execute(sql, sql_params, _Fetch.LAST_ID)
    #
    #
    # def _insert_ocn_state(self, ocn_state):
    #
    #     if not ocn_state:
    #         return 0
    #
    #     sql = '''INSERT INTO
    #                 ocn_states
    #                    (multiplier,
    #                     met_max_iterations,
    #                     met_residual_error,
    #                     met_differential_error,
    #                     met_max_time,
    #                     met_increasing_error,
    #                     run_time,
    #                     iterations,
    #                     decreases
    #                    )
    #              VALUES
    #                 (?, ?, ?, ?, ?, ?, ?, ?, ?)
    #           '''
    #
    #     params = (ocn_state.multiplier,
    #               ocn_state.met_max_iterations,
    #               ocn_state.met_residual_error,
    #               ocn_state.met_differential_error,
    #               ocn_state.met_max_time,
    #               ocn_state.met_increasing_error,
    #               ocn_state.run_time,
    #               ocn_state.iterations,
    #               ocn_state.decreases
    #              )
    #
    #     ocn_state_id = self._execute(sql, params, _Fetch.LAST_ID)
    #
    #     # Now insert the residual error history
    #     sql = '''INSERT INTO
    #                 ocn_residual_errors
    #                     (ocn_state_id, value, progression)
    #              VALUES
    #                 (?, ?, ?)
    #           '''
    #     for i, residual_error in enumerate(ocn_state.residual_errors):
    #         self._execute(sql, (ocn_state_id, residual_error, i) , _Fetch.NONE)
    #
    #     # Now insert all the individual data points for deltab1.
    #     # deltab1 is a numpy array of complex numbers
    #     sql = '''INSERT INTO
    #                 deltab1_points
    #                     (ocn_state_id, progression,
    #                      real_amplitude, imaginary_amplitude)
    #              VALUES
    #                 (?, ?, ?, ?)
    #           '''
    #     for progression, z in enumerate(ocn_state.deltab1):
    #         params = (ocn_state_id, progression, z.real, z.imag)
    #         self._execute(sql, params, _Fetch.NONE)
    #
    #     return ocn_state_id
    #
    #
    # def _insert_ocn(self, parameters):
    #
    #     sql = '''INSERT INTO
    #                 ocn_parameters
    #                     (pulse_type,
    #                      phase_type,
    #                      gradient_refocusing_value,
    #                      tip_angle,
    #                      bandwidth,
    #                      step_size_multiplier,
    #                      step_size_modification,
    #                      excite_band_points,
    #                      b1_immunity_range,
    #                      steps,
    #                      b1_maximum,
    #                      limit_sar,
    #                      sar_factor,
    #                      error_increase_tolerance,
    #                      max_iteration_check,
    #                      max_iterations,
    #                      residual_error_check,
    #                      residual_error_tolerance,
    #                      differential_error_check,
    #                      differential_error_tolerance,
    #                      halt_if_error_increasing,
    #                      halt_on_max_time,
    #                      max_time,
    #                      enforce_symmetry)
    #              VALUES
    #                 (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    #           '''
    #
    #     sql_params = (parameters.pulse_type['db'],
    #                   parameters.phase_type['db'],
    #                   parameters.gradient_refocusing_value,
    #                   parameters.tip_angle,
    #                   parameters.bandwidth,
    #                   parameters.step_size_multiplier,
    #                   parameters.step_size_modification['db'],
    #                   parameters.excite_band_points,
    #                   parameters.b1_immunity_range,
    #                   parameters.steps,
    #                   parameters.b1_maximum,
    #                   parameters.limit_sar,
    #                   parameters.sar_factor,
    #                   parameters.error_increase_tolerance,
    #                   parameters.max_iteration_check,
    #                   parameters.max_iterations,
    #                   parameters.residual_error_check,
    #                   parameters.residual_error_tolerance,
    #                   parameters.differential_error_check,
    #                   parameters.differential_error_tolerance,
    #                   parameters.halt_if_error_increasing,
    #                   parameters.halt_on_max_time,
    #                   parameters.max_time,
    #                   parameters.enforce_symmetry)
    #
    #     return self._execute(sql, sql_params, _Fetch.LAST_ID)
    #
    #
    # def _insert_parameters(self, trans_type, parameters):
    #
    #     if trans_type == constants.TransformationType.CREATE_SLR:
    #         return self._insert_slr(parameters)
    #
    #     elif trans_type == constants.TransformationType.CREATE_GAUSSIAN:
    #         return self._insert_gaussian(parameters)
    #
    #     elif trans_type == constants.TransformationType.CREATE_RANDOMIZED:
    #         return self._insert_randomized(parameters)
    #
    #     elif trans_type == constants.TransformationType.CREATE_IMPORT:
    #         return self._insert_import(parameters)
    #
    #     elif trans_type == constants.TransformationType.CREATE_HYPERBOLIC_SECANT:
    #         return self._insert_hyperbolic_secant(parameters)
    #
    #     elif trans_type == constants.TransformationType.INTERPOLATE_RESCALE:
    #         return self._insert_interpolate_rescale(parameters)
    #
    #     elif trans_type == constants.TransformationType.OCN:
    #         return self._insert_ocn(parameters)
    #
    #     elif trans_type == constants.TransformationType.ROOT_REFLECTION:
    #         return self._insert_root_reflect(parameters)
    #
    #     else:
    #         error_string = "Unknown transformation/params type: Unable to Insert"
    #         raise NotImplementedError(error_string)
    #
    #
    # def _insert_randomized(self, parameters):
    #
    #     sql = '''INSERT INTO
    #                 randomized_pulse_parameters
    #                     (time_steps,
    #                      duration
    #                     )
    #              VALUES
    #                 (?, ?)
    #           '''
    #
    #     params = (parameters.time_points, parameters.duration)
    #
    #     return self._execute(sql, params, _Fetch.LAST_ID)
    #
    #
    # def _insert_result(self, result):
    #     # Currently this code is only called from within a transaction, so
    #     # it doesn't need to open a transaction itself.
    #     gradient_id = self._insert_gradient(result.gradient)
    #     ocn_state_id = self._insert_ocn_state(result.ocn_state)
    #
    #     # Insert information about the container object.
    #     sql = '''INSERT INTO
    #                 results
    #                     (created,
    #                      gradient_id,
    #                      ocn_state_id)
    #              VALUES
    #                 (?, ?, ?)
    #           '''
    #
    #     sql_params = (result.created, gradient_id, ocn_state_id)
    #
    #     result_id = self._execute(sql, sql_params, _Fetch.LAST_ID)
    #
    #     # Now insert all the individual data points for the rf waveforms.
    #     sql = '''INSERT INTO
    #                 rf_waveforms
    #                     (result_id, time_point,
    #                      real_amplitude, imaginary_amplitude)
    #              VALUES
    #                 (?, ?, ?, ?)
    #           '''
    #
    #     rlen = len(result.rf.waveform_x_axis)
    #     for ii in range(rlen):
    #         rtime = result.rf.waveform_x_axis[ii]
    #         z = result.rf.waveform[ii]
    #         zreal = z.real
    #         zimaginary = z.imag
    #         sql_params = (result_id, rtime, zreal, zimaginary)
    #         self._execute(sql, sql_params, _Fetch.NONE)
    #
    #
    #     return result_id
    #
    #
    # def _insert_root_reflect(self, parameters):
    #
    #     sql = '''INSERT INTO
    #                 root_reflect_parameters
    #                     (a_roots_only,
    #                      graph_angle,
    #                      x_axis_start,
    #                      anorm_real,
    #                      anorm_imaginary,
    #                      bnorm_real,
    #                      bnorm_imaginary,
    #                      leading_zeros,
    #                      trailing_zeros)
    #              VALUES
    #                 (?, ?, ?, ?, ?, ?, ?, ?, ?)
    #           '''
    #
    #     sql_params = (parameters.a_roots_only,
    #                   parameters.graph_angle,
    #                   parameters.x_axis_start,
    #                   parameters.anorm.real,
    #                   parameters.anorm.imag,
    #                   parameters.bnorm.real,
    #                   parameters.bnorm.imag,
    #                   parameters.leading_zeros,
    #                   parameters.trailing_zeros)
    #
    #     root_reflect_id = self._execute(sql, sql_params, _Fetch.LAST_ID)
    #
    #     sql = '''INSERT INTO
    #                  a_roots
    #                     (root_reflect_id,
    #                      aroot_real,
    #                      aroot_imaginary,
    #                      was_flipped)
    #              VALUES
    #                 (?, ?, ?, ?)
    #           '''
    #
    #     for i in range(len(parameters.aroots)):
    #         sql_params = (root_reflect_id,
    #                        parameters.aroots[i].real,
    #                        parameters.aroots[i].imag,
    #                        parameters.aroots_flipped[i])
    #         self._execute(sql, sql_params, _Fetch.NONE)
    #
    #     sql = '''INSERT INTO
    #                  b_roots
    #                     (root_reflect_id,
    #                      broot_real,
    #                      broot_imaginary,
    #                      was_flipped)
    #              VALUES
    #                 (?, ?, ?, ?)
    #           '''
    #
    #     for i in range(len(parameters.broots)):
    #         sql_params = (root_reflect_id,
    #                        parameters.broots[i].real,
    #                        parameters.broots[i].imag,
    #                        parameters.broots_flipped[i])
    #         self._execute(sql, sql_params, _Fetch.NONE)
    #
    #     return root_reflect_id
    #
    #
    # def _insert_slr(self, parameters):
    #
    #     sql = '''INSERT INTO
    #                 slr_pulse_parameters
    #                     (tip_angle,
    #                      time_steps,
    #                      duration,
    #                      bandwidth,
    #                      separation,
    #                      is_single_band,
    #                      nc_phase_subtype,
    #                      slr_filter_type,
    #                      pass_ripple,
    #                      reject_ripple)
    #              VALUES
    #                 (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    #           '''
    #
    #     if parameters.use_remez:
    #         slr_filter_type = constants.SLRFilterType.REMEZ['db']
    #     else:
    #         slr_filter_type = constants.SLRFilterType.LEAST_SQUARES['db']
    #
    #     sql_params = (parameters.tip_angle,
    #                   parameters.time_points,
    #                   parameters.duration,
    #                   parameters.bandwidth,
    #                   parameters.separation,
    #                   parameters.is_single_band,
    #                   parameters.nc_phase_subtype['db'],
    #                   slr_filter_type,
    #                   parameters.pass_ripple,
    #                   parameters.reject_ripple)
    #
    #     return self._execute(sql, sql_params, _Fetch.LAST_ID)
    #
    #
    # def _insert_transformations(self, ppid, transformations,
    #                             results_to_skip=[ ]):
    #     # Currently this code is only called from within a transaction, so
    #     # it doesn't need to open a transaction itself.
    #
    #     # See insert_pulse_project() for an explanation of results_to_skip
    #
    #     progression = 1
    #
    #     for transformation in transformations:
    #
    #         parameters_id = self._insert_parameters(transformation.type,
    #                                                 transformation.parameters)
    #
    #         if transformation.result and \
    #            (transformation.result not in results_to_skip):
    #             result_id = self._insert_result(transformation.result)
    #         else:
    #             result_id = 0
    #
    #         sql = '''INSERT INTO
    #                     transformations
    #                         (pulse_project_id,
    #                          progression,
    #                          transformation_type,
    #                          parameters_id,
    #                          result_id)
    #                  VALUES
    #                     (?, ?, ?, ?, ?)
    #               '''
    #
    #
    #         sql_params = (ppid, progression,
    #                       transformation.type['db'],
    #                       parameters_id, result_id)
    #
    #         self._execute(sql, sql_params, _Fetch.NONE)
    #
    #         progression += 1
    #
    #
    #
    #




    #################    Pulse  Public  Methods          #################
    ##################   (in alphabetic order, no less!)   #################

    def count_machine_specs_templates(self, id_=None, name=None):
        """Given an id or a name, returns the number of machine specs
        templates  matching that criterion. Exactly one of id_ or name must be
        supplied.
        """
        if id_:
            column_name = "id"
            parameter = id_
        else:
            column_name = "name"
            parameter = name

        sql = """SELECT
                    count(*)
                 FROM
                    machine_specs
                 WHERE
                    %s = ?
              """ % column_name
        return self._execute(sql, parameter, fetch=_Fetch.SINGLETON)


    def insert_machine_specs(self, machine_specs):
        """ 
        Templates are a little different from the specs associated with
        a pulse design.
        
        """
        is_template = isinstance(machine_specs, rfp_machine_specs.MachineSpecsTemplate)
        
        sql = '''INSERT INTO
                    machine_specs
                        (name, 
                         is_template, 
                         is_default,
                         machine_type,
                         field_strength, 
                         max_b1_field,
                         zero_padding, 
                         min_dwell_time,
                         dwell_time_increment,
                         gradient_raster_time,
                         gradient_slew_rate,
                         gradient_maximum)
                 VALUES
                    (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
              '''
              
        name = machine_specs.name if is_template else None
        is_default = machine_specs.is_default if is_template else False
        
        machine_type = machine_specs.machine_type
        # If this is one of our predefined machine types, use its associated
        # constant. Otherwise, this is freeform text entered by the user so
        # we leave it as-is.
        if machine_type in constants.MachineType.ALL:
            machine_type = machine_type["db"]
        
        sql_params = (name, 
                      is_template, 
                      is_default, 
                      machine_type,
                      machine_specs.field_strength,
                      machine_specs.max_b1_field,
                      machine_specs.zero_padding,
                      machine_specs.min_dwell_time,
                      machine_specs.dwell_time_increment,
                      machine_specs.gradient_raster_time,
                      machine_specs.gradient_slew_rate,
                      machine_specs.gradient_maximum)

        return self._execute(sql, sql_params, _Fetch.LAST_ID)


    def mark_machine_specs_template_as_default(self, id_):
        self.begin_transaction()
        
        # Ensure no template is marked as default
        sql = """UPDATE
                    machine_specs
                 SET
                    is_default = 0
                 """
        self._execute(sql, None, _Fetch.NONE)
        
        # Now mark one of them as default.
        sql = """UPDATE
                    machine_specs
                 SET
                    is_default = 1
                 WHERE
                    id = ?
                 """
        self._execute(sql, id_, _Fetch.NONE)
        
        self.commit()   


    def fetch_default_machine_specs_template(self):
        """
        Returns the default machine specs template. (There's always one).
        
        This is in main db.py because db_upgrader.py needs it
        """
        sql = """SELECT
                    *
                 FROM
                    machine_specs
                 WHERE
                    is_template = 1 AND
                    is_default  = 1
              """                                      
        row = self._execute(sql, None, _Fetch.ONE)
        
        return rfp_machine_specs.MachineSpecsTemplate(row)


    def fetch_transform_kernel(self, id_):
        """
        Given a transform kernel id, returns the associated 
        transform_kernel object.
        
        """
        sql = """SELECT
                    *
                 FROM
                    transform_kernels
                 WHERE
                    id = ?
              """
        row = self._execute(sql, id_, fetch=_Fetch.ONE)

        # Will inflate row into a transform_kernel.
        transform_kernel = rfp_transform_kernel.TransformKernel(row)

        # Fetch pulse sequence params. (There might not be any.)
        sql = """SELECT
                    *
                 FROM
                    transform_kernel_controls
                 WHERE
                    transform_kernel_id = ?
                 ORDER BY
                    display_order
                 """
        rows = self._execute(sql, transform_kernel.id)

        # Some special attention is required here. Unlike most of our
        # tables/objects, there's not a 1:1 map between the column and 
        # attribute names for pulse sequence params. The table 
        # transform_kernel_controls contains a column called 
        # default_ while the TransformKernelControl class has an 
        # attribute called default. Similarly for name_, type_ and
        # variable_ columns. We can't use the latter as a column 
        # name because it's a reserved word in SQL, and neither can we 
        # use it as a column name alias. 
        # So here I rename the column in the returned row. Remember 
        # that the rows returned by self._execute() are _BetterRow objects
        # and need to be turned into proper Python dicts before I 
        # manipulate them.
        rows = [dict(row) for row in rows]

        for row in rows:
            row["name"]     = row["name_"]
            row["type"]     = row["type_"]
            row["default"]  = row["default_"]
            row["variable"] = row["variable_"]
            del row["name_"]
            del row["type_"]
            del row["default_"]
            del row["variable_"]

        transform_kernel.transform_kernel_controls = [rfp_transform_kernel.TransformKernelControl(row) for row in rows]

        transform_kernel.referrers = self.fetch_transform_kernel_referrers(transform_kernel.id)

        return transform_kernel


    def fetch_transform_kernel_referrers(self, id_):
        """
        Given a transform_kernel id, returns a (possibly empty) list of the 
        ids of pulse designs that use this transform_kernel. 
        """
        return None
        
        sql = """SELECT
                    id, name
                 FROM
                    pulse_designs
                 WHERE
                    id IN
                        (SELECT
                            pulse_design_id
                         FROM
                            transforms
                         WHERE
                            transform_kernel_id = ?
                        )
              """
        return [(row[0], row[1]) for row in self._execute(sql, id_)]   
        

    def insert_transform_kernel(self, transform_kernel):
        """
        Saves the transform_kernel described by the param. The
        transform_kernel must not exist in the database.

        If the id is not populated, an id will be assigned.

        It is assumed that all fields in the pulse sequence are valid and
        that database constraints (e.g. unique name) will not be violated.

        When this function completes, the id fields of the parameters
        are populated.
        """
        self.begin_transaction()

        if not transform_kernel.id:
            transform_kernel.id = util_misc.uuid()

        if not transform_kernel.created:
            transform_kernel.created = util_time.now()

        # Insert the transform_kernels
        sql = """INSERT INTO
                    transform_kernels
                        (id, 
                         type, 
                         name, 
                         menu_label, 
                         is_public, 
                         created, 
                         creator, 
                         comment, 
                         algorithm_code, 
                         time_steps,
                         duration, 
                         hide_file1, 
                         hide_file2, 
                         file1_label,
                         file2_label, 
                         hide_time_steps,
                         hide_duration,
                         hide_tip_angle, 
                         hide_bandwidth,
                         tip_angle, 
                         bandwidth)
                 VALUES
                    (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
              """
        params = (transform_kernel.id, 
                  transform_kernel.type,
                  transform_kernel.name,
                  transform_kernel.menu_label,
                  transform_kernel.is_public,
                  transform_kernel.created, 
                  transform_kernel.creator,
                  transform_kernel.comment,
                  transform_kernel.algorithm_code,
                  transform_kernel.time_steps,
                  transform_kernel.duration,
                  transform_kernel.hide_file1,
                  transform_kernel.hide_file2,
                  transform_kernel.file1_label,
                  transform_kernel.file2_label,
                  transform_kernel.hide_time_steps,
                  transform_kernel.hide_duration,
                  transform_kernel.hide_tip_angle,
                  transform_kernel.hide_bandwidth,
                  transform_kernel.tip_angle,
                  transform_kernel.bandwidth)
        self._execute(sql, params, _Fetch.NONE)

        # Insert the params
        sql = """INSERT INTO
                    transform_kernel_controls
                        (transform_kernel_id, 
                        name_, 
                        type_, 
                        default_, 
                        variable_, 
                        display_order)
                 VALUES
                    (?, ?, ?, ?, ?, ?)
              """
        for i, control in enumerate(transform_kernel.transform_kernel_controls):
            params = (transform_kernel.id, 
                      control.name,
                      control.type, 
                      control.default, 
                      control.variable,  i)

            self._execute(sql, params, _Fetch.NONE)


        self.commit()        


    def count_transform_kernels(self, id_=None, name=None):
        """Given an id or a name, returns the number of transform_kernels 
        matching that criterion. Exactly one of id_ or name must be supplied.
        """
        if id_:
            column_name = "id"
            parameter = id_
        else:
            column_name = "name"
            parameter = name

        sql = """SELECT
                    count(*)
                 FROM
                    transform_kernels
                 WHERE
                    %s = ?
              """ % column_name
        return self._execute(sql, parameter, fetch=_Fetch.SINGLETON)


    def count_pulse_designs(self, id_=None, name=None):
        """
        Given an id or a name, returns the number of pulse designs 
        matching that criterion. Exactly one of id_ or name must be supplied.
        """
        if id_:
            column_name = "id"
            parameter = id_
        else:
            column_name = "name"
            parameter = name

        sql = """SELECT
                    count(*)
                 FROM
                    pulse_designs
                 WHERE
                    %s = ?
              """ % column_name
        return self._execute(sql, parameter, fetch=_Fetch.SINGLETON)


    def fetch_machine_specs(self, id_):
        """
        Given the id of a machine specs object, returns that object. It
        doesn't matter if it's a template or not.
        """
        
        sql = """SELECT
                    *
                 FROM
                    machine_specs
                 WHERE
                    id = ?
              """                                      
        row = self._execute(sql, id_, _Fetch.ONE)
        
        if row["is_template"]:
            return rfp_machine_specs.MachineSpecsTemplate(row)
        else:
            return rfp_machine_specs.MachineSpecs(row)
            

    def insert_pulse_design(self, pulse_design, use_transactions=True, results_to_skip=[ ]):
        """
        Given a populated pulse_design object, creates database rows
        representing the new pulse_design and its associated transforms
        (and results - a.k.a. waveforms and gradients).
        
        results_to_skip is a list of result objects. Any result present in 
        this list will be skipped; i.e. not written to the database. This 
        allows Pulse to ignore out-of-sync results without destroying them.
        """        
        if use_transactions:
            self.begin_transaction()

        # These next two if's should not be required
        if not pulse_design.id:
            pulse_design.id = util_misc.uuid()
        if not pulse_design.created:
            pulse_design.created = util_time.now()

        machine_specs_id  = self.insert_machine_specs(pulse_design.machine_specs)
        
        sql = '''INSERT INTO
                    pulse_designs
                        (id, 
                        is_public, 
                        name, 
                        creator, 
                        created,  
                        comment, 
                        calc_resolution, 
                        pulse_bandwidth_type, 
                        machine_specs_id,
                        gyromagnetic_nuclei)
                 VALUES
                    (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
              '''                        
              
        sql_params = (pulse_design.id, 
                      pulse_design.is_public , 
                      pulse_design.name,
                      pulse_design.creator, 
                      pulse_design.created, 
                      pulse_design.comment,
                      pulse_design.calc_resolution, 
                      pulse_design.pulse_bandwidth_type,
                      machine_specs_id,
                      pulse_design.gyromagnetic_nuclei)

        self._execute(sql, sql_params, _Fetch.NONE)

        # Add Transforms
        self._insert_transforms(pulse_design.id, 
                                pulse_design.transforms, 
                                results_to_skip)

        if use_transactions:
            self.commit()

        return pulse_design.id        


    def insert_rfpulse_convert_pulse_design(self, pulse_design, use_transactions=True, results_to_skip=[ ]):
        """
        This is very similar to insert_pulse_design() except that is is only
        called from db_upgrader.DatabaseUpgrader._upgrade_8() where we upgrade 
        from database version 8 to 9 and we need to convert RFPulse pulse_project 
        objects to Pulse pulse_design objects. In this step we don't have the
        gyromagneti_nuclei column in the pulse_design table yet. So, here we
        use the older call to insert the converted pulse_designs.
        
        Given a populated pulse_design object, creates database rows
        representing the new pulse_design and its associated transforms
        (and results - a.k.a. waveforms and gradients).
        
        results_to_skip is a list of result objects. Any result present in 
        this list will be skipped; i.e. not written to the database. This 
        allows Pulse to ignore out-of-sync results without destroying them.
        """        
        if use_transactions:
            self.begin_transaction()

        # These next two if's should not be required
        if not pulse_design.id:
            pulse_design.id = util_misc.uuid()
        if not pulse_design.created:
            pulse_design.created = util_time.now()

        machine_specs_id  = self.insert_machine_specs(pulse_design.machine_specs)
        
        sql = '''INSERT INTO
                    pulse_designs
                        (id, 
                        is_public, 
                        name, 
                        creator, 
                        created,  
                        comment, 
                        calc_resolution, 
                        pulse_bandwidth_type, 
                        machine_specs_id)
                 VALUES
                    (?, ?, ?, ?, ?, ?, ?, ?, ?)
              '''                        
              
        sql_params = (pulse_design.id, 
                      pulse_design.is_public , 
                      pulse_design.name,
                      pulse_design.creator, 
                      pulse_design.created, 
                      pulse_design.comment,
                      pulse_design.calc_resolution, 
                      pulse_design.pulse_bandwidth_type,
                      machine_specs_id)

        self._execute(sql, sql_params, _Fetch.NONE)

        # Add Transforms
        self._insert_transforms(pulse_design.id, 
                                pulse_design.transforms, 
                                results_to_skip)

        if use_transactions:
            self.commit()

        return pulse_design.id     


    def fetch_pulse_design(self, pulse_design_id):
        """Given a pulse design id, returns the associated pulse_design object."""
        sql = """SELECT
                    *
                 FROM
                    pulse_designs
                 WHERE
                    id = ?
              """
        row = self._execute(sql, pulse_design_id, fetch=_Fetch.ONE)

        # Will inflate row into a pulse_design.
        pulse_design = rfp_pulse_design.PulseDesign(row)

        if row["machine_specs_id"]:                        
            machine_specs_id = row["machine_specs_id"]
            pulse_design.machine_specs = self.fetch_machine_specs(machine_specs_id)

        # Fetch transforms for this pulse_design
        sql =   """SELECT
                   id, progression, transform_kernel_id, rf_result_id
                 FROM
                    transforms
                 WHERE
                    pulse_design_id = ?
                 ORDER BY
                    progression
              """
        trows = self._execute(sql, pulse_design_id, fetch=_Fetch.ALL)
        transforms = [ ]

        for row in trows:
            transform_kernel = self.fetch_transform_kernel(row["transform_kernel_id"])
            
            parameters = self._fetch_transform_parameters(row['id'])
            
            rf_result_id = row["rf_result_id"]
            if rf_result_id:
                rf_result = self._fetch_rf_result(rf_result_id) 
            
            transform = rfp_transform.Transform({'transform_kernel' : transform_kernel,
                                                 'parameters' : parameters,
                                                 'result' : rf_result })
            transforms.append(transform)

        pulse_design.transforms = transforms  
        pulse_design.referrers = self.fetch_pulse_design_referrers(pulse_design_id)

        return pulse_design  
        

    def fetch_pulse_design_preview(self, id_):
        """Returns a preview of the pulse design matching the given UUID"""
        previews = self.fetch_pulse_design_previews(id_)
        
        return (previews[0] if previews else None)
            

    def fetch_pulse_design_previews(self, id_=None):
        """
        Returns a list of pulse design previews. If the id_ parameter is 
        populated, only designs with that UUID are returned. Of course there
        will be at most one matching design, but it's still returned in a 
        list.
        """
        parameters = [ ]
        where_clause = ""
        if id_:
            where_clause = " AND id = ? "
            parameters = (id_, )

        sql = """SELECT
                    id, name, creator, created, is_public, comment
                 FROM
                    pulse_designs
                 WHERE
                    1 = 1  %s
                 ORDER BY
                    name
              """ % where_clause

        rows = self._execute(sql, parameters)
        
        previews = [rfp_pulse_design_preview.PulseDesignPreview(row) for row in rows]

        for preview in previews:
            preview.referrers = self.fetch_pulse_design_referrers(preview.id)
            
        return previews


    def fetch_pulse_design_referrers(self, id_):
        """
        Given a pulse design id, returns a (possibly empty) list of the 
        ids of pulse sequences that use this pulse design. 
        """
        sql = """SELECT
                    id, name
                 FROM
                    pulse_sequences
                 WHERE
                    id IN
                        (SELECT
                            pulse_sequence_id
                         FROM
                            pulse_sequence_pulse_designs
                         WHERE
                            pulse_design_id = ?
                        )
              """
        return [(row[0], row[1]) for row in self._execute(sql, id_)]    


    def _fetch_transform_parameters(self, transform_id):
        """
        Currently this code is only called from within a transaction, so
        it doesn't need to open a transaction itself.

        Given a single transform id, delete the transform parameters 
        associated with it. 
        
        """
        # Fetch transform parameters
        sql = """SELECT
                    *
                 FROM
                    transform_parameters
                 WHERE
                    transform_id = ?
                 ORDER BY
                    sort_order
                 """
        rows = self._execute(sql, transform_id)
        rows = [dict(row) for row in rows]
        parameters = [rfp_transform.TransformParameter(row) for row in rows]
        
        return parameters
        

    def _insert_transforms(self, pdid, transforms, results_to_skip=[ ]):
        # Currently this code is only called from within a transaction, so
        # it doesn't need to open a transaction itself.
        
        # See insert_pulse_design() for an explanation of results_to_skip
        
        progression = 1
            
        for transform in transforms:

            transform_kernel_id = transform.transform_kernel.id
            
            if transform.result and (transform.result not in results_to_skip):
                rf_result_id = self._insert_rf_result(transform.result)
            else:
                rf_result_id = 0
                
            sql = '''INSERT INTO
                        transforms
                            (pulse_design_id,
                             progression,
                             transform_kernel_id,
                             rf_result_id)
                     VALUES
                        (?, ?, ?, ?)                                      
                  '''
            sql_params = (pdid, 
                          progression,
                          transform_kernel_id,
                          rf_result_id)
    
            transform_id = self._execute(sql, sql_params, _Fetch.LAST_ID)

            # Now insert the transform parameter values
            sql = '''INSERT INTO
                        transform_parameters
                            (transform_id, 
                             variable,
                             type,
                             value, 
                             sort_order)
                     VALUES
                        (?, ?, ?, ?, ?)
                  '''
            for i, item in enumerate(transform.parameters):
                sql_params = (transform_id, 
                              item.variable,
                              item.type,
                              str(item.value),  i)
                
                self._execute(sql, sql_params, _Fetch.NONE)
            
            progression += 1     


    def _fetch_rf_result(self, rf_result_id):    
        
        sql = '''SELECT 
                    *
                 FROM   
                    rf_results
                 Where 
                    id = ?
              '''

        row = self._execute(sql, rf_result_id, fetch=_Fetch.ONE)           

        created = row["created"]  

        if row["opcon_state_id"]:
            opcon_state = self._fetch_opcon_state(row["opcon_state_id"])
        else:
            opcon_state = None
 
        rf_waveform = _blob_to_numpy_array(row["rf_waveform"],constants.DataTypes.COMPLEX128)
        rf_xaxis    = _blob_to_numpy_array(row["rf_xaxis"],constants.DataTypes.FLOAT64)
        gradient    = None
        grad_xaxis  = None
        
        if row["gradient"]:
            gradient  = _blob_to_numpy_array(row["gradient"],constants.DataTypes.FLOAT64)
        if row["grad_xaxis"]:
            grad_xaxis = _blob_to_numpy_array(row["grad_xaxis"],constants.DataTypes.FLOAT64)

        return rfp_rf_result.RfResults({ "rf_waveform" : rf_waveform, 
                                        "rf_xaxis"    : rf_xaxis, 
                                        "gradient"    : gradient,
                                        "grad_xaxis"  : grad_xaxis,
                                        "opcon_state" : opcon_state,
                                        "created"     : created})    
                                        

    def _insert_rf_result(self, rf_result):    
        # Currently this code is only called from within a transaction, so
        # it doesn't need to open a transaction itself.
        
        rf_waveform = _numpy_array_to_blob(rf_result.rf_waveform)
        rf_xaxis    = _numpy_array_to_blob(rf_result.rf_xaxis)
        gradient    = None
        grad_xaxis  = None
        
        if rf_result.gradient is not None:
            gradient = _numpy_array_to_blob(rf_result.gradient)
        if rf_result.grad_xaxis is not None:
            grad_xaxis = _numpy_array_to_blob(rf_result.grad_xaxis)
        
        opcon_state_id = self._insert_opcon_state(rf_result.opcon_state)
 
        # Insert information about the container object.      
        sql = '''INSERT INTO
                    rf_results
                        (created,
                         rf_waveform,
                         rf_xaxis,
                         gradient,
                         grad_xaxis,
                         opcon_state_id)
                 VALUES
                    (?, ?, ?, ?, ?, ?)                                      
              '''
      
        sql_params = (rf_result.created, rf_waveform, rf_xaxis, gradient, grad_xaxis, opcon_state_id)      

        rf_result_id = self._execute(sql, sql_params, _Fetch.LAST_ID)     
        
        return rf_result_id    


    def _insert_opcon_state(self, opcon_state):

        if not opcon_state:
            return 0

        sql = '''INSERT INTO
                    opcon_states
                       (multiplier,
                        met_max_iterations, 
                        met_residual_error,
                        met_differential_error,
                        met_max_time, 
                        met_increasing_error,
                        run_time, 
                        iterations, 
                        decreases
                       )
                 VALUES
                    (?, ?, ?, ?, ?, ?, ?, ?, ?)
              '''

        params = (opcon_state.multiplier,
                  opcon_state.met_max_iterations, 
                  opcon_state.met_residual_error,
                  opcon_state.met_differential_error,
                  opcon_state.met_max_time, 
                  opcon_state.met_increasing_error,
                  opcon_state.run_time, 
                  opcon_state.iterations, 
                  opcon_state.decreases
                 )      

        opcon_state_id = self._execute(sql, params, _Fetch.LAST_ID)

        # Now insert the residual error history
        sql = '''INSERT INTO
                    opcon_residual_errors
                        (opcon_state_id, value, progression)
                 VALUES
                    (?, ?, ?)
              '''
        for i, residual_error in enumerate(opcon_state.residual_errors):
            self._execute(sql, (opcon_state_id, residual_error, i) , _Fetch.NONE)
        
        # Now insert all the individual data points for deltab1.
        # deltab1 is a numpy array of complex numbers
        sql = '''INSERT INTO
                    opcon_deltab1_points
                        (opcon_state_id, progression, 
                         real_amplitude, imaginary_amplitude)
                 VALUES
                    (?, ?, ?, ?)                                      
              '''
        for progression, z in enumerate(opcon_state.deltab1):
            params = (opcon_state_id, progression, z.real, z.imag)
            self._execute(sql, params, _Fetch.NONE)
        
        return opcon_state_id    
            

##################    Module-level  Internal  Stuff    ##################


def _build_placeholders(length):
    """A convenience function for building a list of parameter placeholders
    for insertion into SQL."""
    return ",".join( ["?"] * length)

def _adapt_boolean(b):
    return 1 if b else 0        # bjs, This is OK Py2 to Py3

def _convert_boolean(s):
    return (s in ['1', b'1'])   # bjs, change for Py2 to Py3
    # Py2 bjs return (s == '1')


class _BetterRow(sqlite.Row):
    """An enhanced version of the sqlite3.Row object which is documented here:
    http://docs.python.org/release/2.5.4/lib/node352.html
    
    The sqlite3.Row class is a tuple/dict hybrid with some of the features
    of both and all the convenience of neither. Tuplishly, it supports
    index access (e.g. row[3]). Dictishly, it also supports access by 
    column name (e.g. row["foo"]). However, it doesn't support iteration over 
    its members (which both tuples and dicts support) nor does it offer handy 
    dict methods like .keys(), .get(), etc. 

    This subclass of sqlite3.Row adds support for tuple-style iteration over
    its members. It also adds dict-like membership testing, a .get() method,
    and a .keys() method.
    
    The row is read-only. It can be converted into a standard Python 
    dictionary by simply passing to the Python's dict constructor, e.g.
    row = dict(row).
    
    Note that when iterating over a tuple, values are returned. When iterating
    over a dict, keys are returned. This tuple/dict hybrid behaves like a 
    tuple during iteration.
    """
    def __init__(self, *args):
        # bjs-2to3 sqlite.Row.__init__(self, *args)
        sqlite.Row.__init__(self)
        # _column_names is populated by the Database._execute() method.
        self._column_names = ( )


    def __contains__(self, key):
        """Support membership testing, e.g. if "foo" in row """
        rc = True
        try:
            self[key]
        except IndexError:
            rc = False

        return rc


    def __iter__(self):
        """Support iteration, e.g. for column_value in row"""
        for i in range(len(self)):
            yield self[i]


    def __str__(self):
        # Printing the data as a tuple mirrors the behavior of the
        # sqlite.Row object.
        columns = [ ]
        for column in self:
            if isinstance(column, str):
                column = column.encode("utf-8")
                
            columns.append(column)
        
        return str(tuple(columns))            


    def __unicode__(self):
        # Printing the data as a tuple mirrors the behavior of the
        # sqlite.Row object. 
        # Note that __unicode__() must return a Unicode object. It's not 
        # necessary to pass an encoding to unicode() here because the 
        # column values are already unicode and won't need to be converted.
        # The only part that's being unicodified are the parens and commas
        # in the string representation of the tuple.
        return str(tuple([column_value for column_value in self]))
            

    def keys(self):
        """Returns a tuple or list of column names."""
        return self._column_names


    def get(self, key, default=None):
        """Implements dict-style .get()"""
        return self[key] if (key in self) else default
