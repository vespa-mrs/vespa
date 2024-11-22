# Python modules

import sys
import tempfile
import importlib
import os
import traceback
import inspect
from codecs import BOM_UTF8

# third party modules

# Our modules
import vespa.common.transform_run_exception as transform_run_exception
import vespa.common.util.time_ as util_time
import vespa.common.util.misc as util_misc
import vespa.public.transform_description as transform_description


# A proper result tuple must have 4 entries.
RESULT_LENGTH = 5

# Note that bool(EMPTY_RESULT) is True because the outer tuple is not empty.
EMPTY_RESULT = tuple([tuple()] * RESULT_LENGTH)

# "run" is the name of the function that we expect in the user's algorithm
FUNCTION_NAME = "run"

# Aliases are the strings that appear in error tracebacks in place of the
# ugly names of the temp files to which we save the user's algorithm code.
ALGORITHM_CODE_ALIAS = "~ your algorithm code ~"
GENERIC_CODE_ALIAS = "~ your transform code ~"


# These three global variables allow this code to make values persistent
# across function calls.
_algorithm_function = None


#    = About Error Handling (unrelated to distributed processing) =
#
# Since run_transform() runs arbitrary code, it can raise any error
# imaginable. Rather than forcing the caller to surround the whole call with
# a try/except that catches all errors, run_transform() traps any errors that
# arise, reformats them a little and passes them back to the caller. See
# _repackage_exception() for details about how & why we reformat the
# exceptions.

# We call sys.exc_info() in lots of places in this code, and that returns
# a tuple that includes a traceback. The doc for sys.exc_info() includes a
# warning --
#    Assigning the traceback return value to a local variable in a
#    function that is handling an exception will cause a circular
#    reference...delete [the traceback] after use...or call exc_info()
#    in a function that does not itself handle an exception.
# I find it easier to do the latter, so if the try/excepts in this code
# look at little awkward, that's probably why.

# If you're modifying this code, here's some testing notes.

# There are eight cases in which this code traps exceptions. Five of those
# cases can occur when executing algorithm code. Since the error
# returned specifies which module was at fault, each of those five cases must
# be tested for algorithm code. That's a total of
# 5 + 3 = 8 cases that need to be unit tested for this code.
# They need to be tested in the transform kernel editor dialog so that you 
# can verify that the proper error message is presented to the user.

#              = Test Cases for Exception Handling =
#
#             == Tests Related to Running the Code ==
# These need to be tested for the algorithm code.
# Exceptions 2, 3 and 4 are raised from within our code.
#
# 1) Exception: code can't be imported.
# Suggested way to create this exception:
# Add this just after 'import pygamma':
# this_garbage_raises_an_error
#
# 2) Exception: run function not present
# Suggested way to create this exception:
# rename run to xrun
#
# 3) Exception: run function must take just one param
# Suggested way to create this exception:
# change this:
# def run(trans_desc):
# to this:
# def run(x, trans_desc):
#
# 4) Exception: run function isn't callable
# Suggested way to create this exception:
# Add this at the end of the code:
# run=42
#
# 5) Exception: error in code during execution (as opposed to during import)
# Suggested way to create this exception:
# Add this inside the run function to generate a division by zero error:
# 42/0
#

#                == Tests Related to the Results ==
#
# 1) Exception: results must be a 3-tuple (or list)
# Suggested way to create this exception:
# Change binning code return to this:
# return None
#
# 2) Exception: result elements must be of the same length
# Suggested way to create this exception:
# Change binning code return to this:
# return ( range(1), range(2), range(2) )
#
# 3) Exception: Results must contain only Python native numeric types
# Suggested way to create this exception:
# Change binning code return to this:
# return ( ["xyz"], [0], [0] )

class NotRunnableError(Exception):
    """A custom error raised if any of our three criteria for the callable
    interface to the sequence and binning code are violated:
    - If the algorithm code doesn't contain a run() function
    - If run() exists but isn't callable
    - If run requires more or fewer than 1 parameter.
    """
    pass


##############    Public  Functions     ############


def exception_to_message(exception):
    """
    Given an exception 3-tuple as returned by run_transform(), returns
    a string formatted as Python would have formatted it.

    This is a convenience function for handling the somewhat odd
    exception info returned by run_transform().
    """
    # trace is a list produced by traceback.extract_tb().
    type_, value, trace = exception
    trace = traceback.format_list(trace)

    # We format the exception just like Python would.
    trace = "Traceback (most recent call last):\n" + "\n".join(trace)
    message = traceback.format_exception_only(type_, value)
    message = trace + "\n".join(message)

    return message


def run_transform(transform, previous=None, extra=None):
    """
    Given a transform, runs the algorithm code. The transform is not
    modified.

    When the previous parameter is populated (non-empty), the previous
    transforms result (waveform and gradient) are also passed into the
    algorithm. This is necessary for "Modify" type transforms but not
    for "Create" type transforms.

    This function is written to trap all exceptions and return them
    repackaged (see below). Despite the fact that it's running arbitrary user
    code, it's designed to never raise an exception and so the caller doesn't
    need to surround it with try/except.

    The function returns a two-tuple of (results, exception). The results
    are a list of dicts, one for each simulation. The dicts are described
    in detail below.

    The second tuple item returned ("exception") is a little complicated. It's
    None if nothing failed. Otherwise it's a 2-tuple of (exception, None) or
    (exception, result dict). It's the former if the exception occurred
    *outside* of the simulations loop (e.g. if the user's code failed to
    compile). It's the latter if the exception occurred *inside* the
    simulations loop. In that case, the second item is a reference to the
    result dict from the algorithm that failed. This allows the calling code
    to report to the user the inputs under which the exception occurred.

    Also note that the exception is a repackaged version of the values
    returned by sys.exc_info(). The convenience function
    exception_to_message() can turn that exception into a typical Python error
    message, which is probably what you want.

    If you're interested in the details of how the exeception info is
    repackaged, see the function _repackage_exception() in this module.

    To summarize, this function can return three things, with the most common
    case (success) first:
    results list, None
    results list, ((exc. type, exc. value, traceback list), None)
    results list, ((exc. type, exc. value, traceback list), one result dict)

    Each dict in the results has these keys: started, completed, rf_waveform,
    gradient. They point to these values:
    - started and completed are Python datetime.datetime() objects
    - rf_waveform and gradient are tuples or lists of numeric Python types
      (floats, ints, or longs)

    A result dict might also have an "exception" key. It's for the internal
    use of this function. Callers shouldn't rely on it nor infer anything
    from it.
    """

    results = [ ]

    # There are two possible sources of exceptions in this code, and I must
    # track them separately. main_exception tracks whether or not an
    # exception has occurred here in the main process. The workers never
    # see main_exception.
    #
    # In contrast, _worker_exception is shared between this (the main)
    # process and the workers. This process creates the shared variable with
    # an initial value of 0. Thereafter, the main process monitors it but
    # never writes to it again. If a worker detects an exception, it writes
    # a non-zero value to the flag.
    main_exception = None
    _worker_exception = False

    # Here we create the algorithm module. If the user's code contains a
    # syntax error, the import will fail. If it fails, then there's no point
    # to the rest of this function.
    algorithm_module = None

    # We have to be prepared for any error here, hence the naked except.
    module_name = "transform_" + transform.transform_kernel.id
    try:
        algorithm_code   = transform.transform_kernel.algorithm_code
        algorithm_module = _create_temporary_module(module_name, "transform_", algorithm_code)

        spec = importlib.util.spec_from_file_location(algorithm_module.__name__, 
                                                      algorithm_module.__file__)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        sys.modules[algorithm_module.__name__] = module

        _algorithm_function = getattr(module, "run")
    except:
        main_exception = (_repackage_exception(ALGORITHM_CODE_ALIAS), None)


    if not main_exception:
        # OK, the algorithm module is constructed. Prepare to run it.
        trans_desc = _build_trans_descs(transform, previous, extra)

        # By the time this code is invoked, I am guaranteed that the algorithm
        # code is importable & callable, so I don't check for errors here.
        results, _worker_exception = _run_transform(trans_desc, _algorithm_function)
        results = [results]     # needed because transforms return only one result 

    #else:
        # exception occurred above

    if algorithm_module:
        _clean_up_temporary_module(algorithm_module)

    # return_exception holds the exception info we return to the caller.
    return_exception = None
    if main_exception:
        return_exception = main_exception
    else:
        # OK, the main process had no errors. Check to see if the worker
        # processes had any.
        if _worker_exception:
            # At least one of the results has an exception attached to it.
            # At present we only report one exception so if there's more than
            # one in the list of results, I ignore all except the first.
            for result in results:
                if "exception" in result:
                    # I found what I need.
                    # Package the exception in the way the caller expects.
                    return_exception = (result["exception"], result)
                    # Storing the exception in the result dict is a bit of a
                    # hack that I use to get the exception from the reaper back
                    # to the main thread. I clean it up here.
                    del result["exception"]
                    break
        #else:
            # All workers completed happily.
            
    results = [result for result in results]

    return results, return_exception


##############    Main  Process  Functions     ############
##############       (internal use only)       ############


def _build_trans_descs(transform, previous, extra):
    """
    A generator function that returns a description of each algorithm.

    See here for more info on generators:
    http://docs.python.org/tutorial/classes.html#generators

    """
    vespa_version = util_misc.get_vespa_version()       #.encode("utf-8")
    
    # I now keep each parameter in its own object, so need to convert to dict
    parameters = transform.parameters_to_dict()

    if previous is not None:
        previous_rf         = previous.rf_waveform
        previous_rf_xaxis   = previous.rf_xaxis
        previous_gradient   = previous.gradient
        previous_grad_xaxis = previous.grad_xaxis
    else:
        previous_rf = None
        previous_rf_xaxis = None
        previous_gradient = None
        previous_grad_xaxis = None

    # Build the trans_desc a.k.a. transform description
    # object that gets passed to the user's algorithm code
    trans_desc = transform_description.TransformDescription()
    trans_desc.vespa_version        = vespa_version
    trans_desc.extra                = extra
    trans_desc.parameters           = parameters
    trans_desc.previous_rf          = previous_rf
    trans_desc.previous_rf_xaxis    = previous_rf_xaxis
    trans_desc.previous_gradient    = previous_gradient
    trans_desc.previous_grad_xaxis  = previous_grad_xaxis
        
    return trans_desc


def _clean_up_temporary_module(module):
    # Nothing fancy here, just encapsulates the fact that two things have
    # to happen during module cleanup. (1) delete the temp file created to
    # hold the module's source code and (2) delete the module itself.
    filename = module.__file__
    os.remove(filename)

    del sys.modules[module.__name__]


def _create_temporary_module(module_name, description, code):
    """
    Given a module name (e.g. "my_sequence_code"), code description
    ("transform_"), and a Python code string that contains a
    function called "run" (a.k.a. FUNCTION_NAME), saves the code to a temp
    file and returns the module object.

    May raise NotRunnableError (qv). In addition, any error imaginable
    might be raised during the import step (since importing a module
    executes it).

    The caller is responsible for cleaning up the temp file.
    """
    # Convert the code from Unicode in order to make it bytestream
    # compatible and then write it to a temp file.
    # If there are non-ASCII characters in a .py file, one must make Python
    # aware of the encoding in which the file is written. There's two ways to
    # do that. The first is to put a comment in a specific format in the
    # first line or two of the file. That's the most clear way of doing it,
    # but if we prepend such a line to the user's code here, it means that
    # the line numbers associated with errors will be off by one.
    # Luckily there's an alternative which is to prefix the file with the
    # 3-byte UTF-8 BOM, so that's what we do.
    # ref: http://docs.python.org/reference/lexical_analysis.html#encoding-declarations
    code = BOM_UTF8 + code.encode("utf-8")
    fd, filename = tempfile.mkstemp(".py")
    os.write(fd, code)
    os.close(fd)

    module = None
    run_function = None

    # Import that temp file as a module. Of course, when Python imports a
    # module, the module is compiled and everything at the top level is
    # executed. Since arbitrary code is being executed, any Python error
    # imaginable is possible.
    # However, even though these errors are sure to occur now and again,
    # I don't use a try/except here. I just allow those errors to percolate
    # up to the caller so they can be reported to the user in a
    # context-appropriate fashion.
    # Note that importing the module here has the happy side effect of
    # compiling the module to a .pyc file. This speeds up subsequent imports.
    
    # module = imp.load_source(module_name, filename)
    spec = importlib.util.spec_from_file_location(module_name, filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    sys.modules[module_name] = module
    

    # If we've survived this far, I also evaluate the module to see that it
    # meets our criteria. Namely, if must offer a function called "run"
    # that accepts one param.
    # If it fails to meet these criteria, I raise a custom error. Using a
    # custom error allows the caller to distinguish this class of errors
    # from an error that occurred during the import step.
    if hasattr(module, FUNCTION_NAME):
        run_function = getattr(module, FUNCTION_NAME)
    else:
        msg = "The '%s' function is not present in your %s code." \
                                                % (FUNCTION_NAME, description)
        raise NotRunnableError(msg)

    if not callable(run_function):
        msg = "'%s' in your %s code is not callable." % (FUNCTION_NAME, description)
        raise NotRunnableError(msg)

    # Check the run function's args. It's allowed to have multiple args, but
    # only one can be required (i.e. all the other must have default
    # values).
    if sys.version_info[1] < 11:
        args, _, _, defaults = inspect.getargspec(run_function)
    else:
        args, _, _, defaults, _, _, _ = inspect.getfullargspec(run_function)

    # defaults is None if there's no default args.
    if not defaults:
        defaults = [ ]

    if len(args) - len(defaults) != 1:
        msg = "The '%s' function in your %s code must have exactly one required argument." \
                                                % (FUNCTION_NAME, description)
        raise NotRunnableError(msg)

    return module


##############    Main & Worker Process Functions     ############
##############       (internal use only)       ############


def _repackage_exception(code_name):
    """
    Given the name of some code (e.g. "~ your sequence code ~"), returns
    a 3-tuple of (exception type, exception value, traceback).

    The type and value are as described in the doc for sys.exc_info().
    However, the traceback has been manipulated a little. First, it's been
    turned from a traceback object into a list. Second, all instances of the
    temp file names to which the user's code is saved has been replaced by the
    code_name param.

    See comments inline to learn why we repackage the exceptions.
    """
    # We repackage exceptions for two reasons, both related to the traceback.

    # First, because user code is saved into a temp file, the error message
    # generated from an unmodified traceback looks something like this:
    #
    #    File "/var/folders/bB/bB1ZYFLPHVmSX3sKQzEQ8k+++TI/-Tmp-/tmpfVM73r.py", line 25, in run
    #      oops_its_an_error
    #    NameError: global name 'oops_its_an_error' is not defined
    # Needless to say, that file name is not too informative. Modifying the
    # file names in the traceback gives us the opportunity to replace the temp
    # file name with a more informative string.

    # Second, some of these exceptions (and therefore tracebacks) are
    # generated by code in a process other than the main one and are returned
    # to the main process as results. That means they must be pickleable.
    # Traceback objects can't be pickled; doing so raises "TypeError: can't
    # pickle traceback objects".
    # Therefore, if we try to return a raw traceback object, pickling will
    # fail. Repackaging them as a list retains the traceback info without
    # offending the pickle gods.
    type_, value, trace = sys.exc_info()

    # Turn the trace into a list.
    trace = traceback.extract_tb(trace)

    # Each entry in the traceback list is a 4-tuple --
    #    (filename, line number, function name, text)
    # We replace the filename as appopriate.
    trace = [_scrub_traceback_entry(code_name, entry) for entry in trace]

    return (type_, value, trace)


def _scrub_traceback_entry(code_name, entry):
    """
    Given the name of some code (e.g. "~ your sequence code ~") and
    an entry from the traceback list, scrubs any references to the names of
    the files that contain the user code and replaces them with code_name.
    """
    # The filenames in the traceback will point to either Pulse/Vespa
    # files or the user's pulse seq code. The latter is saved to the temp dir,
    # so if the temp dir is present in the filename, it must refer to the
    # user's code.
    if tempfile.gettempdir() in entry[0]:
        entry = tuple([code_name] + list(entry[1:]))

    return entry

##############    Worker Process Functions     ############
##############       (internal use only)       ############

def _run_transform(trans_desc, _algorithm_function):
    # Since _run_simulation() is always run in a separate process, its input
    # and output params must be pickle-friendly. Keep that in mind when
    # making changes.
    """
    This is what each worker executes.

    Given a SimulationDescription object, calls the sequence & binning
    code, traps any errors that arise and grabs results. Also verfies that
    the results meet our criteria (e.g. converts to tuples/lists if necessary,
    raises an exception if the ppms, areas and phases arrays are not all the
    same length, etc.)

    If an exception is raised at any point, it sets _worker_exception.

    Returns a result dict. If an exception occurred, the repackaged
    exception is in result["exception"].
    """
    started = util_time.now()

    exception = False
    user_exception = False
    # Execute the user's sequence code
    try:
        result = _algorithm_function(trans_desc)
    except transform_run_exception.TransformRunException as tre:
        user_exception = True
        exception = [tre, tre.message, '']
    except:
        # in Py3 sys.exc_info() resets when it leaves this try-except scope, need Vespa exception process call here
        exception = _repackage_exception(ALGORITHM_CODE_ALIAS)

    if not exception:
        # Sequence code completed OK.
        if result:
            # Sequence code returned the result. There's no need to
            # execute the binning code.
            pass

    if exception:
        result = EMPTY_RESULT
    else:
        # Execution completed with no errors. Let's see if what was returned
        # meets our criteria. First, the result must be an N-tuple, where
        # N == RESULT_LENGTH. As of this writing, RESULT_LENGTH == 2.
        result_length = _safe_len(result)
        if result_length != RESULT_LENGTH:
            result = EMPTY_RESULT
            # I force an error here so I can get the exception including a traceback.
            try:
                raise ValueError("Result returned from your code must be a %d-tuple, but has length %d" % \
                                   (RESULT_LENGTH, result_length))
            except ValueError:
                exception = _repackage_exception(GENERIC_CODE_ALIAS)

        # Our second criteria is that first two elements of the 5-tuple must be the
        # same length.
        lengths = [_safe_len(element) for element in result[0:2]]

        for length in lengths:
            if length != lengths[0]:
                result = EMPTY_RESULT
                # I force an error here so I can get the exception including a traceback.
                try:
                    raise ValueError("Result elements differ in length: %s" % lengths)
                except ValueError:
                    exception = _repackage_exception(GENERIC_CODE_ALIAS)

# For Kernels, we need Numpy arrays in first two results, then None can be in
# any of the remaining. If there are gradient results (indices 2 and 3) then
# they have to be Numpy arrays of the same length as the first two.
#
# FIXME bjs, Not ready to code this check in here yet
#
#     # The user's code is required to return a tuple of iterables. Those
#     # iterables might be lists, numpy arrays, PyGAMMA.DoubleVectors or any
#     # number of other things. PyGAMMA objects in particular don't pickle, and
#     # this function's result needs to be pickleable.
#     # So here we ensure that the result contains only ordinary Python objects.
#     result = map(_tuplify, result)
    
    
# FIXME bjs, ALSO Not ready to code this check in here yet
#
#     # Last but not least, ensure that each value is numeric and a native
#     # Python type.
#     f = lambda an_object: isinstance(an_object, (float, int, long, complex))
#     # Loop through ppms, areas & phases lists
#     for result_chunk in result:
#         # map() allows me to test all the items in the list in one shot.
#         if not all(map(f, result_chunk)):
#             # Ooops, at least one of the results in this list isn't a float,
#             # int, or long.
#             # I force an error here so I can get the exception including a traceback.
#             try:
#                 raise ValueError("Results must contain only floats, ints or longs")
#             except ValueError:
#                 exception = _repackage_exception(GENERIC_CODE_ALIAS)

    # The result (good or bad) is returned as a dict.
    result = dict(list(zip(("rf_waveform", "rf_xaxis", "gradient", "grad_xaxis", "outputs"), result)))
    result["started"] = started
    result["completed"] = util_time.now()

    if exception or user_exception:
        _worker_exception = True
        result["exception"] = exception
    else:
        _worker_exception = False

    return result, _worker_exception


def _safe_len(an_object):
    """Returns len(an_object), or 0 if the object doesn't support len()."""
    return len(an_object) if hasattr(an_object, "__len__") else 0


def _tuplify(an_iterable):
    """Given an iterable (list, tuple, numpy array, pygamma.DoubleVector,
    etc.), returns a native Python type (tuple or list) containing the same
    values.

    The iterable is converted to a tuple if it isn't already a tuple or list.

    If it is a tuple or list, it's returned unchanged.
    """
    if isinstance(an_iterable, (tuple, list)):
        # It's already a native Python type
        return an_iterable
    else:
        return tuple(an_iterable)

