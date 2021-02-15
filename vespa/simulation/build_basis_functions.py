# Python modules

import multiprocessing
import math

# 3rd party modules
import numpy as np

# Our modules
import vespa.common.constants as common_constants
import vespa.common.util.config as util_config
import vespa.common.util.misc as util_misc


"""Code for building basis functions. This code is independent of the GUI
and can be called from the command line. There's only one public function
which is build_basis_functions().

This code iterates over each simulation in the experiment. For each 
simulation, it manipulates the numbers that represent the spectrum lines.

It uses multiprocessing.

"""

# The naive algorithm for this code would be to process one simulation at
# a time. We found that experiments are processed faster if we combine the
# simulation's  lines into chunks prior to processing. 
#
# MAX_CHUNK_SIZE is the largest # of lines that will be grouped together
# in one chunk. This is a hard limit; it's never exceeded except in the
# case where one simulation contains > MAX_CHUNK_SIZE lines.
#
# Choosing MAX_CHUNK_SIZE was the result of a lot of experimentation. 
# Obviously we can't test all possibilities but 75 gave good average 
# performance. Our original chunk size of 5000 caused out-of-memory errors 
# in some situations, so be careful with using big numbers.
#
# To force this code to use the naive algorithm (each simulation = one chunk), 
# set MAX_CHUNK_SIZE = 1.
MAX_CHUNK_SIZE = 75

# CPU_COUNT is the number of processing cores (virtual CPUs) available on 
# this machine. We ask multiprocessing.Pool() to create CPU_COUNT workers.
# CPU_COUNT can be determined from a variety of sources. Users can specify
# it in vespa.ini file, and if they do so that trumps all other sources. We
# accept any int > 0.
CPU_COUNT = 0
config = util_config.VespaConfig()
if "general" in config:
    if "cpu_limit" in config["general"]:
        cpu_count = config["general"]["cpu_limit"]
        
        if util_misc.is_intable(cpu_count):
            cpu_count = int(cpu_count)
            if cpu_count > 0:
                CPU_COUNT = cpu_count

if not CPU_COUNT:
    # OK, the user didn't specify so we ask multiprocessing. Note that 
    # multiprocessing.cpu_count() isn't implemented on all platforms.
    # Where it's not implemented we default to 2 for no really strong reasons.
    try:
        CPU_COUNT = multiprocessing.cpu_count()
    except NotImplementedError:
        CPU_COUNT = 2


PI = math.pi

DEGREES_TO_RADIANS_C = common_constants.DEGREES_TO_RADIANS * 1j

# These are some values that need to be global to this module in order to make
# them accessible to child processes. 
global_resppm = None
global_npts = 0
global_sw = 0
global_const1 = 0
global_const2 = 0
global_nlines_total = 1
# global_progress is a 2 tuple of (progress as a float, progress as string).
# The second part of the tuple is just a formatted version of the first part
# provided for convenience.
# This global is updated as the code works.
global_progress = (0, "0%")



class _Chunk(object):
    """A chunk of lines to process. nlines is always < MAX_CHUNK_SIZE 
    unless a single simulation contains more lines 
    ."""
    def __init__(self, start=0):
        # Each chunk has start, end, ilines, ppms, areas, phases, nlines.
        # start and end track what portion of the total simulation array
        # this chunk contains. This is always true for both start and end:
        #    0 <= N < len(experiment.simulations)
        # Also:
        #    sum(ilines) == nlines
        # Also:
        #    nlines == self.ppms.size
        # Also:
        #    self.ppms.size == self.areas.size == self.phases.size
        self.start = start
        self.end = 0
        self.ilines = [0]
        self.ppms   = np.ndarray(0)
        self.areas  = np.ndarray(0)
        self.phases = np.ndarray(0)

    @property
    def nlines(self):
        return self.ppms.size
        
    def __str__(self):
        # __str__ is useful for debugging
        lines = [ ]
        lines.append("---- chunk ----")
        lines.append("range: %d - %d" % (self.start, self.end))
        lines.append("ilines: %s" % self.ilines)
        lines.append("nlines: %d" % self.nlines)
        lines.append("ppms: %d" % self.ppms.size)
        lines.append("areas: %d" % self.areas.size)
        lines.append("phases: %d" % self.phases.size)
        
        return "\n".join(lines)
        


def build_basis_functions(experiment, dim0, sw, resppm, progress_callback=None, local_cpu=None):
    """
    Builds basis data for the experiment and freq array passed. Uses
    multiprocessing. Returns a numpy array appropriate for the data element 
    of an DataBasis object.
    
    progress_callback should be a function that accepts no args. It's called
    by this code when it begins processing a new chunk.
    
    Note that under Windows, the callback function is so limited (e.g. calls
    to wx cause hangs) so as to be useless. It's recommended that Windows
    callers pass None for progress_callback. 
    """
    global global_nlines_total

    # allow for single cpu in sequence editor dialog
    if local_cpu is not None:
        n_cpu = local_cpu
    else:
        n_cpu = CPU_COUNT


    npts = dim0
    hpp  = sw/npts

    # Calculate number of simulations
    nsims = len(experiment.metabolites)
    for dim in experiment.dims:
        nsims *= len(dim)
        
    # I calculate the total number of lines so that I can report progress.
    # Even on a simulation with ~2.6m lines, this calculation takes just 
    # one tenth of a second on my laptop so it's definitely worthwhile.
    nlines_total = sum([len(simulation.ppms) for simulation 
                                                 in experiment.simulations])
                                                     
    global_nlines_total = nlines_total              

    # const1 and const2 are loop invariants that we calculate here.
    const1 = experiment.b0 / hpp
    const2 = PI * 2 * hpp * 1j

    # PS - If you want to run this code without using multiprocessing (e.g.
    # in order to profile execution), use the 3 lines below in place of
    # the use of multiprocessing.Pool.
    # _initializer(freq1d.resppm, npts, freq1d.sw, const1, const2, nlines_total)
    # chunks = _build_chunks(experiment, progress_callback)
    # results = [_process_chunk(chunk) for chunk in chunks]

    pool = multiprocessing.Pool(n_cpu, _initializer, [resppm, npts, sw, const1, const2, nlines_total])

    
    # The 3rd param to imap_unordered() is a chunksize. These chunks are not
    # to be confused with the chunks returned by _build_chunks()! chunksize
    # just determines how many values will be grabbed from the iterator 
    # at once. Using a chunksize > 1 gives slightly better performance, but
    # only slightly. The advantage of using a chunksize == 1 is that 
    # _build_chunks() is called every time a worker needs new work, so we
    # can use it as a cheap callback/progress indicator.
    results = pool.imap_unordered(_process_chunk, 
                                  _build_chunks(experiment, progress_callback), 
                                  1)
    
    pool.close()
    pool.join()
    
    # The lines from each simulation are combined into one uber-FID
    basis_data = np.empty((nsims, npts), dtype='complex64')
    for result in results:
        # Each result is a list of 2-tuples. Each 2-tuple is (index, sum)
        for index, xx_sum in result:
            basis_data[index] = xx_sum
            
    shape = [len(dim) for dim in experiment.dims[::-1]]
    shape += [len(experiment.metabolites)]
    shape += [npts]

    basis_data.shape = shape
    
    return basis_data


###############       Internal use only below this line     ###############

def _build_chunks(experiment, progress_callback):
    """A generator function. Given an experiment, iterates over the 
    experiment's simulations and returns lines within those simulations
    chunked according to MAX_CHUNK_SIZE.
    
    See here for more info on generators:
    http://docs.python.org/tutorial/classes.html#generators
    """
    current = _Chunk(0)
    nlines_processed = 0
    for isim, simulation in enumerate(experiment.simulations):
        nlines = len(simulation.ppms)
        
        if current.nlines and ((current.nlines + nlines) > MAX_CHUNK_SIZE):
            # This chunk already has stuff in it, and adding more would make 
            # it exceed the max.
            nlines_processed += current.nlines
            _set_progress(nlines_processed)
            if progress_callback:
                progress_callback()
            
            yield current
            current = _Chunk(isim)
        #else:
            # The current chunk is empty or there's still room in the current
            # chunk for the next collection of lines. 
            # One might think that if the former is true then the latter would
            # automatically be true too. But that's not the case when
            # len(simulation.ppms) > MAX_CHUNK_SIZE. It makes no sense to 
            # create an empty chunk, so we just create a chunk that's 
            # larger than our preferred size.

        # Append the contents of this simulation to the current chunk.
        current.ilines.append(nlines)
        current.end = isim
        current.ppms = np.concatenate((current.ppms, simulation.ppms))
        current.areas = np.concatenate((current.areas, simulation.areas))
        current.phases = np.concatenate((current.phases, simulation.phases))

    # Return the final set of lines.
    yield current


def _initializer(resppm, npts, sw, const1, const2, nlines_total):
    # This function is subtle...it's called by each worker process, and is
    # passed the values of the global constants that I need in 
    # _process_chunk(). Under *nix, I can just declare them global and 
    # (thanks to the magic of fork()) the variables and their values will be 
    # copied to the worker processes'. Under Windows, this module is 
    # re-imported once for each worker, and as a result these globals are
    # recreated and re-initialized to 0 in each worker. This function sets 
    # them back to the values they need to be, and that's the only reason
    # it exists.
    global global_resppm
    global global_npts
    global global_sw
    global global_const1
    global global_const2
    global global_nlines_total
    
    global_resppm = resppm
    global_npts = npts
    global_sw = sw
    global_const1 = const1
    global_const2 = const2
    global_nlines_total = nlines_total

    
def _process_chunk(chunk):
    # This is what each worker executes
    resppm = global_resppm
    npts = global_npts
    sw = global_sw
    const1 = global_const1
    const2 = global_const2
    
    # convert freq(ppm)/phase(deg) to freq(hz)/phase(rad)
    nlines = chunk.nlines

    freq = ((npts / 2) - (chunk.ppms - resppm) * const1)  # ppm2pts
    freq = freq * const2
    freq = np.repeat(freq, npts).reshape(nlines, npts)
    phas = chunk.areas * np.exp(chunk.phases * DEGREES_TO_RADIANS_C)
    phas = np.repeat(phas, npts).reshape(nlines, npts)

    # calculate the FIDs for each line
    # xx stands for "X axis". It's the dwell time in seconds.
    xx = (np.arange(npts, dtype='float32') % npts) / sw
    xx = np.tile(xx, nlines)
    xx.shape = nlines, npts

    xx = freq * xx

    xx = phas * np.exp(xx)

    # The results are a list of 2-tuples (index, xx). index is an index 
    # into the basis array -- it's where this xx will reside in the
    # uber-FID.
    result = [ ]
    for i in range((chunk.end - chunk.start) + 1):
        start = sum(chunk.ilines[:i + 1])
        end   = sum(chunk.ilines[:i + 2])
    
        result.append( (chunk.start + i, xx[start:end,:].sum(axis=0)) )
    
    return result


def _set_progress(processed):
    global global_progress

    # When there will be < 100 chunks, we can use an integer to indicate
    # % progress. Otherwise we need to use a float.
    formatter = "%.1f" if (global_nlines_total > MAX_CHUNK_SIZE * 100) else "%d"
    formatter += "%%"

    local_progress = (processed / global_nlines_total) * 100
    global_progress = (local_progress, formatter % local_progress)

