"""
Helper routines for reading & writing MRS data to/from disk.
"""

# Python modules

import struct
import xdrlib
import itertools

# 3rd party modules

# Our modules
import vespa.common.constants as constants
import vespa.common.util.misc as util_misc



class UnreadableDataError(Exception):
    """ Raised when this module can't make sense of data (raw or XDR) due
    to an unexpected format, buffer underrun (less data than expected),
    buffer overrun (more data than expected), etc.
    """
    pass


class IncompleteMetadataError(Exception):
    """ Raised when the metadata associated with a dataset doesn't contain
    required information like the data's type or format.
    """
    pass



def dump_iterable(filename, an_iterable):
    """Given a filename and an iterable (tuple, list, numpy array), dumps the
    contents to the filename in XDR format. Multi-dimensional numpy arrays
    are flattened using ravel().

    This does not support iterator objects (like the object returned from
    reversed()).

    If filename exists, its contents are overwritten.

    This is meant to be a crude method for dumping raw numbers to disk. The 
    file that's written doesn't contain any hints about data type or the
    number of elements in the file, it's pure data.
    """
    if hasattr(an_iterable, "ravel"):
        # It's a numpy array
        an_iterable = an_iterable.ravel().tolist()

    if an_iterable:
        data_type = get_data_type(an_iterable)
        data = encode_xdr(an_iterable, data_type)
    else:
        data = ""

    open(filename, "wb").write(data)


def get_data_type(something, default=None):
    """Given an iterable (list, numpy array, etc.) or a scalar (int, float, 
    complex, etc.), determines the type and returns a constants.DataTypes
    value representing the data type.

    Callers should be aware that passing an empty list will raise a TypeError
    because an empty list contains no clues as to the type of data it can
    contain. This weakness is specific to standard Python types. Numpy arrays 
    are immune from this problem because they specify a data type even when the
    array is empty.

    If a default is specified, this code will return that instead of raising
    an error when it can't determine the type.
    """
    data_type = constants.DataTypes.NONE

    if default:
        # In case I use the default, I map it to a DataTypes value
        default = constants.DataTypes.any_type_to_internal(default)

    if hasattr(something, "dtype"):
        # It's a numpy array or a numpy elemental value (e.g. a single float)
        data_type = constants.DataTypes.any_type_to_internal(str(something.dtype))
    elif util_misc.is_iterable(something, False):
        # It's not a numpy array but it is iterable so it's something like a 
        # list or tuple.
        if len(something):
            something = something[0]
            # It's unusual but possible to encounter a regular Python list 
            # containing numpy values, so here we check again for a numpy 
            # object.
            if hasattr(something, "dtype"):
                data_type = str(something.dtype)
            else:
                # It's a non-numpy elemental type like a Python int or float.
                data_type = type(something)

            data_type = constants.DataTypes.any_type_to_internal(data_type)
        else:
            # empty list
            if default:
                data_type = default
            else:
                raise TypeError("Can't determine the data type of an empty list")
    else:
        # It's a non-numpy elemental type like a Python int or float.
        data_type = constants.DataTypes.any_type_to_internal(type(something))

    if not data_type:
        if default:
            data_type = default
        else:
            raise ValueError("Unable to determine the data type")

    return data_type




###################       collapse/expand complexes       ###################

def collapse_complexes(data, conjugate_flag=False):
    """Given a list or other iterable that's a series of (real, imaginary)
    pairs, returns a list of complex numbers. For instance, given this list --
       [a, b, c, d, e, f]
    this function returns --
       [complex(a, b), complex(c, d), complex(e, f)]

    The returned list is a new list; the original is unchanged.
    """
    # This code was chosen for speed and efficiency. It creates an iterator
    # over the original list which gets called by izip. (izip() is the same
    # as the builtin zip() except that it returns elements one by one instead
    # of creating the whole list in memory.)
    # It's the fastest method of the 5 or 6 I tried, and I think it is also
    # very memory-efficient. 
    # I stole it from here:
    # http://stackoverflow.com/questions/4628290/pairs-from-single-list
    data_iter = iter(data)

    if not conjugate_flag:
        tmp = [complex(r, i) for r, i in zip(data_iter, data_iter)]
    else:
        tmp = [complex(r, -1*i) for r, i in zip(data_iter, data_iter)]

    return tmp
    

def expand_complexes(data):
    """Expands a list or tuple of complex numbers into a list of
    (real, imaginary) pairs. For instance, given this list of complex
    numbers --
       [za, zb, zc]
    this function returns --
       [za.real, za.imag, zb.real, zb.imag, zc.real, zc.imag]

    The returned list is a new list; the original is unchanged.
    
    You can also pass a numpy array as long as it is one-dimensional. 
    """
    # First I double the length of the list by adding empty elements to
    # the front.
    data = ([None] * len(data)) + list(data)

    # Now I overwrite the items in the list with item N being split
    # into real and imag and list[N/2] = real and list[N/2 + 1] = imag.
    j = 0
    for i in range(len(data) // 2, len(data)):
        data[j] = data[i].real
        j += 1
        data[j] = data[i].imag
        j += 1

    return data


#####################       decode/encode XDR       #####################

def decode_xdr(data, data_type, element_count):
    """Given a string of data in XDR format and a data type, returns
    an iterable (tuple or list) of Python objects representing the decoded
    data. data_type must be one of the constants from
    vespa.common.constants.DataTypes.ALL.
    
    element_count is the number of elements expected in the data.
    
    If the data is complex, the element count should be the number of complex
    numbers expected (despite the fact that XDR doesn't understand complex
    numbers).

    The companion function to this is encode_xdr().
    """
    p = xdrlib.Unpacker(data)

    if data_type in (constants.DataTypes.FLOAT32, constants.DataTypes.COMPLEX64):
        unpack_function = p.unpack_float
    elif data_type in (constants.DataTypes.FLOAT64, constants.DataTypes.COMPLEX128):
        unpack_function = p.unpack_double
    elif data_type in (constants.DataTypes.BYTE, constants.DataTypes.INT32):
        unpack_function = p.unpack_int
    elif data_type in (constants.DataTypes.BOOL, ):
        unpack_function = p.unpack_bool
    elif data_type in (constants.DataTypes.INT64,):
        unpack_function = p.unpack_hyper
    else:
        raise ValueError("Unknown data type '%s'" % data_type)

    if constants.DataTypes.is_complex(data_type):
        # XDR doesn't explicitly support complex numbers, so they're written
        # as pairs of floats (or doubles).
        element_count *= 2

    try:
        data = p.unpack_farray(element_count, unpack_function)
    except (xdrlib.Error, xdrlib.ConversionError) as instance:
        raise UnreadableDataError(instance.msg)

    # Calling p.done() here will raise an xdrlib.Error if unextracted
    # data remains (i.e. the code above is buggy or element_count is wrong)
    try:
        p.done()
    except xdrlib.Error:
        raise UnreadableDataError("More data in file than expected (XDR overrun)")

    if constants.DataTypes.is_complex(data_type):
        # Knit the (real, imaginary) pairs back together into complex numbers.
        data = collapse_complexes(data)

    return data
    
    
def encode_xdr(data, data_type):
    """Given a tuple or list of numbers (ints, floats, complex), returns
    an string representing the data in XDR format. You can also pass 
    single-dimension numpy arrays. Multi-dimension arrays raise an error.

    data_type must be one of the values in 
    vespa.common.constants.DataTypes.ALL.
    
    To reverse the encoding, use decode_xdr().
    """
    if constants.DataTypes.is_complex(data_type):
        # XDR doesn't understand complex numbers so I expand these into 
        # (real, imaginary) pairs.
        data = expand_complexes(data)

    p = xdrlib.Packer()
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

    p.pack_farray(len(data), data, function)
    
    return p.get_buffer()


#####################       Unit tests       #####################

def _test_collapse_expand_complexes():
    import random
    import numpy

    random.seed()
    LIST_SIZE = random.randint(0, 1000)

    collapsed = [ ]
    raw = [ ]
    # Generate a bunch of random floats
    random_float = lambda: random.randint(-1000, 1000) + random.random()
    for i in range(LIST_SIZE):
        real = random_float()
        imaginary = random_float()

        raw.append(real)
        raw.append(imaginary)
        collapsed.append(complex(real, imaginary))

    assert(collapse_complexes(raw) == collapsed)
    assert(expand_complexes(collapsed) == raw)

    # Ensure the functions work with numpy arrays
    raw = numpy.array(raw)
    collapsed = numpy.array(collapsed)
    assert((numpy.array(collapse_complexes(raw)) == collapsed).all())
    assert((numpy.array(expand_complexes(collapsed)) == raw).all())


def _test_dump_iterable():
    import tempfile
    import os
    import numpy

    fd, filename = tempfile.mkstemp()
    os.close(fd)

    dump_iterable(filename, tuple(range(20)))       # tuple of ints
    dump_iterable(filename, list(range(20)))        # list of ints
    dump_iterable(filename, numpy.zeros(10))        # numpy 1D array of floats
    dump_iterable(filename, numpy.zeros( (5,5) ))   # numpy 2D array of floats

    os.remove(filename)


def _test_encode_decode_xdr():
    import random
    import numpy
    import vespa.common.util.math_ as util_math

    random.seed()
    LIST_SIZE = random.randint(0, 1000)

    # Generate a bunch of random floats
    random_float = lambda: random.randint(-1000, 1000) + random.random()
    original = [random_float() for i in range(LIST_SIZE)]

    s = encode_xdr(original, constants.DataTypes.FLOAT32)
    
    round_tripped = decode_xdr(s, constants.DataTypes.FLOAT32, LIST_SIZE)
        
    assert(all([util_math.eq(a, b) for (a, b) in zip(original, round_tripped)]))

    # Ensure the functions work with numpy arrays
    original = numpy.array(original)
    s = encode_xdr(original, constants.DataTypes.FLOAT32)
    
    round_tripped = decode_xdr(s, constants.DataTypes.FLOAT32, LIST_SIZE)
    
    round_tripped = numpy.array(round_tripped)
    
    assert(numpy.allclose(original, round_tripped))


def _test_get_data_type():
    # test ability to translate basic Python data types
    assert(get_data_type(True) == constants.DataTypes.BOOL)
    assert(get_data_type(42) == constants.DataTypes.INT32)
    assert(get_data_type(int(42)) == constants.DataTypes.INT64)
    assert(get_data_type(2.2) == constants.DataTypes.FLOAT64)
    assert(get_data_type(complex(1.1, 3.3)) == constants.DataTypes.COMPLEX128)

    # test ability to detect basic Python types inside iterables
    assert(get_data_type( (42, ) ) == constants.DataTypes.INT32)
    assert(get_data_type( [int(42), ] ) == constants.DataTypes.INT64)
    assert(get_data_type( (2.2, ) ) == constants.DataTypes.FLOAT64)
    assert(get_data_type( [complex(1.1, 3.3), ] ) == constants.DataTypes.COMPLEX128)

    # Test use of default
    assert(get_data_type([ ], float) == constants.DataTypes.FLOAT64)

    import numpy

    # Test numpy basic types
    assert(get_data_type(numpy.int8(42)) == constants.DataTypes.INT32)
    assert(get_data_type(numpy.int16(42)) == constants.DataTypes.INT32)
    assert(get_data_type(numpy.int32(42)) == constants.DataTypes.INT32)
    assert(get_data_type(numpy.int64(42)) == constants.DataTypes.INT64)

    assert(get_data_type(numpy.float32(2.2)) == constants.DataTypes.FLOAT32)
    assert(get_data_type(numpy.float64(2.2)) == constants.DataTypes.FLOAT64)

    assert(get_data_type(numpy.complex64(1+1j)) == constants.DataTypes.COMPLEX64)
    assert(get_data_type(numpy.complex128(1+1j)) == constants.DataTypes.COMPLEX128)

    # Test lists of numpy basic types
    assert(get_data_type( [numpy.int8(42)] ) == constants.DataTypes.INT32)
    assert(get_data_type( [numpy.int16(42)] ) == constants.DataTypes.INT32)
    assert(get_data_type( [numpy.int32(42)] ) == constants.DataTypes.INT32)
    assert(get_data_type( [numpy.int64(42)] ) == constants.DataTypes.INT64)

    assert(get_data_type( [numpy.float32(2.2)] ) == constants.DataTypes.FLOAT32)
    assert(get_data_type( [numpy.float64(2.2)] ) == constants.DataTypes.FLOAT64)

    assert(get_data_type( [numpy.complex64(1+1j)] ) == constants.DataTypes.COMPLEX64)
    assert(get_data_type( [numpy.complex128(1+1j)] ) == constants.DataTypes.COMPLEX128)

    # Test numpy arrays
    assert(get_data_type(numpy.zeros(5, numpy.int8)) == constants.DataTypes.INT32)
    assert(get_data_type(numpy.zeros(5, numpy.int16)) == constants.DataTypes.INT32)
    assert(get_data_type(numpy.zeros(5, numpy.int32)) == constants.DataTypes.INT32)
    assert(get_data_type(numpy.zeros(5, numpy.int64)) == constants.DataTypes.INT64)

    assert(get_data_type(numpy.zeros(5, numpy.float32)) == constants.DataTypes.FLOAT32)
    assert(get_data_type(numpy.zeros(5, numpy.float64)) == constants.DataTypes.FLOAT64)

    assert(get_data_type(numpy.zeros(5, numpy.complex64)) == constants.DataTypes.COMPLEX64)
    assert(get_data_type(numpy.zeros(5, numpy.complex128)) == constants.DataTypes.COMPLEX128)



if __name__ == '__main__':
    _test_get_data_type()
    _test_collapse_expand_complexes()
    _test_encode_decode_xdr()
    _test_dump_iterable()


