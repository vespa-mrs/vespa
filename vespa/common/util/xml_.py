# Python modules

import xml.etree.cElementTree as ElementTree
import zlib
import base64
import datetime
import sys
import io

# 3rd party modules
import numpy as np

# Our modules
import vespa.common.constants as constants
import vespa.common.util.fileio as util_fileio
from functools import reduce


# ENCODING_ATTR is a convenience constant. It's the encoding as dict suitable
# for passing to ElementTree as the description of an element's attributes.
ENCODING_ATTR = { "encoding" : constants.NUMERIC_LIST_ENCODING }

# BOOLEANS maps Python booleans to XML booleans and vice versa. 
# Per the W3C --
#     An instance of a datatype that is defined as boolean can 
#     have the following legal literals {true, false, 1, 0}.
# ref: http://www.w3.org/TR/xmlschema-2/#boolean
# In the spirit of "Be generous in what you accept", we also accept 
# "True" and "False". 
# Note that when creating XML (mapping Python ==> XML), we only map to 
# "true" and "false".
BOOLEANS = { True   : "true", False   : "false",
             "true" : True  , "false" : False  ,
             "1"    : True  , "0"     : False  ,
             "True" : True  , "False" : False  ,
           }  


# _SIMPLE_TYPE_NAMES maps simple types to their string counterparts, and 
# vice versa. They're for the functions dict_to_elements() and 
# element_to_dict().
_SIMPLE_TYPE_NAMES = { bool  : "bool",   int : "int", int : "long", 
                       float : "float",  complex : "complex",
                       str   : "string", str : "string", 
                     }
_SIMPLE_TYPE_NAMES.update(dict(list(zip(list(_SIMPLE_TYPE_NAMES.values()), 
                                   list(_SIMPLE_TYPE_NAMES.keys())))))
# I make one small alteration to force string elements to be extracted
# as unicode.
_SIMPLE_TYPE_NAMES["string"] = str

# bjs - bugfix june 8, 2020
_SIMPLE_TYPE_NAMES["int"] = int

# This file offers a number of utility functions for translating Python
# objects into XML (actually, into ElementTree objects) and vice versa.
# Many of the functions are reciprocal pairs, e.g. element_to_numpy_array() 
# and numpy_array_to_element().

################################################################ 
# Some general purpose internal functions are here
################################################################ 
def _is_homogenous(an_iterable):
    # Given an iterable (list or tuple), returns True if all elements are
    # of the same type or the list is empty, False otherwise. 
    if an_iterable:
        # Get the type of the first element
        type_ = type(an_iterable[0])

        # Create a list of bools that state whether or not the type of the
        # remaining elements matches that of the first.
        matches = [(type(item) == type_) for item in an_iterable[1:]]

        return all(matches)
    else:
        return True


def _is_numeric(x):
    # Given an object, returns True if its type is numeric, False otherwise.
    # Python >= 2.6 has the numbers module which would make this slightly
    # simpler. This code is compatible with Python 2.5. 
    return isinstance(x, (int, complex, float))



################################################################ 
# Some public functions
################################################################ 
def TextElement(tag, content, attributes={}, **extra):
    """A factory function identical to ElementTree.Element() except that
    this function accepts a parameter ("content") which is assigned to the 
    element's .text attribute.
    
    In addition, this function offers a small convenience to the caller. The
    caller can pass values of types int, long, float, boolean and
    datetime.datetime in the content param, and they'll automatically be
    converted to string (via repr(), or .isoformat() in the case of datetime
    objects).
    
    Note that only scalars are offered this convenience. Lists, arrays, etc.
    must still be converted to text by the caller.
    """
    e = ElementTree.Element(tag, attributes, **extra)
    
    # Convert to text if necessary. Note that bool is a type of int, so the
    # test to see if content is of boolean type must come before the test 
    # for int.
    if isinstance(content, bool):
        content = BOOLEANS[content]
    elif isinstance(content, (int, float, complex)):
        content = repr(content)
    elif isinstance(content, datetime.datetime):
        content = content.isoformat()
    #else:
        # It had better be a string already or there's gonna be trouble...

    if content:
        e.text = content
        
    return e
    

def TextSubElement(parent, tag, content, attributes={}, **extra):
    """A factory function identical to ElementTree.SubElement() except that
    this function accepts a parameter which is assigned to the element's 
    .text attribute.
    
    Also see the docstring for TextElement().
    """
    e = TextElement(tag, content, attributes, **extra)
    parent.append(e)
        
    return e


def find_settings(element, alternate_name):
    """Given an ElementTree.Element and a string, searches the element for 
    a subelement called "settings". If it doesn't find it, searches for
    a subelement using the alternate_name.

    Prior to Vespa 0.7.0, the settings in Analysis blocks had unique class 
    and XML tag names, e.g. BlockSpectralSettings/block_spectral_settings.
    As of 0.7.0, they're all just called "settings". This function finds
    whichever of the two is present.
    """
    # First search under the preferred name.
    settings_element = element.find("settings")

    if settings_element is None:
        settings_element = element.find(alternate_name)

    return settings_element
    

def indent(element, level=0):
    """Given an ElementTree.Element, indents it for pretty printing"""
    # Based on Fredrik Lundh's code from effbot.org 
    # Copyright (c) 1999-2006 by Secret Labs AB
    # Copyright (c) 1999-2006 by Fredrik Lundh
    # 
    # By obtaining, using, and/or copying this software and/or its
    # associated documentation, you agree that you have read, understood,
    # and will comply with the following terms and conditions:
    # 
    # Permission to use, copy, modify, and distribute this software and its
    # associated documentation for any purpose and without fee is hereby
    # granted, provided that the above copyright notice appears in all
    # copies, and that both that copyright notice and this permission notice
    # appear in supporting documentation, and that the name of Secret Labs
    # AB or the author not be used in advertising or publicity pertaining to
    # distribution of the software without specific, written prior
    # permission.
    # 
    # SECRET LABS AB AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO
    # THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
    # FITNESS.  IN NO EVENT SHALL SECRET LABS AB OR THE AUTHOR BE LIABLE FOR
    # ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
    # WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
    # ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT
    # OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
    INDENT = "\t"
    i = "\n" + (level * INDENT)
    if len(element):
        if not element.text or not element.text.strip():
            element.text = i + INDENT
        for e in element:
            indent(e, level + 1)
        if not e.tail or not e.tail.strip():
            e.tail = i
    if level and (not element.tail or not element.tail.strip()):
        element.tail = i


#############   Function pairs start here (to/from ElementTree.Elements)

# _dict_to_file() and _file_to_dict() are public convenience functions. 
# They're good enough for Vespa developers to use while debugging, but too 
# crude for use in Vespa production code. 
def _dict_to_file(the_dict, filename, root_element_name="root"):
    # Given a dict and a filename, converts the dict to XML and writes it.
    # Uses dict_to_elements()  (q.v.)
    root = ElementTree.Element(root_element_name)
    dict_to_elements(the_dict, root)
    indent(root)
    ElementTree.ElementTree(root).write(filename, "utf-8")

def _file_to_dict(filename):
    # Given the name of a file created by _dict_to_file(), returns the 
    # appropriate dict. Uses element_to_dict() (q.v.)
    s = open(filename, "rb").read()
    return element_to_dict(ElementTree.fromstring(s))


def _any_type_to_element(tag, value):
    # A helper function for dict_to_elements().
    #
    # Given a tag name and a value, creates an ElementTree.Element 
    # containing that value. The value must be one of the types valid for
    # dict_to_elements().
    # Note that this function is recursive.
    
    # We perform the trick of converting arbitrary objects to XML by writing
    # type info as an attribute associated with the element. So if this
    # function is invoked with params of ("answer", 42), it will create a
    # tag like this:
    #    <answer type="int">42</answer>
    
    if value is None:
        # I create an empty element
        element = ElementTree.Element(tag, {"type" : "none"})
    elif isinstance(value, (bool, int, float, str, complex)):
        # It's a simple type, I convert it to text
        attributes = { "type" : _SIMPLE_TYPE_NAMES[type(value)] }
        element = TextElement(tag, value, attributes)
    elif isinstance(value, np.ndarray):
        # It's a (relatively) simple type
        element = numpy_array_to_element(value, tag)
    elif isinstance(value, (list, tuple)):
        # It's a container type. 
        if value and _is_numeric(value[0]) and _is_homogenous(value):
            # The list is all numbers, so we can take advantage of
            # numeric_list_to_element() and get a compact representation.
            element = numeric_list_to_element(value, type(value[0]), tag)
        else:
            # Add an element as the root of this fragment and recurse.
            element = ElementTree.Element(tag, {"type" : "list"})
        
            for item in value:
                element.append(_any_type_to_element("value", item))
    elif isinstance(value, dict):
        # It's a container type. Add an element as the root of this 
        # fragment and recurse.
        element = ElementTree.Element(tag, {"type" : "dict"})
        # This is a dict inside a dict; inner_xxx refer to the inner dict.
        for inner_key, inner_value in value.items():
            element.append(_any_type_to_element(inner_key, inner_value))
    else:
        raise TypeError("Can't convert the type '%s'" % str(type(value)))
        
    return element


def _element_to_any_type(element):
    # A helper function for element_to_dict().
    #
    # Given an element created by _any_type_to_element(), returns the 
    # Python object that created the element.
    # Note that this function is recursive.
    
    # All elements have a type attribute except numpy arrays and numeric 
    # lists. Numpy arrays are easy to spot because they have shape and
    # data_type attributes. Numeric lists don't have a shape attribute but
    # they do have data_type. 
    type_ = element.get("type")
    
    if element.get("shape"):
        # It's a numpy array
        value = element_to_numpy_array(element)
    elif element.get("data_type"):
        # It's a numeric list
        value = element_to_numeric_list(element)
    elif type_ == "none":
        value = None
    elif type_ == "bool":
        value = BOOLEANS[element.text]
    elif type_ in _SIMPLE_TYPE_NAMES:
        # It's a simple type (int, string, etc.) so the conversion is
        # straightforward. 
        value = element.text
        if type_ == "complex":
            # In the case of complex numbers, there's a hitch. In some (most?)
            # versions of Python 2.5, complex() doesn't understand repr()'s
            # version of a complex number because it doesn't like the 
            # surrounding parens. I don't think an actual bug was filed for
            # this, but it was discussed on the Python mailing list and
            # a patch was proposed. See this conversation:
            # http://code.activestate.com/lists/python-list/515842/
            # Or Google for: "Complex evaluation bug" python 2006
            major, minor = sys.version_info[0:2]
            if (major == 2) and (minor <= 5):
                value = value.strip("()")
            #else:
                # Be conservative, don't alter it.

        value = _SIMPLE_TYPE_NAMES[type_](value)
    elif type_ == "list":
        # It's a container type. Create the container and populate it with
        # the element's children via recursion.
        value = [ ]
        for child in element.getchildren():
            value.append(_element_to_any_type(child))
    elif type_ == "dict":
        # It's a container type. Create the container and populate it with
        # the element's children via recursion.
        value = { }
        for child in element.getchildren():
            value[child.tag] = _element_to_any_type(child)
    else:
        raise TypeError("Can't convert the type '%s'" % type_)

    return value


def dict_to_elements(d, parent):
    """Given an arbitrary dict (with some restrictions), converts it to a
    set of ElementTree.Elements which become children of the parent param.
    There is one child for each key in the dictionary, hence the plural
    "elements" in the function name.
    
    The restrictions on the input dictionary are as follows:
    1. All keys must be strings that are valid XML entity names. If you 
       restrict your strings to upper and lower case ASCII alphanumerics and 
       underscore, you'll be OK. Spaces are NOT permitted in key names.

       If you want the precise details on what's permissible, see here:
       http://www.w3.org/TR/2008/REC-xml-20081126/#sec-starttags
    
    2. All values must be one of the following types:
       - bool, int, long, float, complex, string (unicode or str)
       - a numpy array
       - list, tuple, or dict
       You can nest containers (dicts, lists and tuples) as deeply as you 
       like, as long as the leaf nodes are one of the non-container types
       listed above.
       
    3. Unicode strings and pure ASCII strings are trouble-free. All other
       strings must be UTF-8 encoded.
       
    When a dict is round-tripped through this function and element_to_dict(),
    the resulting dict can differ from the original in two ways --
    - Tuples become lists.
    - All string values become unicode.
    """
    for key, value in d.items():
        parent.append(_any_type_to_element(key, value))


def element_to_dict(root):
    """Given the parent element passed to dict_to_elements(), returns 
    the corresponding dict.
    """
    d = { }
    for element in root.getchildren():
        d[element.tag] = _element_to_any_type(element)

    return d


def list_to_elements(tag, value):
    """Given an arbitrary list (with some restrictions), converts it to a
    set of ElementTree.Elements which become children of the parent param.
    There is one child for each entry in the list, hence the plural
    "elements" in the function name.
    
    The restrictions on the input list are as follows:
    1. All values must be one of the following types:
       - bool, int, long, float, complex, string (unicode or str)
       - a numpy array
       - list, tuple, or dict
       You can nest containers (dicts, lists and tuples) as deeply as you 
       like, as long as the leaf nodes are one of the non-container types
       listed above.
       
    2. Unicode strings and pure ASCII strings are trouble-free. All other
       strings must be UTF-8 encoded.
       
    When a list is round-tripped through this function and element_to_list(),
    the resulting list can differ from the original in two ways --
    - Tuples become lists.
    - All string values become unicode.
    """
    if not isinstance(value, list):
        return None
    
    return _any_type_to_element(tag, value)

def element_to_list(root):
    """Given the parent element passed to list_to_elements(), returns 
    the corresponding list.
    """
    if root.get("type") != "list":
        return None
    
    d = []
    for element in root.getchildren():
        d.append(_element_to_any_type(element))

    return d


def numeric_list_to_element(an_iterable, data_type, tag_name):
    """Given a list (or other iterable, e.g. a tuple) of numbers (int, float,
    complex), returns an ElementTree.Element containing a string 
    representation of the list containing enough info to reconstitute the
    original list.
    
    If your "list" is a numpy array, the function numpy_array_to_element() is
    probably a better choice than this one since it preserves shape info.
    """
    # We need a string representation of the datatype; the numpy strings
    # are convenient so we use them.
    data_type = constants.DataTypes.any_type_to_numpy(data_type)
    
    attrs = {   "encoding"  : constants.NUMERIC_LIST_ENCODING,
                "data_type" : data_type
            }
    e = ElementTree.Element(tag_name, attrs)
    
    data_type = constants.DataTypes.any_type_to_internal(data_type)
    e.text = encode_numeric_list(an_iterable, data_type)
    
    return e
    
    
def element_to_numeric_list(e):
    """Given an ElementTree.Element written by numeric_list_to_element(), 
    extracts the data therein and returns a list containing the data.
    """
    data_type = e.get("data_type")
    
    data_type = constants.DataTypes.any_type_to_internal(data_type)
    
    return decode_numeric_list(e.text, e.get("encoding"), data_type)
    

def array3d_to_element(array3d, tag_name):
    """Given a 3D array of numpy arrays (a numpy array where each element
    is also a numpy array), returns an ElementTree.Element containing the shape 
    subelements representing each numpy array.
    """
    shape = ','.join(map(str, array3d.shape))

    e = ElementTree.Element(tag_name, {"shape" : shape})

    for item in array3d.flat:
        e.append(numpy_array_to_element(item, "array"))
    
    return e
    
    
def element_to_array3d(e):
    """Given an ElementTree.Element written by array3d_to_element(), 
    extracts the data therein and returns the 3D array.
    """
    # An earlier version of this function contained a bug that wrote a
    # four element shape, i.e. it recorded 4 dimensions for a 3D list. The
    # dimensions were from Analysis data; the first dimension was the number
    # of data points (e.g. 4096) and was irrelevant. That's the one we need to
    # ignore if we see it.
    shape = e.get("shape")
    shape = shape.split(',')
    if len(shape) == 4:
        # This is old style; ignore the first dimension.
        shape = shape[1:]

    shape = [int(dim) for dim in shape]

    # Create an empty, flat numpy array, read each XML element, convert it to a
    # numpy array, and store it in the flat numpy array. 
    # This should be possible to do with one line of code which is faster
    # and more readable, IMO:
    #    array3d = np.array(map(element_to_numpy_array, e.getiterator("array")))
    # However when numpy sees the arrays inside the list that map() generates,
    # it combines them into one single 2D array. I don't want a 2D array, I 
    # want a 1D array that happens to have numpy arrays in each element.
    size = reduce(lambda x, y: x * y, shape)

    array3d = np.empty(size, object)
    i = 0
    for array_element in e.getiterator("array"):
        array3d[i] = element_to_numpy_array(array_element)
        i += 1

    array3d.shape = shape

    return array3d
    

def element_to_numpy_array(e):
    """
    Given an ElementTree.Element written by numpy_array_to_element(), 
    extracts the data therein and returns a properly shaped numpy array.
    
    """
    encoding = e.get("encoding")
    
    data_type = constants.DataTypes.any_type_to_internal(e.get("data_type"))
    
    data = decode_numeric_list(e.text, encoding, data_type)

    if 'xdr' in encoding:
        data_type = constants.DataTypes.any_type_to_numpy(data_type)
        ndarray = np.fromiter(data, data_type)
        shape = e.get("shape")
        shape = [int(dim) for dim in shape.split(',')]
        data = ndarray.reshape(shape)

    return data


def numpy_array_to_element(array, tag_name, encoding=''):
    """
    Given a numpy array, returns an ElementTree.Element containing a
    string representation of the array containing enough info to reconstitute
    the array and its shape.

    This function only supports numeric and boolean arrays. It can't handle
    arrays of objects.
    
    """
    encode = constants.NUMERIC_LIST_ENCODING if encoding=='' else encoding

    attrs = {   "encoding"  : encode,
                "shape"     : ','.join([str(dim) for dim in array.shape]),
                "data_type" : str(array.dtype)
            }
    e = ElementTree.Element(tag_name, attrs)
    
    data_type = constants.DataTypes.any_type_to_internal(str(array.dtype))
    
    # Recording the shape info and encode the array 
    e.text = encode_numeric_list(array, data_type)
    
    return e

   
def decode_numeric_list(data, encoding, data_type):
    """
    2019-12-2 We have added the Numpy 'load()/save()' cross-platform format
        as an optional 'encoding' option. So this method may return a list 
        of numbers, or an actual Numpy array, depending on the encoding 
        list of strings. (e.g. npy zlib base64)
    
    Given a string (data), returns it decoded into a Python list. The 
    encoding parameter describes how the data was encoded from a list to
    a string (e.g. "xdr zlib base64").
    
    data_type must be one of the values in 
    vespa.common.constants.DataTypes.ALL.
    
    Typically the data string comes from XML and was created by 
    encode_numeric_list().
    
    """
    encoding = encoding.strip().split()
    encoding.reverse()
    
    for transform in encoding:
        if transform == "xdr":
            element_count = len(data) // constants.DataTypes.XDR_TYPE_SIZES[data_type]

            # Note that if the data type is complex, decode_xdr will return
            # a list of complex numbers for me. I don't have to call
            # util_fileio.collapse_complexes() here.
            data = util_fileio.decode_xdr(data, data_type, element_count)
        elif transform == "npy":
            buf = io.BytesIO(data)              # was OK taking a 'str' in Py27
            data = np.load(buf)
        elif transform == "zlib":
            data = zlib.decompress(data)        # 'str' in returns 'str' in Py27
        elif transform == "base64":
            data = base64.b64decode(data)       # 'str' in from Py27 created data -> 'bytes' out
        else:
            raise ValueError("Unrecognized transformation '%s'" % transform)

    return data


def encode_numeric_list(data, data_type):
    """
    2019-12-2 We have added the Numpy 'load()/save()' cross-platform format
        as an optional 'encoding' option. Had to modify input and 'xdr' option 
        to accomodate the 'npy' format. Now we can pass in a numpy array, and 
        in the 'xdr' step it will have ravel() and tolist() applied.
              
        So, generally we now assume that a numpy array is being passed in
        the 'data' parameter. And that will be either converted into a
        Numpy save byte-string, or a ravel().tolist() list of values. Thus,
        the name of the function is only loosely requiring an actual list
        of numbers.         

    Original documentation:

    Given a list of numbers of the same type (int, float, or complex), 
    encodes it into a string suitable for writing to XML.
    
    data_type must be one of the values in 
    vespa.common.constants.DataTypes.ALL.
    
    To turn the data back into a list, use decode_numeric_list().
    
    """
    for transform in constants.NUMERIC_LIST_ENCODING.strip().split():
        if transform == "xdr":
            data = util_fileio.encode_xdr(data.ravel().tolist(), data_type)
        elif transform == "npy":
            buf = io.BytesIO()
            np.save(buf, data)      # data still numpy array here
            data = buf.getvalue()   # returns a 'byte' representation of buffer
        elif transform == "zlib":
            data = zlib.compress(data, 9)   # in py27 'str' input returns type(data)='str'
        elif transform == "base64":
            data = base64.b64encode(data)   # in py27 'str' input returns type(data)='str'
        else:
            raise ValueError("Unrecognized data format '%s'" % transform)

    data = data.decode('utf-8')     # new with Py3 since some steps above now return 'byte' instead of 'str'

    return data

    

if __name__ == "__main__":
    # test _is_homogenous()
    assert(_is_homogenous( [ ] ))
    assert(_is_homogenous( [1, 2, 3] ))
    assert(_is_homogenous( [1.0, 2.0, 3.0] ))
    assert(_is_homogenous( ['a', 'b', 'c'] ))
    assert(_is_homogenous( [None, None, None] ))
    assert(_is_homogenous( [ElementTree.Element('x'), ElementTree.Element('y')] ))
    assert(not _is_homogenous( [1, 2, 3.14159] ))
    assert(not _is_homogenous( [1, 'a'] ))
    assert(not _is_homogenous( [1.0, complex(1.0, 0.0)] ))
    
    # Test dict_to_elements() & element_to_dict()
    
    # This dict contains at least one instance of each type accepted by
    # dict_to_elements() except for a tuple. That's because when round
    # tripped, a tuple becomes a list, so if I include a tuple in the 
    # input the final assert() will always fail.
    original = {  "blah" : True,                    # bool
                  "foo" : 42,                       # int
                  "very_big" : -999999999999999,    # long
                  "zzz" : 3.14,                     # float
                  "z" : complex(5.6, -9),           # complex
                  "greeting" : "hello world",       # str
                  "u_greeting" : "hello world",    # unicode
                  "stuff1" : [6, 7, 8],             # list (numeric, homogenous)
                  "stuff2" : [6, 7, 8.0],           # list (numeric, mixed)
                  "stuff3" : [6, "dog", 8.0],       # list (mixed)
                  "stuff_in_stuff" : [ list(range(4)), list(range(5)) ], # nested lists
                  "monkey" : {                      # dict
                               "aaa" : 42,          # int in dict
                               "bbb" : [ 6,7 ],     # list in dict
                               "ccc" : { "xyz" : [ ] }, # dict in dict
                               "nowt" : None        # None in dict
                             },
                  "some_data" : np.random.randint(1000, size=(5, 5)),
               }

    root = ElementTree.Element("root")
    
    dict_to_elements(original, root)
    
    tree = ElementTree.ElementTree(root)
    # Here you can print the XML if you want to see what it looks like.
    # indent(root)
    # import sys
    # tree.write(sys.stdout, "utf-8")
    
    round_tripped = element_to_dict(root)

    # OK, we're ready to compare the orginal & round tripped dicts.
    # The numpy arrays need special treatment when being compared.
    assert(all(np.equal(original["some_data"], round_tripped["some_data"]).tolist()))
    
    # Get rid of those fussy numpy arrays.
    del original["some_data"]
    del round_tripped["some_data"]

    assert(original == round_tripped)
