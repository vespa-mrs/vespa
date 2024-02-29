from __future__ import division
from __future__ import unicode_literals

import sys
import os
import zlib
import base64
import xdrlib

import test1_encoded

# Our modules
import vespa.common.constants as constants
import vespa.common.util.xml_ as util_xml
import vespa.common.util.fileio as util_fileio




# We encode numeric lists (also numpy arrays) in a three step process.
# First is XDR (http://en.wikipedia.org/wiki/External_Data_Representation),
# second is zlib to save space, third is base64 to make the output of
# zlib palatable to XML.
NUMERIC_LIST_ENCODING = "xdr zlib base64"


def read_fake_bitmap(fname='test1.bin'):

    bitmap = []

    try:
        bytes_read = open(fname, "rb").read()
        for item in bytes_read:
            bitmap.append(ord(item))    
    except:
        e = "bitmap error" + sys.exc_info()[0]
        print "======> bjs sample.py      python got here 3 - except = "+str(e)

    return bitmap


def encode_fake_bitmap(bitmap):
    
    data = bitmap
    p = xdrlib.Packer()
    p.pack_farray(len(data), data, p.pack_int)
    data = p.get_buffer()
    data = zlib.compress(data, 9)
    data = base64.b64encode(data)

    return data

 
def write_fake_bitmap_to_file(bitmap, fname='temp.txt'):

    textfile = open(fname, 'w')
    textfile.write(bitmap)
    textfile.close()


def decode_fake_bitmap(bitmap):
    
    data = bitmap
    
    data = base64.b64decode(data)
    data = zlib.decompress(data)
    p = xdrlib.Unpacker(data)
    data = p.unpack_farray(len(data) // 4, p.unpack_int)     # 4 since we expect INT32 here
       
    return data







if __name__ == '__main__':
    
    bitmap0 = read_fake_bitmap()
    
    bitmap1 = encode_fake_bitmap(bitmap0)
    
    bitmap2 = decode_fake_bitmap(bitmap1)
    
    bitmap3 = decode_fake_bitmap(test1_encoded.test1)
    
    if bitmap2 == bitmap0:
        print 'bitmap2 == bitmap0'
    if bitmap3 == bitmap0:
        print 'bitmap3 == bitmap0'
    
    bob = 10
    bob += 1