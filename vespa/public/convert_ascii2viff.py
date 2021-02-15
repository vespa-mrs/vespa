

import os
import sys
import argparse

import numpy as np

import vespa.common.util.export as util_export
import vespa.common.util.time_ as util_time
import vespa.common.wx_gravy.common_dialogs as common_dialogs
import vespa.common.mrs_data_raw as mrs_data_raw

DESC =  \
"""This utility converts single voxel MRS data from text file format 
into the Vespa VIFF file format readable by Vespa-Analysis.

The text file should have two columns of numbers separated
by whitespace. Column one will contain the real data values
and column two will contain the imaginary data values. One
data point will be on each line. The number of lines in the
text file (that are not comments) that have data points will 
be assumed to be the number of points in SVS MRS data.

Use the command line flags to set SW, FREQ, RESPPM, MIDPPM,
and NoWater settings for the VIFF file.

"""



def convert_ascii_to_viff(fin, sw=6002.4, freq=123.7, resppm=4.7, scale=1.0, fout='', haswater=False, normalize=False, verbose=False):

        
    if not fout:
        fout = fin
        
    fmetab = os.path.splitext(fout)[0] + "_metab.xml"
    fwater = os.path.splitext(fout)[0] + "_water.xml"
        
        
    # Read in entire file 
    with open(fin) as f:
        all_lines = f.read().splitlines()
    
    data = []
    hdr  = []
    for line in all_lines:
        # get rid of whitespace and check for 'comment' character in front
        if line.lstrip().startswith('#'):
            hdr.append(line)
        else:
            data.append(line)

    dim0 = len(data)

    metab = np.ndarray((dim0,),complex)
    water = np.ndarray((dim0,),complex)
    for i,line in enumerate(data):
        vals = [float(item) for item in line.split()]
        metab[i] = vals[0] + 1j*vals[1]
        if haswater:
            water[i] = vals[2] + 1j*vals[3]

    water[0] *= 2.0   # need this because Analysis multiplies all FID data first
    metab[0] *= 2.0   # data points by 0.5 and they made this data up funny.

    if normalize:
        if np.round(np.abs(water[0])) != 0.0:
            water = water / np.abs(water[0])
        metab = metab / np.abs(metab[0])

    water = water * scale
    metab = metab * scale

    metab.shape = 1,1,1,dim0
    water.shape = 1,1,1,dim0

    # Save water and metabolite data to files
        
    stamp = util_time.now(util_time.ISO_TIMESTAMP_FORMAT).split('T')

    lines     = ['Convert_ASCII_to_VIFF ']
    lines.append('------------------------------------------------')
    lines.append('The following information is a summary of the enclosed MRS data.')
    lines.append(' ')
    lines.append('Creation_date    - '+stamp[0])
    lines.append('Creation_time    - '+stamp[1])
    lines.append(' ')
    lines.append('ASCII File       - '+fin)
    lines.append('VIFF File base   - '+fout)
    lines.append(' ')
    lines.append('Frequency [MHz]    '+str(freq))
    lines.append('Sweep width [Hz]   '+str(sw))
    lines.append('Number of points   '+str(dim0))
    lines.append('Resonance PPM      '+str(resppm))
    lines.append(' ')
    lines.append('------------------------------------------------')
    lines.append('ASCII Comment Lines')
    lines.append(' ')
    lines = lines + hdr
    lines = "\n".join(lines)
    if (sys.platform == "win32"):
        lines = lines.replace("\n", "\r\n")

    msg = ''

    if haswater:
        wat = mrs_data_raw.DataRaw() 
        wat.data_sources    = [fwater]
        wat.headers         = [lines]
        wat.sw              = sw
        wat.frequency       = freq
        wat.resppm          = resppm
        wat.data            = water
        filename            = fwater
        try:
            util_export.export(filename, [wat], None, lines, False)
        except IOError:
            msg = """I can't write the file "%s".""" % filename
    
        if msg:
            common_dialogs.message(msg, style=common_dialogs.E_OK)
            return


    met = mrs_data_raw.DataRaw() 
    met.data_sources    = [fmetab]
    met.headers         = [lines]
    met.sw              = sw
    met.frequency       = freq
    met.resppm          = resppm
    met.data            = metab
    filename            = fmetab
    try:
        util_export.export(filename, [met], None, lines, False)
    except IOError:
        msg = """I can't write the file "%s".""" % filename

    if msg:
        common_dialogs.message(msg, style=common_dialogs.E_OK)
        return


#------------------------------------------------------------------------------
# Test routines

def bjs_float(x):
    x = float(x)
    return x

def create_parser():

    parser = argparse.ArgumentParser(prog='convert_ascii2viff', 
                                     usage='%(prog)s [options]',
                                     description=DESC)

    parser.add_argument('infile', nargs='?', default=None,
                                   help='name of ASCII file to be processed')
    
    parser.add_argument('-s', '--sw',     type=float, dest='sw', default=6002.4, help='float, sweep width in Hz ')

    parser.add_argument('-f', '--freq',   type=float, dest='freq', default=123.7, help='float, center frequency in MHz ')

    parser.add_argument('-r', '--resppm', type=float, dest='resppm', default=4.7,help='float, on-resonance PPM value ')

    parser.add_argument('-c', '--scale', type=float, dest='scale', default=1.0,help='float, scaling factor applied to FID data')

    parser.add_argument('-w', '--haswater',  dest='haswater', action="store_true", 
                                help='flag, whether columns 3 and 4 contain water unsuppressed FID data on each line ')

    parser.add_argument('-n', '--normalize',  dest='normalize', action="store_true", 
                                help='flag, normalize by first point of FID ')

    parser.add_argument('-d', '--indir',  dest='indir', action="store_true", 
                                help='flag, whether infile should be treated as a directory to be parsed ')

    parser.add_argument('-v', '--verbose', dest='verbose', action="store_true", 
                                help='increase output verbosity')
    return parser    


def main():

    import glob

    parser = create_parser()
    args = parser.parse_args()

    msg = ''
    if args.infile is None:
        msg = "The 'infile' file name argument is always required. "

    if msg:
        raise parser.error(msg)

    if args.indir:
        os.chdir(args.infile)
        filelist = glob.glob('./*.txt')
    else:
        filelist = [args.infile,]

    for item in filelist:
        convert_ascii_to_viff(item, sw=args.sw, 
                                    freq=args.freq, 
                                    resppm=args.resppm, 
                                    scale=args.scale, 
                                    haswater=args.haswater, 
                                    normalize=args.normalize,
                                    verbose=args.verbose)


if __name__ == '__main__':
    main()

                    
                    