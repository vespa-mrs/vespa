# Python modules
import argparse
import os
import sys
import xml.etree.cElementTree as ElementTree
import collections

# 3rd party modules
import numpy as np

# Our modules
import vespa.analysis.util_import as util_import
import vespa.analysis.fileio.raw_reader as raw_reader 
import vespa.common.util.export as util_export
import vespa.common.twix_parser_multi_raid as twix_multi_raid
import vespa.common.util.misc as util_misc

from vespa.common.mrs_data_raw import DataRawFidsum


# Change to True to enable the assert() statements sprinkled through the code
ASSERTIONS_ENABLED = False

# data is complex64 per Siemens documentation for Twix
NUMPY_DATA_TYPE = np.complex64
# BYTES_PER_ELEMENT expresses how many bytes each element occupies. You
# shouldn't need to change this definition.
BYTES_PER_ELEMENT = np.zeros(1, dtype=NUMPY_DATA_TYPE).nbytes

# RAWDATA_SCALE is based on ICE_RAWDATA_SCALE in SpecRoFtFunctor.cpp
RAWDATA_SCALE = 131072.0 * 256.0


DESC =  \
"""Command line interface to process MRS data in Vespa-Analysis. 
 Data filename, preset file name, data type string and CSV output 
 file name values are all required for this command to function 
 properly.
  
 Note. You may have to enclose data/preset/output strings in double 
 quotation marks for them to process properly if they have  
 spaces or other special characters embedded in them.
"""

def clean_header(header):
    """ converts all values in ICE dict into a long string"""
    return "need to write"


def analysis_cli_gulin(data, preset, verbose=False, debug=False):
    
    
    # Load PRESET data ----------------------------------------------
    
    if verbose: print("Load Preset into the Dataset object")    
    try:
        msg = ""
        
        # update dataset object with preset blocks and chains
        dataset.apply_preset(preset, voxel=(0,0,0))
    except:
        msg = """Unknown exception reading Preset file "%s".""" % filename 
        print(msg, file=sys.stderr)
        sys.exit(-1)


    # Process and fit data ------------------------------------------
    
    if verbose: print("Running dataset chain objects")

    #FIXME-bjs  Need a routine to 'get all voxels' here
    
    voxel = [(0,0,0),]
    for key in list(dataset.blocks.keys()):
        block = dataset.blocks[key]
        block.chain.run(voxel, entry='all')

    # Create unique name ID for this dataset ------------------------
    #
    # Some data files use the same base name over and over again (eg. twix 
    # uses meas.dat for all data files), so we mix either one or two
    # parent directory names to get a unique name to report in the csv file
    # and to use in the output XML file.
    
    base = os.path.basename(datafile)
    base,  _     = os.path.splitext(base)
    path1, _     = os.path.split(datafile)  
    path3, path2 = os.path.split(path1)
    if not twodir:
        # unique output name using one parent directory
        outxml = os.path.join(path1, path2+'_'+base+'.xml')
    else: 
        # unique output name using two parent directories
        _ , path4 = os.path.split(path3)
        outxml = os.path.join(path1, path4+'_'+path2+'_'+base+'.xml')
    
    dataset.dataset_filename = outxml

    # Save results to CSV file --------------------------------------

    if verbose: print("""Saving results to CSV file "%s". """ % csvfile)
    
    fit = dataset.blocks["fit"]
    data_source = dataset.blocks["raw"].get_data_source(voxel)
    
    val, hdr = fit.results_as_csv(voxel[0], fit.chain.fitted_lw,
                                            fit.chain.minmaxlw[0],
                                            fit.chain.minmaxlw[1], 
                                            data_source, outxml)
    nhdr = len(hdr)
    val = ",".join(val)
    hdr = ",".join(hdr)
    val += "\n"
    hdr += "\n"
     
    hdr_flag = True
    if os.path.isfile(csvfile):
        with open(csvfile, 'r+') as f:
            data = f.readlines()
            if len(data)>1:
                last = data[-1]
                nlast = len(last.split(','))
                if nlast == nhdr:
                    hdr_flag = False
                
    with open(csvfile, 'a') as f:
        if hdr_flag:
            f.write(hdr)
        f.write(val)

    # Save results to XML if flag set -------------------------------
    
    if savexml:
        if verbose: print("""Saving dataset to XML file "%s". """ % outxml)
        
        try:
            util_export.export(outxml, [dataset], None, None, False)
        except:
            msg = """I can't write the file "%s".""" % outxml
            print(msg, file=sys.stderr)
            sys.exit(-1)
            
               








#------------------------------------------------------------------------------

class RawReaderSiemensTwixSlaserCmrrVe(raw_reader.RawReader):
    # This inherits from raw_reader.RawReader (q.v.). The only methods you
    # need to implement are __init__() and read_raw(). You *may* want to 
    # override or supplement some of raw_reader.RawReader's other methods.

    def __init__(self):
        raw_reader.RawReader.__init__(self)
        
        # The twix files all have an extension of *.dat, thus the filter below.
        self.filetype_filter = "Twix File (*.dat)|*.dat;*.DAT;"
        self.multiple = False
        
        
    def read_raw(self, filename, ignore_data=False, open_dataset=None):
        """
        Given the name of a twix file, returns a DataRawFidsum object
        populated with the parameters and data represented by the file. 
        The 'Fidsum' flavor of raw data object has individual FIDs returned in
        it rather than one summed FID.

        The ignore_data option is not implemented for this reader. 
        
        The open_dataset attribute is not used in this reader. 
        
        """
        twix = twix_multi_raid.TwixMultiRaid()
        
        twix.populate_from_file(filename)
        
        raws = []
        
        meas = twix.measurements[-1]
        
        scans, evps = meas.scans, meas.evps

        # Fill a dictionary with the minimum necessary parameters needed from
        # the twix file in order to create the DataRawFidsum objects.
        d = _extract_parameters(evps)

        if d["remove_os"]: d["sw"]   = d["sw"] / 2.0

        #----------------------------------------------------------------------
        # dataset1 - scan 0, water unsuppressed for coil combine
        
        d["filename"] = filename+'.combine'

        # Parse scan data into array
        tmp = scans[0]
        dat = _read_data(tmp, d)

        d["data"] = dat
        d["data_source"] = filename+'.combine'
        raw = DataRawFidsum(d)

        raws.append(raw)

        #----------------------------------------------------------------------
        # dataset1 - scan 1-2, water unsuppressed for ECC
        
        # Parse scan data into array
        tmp = scans[1:3]
        dat = _read_data(tmp, d)

        d["data"] = dat
        d["data_source"] = filename+'.ecc1'
        raw = DataRawFidsum(d)

        raws.append(raw)

        #----------------------------------------------------------------------
        # dataset3 - scans 3-4, water unsuppressed for water scale
        
        # Parse scan data into array
        tmp = scans[3:5]
        dat = _read_data(tmp, d)

        d["data"] = dat
        d["data_source"] = filename+'.water1'
        raw = DataRawFidsum(d)

        raws.append(raw)

        #----------------------------------------------------------------------
        # dataset4 - scans 5-68 (64 total), metabolite data with WS

        # Parse scan data into array
        tmp = scans[5:69]
        dat = _read_data(tmp, d)

        d["data"] = dat
        d["data_source"] = filename+'.metab64'
        raw = DataRawFidsum(d)

        raws.append(raw)

        #----------------------------------------------------------------------
        # dataset5 - scans 69-70 (2 total), water unsuppressed for ecc

        # Parse scan data into array
        tmp = scans[69:71]
        dat = _read_data(tmp, d)

        d["data"] = dat
        d["data_source"] = filename+'.ecc2'
        raw = DataRawFidsum(d)

        raws.append(raw)

        #----------------------------------------------------------------------
        # dataset6 - scans 71-72 (2 total), water unsuppressed for water scale

        # Parse scan data into array
        tmp = scans[71:73]
        dat = _read_data(tmp, d)

        d["data"] = dat
        d["data_source"] = filename+'.water2'
        raw = DataRawFidsum(d)

        raws.append(raw)
        # TODO - bjs
        # - check scaling to get these numbers up
        # - check if ECC associate works in concert with other associations 

        return raws




####################    Internal functions start here     ###############

def _extract_parameters(evps):
    """
    Given the contents of an SPAR file as a string, extracts a few specific
    parameters and returns a flat dict containing those parameters and their 
    value. 
    The returned dict is appropriate for passing to DataRaw.inflate().
    
    """
    d = { }

    header, clean_header = _parse_protocol_data(evps[3][1])

    # A copy of this goes into the dict.
    d["header"] = clean_header

    remove_oversample_flag = header.get("sSpecPara.ucRemoveOversampling", "0x0")
    remove_oversample_flag = (remove_oversample_flag.strip() == "0x1")

    d["sw"]             = 1.0 / (float(header.get("sRXSPEC.alDwellTime[0]", 1.0)) * 1e-9)
    d["remove_os"]      = remove_oversample_flag
    d["readout_os"]     = float(_get_siemens_xprotocol(evps[0][1], "ReadoutOS", 1.0))
    d["sequence_type"]  = header.get("tSequenceFileName", "wbnaa")
    d["frequency"]      = float(header["sTXSPEC.asNucleusInfo[0].lFrequency"])/1000000.0
    d["dims"]           = [1,1,1,2048]
    d["dims"][0]        = int(header["sSpecPara.lVectorSize"]) 
    d["dims"][1]        = 1 # concat will take care of header["lAverages"]
    d["dims"][2]        = 1
    d["dims"][3]        = 1 
    d["seqte"]          = float(header["alTE[0]"])*1e-6
    
    d["nucleus"] = header["sTXSPEC.asNucleusInfo[0].tNucleus"].replace('"',' ').strip()

    return d  
    
    
    
def _read_data(scans, d):
    """
    This Import method was originally developed to read/process from Siemens
    IceSpectroEdit, so I leave in the documentation below describing where the
    processing code came from.
    
    Given a list of scans and the extracted parameter dictionary we process
    the data in as similar a way to Siemens ICE program IceSpectroEdit as
    we can. My way of removing oversampling seems to not match 100% but the
    other steps are in the order and perform very similarly to the ICE program.

    When we are done parsing the ON and OFF state FIDs into arrays, we create
    SUM and DIF arrays to also send back.
    
    """
    dims    = d["dims"]
    dim0    = dims[0]
    acqdim0 = dims[0] * int(d["readout_os"])
    remove_os = d["remove_os"]
        
    ncoils2  = scans[0].scan_header.used_channels
    coil_ids = sorted(set([chan[0].channel_id for scan in scans for chan in scan.channels ]))
    nscans   = len(scans)
    ncoils   = len(coil_ids)

    data = np.ndarray([ncoils,nscans,dim0],dtype=complex)

    # determine if there are any extra points at the beginning or end of FID
    start_point = scans[0].scan_header.free_parameters[0]
    end_point   = start_point + acqdim0

    # lead to typical integral values of 1 or larger which are nicely displayed
    scale = RAWDATA_SCALE / float(nscans) 

    # Parse each group of FIDs for N channels for each average as separate
    # from all other scans. Perform the following steps:
    # - scale, all channels
    # - if needed, we remove oversampling, all channels
    # - convert to complex conjugate for proper plotting
    dat = []

    for iscan, scan in enumerate(scans):
        
        # for each average, remove pre-/post- points and apply phase to  
        # correct for coil geometry
        chans = []
                
        for ichan, item in enumerate(scan.channels):
            chan = item[1]
            chan = np.array(chan[start_point:end_point])*scale

            if remove_os:
                #
                # Note that this is a simple OS removal, not the full Siemens
                # algorithm. I have not had time to implement and test that.
                #
                # There are five steps to remove oversampling:
                # 1. forward transform data at acqdim0 resolution
                # 2. shift data by dim0/2 so we have a FID centered at dim0/2 before snip
                # 3. snip out first dim0 points and roll back by -dim0/2
                # 4. inverse transform data at dim0 resolution
                # 5. (omitted) multiply by 0.5 to maintain same scale as original data
                #    (based on the results in the actual IceSpecEdit program,
                #     but still these results are a few percent off of Siemens 
                #     results, ah well) 
                chan = np.fft.fft(chan)
                chan = np.roll(chan, int(dim0/2))
                chan = np.fft.ifft(np.roll(chan[:dim0], int(-dim0/2))) 

            # apply complex conjugate to swap x-axis for proper display of data
            chan = np.conjugate(chan)

            data[ichan,iscan,:] = chan 

    return data

    
    
def _parse_protocol_data(protocol_data):
    """
    Returns a dictionary containing the name/value pairs inside the
    "ASCCONV" section of the MrProtocol or MrPhoenixProtocol elements
    of a Siemens CSA Header tag.
    """
    # Protocol_data is a large string (e.g. 32k) that lists a lot of
    # variables in a JSONish format with which I'm not familiar. Following
    # that there's another chunk of data delimited by the strings you see
    # below.
    # That chunk is a list of name=value pairs, INI file style. We
    # ignore everything outside of the ASCCONV delimiters. Everything inside
    # we parse and return as a dictionary. 
    #
    # As of the Siemens VD scanner software version the starting string is
    # no longer ### ASCCONV BEGIN ### rather it seems to have some other
    # info about what was converted inserted after the BEGIN and before 
    # the ### delimiter. To get around this for now, we search just for the  clean_header[-1] len(protocol_data)
    # beginning of the string ### ASCONV BEGIN, and then throw away the
    # first line after we split the string into lines.
    #
    start = protocol_data.find("### ASCCONV BEGIN")
    end = protocol_data.find("### ASCCONV END ###")

    _my_assert(start != -1)
    _my_assert(end != -1)

    clean_start = start
    clean_end   = end + len("### ASCCONV END ###")
    clean_header = protocol_data[clean_start:clean_end]

    start += len("### ASCCONV BEGIN ###")
    protocol_data = protocol_data[start:end]

    lines = protocol_data.split('\n')
    lines = lines[1:]

    # The two lines of code below turn the 'lines' list into a list of
    # (name, value) tuples in which name & value have been stripped and
    # all blank lines have been discarded.
    f = lambda pair: (pair[0].strip(), pair[1].strip())
    lines = [f(line.split('=')) for line in lines if line]

    return dict(lines), clean_header


def _my_assert(expression):
    if ASSERTIONS_ENABLED:
        assert(expression)
        
 
def _get_siemens_xprotocol(head_only, key, default):
    
    head = util_misc.normalize_newlines(head_only)   
    head = head.split("\n")
    
    # find substring 'key' in list even if item in list is not iterable
    items = [el for el in head if isinstance(el, collections.Iterable) and (key in el)]
    
    for item in items:
        start = item.find("{")
        end = item.find("}")
        if start != -1 and end != -1:
            temp = item[start+1:end]
            # remove duplicate white space
            temp = " ".join(temp.split())
            temp = temp.split()
            
            if len(temp) == 1:
                return temp[0]
            elif len(temp) == 3:
                return temp[2]
            else:
                return temp
    
    return default



def _remove_oversampling(scan, d):   
    """ not currently in use as incomplete """
    
    vector_size = d["vector_size"]
    remove_os   = d["remove_os"]
    left_points = scan.pre
    right_points = scan.post
    reduced_points = scan.samples_in_scan - scan.pre - scan.post
    half_vector_size = int(vector_size / 2)
    
    if (reduced_points % vector_size) != 0:
        raise ValueError('remove_oversampling: final data size not multiple of vector size.')


    if not remove_os:
        # keep oversampled points but remove extra points
        start_point = scan.pre
        scan.data = scan.data[start_point:start_point+reduced_points]
    else:
        # remove oversampled data
        shift_points = scan.post if scan.post < scan.pre else scan.pre
        
        if shift_points == 0:
            # no extra pts available, final data will show truncation artifact
            start_point = scan.pre
            data = np.array(scan.data[start_point:start_point+reduced_points])
            data = np.fft.fft(data, n=vector_size)
            data = np.fft.ifft(data) * 0.5
            scan.data = data.tolist()
            
        else:
            # Extra pts available to use for removing truncation artifact.
            # Process data twice, centering signal to the left and right of kSpaceCentreColumn (TE)
            # Retrieve half of final signal from each set.
            pass


#------------------------------------------------------------------------------

def process_ice(scans, sw, freq, dim0, seqte, ncoils, nscans, npts, preset_string):
    """
    Assume that data is a list of FID data arrays in (scan, chan, npts) order.
    I need to process each FID to remove oversampling and then store into 
    a numpy array.
    
    
    """
    d = []
    
    filename = "Siemens_VB19_ICE_Transfer_semi-LASER"
    
    d["filename"]       = filename
    d["data_source"]    = filename
    d["sw"]             = 1.0 / float(sw * 1e-9)
    d["readout_os"]     = 2
    d["remove_os"]      = True
    d["sequence_type"]  = 'slaser'
    d["frequency"]      = float(freq)/1000000.0
    d["dims"]           = [1,1,1,2048]
    d["dims"][0]        = int(dim0) * int(d["readout_os"])
    d["dims"][1]        = 1 # this updates as raw.concatenate is applied
    d["dims"][2]        = 1
    d["dims"][3]        = 1 
    d["seqte"]          = float(seqte)*1e-6
    d["nucleus"]        = '1H'
    d["preset_string"]  = preset_string

    # load DATA into DATASETS ---------------------------------------
    
    data = np.ndarray([ncoils,nscans,dim0],dtype=complex)

    # determine if there are any extra points at the beginning or end of FID
    start_point = 0
    end_point   = start_point + (dim0 * 2)

    # lead to typical integral values of 1 or larger which are nicely displayed
    scale = RAWDATA_SCALE / float(nscans) 
    
    for iscan in range(nscans):
        
        # for each average, remove pre-/post- points and apply phase to  
        # correct for coil geometry
        chans = []
                
        for ichan in range(nchan):
            
            indx = ichan + iscan * nchan
            
            chan = scans[indx]
            chan = np.array(chan[start_point:end_point])*scale

            if remove_os:
                chan = np.fft.fft(chan)
                chan = np.roll(chan, int(dim0/2))
                chan = np.fft.ifft(np.roll(chan[:dim0], int(-dim0/2))) 

            # apply complex conjugate to swap x-axis for proper display of data
            chan = np.conjugate(chan)

            data[ichan,iscan,:] = chan 

        if d["remove_os"]: d["sw"]   = d["sw"] / 2.0

        datasets = []

        #----------------------------------------------------------------------
        # dataset1 - scans 1-4, water unsuppressed
        
        d["filename"]    = filename+'.water1'
        d["data"]        = data[:,1:5,:]
        raw = DataRawFidsum(d)

        datasets.append(raw)

        #----------------------------------------------------------------------
        # dataset2 - scans 5-68 (64 total), metabolite data with WS

        d["filename"]    = filename+'.metab64'
        d["data"]        = data[:,5:69,:]
        raw = DataRawFidsum(d)

        datasets.append(raw)

        #----------------------------------------------------------------------
        # dataset3 - scans 69-72 (4 total), water unsuppressed

        d["filename"]    = filename+'.water2'
        d["data"]        = data[:,69:73,:]
        raw = DataRawFidsum(d)

        datasets.append(raw)
        

    # Load PRESET data ----------------------------------------------
    
    if verbose: print("Load Preset into the Dataset object")    
    try:
        msg = ""
        
        preset_string  = header["preset_string"]
        preset_element = ElementTree.fromstring(preset_string)
        
        try:
            importer = util_import.DatasetImporter(preset_string)
        except IOError:
            msg = """I can't read the preset string "%s".""" % preset_string
        except SyntaxError:
            msg = """The preset string "%s" isn't valid Vespa Interchange File Format.""" % preset_string

        if msg:
            pass
            # TODO need to exit gracefully from this
        else:
            # Time to rock and roll!
            presets = importer.go()
            preset  = presets[0]
    except:
        msg = """Unknown exception reading Preset string "%s".""" % preset_string 
        # TODO need to exit gracefully from this


    img, outxml = analysis_cli_gulin(datasets, 
                                     preset, 
                                     verbose=args.verbose, 
                                     debug=args.debug)

    
    


def process_twix(datafile, presetfile, verbose=False, debug=False):
    
    # Load DATASET object ----------------------------------------------
    
    if verbose: print("Load Datafile into the Dataset object")    
    try:
        data_parser = RawReaderSiemensTwixSlaserCmrrVe()
        datasets = data_parser.read_raw(datafile)
    except:
        msg = """Unknown exception parsing/reading Dataset file "%s".""" % datafile 
        print(msg, file=sys.stderr)
        sys.exit(-1)

    # Load PRESET object ----------------------------------------------
    
    if verbose: print("Load Presetfile into the Preset object")    
    try:
        msg = ""
        try:
            importer = util_import.DatasetImporter(presetfile)
        except IOError:
            msg = """I can't read the preset file "%s".""" % presetfile
        except SyntaxError:
            msg = """The preset file "%s" isn't valid Vespa Interchange File Format.""" % presetfile

        if msg:
            print(msg, file=sys.stderr)
            sys.exit(-1)
        else:
            # Time to rock and roll!
            presets = importer.go()
            preset  = presets[0]
    except:
        msg = """Unknown exception reading Preset file "%s".""" % presetfile 
        print(msg, file=sys.stderr)
        sys.exit(-1)

    
    return datasets, preset

    
    
def create_parser():

    parser = argparse.ArgumentParser(prog='analysis_cli', 
                                     usage='%(prog)s [options] datafile presetfile',
                                     description=DESC)
    
    parser.add_argument("datafile",   help='file name of data to be processed')
    parser.add_argument("presetfile", help='xml file of preset values to use in processing')
    
    parser.add_argument('-x', '--savexml', dest='savexml',
                                action="store_true", 
                                help='save dataset and results to VIFF XML format')
    parser.add_argument('-d', '--debug',   dest='debug',
                                action="store_true", 
                                help='print out command line to console only')
    parser.add_argument('-v', '--verbose', dest='verbose', 
                                action="store_true", 
                                help='increase output verbosity')
    return parser    



def main():

    parser = create_parser()
    args = parser.parse_args()
      
    datasets, preset = process_twix( args.datafile, 
                                     args.presetfile, 
                                     verbose=args.verbose, 
                                     debug=args.debug)

    img, outxml = analysis_cli_gulin(datasets, preset, 
                                     verbose=args.verbose, 
                                     debug=args.debug)

    
#     # generalize the call to the local testdata directory, wherever it is
#     fpath = os.path.dirname(os.path.realpath(__file__))
#  
#     datafile   = fpath+r'/testdata/1/meas.dat'
#     presetfile = fpath+r'/testdata/preset_analysis_wbnaa_from_paper.xml'
#     csvfile    = fpath+r'/testdata/outcsv1.txt'
#    
#     analysis_cli(datafile, 'wbnaa', presetfile, csvfile, twodir=True, verbose=True)
        

#     # this gets all files *.dat in all subdirectories of STARTDIR
#     STARTDIR = 'D:\\bsoher\\Dropbox\\23_VB_Datasets4BrianSoher'
#     datfiles = []
#     for dirpath, dirnames, filenames in os.walk(STARTDIR):
#         for filename in [f for f in filenames if f.endswith(".dat")]:
#             datfiles.append(os.path.join(dirpath, filename))
#             print os.path.join(dirpath, filename)
#      
#     # this demonstrates how to use one or two levels of parent directories
#     # to create a uniques filename for a twix meas.dat        
#     twodir = True
#     for datafile in datfiles:
#         base = os.path.basename(datafile)
#         base,  _     = os.path.splitext(base)
#         path1, _     = os.path.split(datafile)  # unique output name for XML save file if needed
#         path3, path2 = os.path.split(path1)
#         if not twodir:
#             outxml = os.path.join(path1, path2+'_'+base+'.xml')
#         else: 
#             _ , path4 = os.path.split(path3)
#             outxml = os.path.join(path1, path4+'_'+path2+'_'+base+'.xml')   
#         print outxml 

#     # this gets all files *.dat in all subdirectories of STARTDIR
#     STARTDIR = 'D:\\bsoher\\3xWBNAA_Young_Controls'
#     datafiles = []
#     for dirpath, dirnames, filenames in os.walk(STARTDIR):
#         for filename in [f for f in filenames if f.endswith(".dat")]:
#             datafiles.append(os.path.join(dirpath, filename))
#             print os.path.join(dirpath, filename)
# 
#     presetfile = STARTDIR+'\\preset_wbnaa_analysis_v4.xml'
#     csvfile    = STARTDIR+'\\outcsv1.txt'
# 
#     for datafile in datafiles:
#         analysis_cli(datafile, 'wbnaa', presetfile, csvfile, twodir=True, verbose=True)

        
if __name__ == '__main__':
    main()        
        