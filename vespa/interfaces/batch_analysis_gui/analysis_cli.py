# Python modules

import argparse
import os
import sys
import platform
import importlib

# 3rd party modules


# Our modules
import vespa.analysis.mrs_dataset as mrs_dataset
import vespa.analysis.block_prep_wbnaa as block_prep_wbnaa
import vespa.analysis.util_import as util_import

import vespa.common.configobj as configobj
import vespa.common.util.misc as util_misc
import vespa.common.util.export as util_export


SUPPORTED = ['wbnaa',]

DESC =  \
"""Command line interface to process MRS data in Vespa-Analysis. 
 Data filename, preset file name, data type string and CSV output 
 file name values are all required for this command to function 
 properly.
  
 Note. You may have to enclose data/preset/output strings in double 
 quotation marks for them to process properly if they have  
 spaces or other special characters embedded in them.
"""

def analysis_cli(datafile, datatype, presetfile, csvfile, savexml=True,
                                                          twodir=False,
                                                          verbose=False, 
                                                          debug=False):
    
    # Test input arguments for consistency --------------------------
    
    msg = ''
    if not os.path.isfile(datafile):
        msg = """DATAFILE does not exist "%s".""" % datafile 

    if not os.path.isfile(presetfile):
        msg = """PRESETFILE does not exist "%s".""" % presetfile 
        
    if not str(datatype).lower() in SUPPORTED:
        msg = """DATATYPE not supported - "%s".""" % datatype 

    csvpath, _ = os.path.split(csvfile)
    if not os.path.exists(csvpath):
        msg = """CSVFILE path does not exist "%s".""" % csvpath 
        
    if msg:        
        print(msg, file=sys.stderr)
        sys.exit(-1)
        
    if not os.path.isfile(csvfile):
        if verbose:
            print("""Output CSV file will be created - %s""" % csvfile)
    
    
    # Load DATA data into a Dataset object --------------------------
    
    if verbose: print("Load Data data into a Dataset object")    
    try:
        # some datatypes are user derived and describe in this ini file
        filename = "analysis_import_menu_additions.ini"
        filename = os.path.join(util_misc.get_data_dir(), filename)
        menu_config = configobj.ConfigObj(filename, encoding="utf-8")
        
        if datatype == 'wbnaa':

            if 'siemens_twix_wbnaa' not in list(menu_config.keys()):
                msg = """WBNAA VB reader not available to Analysis - please locate."""
                print(msg, file=sys.stderr)
                sys.exit(-1)
            
            module_name = "siemens_twix_wbnaa"
            section = menu_config[module_name]
            spath   = section["path"]
            if os.path.exists(spath):
                # module = imp.load_source(module_name, spath)
                spec = importlib.util.spec_from_file_location(module_name, spath)
                module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(module)
                sys.modules[module_name] = module
            
            klass = getattr(module, section["class_name"])
            reader = klass()
            reader.filenames = [datafile,]
            
            raws = reader.read_raws()   # no open dataset

            d = {}
            d["prep"] = block_prep_wbnaa.BlockPrepWbnaa
            dataset = mrs_dataset.dataset_from_raw(raws[0], d, 0)
            
        #else:
            # should not get here with consistency check above

    except:
        msg = """Unknown exception reading Data file "%s".""" % filename 
        print(msg, file=sys.stderr)
        sys.exit(-1)
    
    
    # Load PRESET data ----------------------------------------------
    
    if verbose: print("Load Preset into the Dataset object")    
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
            
               






    
    
def create_parser():

    parser = argparse.ArgumentParser(prog='analysis_cli', 
                                     usage='%(prog)s [options] datafile datatype presetfile csvfile',
                                     description=DESC)
    
    parser.add_argument("datafile",   help='file name of data to be processed')
    parser.add_argument("datatype",   help='string name of type of data in file')
    parser.add_argument("presetfile", help='xml file of preset values to use in processing')
    parser.add_argument("csvfile",    help='name of file where CSV result values are saved')
    
    parser.add_argument('-x', '--savexml', dest='savexml',
                                action="store_true", 
                                help='save dataset and results to VIFF XML format')
    parser.add_argument('-t', '--twodir', dest='twodir',
                                action="store_true", 
                                help='use two directory names in unique data id')
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
      
    analysis_cli(args.datafile, 
                 args.datatype, 
                 args.presetfile, 
                 args.csvfile, 
                 savexml=args.savexml,
                 twodir=args.twodir,
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
        