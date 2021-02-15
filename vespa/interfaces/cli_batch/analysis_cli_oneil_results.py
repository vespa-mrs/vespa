# Python modules
import os
import sys

# 3rd party modules
import matplotlib.pyplot as plt

# Our modules
import vespa.analysis.util_import as util_import
import vespa.analysis.figure_layouts as figure_layouts 
import vespa.common.util.misc as util_misc


SUPPORTED = ['wbnaa', 'siemens dicom']

DESC =  \
"""Command line interface to process MRS data in Vespa-Analysis. 
 Data filename, preset file name, data type string and CSV output 
 file name values are all required for this command to function 
 properly.
  
 Note. You may have to enclose data/preset/output strings in double 
 quotation marks for them to process properly if they have  
 spaces or other special characters embedded in them.
"""


def analysis_cli_oneil_results(dataset, csvfile, outbase='D:\\Users\\bsoher\\myplot',
                                                 viffpath='', 
                                                 vespa_version='',
                                                 timestamp='',
                                                 fixphase=False,
                                                 savetype='pdf',
                                                 fontname='Courier New',
                                                 minplot=0.1,
                                                 maxplot=4.9,
                                                 verbose=False, debug=False):
    """
    Some typical save type formats = 'svg' 'eps' 'pdf' 'png' 'raw' 'rgba' 'ps' 'pgf' etc.
    Typical fontnames  'Consolas' 'Calibri' 'Courier New' 'Times New Roman'
    
    minplot and maxplot are in PPM
    
    """
    fsupported = plt.gcf().canvas.get_supported_filetypes()
    if savetype not in list(fsupported.keys()):
        msg = r"Output file format '%s' not supported by current Matplotlib backend, Returning." % savetype
        raise ValueError(msg)
    
    
    # Print Control Setting ----------------------

    fig = figure_layouts.lcm_like(dataset, 
                                  viffpath=viffpath, 
                                  vespa_version=vespa_version,
                                  timestamp=timestamp,
                                  fontname=fontname,
                                  minplot=minplot,
                                  maxplot=maxplot,
                                  extfig=None,
                                  fixphase=fixphase,
                                  verbose=False, debug=False)
 
    outname = outbase+'.'+savetype
 
    fig.savefig(outname, pad_inches=0.5)

    # Save results to CSV file --------------------------------------
  
    if verbose: print("""Saving results to CSV file "%s". """ % csvfile)
      
    voxel = dataset.all_voxels
    fit = dataset.blocks["fit"]
    data_source = dataset.blocks["raw"].get_data_source(voxel)
      
    val, hdr = fit.results_as_csv(voxel[0], fit.chain.fitted_lw, 0.0, 0.0, data_source, viffpath)
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



def _open_viff(datafile):
    
    datasets = []

    filename = datafile
    timestamp = ''

    msg = ""
    try:
        importer = util_import.DatasetCliImporter(filename)
    except IOError:
        msg = """I can't read the file "%s".""" % filename
    except SyntaxError:
        msg = """The file "%s" isn't valid Vespa Interchange File Format.""" % filename


    if msg:        
        print(msg, file=sys.stderr)
        print(msg, file=sys.stdout)
        sys.exit(-1)
    else:
        # Time to rock and roll!
        dsets, timestamp = importer.go()

        for item in dsets:
            datasets.append(item)


    if datasets:

        for dataset in datasets:
            if dataset.id == datasets[-1].id:
                dataset.dataset_filename = filename
                # dataset.filename is an attribute set only at run-time
                # to maintain the name of the VIFF file that was read in
                # rather than deriving a filename from the raw data
                # filenames with *.xml appended. But we need to set this
                # filename only for the primary dataset, not the associated
                # datasets. Associated datasets will default back to their
                # raw filenames if we go to save them for any reason
            else:
                dataset.dataset_filename = ''

    return datasets, timestamp

    


def main():

    verbose = True

    datatype = 'siemens dicom'

#     # Processing of SVS_EDIT_OFF files
#     STARTDIR = 'D:\\Users\\bsoher\\temp\\dmx'           # \\fitted_pass1'
#     csvfile_metab = STARTDIR+'\\bjs_csv_output_file_off_metab.txt'
#     csvfile_water = STARTDIR+'\\bjs_csv_output_file_off_water.txt'

#     # Processing of SVS_EDIT_OFF files
#     STARTDIR = 'D:\\Users\\bsoher\\projects\\2017_oneil_ucla_lipid_contam_svs\\MRS_RAW_DATA\\_all_svs_edit_off' 
#     csvfile_metab = STARTDIR+'\\bjs_csv_output_file_off_metab.txt'
#     csvfile_water = STARTDIR+'\\bjs_csv_output_file_off_water.txt'

#     # Processing of SVS_EDIT_DIFF files
#     STARTDIR = 'D:\\Users\\bsoher\\projects\\2017_oneil_ucla_lipid_contam_svs\\MRS_RAW_DATA\\_all_svs_edit_diff' 
#     csvfile_metab = STARTDIR+'\\bjs_csv_output_file_diff_metab.txt'
#     csvfile_water = STARTDIR+'\\bjs_csv_output_file_diff_water.txt'
# 
# 
#     # this gets all files *.IMA in all subdirectories of STARTDIR
#     imafiles = []
#     for dirpath, dirnames, filenames in os.walk(STARTDIR):
#         for filename in [f for f in filenames if f.endswith(".xml")]:
#             imafiles.append(os.path.join(dirpath, filename))
#             print os.path.join(dirpath, filename)

    vespa_version = util_misc.get_vespa_version()

    i = 0

    STARTDIR = 'D:\\Users\\bsoher\\projects\\2017_oneil_ucla_lipid_contam_svs\\MRS_RAW_DATA\\_all_svs_edit_diff' 
    imafiles = ['D:\\Users\\bsoher\\projects\\2017_oneil_ucla_lipid_contam_svs\\MRS_RAW_DATA\\_all_svs_edit_diff\\VTP_7017_MRS.MR.LEE_EVP5.0004.0003.2014.11.25.14.38.59.0.26120564.xml',]
    csvfile_metab = STARTDIR+'\\bob_csv_output_file_diff_metab.txt'
    csvfile_water = STARTDIR+'\\bob_csv_output_file_diff_water.txt'

    for datafile in imafiles:

        # Test input arguments for consistency --------------------------
        
        msg = ''
        if not os.path.isfile(datafile):
            msg = """Main DATAFILE does not exist "%s".""" % datafile 
            
        if msg:        
            print(msg, file=sys.stderr)
            print(msg, file=sys.stdout)
            sys.exit(-1)
            
        if not os.path.isfile(csvfile_metab):
            if verbose:
                pass
                print("""Output Metab CSV file will be created - %s""" % csvfile_metab)

        if not os.path.isfile(csvfile_water):
            if verbose:
                pass
                print("""Output Water CSV file will be created - %s""" % csvfile_water)
        
        
        # Load Main Dataset --------------------------
        
#        if verbose: print """%s - Load Data into a Dataset object - %s""" % (str(i), datafile)    
        datasets, timestamp = _open_viff(datafile)
        
               
        dataset = datasets[-1]
        outbase = datafile+'.metab'
        print(str(i)+' : '+' - '+datafile)
        analysis_cli_oneil_results( dataset, 
                                    csvfile_metab, 
                                    outbase=outbase,
                                    viffpath=datafile, 
                                    vespa_version=vespa_version,
                                    timestamp=timestamp,
                                    fixphase=True,
                                    verbose=verbose, 
                                    debug=False)

        dataset = datasets[0]
        outbase = datafile+'.water'
        print(str(i)+' : '+' - '+datafile)
        analysis_cli_oneil_results(dataset, 
                                   csvfile_water, 
                                   outbase=outbase,
                                   viffpath=datafile, 
                                   vespa_version=vespa_version,
                                   timestamp=timestamp,
                                   minplot=-100.0,
                                   maxplot=100.0,
                                   fixphase=False,
                                   verbose=verbose, 
                                   debug=False)
        
        i += 1
#        if i >= 1: break    # debug statement to exit after one file processed
        
    bob = 10
    bob += 1
    
    
  


        
if __name__ == '__main__':
    main()        

