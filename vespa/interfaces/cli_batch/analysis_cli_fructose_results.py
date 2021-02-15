# Python modules
import os
import sys

# 3rd party modules


# Our modules
import vespa.analysis.util_import as util_import


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


def do_results(dataset, datafile, csvfile, verbose=False, debug=False):
    """
    Some typical save type formats = 'svg' 'eps' 'pdf' 'png' 'raw' 'rgba' 'ps' 'pgf' etc.
    Typical fontnames  'Consolas' 'Calibri' 'Courier New' 'Times New Roman'
    
    minplot and maxplot are in PPM
    
    """

    fpath, fname = os.path.split(datafile)

    # Save results to CSV file --------------------------------------
  
    if verbose: print("""  - saving results to CSV file "%s". """ % csvfile)
      
    voxels = dataset.all_voxels
    fit = dataset.blocks["fit"]
    
    # FIXME - this fails if no Prep block, need a property at the dataset level?
    measure_times = dataset.blocks["prep"].measure_time
    
      
    lines = fit.results_as_csv_all_voxels_areas(voxels, fname, measure_times)

    lines = [",".join(line) for line in lines]
    lines = "\n".join(lines)
    lines += '\n'
       
    with open(csvfile, 'a') as f:
        pass
        f.write(lines)



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

    # Processing of SVS_EDIT_DIFF files
    STARTDIR = 'D:\\Users\\bsoher\\temp\\current\\data'
    csvfile  = STARTDIR+'\\_csv_output_file_metabs.txt'

    skip = ['raw_fruct-002','raw_fruct-003','raw_fruct-006','raw_fruct-011','raw_fruct-012',
            'raw_fruct-016','raw_fruct-017','raw_fruct-018','raw_fruct-020','raw_fruct-021','raw_fruct-023',
            'raw_fruct-025','raw_fruct-026','raw_fruct-028','raw_fruct-035','raw_fruct-046',
            'raw_fruct-047','raw_fruct-048','raw_fruct-055','raw_fruct-057','raw_fruct-059',
            'raw_fruct-061','_archive_results_pass1','_archive_results_pass2']

    # get all paths in data directory, remove unusable datasets, convert into filenames
    paths = list(filter(os.path.isdir, [os.path.join(STARTDIR,f) for f in os.listdir(STARTDIR) if f not in skip]))
    
    datafiles = [path+'_all_files.xml' for path in paths]

    for i,datafile in enumerate(datafiles):

        # Load Main Dataset --------------------------
        
        fpath, fname = os.path.split(datafile)
        
        if verbose: print("""%s - Load Data into a Dataset object - %s""" % (str(i), fname))    
        datasets, _ = _open_viff(datafile)
        dataset = datasets[-1]

        do_results( dataset, datafile, csvfile, verbose=verbose, debug=False)
        
        # if i >= 3: break    # debug statement to exit after one file processed
        
    bob = 10
    bob += 1
    
    
  


        
if __name__ == '__main__':
    main()        

