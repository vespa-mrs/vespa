# Python modules
import os
import sys

# 3rd party modules


# Our modules
import vespa.analysis.mrs_dataset as mrs_dataset
import vespa.analysis.block_prep_timeseries as block_prep_timeseries
import vespa.analysis.util_import as util_import

import vespa.analysis.fileio.dicom_siemens_timeseries as fileio_dicom_siemens_timeseries
import vespa.analysis.fileio.util_exceptions as util_exceptions
import vespa.common.util.export as util_export
import vespa.common.mrs_data_raw_timeseries as mrs_data_raw_timeseries


DESC =  \
"""Command line interface to process MRS data in Vespa-Analysis. 
 Data filename, preset file name, data type string and CSV output 
 file name values are all required for this command to function 
 properly.
  
 Note. You may have to enclose data/preset/output strings in double 
 quotation marks for them to process properly if they have  
 spaces or other special characters embedded in them.
"""





def do_analysis(dataset, preset, verbose=False, debug=False):
    

    # Update dataset with preset ------------------------------------
    dataset.apply_preset(preset, voxel=(0,0,0))

    # Process and fit data ------------------------------------------
    if verbose: print("  Running dataset chain objects")
    _process_all_blocks(dataset)

    # Save Dataset with results to VIFF XML if flag set -------------------------------

    dirpath, _     = os.path.split(dataset.blocks["raw"].data_source)
    dirpath, _     = os.path.split(dirpath)
    dirpath, fbase = os.path.split(dirpath)
    filename = fbase + "_all_files.xml"
    filename = os.path.join(dirpath, filename)
    
    if verbose: print("""  Saving dataset to XML file "%s". """ % filename)

    dataset.dataset_filename = filename
    
    try:
        bob = 10
        util_export.export(filename, [dataset], None, None, False)
    except:
        msg = """I can't write the file "%s".""" % outxml
        print(msg, file=sys.stderr)
        print(msg, file=sys.stdout)
        sys.exit(-1)


def _process_all_blocks(dataset):
    """ for all voxels, run chain in all blocks to update """
    
    tmp = dataset.blocks['prep'].chain.run([0,0,0])
    dataset.batch_fit_all()
    


def _import_preset(presetfile):

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
            print(msg, file=sys.stdout)
            sys.exit(-1)
        else:
            # Time to rock and roll!
            presets = importer.go()
            preset  = presets[0]
            
            return preset

    except:
        msg = """Unknown exception reading Preset file "%s".""" % filename 
        print(msg, file=sys.stderr)
        print(msg, file=sys.stdout)
        sys.exit(-1)    
            
               
def _import_siemens_dicom_timeseries(filenames, open_dataset=None):
    """
    Stolen from Analysis main.py module - trimmed for CLI usage
    
    Assumption here is that we are opening one file in the reader. If there is
    an 'open_dataset' sent in, then the current reader is for an associated 
    file. We will associate the current file with the open one at the end of
    the code.
    
    """
    try:

        reader = fileio_dicom_siemens_timeseries.RawReaderDicomSiemensTimeseries()
        reader.filenames = filenames
        
        datasets = [ ]
    
        msg = ""
        try:
            # Step 1
            #
            # Return one or more DataRawXxxx objects that indicate what
            # sort of data was read in by the reader
            raws = reader.read_raws(open_dataset=open_dataset)
        except IOError:
            msg = "One or more of the files couldn't be read due to a disk error."
    
        except util_exceptions.MultifileAttributeMismatchError:
            msg = _MSG_MULTIFILE_ATTRIBUTE_MISMATCH
    
        except util_exceptions.MultifileTypeMismatchError:
            msg = _MSG_MULTIFILE_TYPE_MISMATCH
    
        except util_exceptions.UnsupportedDimensionalityError:
            # Note that this also catches SIDataError which is a
            # subclass of UnsupportedDimensionalityError
            msg = _MSG_UNSUPPORTED_DIMENSIONALITY
    
        except util_exceptions.OpenFileAttributeMismatchError:
            msg = _MSG_OPEN_ATTRIBUTE_MISMATCH
    
        except util_exceptions.OpenFileTypeMismatchError:
            msg = _MSG_OPEN_TYPE_MISMATCH
    
        except util_exceptions.FileNotFoundError as error_instance:
            msg = str(error_instance)
    
        except util_exceptions.OpenFileUserReadRawError as error_instance:
            if not error_instance:
                error_instance = "User read_raw raised OpenFileUserReadRawError"
            msg = str(error_instance)
    
    
        if msg:
            print(msg, file=sys.stdout)
            print(msg, file=sys.stderr)
            sys.exit(-1)
        else:
            # All is well. Convert these raw objects into fully-fledged
            # dataset objects.
            if open_dataset:
                zero_fill_multiplier = open_dataset.zero_fill_multiplier
            else:
                zero_fill_multiplier = 0
    
            # Step 2
            #
            # See if any data types need special classes. We usually only
            # look for raw fidsum classes which trigger a prep fidsum block.
            block_class_specs = [ ]
            for raw in raws:
                d = { }
                if isinstance(raw, mrs_data_raw_timeseries.DataRawTimeseries):
                    d["prep"] = block_prep_timeseries.BlockPrepTimeseries
                    
                block_class_specs.append(d)
    
            f = lambda raw, block_classes: mrs_dataset.dataset_from_raw(raw,
                                                              block_classes,
                                                       zero_fill_multiplier)
            datasets = list(map(f, raws, block_class_specs))
    
        if datasets:
    
            if open_dataset is not None:
                open_dataset.blocks['raw'].set_associated_datasets([datasets[0], ])
    
            return datasets[0], open_dataset
        
        else:
            return None, open_dataset

    except:
        msg = """Unknown exception reading Data file "%s".""" % filenames[0] 
        print(msg, file=sys.stderr)
        print(msg, file=sys.stdout)
        sys.exit(-1)

    


def main():

    verbose = True

    # Processing of SVS_EDIT_DIFF files
    STARTDIR    = 'D:\\Users\\bsoher\\temp\\current\\data'
    presetfile  = STARTDIR+'\\_preset_raw_fruct_v4.xml'
    
    paths = list(filter(os.path.isdir, [os.path.join(STARTDIR,f) for f in os.listdir(STARTDIR)]))
    paths = paths[::-1]

    i = 0

    for path in paths:

        path = path + '\\series_sum'
        filenames = list(filter(os.path.isfile, [os.path.join(path,f) for f in os.listdir(path)]))
        
        print('Processing Subject = '+path+' with '+str(len(filenames))+' files')
        
        # Load Main Dataset --------------------------
        
        dataset = _import_siemens_dicom_timeseries(filenames, open_dataset=None)
        dataset = dataset[0]
          
        # Load Preset data for Main Dataset ----------------------------------------------
        
        if verbose: print("  Read Preset object")    
        preset = _import_preset(presetfile)
          
          
        do_analysis(dataset, preset, verbose = verbose, debug = False)
        
        i += 1
#        if i >= 3: break    # debug statement to exit after one file processed
        
    bob = 10
    bob += 1
    
    
    


        
if __name__ == '__main__':
    main()        
        