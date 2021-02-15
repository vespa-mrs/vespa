# Python modules
import os
import sys

# 3rd party modules


# Our modules
import vespa.analysis.mrs_dataset as mrs_dataset
import vespa.analysis.block_prep_wbnaa as block_prep_wbnaa
import vespa.analysis.util_import as util_import
import vespa.analysis.fileio.dicom_siemens as fileio_dicom_siemens
import vespa.common.util.export as util_export
import vespa.common.mrs_data_raw_wbnaa as mrs_data_raw_wbnaa

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





def analysis_cli_oneil(dataset, preset, ecc_dataset           = None, 
                                        preset_ecc_dataset    = None, 
                                        watref_dataset        = None, 
                                        preset_watref_dataset = None, 
                                        verbose               = False, 
                                        debug                 = False):
    

    # Update dataset with preset ------------------------------------
    dataset.apply_preset(preset, voxel=(0,0,0))

    # Update ECC dataset with ECC preset ------------------------------------
    ecc_dataset.apply_preset(preset_ecc_dataset, voxel=(0,0,0))

    # Process and fit ECC dataset if needed -------------------------
    if verbose: print("Running ECC dataset chain objects")
    _process_all_blocks(ecc_dataset)


    # Attach ECC datasets to the preset/datafile
    if ecc_dataset is not None:
        block    = ecc_dataset.blocks["raw"]
        raw_data = block.data.copy() * 0
        dims     = ecc_dataset.raw_dims
        for k in range(dims[3]):
            for j in range(dims[2]):
                for i in range(dims[1]):
                    dat = block.data[k,j,i,:].copy() / block.data[k,j,i,0]
                    raw_data[k,j,i,:] = dat

        dataset.blocks['spectral'].set.ecc_dataset    = ecc_dataset
        dataset.blocks['spectral'].set.ecc_dataset_id = ecc_dataset.id
        dataset.blocks['spectral'].set.ecc_raw        = raw_data
        dataset.blocks['spectral'].set.ecc_filename   = block.data_source
 
#     dataset.blocks['quant'].set.watref_dataset    = watref_dataset
#     dataset.blocks['quant'].set.watref_dataset_id = watref_dataset.id

    # Process and fit data ------------------------------------------
    if verbose: print("Running dataset chain objects")
    chain_outputs = _process_all_blocks(dataset)

    # Save Dataset with results to VIFF XML if flag set -------------------------------

    dirpath, filename = os.path.split(dataset.blocks["raw"].data_source)
    filename = os.path.splitext(filename)[0] + ".xml"
    filename = os.path.join(dirpath, filename)
    
    if verbose: print("""Saving dataset to XML file "%s". """ % filename)

    dataset.dataset_filename = filename
    
    try:
        util_export.export(filename, [dataset], None, None, False)
    except:
        msg = """I can't write the file "%s".""" % outxml
        print(msg, file=sys.stderr)
        print(msg, file=sys.stdout)
        sys.exit(-1)


def _process_all_blocks(dataset):
    """ for all voxels, run chain in all blocks to update """
    
    chain_outputs = {}
    
    voxel = dataset.all_voxels
    for key in list(dataset.blocks.keys()):
        if key == 'spectral':
            key = 'spectral'
            block = dataset.blocks[key]
            tmp = block.chain.run(voxel, entry='all')
            chain_outputs[key] = tmp
            if 'fit' in list(dataset.blocks.keys()):
                key = 'fit'
                block = dataset.blocks[key]
                block.chain.run(voxel, entry='initial_only')
                key = 'spectral'
                block = dataset.blocks[key]
                block.set_do_fit(True, voxel[0])
                tmp = block.chain.run(voxel, entry='all')
                chain_outputs[key] = tmp
        else:
            block = dataset.blocks[key]
            tmp = block.chain.run(voxel, entry='all')
            chain_outputs[key] = tmp

    return chain_outputs


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
            
               
def _import_siemens_dicom(datafile, open_dataset=None):
    """
    Stolen from Analysis main.py module - trimmed for CLI usage
    
    Assumption here is that we are opening one file in the reader. If there is
    an 'open_dataset' sent in, then the current reader is for an associated 
    file. We will associate the current file with the open one at the end of
    the code.
    
    """
    try:

        reader = fileio_dicom_siemens.RawReaderDicomSiemens()
        reader.filenames = [datafile,]
        
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
                if isinstance(raw, mrs_data_raw_wbnaa.DataRawWbnaa):
                    d["prep"] = block_prep_wbnaa.BlockPrepWbnaa
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
        msg = """Unknown exception reading Data file "%s".""" % filename 
        print(msg, file=sys.stderr)
        print(msg, file=sys.stdout)
        sys.exit(-1)

    


def main():

    verbose = True

    datatype = 'siemens dicom'

#     # Processing of SVS_EDIT_OFF files
#     STARTDIR = 'D:\\Users\\bsoher\\projects\\2017_oneil_ucla_lipid_contam_svs\\MRS_RAW_DATA\\_all_svs_edit_off'
#     presetfile      = STARTDIR+'\\_preset_edit_off_v03.xml'
#     preset_ecc_file = STARTDIR+'\\_preset_edit_off_waterpeak_v03.xml' 

    # Processing of SVS_EDIT_DIFF files
    STARTDIR = 'D:\\Users\\bsoher\\projects\\2017_oneil_ucla_lipid_contam_svs\\MRS_RAW_DATA\\_all_svs_edit_diff'
    presetfile      = STARTDIR+'\\_preset_edit_diff_v01.xml'
    preset_ecc_file = STARTDIR+'\\_preset_edit_diff_waterpeak_v01.xml' 


    # this gets all files *.IMA in all subdirectories of STARTDIR
    imafiles = []
    for dirpath, dirnames, filenames in os.walk(STARTDIR):
        for filename in [f for f in filenames if f.endswith(".IMA")]:
            imafiles.append(os.path.join(dirpath, filename))
            print(os.path.join(dirpath, filename))

    if len(imafiles) % 2:
        msg = """IMA directory must have even number of files in 'Data','Ecc Data' order, returning."""
        print(msg, file=sys.stderr)
        print(msg, file=sys.stdout)
        sys.exit(-1)

        

    i = 0

    for datafile, eccfile in zip(imafiles[0::2], imafiles[1::2]):

        # Test input arguments for consistency --------------------------
        
        msg = ''
        if not os.path.isfile(datafile):
            msg = """Main DATAFILE does not exist "%s".""" % datafile 
        if not os.path.isfile(eccfile):
            msg = """ECC DATAFILE does not exist "%s".""" % eccfile 
        if not os.path.isfile(presetfile):
            msg = """PRESETFILE does not exist "%s".""" % presetfile 
        if not os.path.isfile(preset_ecc_file):
            msg = """PRESET_ECC_FILE does not exist "%s".""" % preset_ecc_file 
        if not str(datatype).lower() in SUPPORTED:
            msg = """DATATYPE not supported - "%s".""" % datatype 
            
        if msg:        
            print(msg, file=sys.stderr)
            print(msg, file=sys.stdout)
            sys.exit(-1)
            
        
        # Load Main Dataset --------------------------
        
        if verbose: print("""%s - Load Data into a Dataset object - %s""" % (str(i), datafile))    
        dataset = _import_siemens_dicom(datafile, open_dataset=None)
        dataset = dataset[0]
          
        # Load Preset data for Main Dataset ----------------------------------------------
        
        if verbose: print("Read Preset object")    
        preset = _import_preset(presetfile)
          
        # Load ECC Dataset --------------------------
        
        if verbose: print("Load ECC into a Dataset object")    
        ecc_dataset = _import_siemens_dicom(eccfile, open_dataset=None)
        ecc_dataset = ecc_dataset[0]

        # Load Preset data for ECC Dataset ----------------------------------------------
        
        if verbose: print("Read Preset ECC object")    
        preset_ecc_dataset = _import_preset(preset_ecc_file)

          
        analysis_cli_oneil(dataset, preset, ecc_dataset           = ecc_dataset, 
                                            preset_ecc_dataset    = preset_ecc_dataset, 
                                            watref_dataset        = None, 
                                            preset_watref_dataset = None,
                                            verbose               = verbose, 
                                            debug                 = False)
        
        i += 1
#        if i >= 2: break    # debug statement to exit after one file processed
        
    bob = 10
    bob += 1
    
    
    


        
if __name__ == '__main__':
    main()        
        