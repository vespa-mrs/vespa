# Python modules
import os
import sys
import importlib

# 3rd party modules
import numpy as np

# Our modules
import vespa.analysis.mrs_dataset as mrs_dataset
import vespa.analysis.util_import as util_import
import vespa.analysis.fileio.util_exceptions as util_exceptions
import vespa.common.util.misc as util_misc
import vespa.common.configobj as configobj

from vespa.common.mrs_data_raw import DataRawFidsum



_MSG_MULTIFILE_ATTRIBUTE_MISMATCH = """
The dimensions and/or sweep widths in the selected data files differ.
 
Multifile selection is only available for single voxel data files.
 
Please select single voxel data files with same dimensions and sweep width."""
 
_MSG_MULTIFILE_TYPE_MISMATCH = """
More that one set of raw data was read from the file(s) selected, however, not all were of the same format type (summed FIDs or not).
 
Where multifile selection is allowed, or multiple datasets are being read out of one file, all files must return the same types of data."""
 

_MSG_OPEN_ATTRIBUTE_MISMATCH = """
The dimensions and/or sweep width of the currently open datasets differ from those of the file(s) you asked to open.
  
You can open these files, but first you have to close all currently open datasets.
"""
  
_MSG_OPEN_TYPE_MISMATCH = """
The currently open datasets differ from those of the file(s) you asked to open in that they are not all of the same format type (summed FIDs or not).
 
You can open these files, but first you have to close all currently open datasets.
"""
 
_MSG_UNSUPPORTED_DIMENSIONALITY = """
The file(s) you opened contains SI datasets. Vespa doesn't support SI at this time.
"""

_MSG_INCORRECT_DIMENSIONALITY = """
The file you selected is likely an SVS file but has a dimension that is not supported by this Import class. One possibility is that this is a Fidsum format.

Please try opening this file with a different Import format.
"""
 
_MSG_PRESET_MISMATCH = """
One or more of the selected VIFF files is an Analysis Preset file. These can not be opened as datasets.
"""
 
_MSG_OPEN_ATTRIBUTE_MISMATCH = """
The dimensions and/or sweep width of the currently open datasets differ from those of the file(s) you asked to open.
 
You can open these files, but first you have to close all currently open datasets.
"""
 
_MSG_OPEN_ZEROFILL_MISMATCH = """
The zerofill factor of the currently open datasets differ from those of the file you asked to open.
 
You can open this file, but the zero fill factor of the open datasets needs to be changed.
"""
 
_MSG_NO_DATASETS_FOUND = """The file "%s" doesn't contain any datasets."""

_MSG_INCOMPLETE_HEADER_PARAMETERS = """Import header missing a necessary parameter -> """




########################################################################
# This is a non-GUI based module that sets up a dictionary that maps
# import data types to Python modules with parser classes/objects. Vespa
# Analysis used to have a few hard coded 'standard' formats (eg. siemens
# DICOM, siemens *.rda, GE probe, etc.) but we have switched over to
# map all formats through this Utility module. It also provides CLI
# (command line interface) usage without GUI imports.
#
# NB. If no 'filename' is provided to the program, eg in CLI mode, then
#     only the default data classes are returned.
#
########################################################################


# these are the default import data classes - used to be hard coded in main()

STANDARD_CLASSES_LOFL = [
'[separator100]',
'[import_bruker]',
    'path=vespa.analysis.fileio.bruker',
    'class_name=RawReaderBruker',
    'menu_item_text=Bruker',
    'ini_file_name=import_bruker',
'[import_ge_probep]',
    'path=vespa.analysis.fileio.ge_pfile',
    'class_name=RawReaderGeProbep',
    'menu_item_text=GE PROBE-P (*.7)',
    'ini_file_name=import_ge_probep',
'[import_ge_oslaser_cmrr]',
    'path=vespa.analysis.fileio.ge_pfile',
    'class_name=RawReaderGeOslaserCmrr',
    'menu_item_text=GE sLASER CMRR (*.7)',
    'ini_file_name=import_ge_oslaser_cmrr',
'[import_nifti_mrs]',
    'path=vespa.analysis.fileio.nifti_mrs',
    'class_name=RawReaderNiftiMrs',
    'menu_item_text=NIfTI-MRS (*.nii/.nii.gz)',
    'ini_file_name=import_nifti_mrs',
'[import_philips_dicom]',
    'path=vespa.analysis.fileio.dicom_philips',
    'class_name=RawReaderDicomPhilips',
    'menu_item_text=Philips DICOM',
    'ini_file_name=import_philips_dicom',
'[import_philips_dicom_fidsum]',
    'path=vespa.analysis.fileio.dicom_philips',
    'class_name=RawReaderDicomPhilipsFidsum',
    'menu_item_text=Philips DICOM Sum FIDs',
    'ini_file_name=import_philips_dicom_fidsum',
'[import_philips_spar]',
    'path=vespa.analysis.fileio.philips_spar',
    'class_name=RawReaderPhilipsSpar',
    'menu_item_text=Philips (*.spar/sdat)',
    'ini_file_name=import_philips_spar',
'[import_philips_fidsum]',
    'path=vespa.analysis.fileio.philips_fidsum',
    'class_name=RawReaderPhilipsFidsum',
    'menu_item_text=Philips Sum FIDs',
    'ini_file_name=import_philips_fidsum',
'[import_siemens_dicom]',
    'path=vespa.analysis.fileio.dicom_siemens',
    'class_name=RawReaderDicomSiemens',
    'menu_item_text=Siemens DICOM',
    'ini_file_name=import_siemens_dicom',
'[import_siemens_dicom_fidsum]',
    'path=vespa.analysis.fileio.dicom_siemens',
    'class_name=RawReaderDicomSiemensFidsum',
    'menu_item_text=Siemens DICOM Sum FIDs',
    'ini_file_name=import_siemens_dicom_fidsum',
'[import_siemens_dicom_timeseries]',
    'path=vespa.analysis.fileio.dicom_siemens_timeseries',
    'class_name=RawReaderDicomSiemensTimeseries',
    'menu_item_text=Siemens DICOM Timeseries',
    'ini_file_name=import_siemens_dicom_timeseries',
'[import_siemens_rda]',
    'path=vespa.analysis.fileio.siemens_rda',
    'class_name=RawReaderSiemensRda',
    'menu_item_text=Siemens Export (*.rda)',
    'ini_file_name=import_siemens_rda',

'[siemens_twix_svs_se]',
    'path=vespa.analysis.fileio.siemens_twix_svs_se',
    'class_name=RawReaderSiemensTwixSvsSe',
    'menu_item_text=Siemens Twix svs_se',
    'ini_file_name=import_siemens_twix_svs_se',
'[siemens_twix_svs_slaser_cmrr_vb]',
    'path=vespa.analysis.fileio.siemens_twix_slaser_cmrr',
    'class_name=RawReaderSiemensTwixSlaserCmrrVb',
    'menu_item_text=Siemens Twix sLASER CMRR VB',
    'ini_file_name=import_siemens_twix_slaser_cmrr_vb',
'[siemens_twix_slaser_cmrr_ve]',
    'path=vespa.analysis.fileio.siemens_twix_slaser_cmrr',
    'class_name=RawReaderSiemensTwixSlaserCmrrVe',
    'menu_item_text=Siemens Twix sLASER CMRR VE',
    'ini_file_name=import_siemens_twix_slaser_cmrr_ve',
'[siemens_twix_svs_edit]',
    'path=vespa.analysis.fileio.siemens_twix_svs_edit',
    'class_name=RawReaderSiemensTwixSvsEdit',
    'menu_item_text=Siemens Twix svs_edit WIP529',
    'ini_file_name=import_siemens_twix_svs_edit',

'[siemens_twix]',
    'path=vespa.analysis.fileio.siemens_twix',
    'class_name=RawReaderSiemensTwix',
    'menu_item_text=Siemens Twix (generic)',
    'ini_file_name=import_siemens_twix_generic',
'[import_varian]',
    'path=vespa.analysis.fileio.varian',
    'class_name=RawReaderVarian',
    'menu_item_text=Varian',
    'ini_file_name=import_varian',
'[import_vasf]',
    'path=vespa.analysis.fileio.vasf',
    'class_name=RawReaderVasf',
    'menu_item_text=VASF (*.rsd/rsp)',
    'ini_file_name=import_vasf',
'[import_vasf_fidsum]',
    'path=vespa.analysis.fileio.vasf',
    'class_name=RawReaderVasfFidsum',
    'menu_item_text=VASF Sum FIDs (*.rsd/rsp)',
    'ini_file_name=import_vasf_fidsum',
'[separator101]',
# '[import_siemens_dicom_cmrr_mpress]',
#     'path=vespa.analysis.fileio.dicom_siemens_cmrr_mpress',
#     'class_name=RawReaderDicomSiemensXaMpress',
#     'menu_item_text=SiemensXA DICOM CMRR MPress',
#     'ini_file_name=import_siemens_dicom_cmrr_mpress',
# '[import_siemens_dicom_cmrr_mpress_fidsum]',
#     'path=vespa.analysis.fileio.dicom_siemens_cmrr_mpress',
#     'class_name=RawReaderDicomSiemensFidsumXaMpress',
#     'menu_item_text=SiemensXA DICOM CMRR MPress Fidsum',
#     'ini_file_name=import_siemens_dicom_cmrr_mpress_fidsum',
'[import_siemens_dicom_eja_svs_mpress]',
    'path=vespa.analysis.fileio.dicom_siemens_eja_svs_mpress',
    'class_name=RawReaderDicomSiemensXaEjaSvsMpress',
    'menu_item_text=SiemensXA DICOM eja_svs_mpress',
    'ini_file_name=import_siemens_dicom_eja_svs_mpress',
'[import_siemens_dicom_eja_svs_mpress_onoff]',
    'path=vespa.analysis.fileio.dicom_siemens_eja_svs_mpress',
    'class_name=RawReaderDicomSiemensXaEjaSvsMpressOnOff',
    'menu_item_text=SiemensXA DICOM eja_svs_mpress - On/Off only',
    'ini_file_name=import_siemens_dicom_eja_svs_mpress_onoff',
#'[import_siemens_dicom_eja_svs_mpress_onoff_indiv]',
#    'path=vespa.analysis.fileio.dicom_siemens_eja_svs_mpress',
#    'class_name=RawReaderDicomSiemensXaEjaSvsMpressOnOffIndiv',
#    'menu_item_text=SiemensXA DICOM eja_svs_mpress (indiv) - On/Off only',
#    'ini_file_name=import_siemens_dicom_eja_svs_mpress_onoff_indiv',
'[import_siemens_dicom_eja_svs_mpress_fidsum]',
    'path=vespa.analysis.fileio.dicom_siemens_eja_svs_mpress',
    'class_name=RawReaderDicomSiemensFidsumXaEjaSvsMpress',
    'menu_item_text=SiemensXA DICOM eja_svs_mpress Fidsum',
    'ini_file_name=import_siemens_dicom_eja_svs_mpress_fidsum',
'[import_siemens_dicom_eja_svs_mpress_fidsum_onoff]',
    'path=vespa.analysis.fileio.dicom_siemens_eja_svs_mpress',
    'class_name=RawReaderDicomSiemensFidsumXaEjaSvsMpressOnOff',
    'menu_item_text=SiemensXA DICOM eja_svs_mpress Fidsum - On/Off only',
    'ini_file_name=import_siemens_dicom_eja_svs_mpress_fidsum_onoff',
'[import_siemens_dicom_eja_svs_mpress_fidsum_onoff_adv]',
    'path=vespa.analysis.fileio.dicom_siemens_eja_svs_mpress',
    'class_name=RawReaderDicomSiemensFidsumXaEjaSvsMpressOnOffAdv',
    'menu_item_text=SiemensXA DICOM Advanced eja_svs_mpress Fidsum - On/Off only',
    'ini_file_name=import_siemens_dicom_eja_svs_mpress_fidsum_onoff_adv',
'[separator102]',
]


def set_import_data_classes(filename=''):
    """
    This module creates a configObj with all available data parsers
    - There are standard data classes we guarantee are provided
    - (optional) There are user provided data classes from a file
    We gather both here and then merge them

    """
    full_cfg = configobj.ConfigObj(STANDARD_CLASSES_LOFL, encoding="utf-8")
    if filename:
        user_cfg = configobj.ConfigObj(filename, encoding="utf-8")
        full_cfg.merge(user_cfg)
    
    # fill dict with reader objects and import paths for reader classes
    items = { }
    msg   = ''
    
    for module_name in list(full_cfg.keys()):
        if not'separator' in module_name.lower():
            section = full_cfg[module_name]
            path    = section["path"]
            if 'vespa.analysis.fileio' in path:
                module = importlib.import_module(path)
            else:
                # to be polite we check to be sure that the path and file exists,
                if os.path.exists(path):
                    spec   = importlib.util.spec_from_file_location(module_name, path)
                    module = importlib.util.module_from_spec(spec)
                    spec.loader.exec_module(module)
                    sys.modules[module_name] = module
                else:
                    msg1 = """\nI couldn't find the file "{0}" referenced in {1}.""".format(path,filename)
                    msg += msg1.format(path)
                    continue

            # Save the reader class & INI file name associated with this menu item.
            klass = getattr(module, section["class_name"])
            items[module_name] = (klass, section["ini_file_name"])
            
    return items, full_cfg, msg


def get_datasets_cli(datafnames, section_name, open_dataset=None):
    """
    This is only called by command line (CLI) scripts. It shortcuts all the GUI
    complexity. It does assume that Vespa is installed and has INI files in 
    their normal locations.
    
    """
    if (sys.version_info > (3, 0)):
        # Python 3 code in this block
        fname = os.path.join(util_misc.get_data_dir(), "analysis_import_menu_additions.ini")
    else:
        fname = os.path.join(util_misc.get_data_dir(), "analysis_import_menu_additions_py2.ini")
    classes, full_cfg, msg = set_import_data_classes(filename=fname)

    if section_name not in list(classes.keys()):
        raise ValueError('Could not find section named %s in FileImportClasses dict' % (section_name,)) 
    reader, _ = classes[section_name]
    reader = reader()
    if isinstance(datafnames,(list,tuple)):
        reader.filenames = datafnames
    else:
        reader.filenames = [datafnames,]
    datasets, msg = get_datasets(reader, open_dataset)

    return datasets
    

def get_datasets(reader, open_dataset=None):
    """ This is a 3 step process - coming in, assumes 'reader' is set to go """
    msg = ''
    datasets = []

    # Step 1 - get DataRawXxxx object(s), shows what type of data we have

    try:
        raws = reader.read_raws(open_dataset=open_dataset)

    except IOError:
        msg = "One or more of the files couldn't be read due to a disk error."
    except util_exceptions.MultifileAttributeMismatchError:
        msg = _MSG_MULTIFILE_ATTRIBUTE_MISMATCH
    except util_exceptions.MultifileTypeMismatchError:
        msg = _MSG_MULTIFILE_TYPE_MISMATCH
    except util_exceptions.UnsupportedDimensionalityError:
        msg = _MSG_UNSUPPORTED_DIMENSIONALITY               # also catches SIDataError
    except util_exceptions.IncorrectDimensionalityError as e:
        msg = _MSG_INCORRECT_DIMENSIONALITY+'\n'+str(e)                 # also catches SIDataError
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
    except util_exceptions.IncompleteHeaderParametersError as e:
        if not e:
            e = 'None'
        msg = str(_MSG_INCOMPLETE_HEADER_PARAMETERS) + str(e)

    # Step 2 - see if any special 'blocks' are required for dataset class

    if not msg:

        block_class_specs = []
        zf_mult = open_dataset.zero_fill_multiplier if open_dataset else 0

        for raw in raws:
            d = { }

            raw_class = type(raw).__name__

            # these associate multiple files in the Raw block (e.g. On/Off/Sum/Dif for edited data)
            if raw_class == 'DataRawProbep':
                d["raw"] = getattr(importlib.import_module('vespa.analysis.block_raw_probep'), 'BlockRawProbep')
            elif raw_class == 'DataRawEdit':
                d["raw"] = getattr(importlib.import_module('vespa.analysis.block_raw_edit'), 'BlockRawEdit')
            elif raw_class == 'DataRawEditFidsum':
                d["raw"] = getattr(importlib.import_module('vespa.analysis.block_raw_edit_fidsum'), 'BlockRawEditFidsum')
            elif raw_class == 'DataRawCmrrSlaser':
                d["raw"] = getattr(importlib.import_module('vespa.analysis.block_raw_cmrr_slaser'), 'BlockRawCmrrSlaser')

            # these require a Preprocess tab to be included in the workflow
            if isinstance(raw, DataRawFidsum):
                d["prep"] = getattr(importlib.import_module('vespa.analysis.block_prep_fidsum'), 'BlockPrepFidsum')

            # these have special Preprocess tabs to include in workflow
            if raw_class == 'DataRawTimeseries':
                d["prep"] = getattr(importlib.import_module('vespa.analysis.block_prep_timeseries'), 'BlockPrepTimeseries')
            elif raw_class == 'DataRawWbnaa':
                d["prep"] = getattr(importlib.import_module('vespa.analysis.block_prep_wbnaa'), 'BlockPrepWbnaa')

            block_class_specs.append(d)

        f = lambda raw, block_classes: mrs_dataset.dataset_from_raw(raw, block_classes, zf_mult)
        datasets = list(map(f, raws, block_class_specs))

    # Step 3 - set any 'associated' datasets

    if datasets:
        for i,dataset in enumerate(datasets):

            # GE PROBE-P - returns both water and metabolite data
            if type(dataset.blocks['raw']).__name__ == 'BlockRawProbep':
                if len(datasets) > 1:
                    dataset.blocks['raw'].set_associated_datasets([datasets[0], datasets[1]])

            # Edited SVS Object - four datasets, acquisition ON/OFF and calculated SUM/DIFF.
            if len(datasets) == 4:
                if type(raws[0]).__name__ == 'DataRawEditFidsum':
                    dataset.blocks['raw'].set_associated_datasets([datasets[0], datasets[1], datasets[2], datasets[3]])
                if type(raws[0]).__name__ == 'DataRawEdit':
                    dataset.blocks['raw'].set_associated_datasets([datasets[0], datasets[1], datasets[2], datasets[3]])
            elif len(datasets) == 2:
                if type(raws[0]).__name__ == 'DataRawEditFidsum':
                    dataset.blocks['raw'].set_associated_datasets([datasets[0], datasets[1]])
                if type(raws[0]).__name__ == 'DataRawEdit':
                    dataset.blocks['raw'].set_associated_datasets([datasets[0], datasets[1]])

            # CMRR sLASER - 3 or 6 datasets:
            #   coil 1 FID, ecc1 2 FIDs, water1 2 FIDs,  metab64 64 FIDs, ecc2 2 FIDs (opt) and water2 2 FIDs (opt)
            if type(raws[0]).__name__ == 'DataRawCmrrSlaser':
                if len(datasets) == 6:
                    dataset.blocks['raw'].set_associated_datasets([datasets[0], datasets[1], datasets[2], datasets[3], datasets[4], datasets[5]])
                elif len(datasets) == 3:
                    dataset.blocks['raw'].set_associated_datasets([datasets[0], datasets[1], datasets[2]])
                elif len(datasets) == 4:
                    dataset.blocks['raw'].set_associated_datasets([datasets[0], datasets[1], datasets[2], datasets[3]])

    return datasets, msg




#------------------------------------------------------------------------------

def open_viff_dataset_file(filenames):
    """
    VIFF - Vespa Interchange File Format - is the XML format for data
    saved from or opened up into the Analysis application. This is the
    only format that is actually 'opened' by Analysis, all other formats
    are considered 'imports'.

    Note that we only allow people to open a single VIFF file as opposed
    to DICOM, VASF or other imported formats where we allow users to
    open multiple files which are then concatenated into one Dataset
    object.

    If the open is successful (and the dimensions match to any existing
    data), the dataset is opened into a new dataset tab.

    The Dataset object is returned (or None if the user
    doesn't choose a file), along with a list of the filenames opened.
    """
    datasets = []
    for filename in filenames:
        msg = ""
        try:
            importer = util_import.DatasetImporter(filename)
        except IOError:
            msg = """I can't read the file "%s".""" % filename
        except SyntaxError:
            msg = """The file "%s" isn't valid Vespa Interchange File Format.""" % filename

        if msg:
            return [None,], msg
        
        # Time to rock and roll!
        dsets = importer.go()

        for dataset in dsets:
            # check to ensure that none of the selected files is
            # actually an Analysis Preset file
            if dataset.behave_as_preset:
                msg = "Analysis - Preset Filetype Mismatch - No data in Preset file, can't load"
                return [None,], msg
            
        for item in dsets:
            datasets.append(item)


    if datasets:

        for dataset in datasets:
            
            # Check that all datasets align to the first. The attributes
            # open_dataset.raw_dims of the currently open dataset(s) must match
            # those of the dataset(s) that we're trying to open. To compare, we
            # grab one of the currently open datasets. It doesn't matter which
            # one since they all have matching attributes.
            #
            # Note. Dimensionality rules also apply to zerofill

            open_dataset = datasets[0]
            if (dataset.raw_dims[0] == open_dataset.raw_dims[0]) and \
               (np.round(dataset.sw,2) == np.round(open_dataset.sw,2)):
                # All is well!
                pass
            else:
                # The dimensions don't match. We can't open these files.
                return [None,], _MSG_OPEN_ATTRIBUTE_MISMATCH

            if (dataset.spectral_dims == open_dataset.spectral_dims):
                # All is well!
                pass
            else:
                # The zerofill factors don't match. We can't open these files.
                return [None,], _MSG_OPEN_ZEROFILL_MISMATCH


        for dataset in datasets:
            dataset.set_associated_datasets(datasets)
            if dataset.id == datasets[-1].id:
                dataset.dataset_filename = filename
                # dataset.filename is only at run-time to save the name of the
                # VIFF file read in. It is not for saving a filename derived
                # from the raw data file(s) used to create the dataset object.
                # We set this filename only for the primary dataset, not the
                # associated datasets.
            else:
                dataset.dataset_filename = ''
                
        return datasets, ''

    else:
        if not datasets:
            return [None,], _MSG_NO_DATASETS_FOUND % filename
