
# Python modules
import os
import io
import base64
import traceback
import datetime
import pathlib

# 3rd party modules
import matplotlib
matplotlib.use('Agg')
import numpy as np
from pydicom import Dataset, FileDataset, dcmread, read_file

# Our modules
import vespa.interfaces.inline.vespa_inline_engine as vie

import vespa.analysis.figure_layouts as figure_layouts
import vespa.analysis.fileio.util_philips as util_philips
import vespa.common.util.time_ as util_time
import vespa.common.util.misc as util_misc

from vespa.interfaces.inline.vespa_inline_engine import VespaInlineError, VespaInlineSettings

VERSION = '0.1.0'

#==============================================================================

def run(settings, verbose=True):
    """
    There are 4 processing steps:
    1. collate all files from specified 'datadir', and sort into 'water', 'metab' etc.
    2. load file names into VIE dataset_filename and preset_filename dicts
    3. run files through the Vespa-Analysis inline engine
        3a. (optional) save provenance XML file and/or PNG/PDF images for debugging.
    4. output 'screenshot image' gets put into a pydicom secondary capture RGB DICOM

    """
    msg = ''

    # these are set here for error checking reasons
    fdatasets = {'metab':None, 'water':None, 'ecc':None, 'coil':None}
    fpresets  = {'metab':None, 'water':None, 'ecc':None, 'coil':None}
    dcm_cur = ''

    try:
        settings.vespa_version = util_misc.get_vespa_version()+'-VIE'  # not really up to the user

        # ---------------------------------------------------------------
        # 1. Get filenames from known DATADIR directory and sort
        #    - may move to separate module in future as formats accrue

        if settings.dataformat == 'philips_press28_dicom':
            settings.import_class = 'import_philips_dicom'
            mrs_files = []
            other_files = []
            for dirpath, dirnames, filenames in os.walk(settings.data_dir):
                for filename in filenames:
                    ftest = os.path.join(dirpath, filename)
                    if vie.is_dicom(ftest):
                        dataset = read_file(ftest, defer_size=1024)
                        if util_philips.is_mrs_dicom(dataset):
                            mrs_files.append(ftest)
                            if verbose: print('Found DICOM MRS file - '+ftest)
                    else:
                        other_files.append(ftest)

            if (len(mrs_files) != 2):
                msg  = 'Exception (do_main): Wrong number of DICOM MRS files found in - '+settings.data_dir
                if verbose: print(msg)
                raise VespaInlineError(msg)

            fname_metab, fname_water, fname_ecc, fname_coil = None, None, None, None
            fname_water = mrs_files[0]
            fname_metab = mrs_files[1]

            fname_metab_preset, fname_water_preset, fname_ecc_preset, fname_coil_preset = None, None, None, None
            fname_metab_preset = os.path.join(settings.preset_dir,'preset_philips_dicom_press28_metab.xml')
            fname_water_preset = os.path.join(settings.preset_dir,'preset_philips_dicom_press28_water.xml')

            dcm_cur = dcmread(fname_metab)

        elif settings.dataformat == 'philips_slaser30_cmrr_spar':
            settings.import_class = 'import_philips_spar'
            mrs_files = []
            other_files = []
            for dirpath, dirnames, filenames in os.walk(settings.data_dir):
                for filename in filenames:
                    ftest = os.path.join(dirpath, filename)
                    if pathlib.Path(ftest).suffix in ['spar','sdat','.SPAR','.SDAT']:
                        mrs_files.append(ftest)
                        if verbose: print('Found Spar/Sdat MRS file - '+ftest)
                    else:
                        other_files.append(ftest)

            if len(mrs_files) != 4:
                msg  = 'Exception (do_main): Wrong number of Spar/Sdat datasets found in - '+settings.data_dir
                if verbose: print(msg)
                raise VespaInlineError(msg)

            fname_metab, fname_water, fname_ecc, fname_coil = None, None, None, None
            for fname in mrs_files:
                if '_act.spar' in fname.lower():
                    fname_metab = fname
                if '_ref.spar' in fname.lower():
                    fname_water = fname

            if fname_metab is None: msg += '\nException (do_main): Metabolite data Spar/Sdat not found in - '+settings.data_dir
            if fname_water is None: msg += '\nException (do_main): Water reference Spar/Sdat not found in - '+settings.data_dir
            if msg:
                if verbose: print(msg)
                raise VespaInlineError(msg)

            fname_metab_preset, fname_water_preset, fname_ecc_preset, fname_coil_preset = None, None, None, None
            fname_metab_preset = os.path.join(settings.preset_dir,'preset_philips_slaser30_cmrr_spar_metab.xml')
            fname_water_preset = os.path.join(settings.preset_dir,'preset_philips_slaser30_cmrr_spar_water.xml')

            dcm_cur = ''

        # ----------------------------------------------------------
        # 2. load filenames into parameter dicts

        fdatasets['metab'] = fname_metab
        fdatasets['water'] = fname_water
        fdatasets['ecc']   = fname_ecc
        fdatasets['coil']  = fname_coil

        fpresets['metab'] = fname_metab_preset
        fpresets['water'] = fname_water_preset
        fpresets['ecc']   = fname_ecc_preset
        fpresets['coil']  = fname_coil_preset

        fbasis_mmol = None       # mmol   fbase+'\\basis_mmol_dataset_seadMM2014_truncat2048pts_normScale100dc015.xml'

        # ----------------------------------------------------------
        # 3. Run the processing

        params = [fdatasets, fpresets, fbasis_mmol, settings]
        png_buf, pdf_buf, _ = vie.analysis_kernel( params, verbose=verbose )
        buf_shape = int(10.24 * settings.png_dpi), int(10.24 * settings.png_dpi), 3

    except VespaInlineError as e:

        if verbose: print('Exception: VespaInlineError - see error report')
        trace = ''      # returned in the VespaInlineError msg
        png_buf = do_error_processing(e, fdatasets, fpresets, trace, settings)
        buf_shape = int(10.24 * settings.err_dpi), int(10.24 * settings.err_dpi), 3

    except Exception as e:

        if verbose: print('Exception: GeneralException - see error report')
        trace = traceback.format_exc()
        png_buf = do_error_processing(e, fdatasets, fpresets, trace, settings)
        buf_shape = int(10.24 * settings.err_dpi), int(10.24 * settings.err_dpi), 3

    # ----------------------------------------------------------
    # 4.  dump dcm_buf to a DICOM RGB file ....

    if settings.save_dcm:
        dcm_out = rgb2dcm(png_buf, dcm_cur, buf_shape)
        dcm_out.save_as(settings.dcm_fname)

    if settings.save_dcm_pdf:
        dcm_pdf_out = pdf2dcm(pdf_buf, dcm_cur)
        dcm_pdf_out.save_as(settings.dcm_pdf_fname)

    if verbose: print('fname_dicom = ' + settings.dcm_fname)

    return


def do_error_processing(e, fdatasets, fpresets, trace, settings):

    fig = figure_layouts.inline_error(e, fdatasets, fpresets, trace,
                                         fontname='Courier New',
                                         dpi=settings.err_dpi)
    dcm_buf = fig[0].canvas.tostring_rgb()
    dcm_buf = np.fromstring(dcm_buf, dtype=np.uint8)

    if settings.save_err:
        fname, _ = os.path.splitext(settings.err_fname)
        if settings.err_fname_unique:
            fname += util_time.filename_timestamp()  # yyyymmdd.hhmmss.usecs

        fig[0].savefig(fname + '.png',
                       dpi=settings.err_dpi,
                       pad_inches=settings.err_pad_inches)
    return dcm_buf


SSC_TEMPLATE = b'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'

def rgb2dcm(png_buf, dcm_cur, buf_shape):

    dcm_ssc = dcmread(io.BytesIO(base64.b64decode(SSC_TEMPLATE)))

    if not isinstance(png_buf, np.ndarray):
        raise VespaInlineError('rgb2dcm(): png_buf was not a Numpy array, returning!')

    png_buf = png_buf.tobytes()

    if dcm_cur == '':
        dcm_cur = dcmread(io.BytesIO(base64.b64decode(SSC_TEMPLATE)))  # just a dummy default

    dt = datetime.datetime.now()                # Create time format strings and ids
    date_str = dt.strftime('%Y%m%d')
    time_str = dt.strftime('%H%M%S.%f')[:-3]  # long format with milliseconds
    unique_ssc_str = '1.3.46.670589.11.71.5.0.10236.' + date_str + time_str.replace('.', '')

    dcm_ssc.file_meta.MediaStorageSOPInstanceUID = unique_ssc_str                   # (0002,0003)
    dcm_ssc.InstanceCreationDate = dt.strftime('%Y%m%d')                            # (0008,0012)
    dcm_ssc.InstanceCreationTime = time_str                                         # (0008,0013)
    dcm_ssc.SOPInstanceUID = dcm_ssc.file_meta.MediaStorageSOPInstanceUID           # (0008,0018)
    dcm_ssc.SeriesDate = date_str                                                   # (0008,0021)
    dcm_ssc.SeriesTime = time_str                                                   # (0008,0031)
    dcm_ssc.AcquisitionDate = date_str                                              # (0008,0022)
    dcm_ssc.AcquisitionTime = time_str                                              # (0008,0032)
    dcm_ssc.ContentDate = date_str                                                  # (0008,0023)
    dcm_ssc.ContentTime = time_str                                                  # (0008,0033)
    dcm_ssc.Manufacturer = dcm_cur.Manufacturer                                     # (0008,0070)
    dcm_ssc.SeriesDescription = "xReport_"                                          # (0008,103E)
    dcm_ssc.PatientName = dcm_cur.PatientName                                       # (0010,0010)
    dcm_ssc.PatientID = dcm_cur.PatientID                                           # (0010,0020)
    dcm_ssc.PatientBirthDate = dcm_cur.PatientBirthDate                             # (0010,0030)
    dcm_ssc.PatientSex = dcm_cur.PatientSex                                         # (0010,0030)
    dcm_ssc.PatientWeight = dcm_cur.PatientWeight                                   # (0010,1030)
    dcm_ssc.DateOfSecondaryCapture = date_str                                       # (0018,1012)
    dcm_ssc.TimeOfSecondaryCapture = time_str                                       # (0018,1014)
    dcm_ssc.ProtocolName = dcm_cur.ProtocolName + time_str                          # (0018,1030)
    dcm_ssc.StudyID = dcm_cur.StudyID                       # do not change         # (0020,0010)
    dcm_ssc.StudyInstanceUID = dcm_cur.StudyInstanceUID     # do not change         # (0020,000D)

    new_id = unique_ssc_str + datetime.datetime.now().strftime('%Y%m%d') + datetime.datetime.now().strftime('%H%M%S%f')

    dcm_ssc.SeriesInstanceUID = new_id                                              # (0020,0010)
    dcm_ssc.SeriesNumber = ''
    dcm_ssc.AcquisitionNumber = ''                                                  # (0020,0012)
    # dcm_ssc.AcquisitionNumber = getattr(dcm_cur,'AcquisitionNumber',1)            # (0020,0012)
    dcm_ssc.SamplesPerPixel = int(3)                                                # (0028,0002)
    dcm_ssc.PhotometricInterpretation = 'RGB'                                       # (0028,0004)
    dcm_ssc.Rows = buf_shape[0]                                                     # (0028,0010)
    dcm_ssc.Columns = buf_shape[1]                                                  # (0028,0011)
    dcm_ssc.PlanarConfiguration = 0               # np shape=[H,W,Ch]               # (0028,0006)

    # dcm_ssc.SecondaryCaptureDeviceManufacturer          = 'Philips'               # (0018,1016)
    # dcm_ssc.SecondaryCaptureDeviceManufacturerModelName = 'ISD v3 CNode'          # (0018,1018)
    # dcm_ssc.SecondaryCaptureDeviceSoftwareVersions      = '5.6.1'                 # (0018,1019)

    dcm_ssc.PixelData = png_buf

    return dcm_ssc


def pdf2dcm(pdf_buf, dcm_cur):

    if not isinstance(pdf_buf, io.BytesIO):
        raise VespaInlineError('pdf2dcm(): pdf_buf was not a byte array, returning!')

    if dcm_cur == '':
        dcm_cur = dcmread(io.BytesIO(base64.b64decode(SSC_TEMPLATE)))  # just a dummy default

    dt = datetime.datetime.now()                 # Create time format strings and ids
    date_str = dt.strftime('%Y%m%d')
    time_str = dt.strftime('%H%M%S.%f')[:-3]     # long format with milliseconds
    unique_ssc_str = '1.3.46.670589.11.71.5.0.10236.' + date_str + time_str.replace('.', '')

    meta = Dataset()
    meta.MediaStorageSOPClassUID    = '1.2.840.10008.5.1.4.1.1.104.1'
    meta.MediaStorageSOPInstanceUID = '2.16.840.1.114430.287196081618142314176776725491661159509.60.1'
    meta.ImplementationClassUID     = '1.3.46.670589.50.1.8.0'

    ds = FileDataset(None, {}, file_meta=meta, preamble=b"\0" * 128)
    ds.is_little_endian = True
    ds.is_implicit_VR   = True
    ds.ContentDate = dt.strftime('%Y%m%d')
    ds.ContentTime = dt.strftime('%H%M%S.%f')
    ds.SOPClassUID = '1.2.840.10008.5.1.4.1.1.104.1'
    ds.MIMETypeOfEncapsulatedDocument = 'application/pdf'

    new_id = unique_ssc_str + datetime.datetime.now().strftime('%Y%m%d') + datetime.datetime.now().strftime('%H%M%S%f')

    ds.Modality = dcm_cur.Modality
    ds.PatientName = dcm_cur.PatientName
    ds.PatientID = dcm_cur.PatientID
    ds.PatientSex = dcm_cur.PatientSex
    ds.PatientWeight = dcm_cur.PatientWeight
    ds.PatientBirthDate = dcm_cur.PatientBirthDate
    ds.StudyInstanceUID = dcm_cur.StudyInstanceUID  # do not change
    ds.SeriesInstanceUID = new_id
    ds.SOPInstanceUID = ds.file_meta.MediaStorageSOPInstanceUID
    ds.MIMETypeOfEncapsulatedDocument = 'application/pdf'
    ds.Manufacturer = dcm_cur.Manufacturer
    ds.SeriesDescription = "xReport_encap_pdf"
    ds.StudyID = dcm_cur.StudyID  # do not change
    ds.SeriesNumber = ''
    ds.AcquisitionNumber = ''
    ds.EncapsulatedDocument = pdf_buf.getvalue()

    return ds



if __name__ == '__main__':

    fname_tstamp = util_time.filename_timestamp()  # yyyymmdd.hhmmss.usecs

    # get defaults
    settings = VespaInlineSettings()

    # reset relative to 'this' filename, not to vespa_inline_engine location
    settings.base_path = os.path.dirname(os.path.abspath(__file__))
    settings.data_dir = os.path.join(settings.base_path, 'datadir')
    settings.preset_dir = os.path.join(settings.base_path, 'presets')
    settings.output_dir = os.path.join(settings.base_path, 'output')
    settings.debug_dir = os.path.join(settings.base_path, 'debug')

    #settings.dataformat = 'philips_press28_dicom'
    settings.dataformat = 'philips_slaser30_cmrr_spar'

    settings.save_err = True
    settings.save_xml = True
    settings.save_pdf = True
    settings.save_png = True
    settings.save_dcm = True
    settings.save_dcm_pdf = True

    settings.err_fname_unique = True
    settings.xml_fname_unique = True
    settings.pdf_fname_unique = True
    settings.png_fname_unique = True
    settings.dcm_fname_unique = True
    settings.dcm_pdf_fname_unique = True

    settings.err_fname = os.path.join(settings.debug_dir, "debug_vespa_viff.png")
    settings.xml_fname = os.path.join(settings.debug_dir, "debug_xml_last_run.xml")
    settings.pdf_fname = os.path.join(settings.debug_dir, "debug_pdf_philips.pdf")
    settings.png_fname = os.path.join(settings.output_dir,     "results_vespa_inline_philips.png")
    settings.dcm_fname = os.path.join(settings.output_dir,     "results_vespa_inline_dicom.dcm")
    settings.dcm_pdf_fname = os.path.join(settings.output_dir, "results_vespa_inline_dicom_pdf.dcm")

    settings.pdf_plotstyle = 'lcm_multi'
    settings.pdf_file_label = 'Analysis- Philips PRIDE Inline'
    settings.pdf_minppm = 0.5
    settings.pdf_maxppm = 4.2
    settings.pdf_apply_phase = False
    settings.pdf_remove_base = False
    settings.pdf_fontname = 'Courier New'
    settings.pdf_dpi = 300
    settings.pdf_pad_inches = 0.5

    settings.png_plotstyle = 'lcm_square'
    settings.png_file_label = 'Analysis- Philips PRIDE Inline'
    settings.png_minppm = 0.5
    settings.png_maxppm = 4.2
    settings.png_apply_phase = False
    settings.png_remove_base = False
    settings.png_fontname = 'Courier New'
    settings.png_dpi = 100  # gives 2048x2048 image
    settings.png_pad_inches = 0.5

    settings.err_dpi = 200
    settings.err_pad_inches = 0.5

    settings.debug = False

    run(settings)


