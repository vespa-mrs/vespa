
Instructions for creating conda environment for Vespa 1.0.0.rc3

1. create an environment with python 3.7.x 

2. Using 'conda install', add these packages 
    numpy
    scipy
    matplotlib
    wxPython=4.0.x  (can not be 4.1.x)
    
3. Using 'conda install -c conda-forge', add these packages
    lmfit
    pydicom
    pypubsub
    
4, Using 'pip', do:
    >pip install packaging
    >pip install -i https://test.pypi.org/simple/ pygamma
    >pip install Vespa_Suite-1.0.0rc3-py37-none-any.whl
    
    
Main module for running Philips inline processing:

    vespa.interfaces.inline.philips.run_inline_vespa_philips
    
The module that loads files, runs pipeline, creates output (but does not always save them to file):

    vespa.interfaces.inline.vespa_inline_engine   ( method = analysis_kernel() )
    
Create a VespaInlineSettings object from vespa_inline_engine to control what processing is done, in which directories, and what output filenames are created.
    



===================================

Issue 01)  with creating "min Py38 install" on OpenSUSE 11.4 for use on MRIR (fingers crossed)

Did a min install on VM: >conda, then >conda install numpy scipy matplotlib configobj, then >conda install -c conda-forge lmfit pydicom

Test ran the Philips run_inline_vespa_philips.py for press28 dicom files and caught an error after presets read: 

    (base) bsoher@linux-rb4d:~/Downloads> python run_inline_vespa_philips.py 
    Found DICOM MRS file - /home/bsoher/Downloads/datadir/XX_0025_met
    Found DICOM MRS file - /home/bsoher/Downloads/datadir/XX_0022_wat
    single - Begin - metab filename = /home/bsoher/Downloads/datadir/XX_0022_wat
    single - load_preset - fname = /home/bsoher/Downloads/presets/preset_philips_dicom_press28_metab.xml
    single - load_preset - fname = /home/bsoher/Downloads/presets/preset_philips_dicom_press28_water.xml
    python: symbol lookup error: /home/bsoher/miniconda3/lib/python3.8/site-packages/mkl/../../../libmkl_intel_thread.so.1: undefined symbol: omp_get_num_procs

Sought help online. Ended up installing (from conda-forge) llvm-openmp package, since it has the lib that was missing, too. And it finished.
https://github.com/ContinuumIO/anaconda-issues/issues/10195



Issue 02) ran the Philips run_inline_vespa_philips.py for press28 dicom files:

XX_0022_wat and XX_0025_met

In browser, these files show up in this order, same as on Windows

In fitting data, the water file was fitted as metab and metab as water.  When I rename XX_0022_wat to XX_1022_wat the fit is GOOD, same as on Win

So, somewhere in reading files in Python the order gets flipped for Philips. This was always a 'magic number' assumption anyway, but for now we ignore since we are trying to get this to work on Siemens
