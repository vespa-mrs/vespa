
import exceptions  - not a module anymore?

need  pip install hlsvdpro worked, but need to test
need  pip install lmfit
need  pip install pypubsub

import dicom --> import pydicom as dicom

change wintypes.create_unicode_buffer() to ctypes.create_unicode_buffer()  (not sure if this was a bug or API change)

change configobj from version 4 to version 5 --> consider using the pip install configobj way of using?

A few raise and print statements that the fixer did not fix

lines in import parsers like: open(parameters_filename, "rb").read() return bytes, need str() around them?

in "block" code need to not encode to bytes in the __str__ call

    def __str__(self):
        return self.__unicode__()  #.encode("utf-8")

Now is where I needed to start re-creating wxGlade code using Python3
- did that, but saw no major impact on performance

Stuck on 'sorting' hlsvd results in Spectral Tab because "numpy does not support '-' subtraction of booleans deprecation warning. Turns out that I was inserting numpy.float64 into the CheckListCtrl from my numpy results. Solution was to cast those to Python float() as they were inserted and numpy issue went away.

Numpy is being pickier. Had code like int(round(util_ppm.ppm2pts())) that returned ndarry from util_ppm call, but the python round() would not work on it. Had to change to np.round() instead

The imp module is deprecated for importlib. Had to find all refs to imp and replace with other calls.
    module = imp.load_source(module_name, filename)
    
        becomes
    
    spec = importlib.util.spec_from_file_location(module_name, filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

string vs byte_array issue in db.py at

    def _convert_boolean(s):
        return (s in ['1', b'1'])
        # Py2 bjs return (s == '1')
        
in Py3 sys.exc_info() resets when it leaves the initial try-except scope, so I needed to be sure to call the Vespa exception processing/translation call right there, and not in a subsequent if-then-else check.  See Puse run_transform_controller._run_transform() call around line 525 for an example.  And similarly in run_experiment_controller() call in Simulation.

In Py3 the struct module (for reading binary data from (?) base64 encoded strings - eg. data from XML files) is a pain in the butt. had to include some strategic .encode() and .decode() processing of byte arrays and strings to make it all happy.  Annoyingly vague in the 2 to 3 docs.

