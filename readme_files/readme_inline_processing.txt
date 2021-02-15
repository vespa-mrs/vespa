
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
	
