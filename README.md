# vespa
Python tools for Magnetic Resonance Spectroscopy - RF Pulses, Spectral Simulation and Data Analysis

- Contact: Brian J. Soher, vespa@briansoher.com
- Licence: BSD, specifically a "three-clause" BSD license
- For install, usage, and technical information see Github Pages: http://github.com/vespa-mrs/vespa/

## Description

Vespa stands for Versatile Simulation, Pulses, and Analysis. It is an integrated, open source, open development platform that contains four magnetic resonance spectroscopy (MRS) software applications, written in Python, called: 
	
1. Pulse      - RF pulse design
2. Simulation - spectral simulation and prototyping
3. Priorset   - application for creating 'fake' MRS data sets 
4. Analysis   - spectral data processing and analysis
	
The Vespa project addresses previous software limitations, of non-standard data access, closed source and multiple language software that complicates algorithm extension, and a lack of integration between programs by porting all the code to Python, and using a shared database design, with high-quality documentation. 

These applications can be run separately, but can also communicate via a shared database of objects/results. Integration allows one application to use the output from another as input. For example, Simulation can make use of an RF pulse designed in Pulse to create a more realistic MR simulation.

## Background

The Vespa package is an extensive redesign of three previous MRS software tools:

- MatPulse   (Pulse) - software for RF pulse design written in Matlab, 
- GAVA/Gamma (Simulation) - software for spectral simulation code in IDL
- IDL_Vespa  (Analysis) - spectral data processing and analysis code in IDL 

Thanks to the NIH (grant number 1R01EB008387-01A1) for funding the maintenance and extension of these separate applications into a combined environment based entirely on the Python language.

## Compatibility

Vespa has been tested and is certified to run on the following systems: Windows, Macintosh OSX, and Linux. However, it should run on any system that supports python and wxpython.
