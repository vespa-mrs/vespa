# Vespa

**Versatile Simulation Pulses and Analysis** - Python tools for Magnetic Resonance Spectroscopy

- Contact: Brian J. Soher, vespa@briansoher.com
- Bug Reports: vespa.bugs@gmail.com 
- Licence: BSD, specifically a "three-clause" BSD license
- Vespa can be installed via **PyPI** (see the Vespa-Suite package). 
- This repository is primarity for contributing towards development.
- For more Installation, User Manuals, and Technical Information see Github Pages: 
  
  https://vespa-mrs.github.io/vespa.io/

## Description

Vespa is an integrated, open source, open development platform that contains four magnetic resonance spectroscopy (MRS) software applications, written in Python, called: 
	
1. Pulse      - for RF pulse design
2. Simulation - allows spectral simulation and prototyping
3. DataSim    - an application for creating 'fake' MRS data sets 
4. Analysis   - interactive spectral data processing and analysis

These applications can be run separately, but can also communicate via a shared database of objects/results. Integration allows one application to use the output from another as input. For example, Simulation can make use of an RF pulse designed in Pulse to create a more realistic MR simulation.

The Vespa project addresses previous software limitations such as: non-standard data access, closed source and multiple language software (that complicate algorithm extension), and a lack of integration between programs. 

## Background

The Vespa package is an extensive redesign of three previous MRS software tools:

- MatPulse   (Pulse) - software for RF pulse design written in Matlab, 
- GAVA/Gamma (Simulation) - software for spectral simulation code in IDL
- IDL_Vespa  (Analysis) - spectral data processing and analysis code in IDL 

Thanks to the NIH (grant number 1R01EB008387-01A1) for funding the maintenance and extension of these separate applications into a combined environment based entirely on the Python language.

## Compatibility

Vespa has been tested and is certified to run on the following systems: Windows, Macintosh OSX, and Linux. However, it should run on any system that supports Python and wxpython. As of version 1.0.0, Vespa now runs under Python 3.7 and later.
