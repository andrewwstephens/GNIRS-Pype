Introduction
============

GNIRS-Pype, for now, uses Python 2.7. Please keep this in mind.

Running from the Command Line
=============================

GBIRS-Pype, for now, is started simply by running the gnirsPipeline script from the command line.

Running from Python
===================

Python Programmers: even though GNIRS-Pype's API has not been defined yet, you can still run GNIRS-Pype Pipeline, Steps, Routines, and Tasks 
from a Python interpeter by importing them.

For example, to run the combine 2D spectra from the Python interpreter:

.. code-block:: python

  This is currently under development.

You might be able to figure out what imports you need by checking the tops of the relevant scripts.

Examples of Running from the Command Line
-----------------------------------------

Starting a Data Reduction from the Beginning


This is currently under development.

Starting a Data Reduction from a Specified Point


You can run each step, one at a time, from the command line like so. *You need a gnirs.cfg file in your current working directory to run individual steps.
Each step requires the general config section and its unique and/or additional config section to be populated*.

You can also run an individual step by turning them on or off in gnirsPipeline config and running the gnirsPipeline.

**gnirsSort:** To only sort and copy the GNIRS raw data or create symbolic links in between the qorking directories use a gnirs.cfg file like this:

.. code-block:: text

  # GNIRS-Pype configuration file.

  [defaults]
  gnirs-pipeVersion = 1.0.0
  manualMode = False
  overwrite = False

  [ScienceDirectories]

  [TelluricDirectories]

  [CalibrationDirectories]

  [gnirsPipeline]
  getData = False
  sort = True
  checkData = False
  calibrationReduction = False
  scienceReduction = False
  telluricReduction = False
  combineSpectra2D = False
  extractSpectra1D = False
  telluricCorrection = False
  fluxCalibration = False
  combineOrdersXD = False
  calculateSpectrumSNR = False
  writeDataSheet = False

  [interactive]
  nsprepareInter = False
  nsflatInter = False
  nscombineInter = False
  nssdistInter = False
  nswavelengthInter = False
  nsfitcoordsInter = False
  nsextractInter = False
  hLineInter = False
  continuumInter = False
  telluricInter = False
  tempInter = False

  [getData]
  rawPath = rawData/
  program = 

  [sort]
  proprietaryCookie = 
  telluricTimeThreshold = 5400.0  ; units [seconds]

  [checkData]
  # No input parameters required

  [calibrationReduction]
  Start = 1
  Stop = 5
  cleanir_IRflats = False
  cleanir_QHflats = False
  cleanir_arcs = False
  cleanir_pinholes = False

  [scienceReduction]
  Start = 1
  Stop = 5
  cleanir = False
  radiationCorrectionMethod = fixpix
  radiationThreshold = 20
  checkPeaksMatch = True

  [telluricReduction]
  Start = 1
  Stop = 5
  cleanir = False
  radiationCorrectionMethod = fixpix
  radiationThreshold = 50

  [combineSpectra2D]
  # No input parameters required

  [extractSpectra1D]
  useApall = True
  extractionApertureRadius = 4
  checkPeaksMatch = True
  toleranceOffset = 5
  extractionFullSlit = False
  extractionStepwise = False
  extractionStepSize = 6  ; units [pixels]

  [telluricCorrection]
  Start = 1
  Stop = 5
  hLineMethod = vega

  [fluxCalibration]
  Start = 1
  Stop = 6
  telluricRA = 
  telluricDEC = 
  telluricSpectralType =  
  telluricMagnitude = 
  telluricTemperature = 

  [combineOrdersXD]
  # This section is currently under development.
  Start = 1
  Stop = 1

  [calculateSpectrumSNR]
  # This section is currently under development.
  Start = 1
  Stop = 1

  [writeDataSheet]
  # This section is currently under development.
  Start = 1
  Stop = 1

  [TelluricRegions]  ; regions to use when tweaking the Telluric line removal
  3 = *
  4 = 14200:18000
  5 = 11200:13400,14200:15200
  6 = 11000:12000
  7 = 9000:10000
  8 = 8500:9500

  [HLineRegions]  ; regions to use when tweaking the H line removal
  3 = 21400:21900
  4 = 15500:18000
  5 = 12750:12900
  6 = 10700:11200
  7 = 9000:10500
  8 = 8500:9300

  [ContinuumRegions]  ; regions to use when fitting the Telluric continuum
  3 = *
  4 = *
  5 = 11200:13300,14300:15200
  6 = 9500:10900,11700:12600
  7 = 8500:9200,9800:10700
  8 = 7500:9500

And run the gnirsPipeline with:

.. code-block:: text

 This section is currently under development.

**gnirsBaselineCalibration:** To only reduce calibrations use a gnirs.cfg file like this:

.. code-block:: text

  # GNIRS-Pype configuration file.

  [defaults]
  gnirs-pipeVersion = 1.0.0
  manualMode = False
  overwrite = False

  [ScienceDirectories]

  [TelluricDirectories]

  [CalibrationDirectories]

  [gnirsPipeline]
  getData = False
  sort = False
  checkData = False
  calibrationReduction = True
  scienceReduction = False
  telluricReduction = False
  combineSpectra2D = False
  extractSpectra1D = False
  telluricCorrection = False
  fluxCalibration = False
  combineOrdersXD = False
  calculateSpectrumSNR = False
  writeDataSheet = False

  [interactive]
  nsprepareInter = False
  nsflatInter = False
  nscombineInter = False
  nssdistInter = False
  nswavelengthInter = False
  nsfitcoordsInter = False
  nsextractInter = False
  hLineInter = False
  continuumInter = False
  telluricInter = False
  tempInter = False

  [getData]
  rawPath = rawData/
  program = 

  [sort]
  proprietaryCookie = 
  telluricTimeThreshold = 5400.0  ; units [seconds]

  [checkData]
  # No input parameters required

  [calibrationReduction]
  Start = 1
  Stop = 5
  cleanir_IRflats = False
  cleanir_QHflats = False
  cleanir_arcs = False
  cleanir_pinholes = False

  [scienceReduction]
  Start = 1
  Stop = 5
  cleanir = False
  radiationCorrectionMethod = fixpix
  radiationThreshold = 20
  checkPeaksMatch = True

  [telluricReduction]
  Start = 1
  Stop = 5
  cleanir = False
  radiationCorrectionMethod = fixpix
  radiationThreshold = 50

  [combineSpectra2D]
  # No input parameters required

  [extractSpectra1D]
  useApall = True
  extractionApertureRadius = 4
  checkPeaksMatch = True
  toleranceOffset = 5
  extractionFullSlit = False
  extractionStepwise = False
  extractionStepSize = 6  ; units [pixels]

  [telluricCorrection]
  Start = 1
  Stop = 5
  hLineMethod = vega

  [fluxCalibration]
  Start = 1
  Stop = 6
  telluricRA = 
  telluricDEC = 
  telluricSpectralType =  
  telluricMagnitude = 
  telluricTemperature = 

  [combineOrdersXD]
  # This section is currently under development.
  Start = 1
  Stop = 1

  [calculateSpectrumSNR]
  # This section is currently under development.
  Start = 1
  Stop = 1

  [writeDataSheet]
  # This section is currently under development.
  Start = 1
  Stop = 1

  [TelluricRegions]  ; regions to use when tweaking the Telluric line removal
  3 = *
  4 = 14200:18000
  5 = 11200:13400,14200:15200
  6 = 11000:12000
  7 = 9000:10000
  8 = 8500:9500

  [HLineRegions]  ; regions to use when tweaking the H line removal
  3 = 21400:21900
  4 = 15500:18000
  5 = 12750:12900
  6 = 10700:11200
  7 = 9000:10500
  8 = 8500:9300

  [ContinuumRegions]  ; regions to use when fitting the Telluric continuum
  3 = *
  4 = *
  5 = 11200:13300,14300:15200
  6 = 9500:10900,11700:12600
  7 = 8500:9200,9800:10700
  8 = 7500:9500

And run the gnirsPipeline with:

.. code-block:: text

 This section is currently under development.

**gnirsReduce Science:** To only reduce science data use a config.cfg file like this:
*Make sure to populate ScienceDirectories, TelluricDirectories and CalibrationDirectories before running!*

.. code-block:: text

  # GNIRS-Pype configuration file.

  [defaults]
  gnirs-pipeVersion = 1.0.0
  manualMode = False
  overwrite = False

  [ScienceDirectories]

  [TelluricDirectories]

  [CalibrationDirectories]

  [gnirsPipeline]
  getData = False
  sort = False
  checkData = False
  calibrationReduction = False
  scienceReduction = True
  telluricReduction = False
  combineSpectra2D = False
  extractSpectra1D = False
  telluricCorrection = False
  fluxCalibration = False
  combineOrdersXD = False
  calculateSpectrumSNR = False
  writeDataSheet = False

  [interactive]
  nsprepareInter = False
  nsflatInter = False
  nscombineInter = False
  nssdistInter = False
  nswavelengthInter = False
  nsfitcoordsInter = False
  nsextractInter = False
  hLineInter = False
  continuumInter = False
  telluricInter = False
  tempInter = False

  [getData]
  rawPath = rawData/
  program = 

  [sort]
  proprietaryCookie = 
  telluricTimeThreshold = 5400.0  ; units [seconds]

  [checkData]
  # No input parameters required

  [calibrationReduction]
  Start = 1
  Stop = 5
  cleanir_IRflats = False
  cleanir_QHflats = False
  cleanir_arcs = False
  cleanir_pinholes = False

  [scienceReduction]
  Start = 1
  Stop = 5
  cleanir = False
  radiationCorrectionMethod = fixpix
  radiationThreshold = 20
  checkPeaksMatch = True

  [telluricReduction]
  Start = 1
  Stop = 5
  cleanir = False
  radiationCorrectionMethod = fixpix
  radiationThreshold = 50

  [combineSpectra2D]
  # No input parameters required

  [extractSpectra1D]
  useApall = True
  extractionApertureRadius = 4
  checkPeaksMatch = True
  toleranceOffset = 5
  extractionFullSlit = False
  extractionStepwise = False
  extractionStepSize = 6  ; units [pixels]

  [telluricCorrection]
  Start = 1
  Stop = 5
  hLineMethod = vega

  [fluxCalibration]
  Start = 1
  Stop = 6
  telluricRA = 
  telluricDEC = 
  telluricSpectralType =  
  telluricMagnitude = 
  telluricTemperature = 

  [combineOrdersXD]
  # This section is currently under development.
  Start = 1
  Stop = 1

  [calculateSpectrumSNR]
  # This section is currently under development.
  Start = 1
  Stop = 1

  [writeDataSheet]
  # This section is currently under development.
  Start = 1
  Stop = 1

  [TelluricRegions]  ; regions to use when tweaking the Telluric line removal
  3 = *
  4 = 14200:18000
  5 = 11200:13400,14200:15200
  6 = 11000:12000
  7 = 9000:10000
  8 = 8500:9500

  [HLineRegions]  ; regions to use when tweaking the H line removal
  3 = 21400:21900
  4 = 15500:18000
  5 = 12750:12900
  6 = 10700:11200
  7 = 9000:10500
  8 = 8500:9300

  [ContinuumRegions]  ; regions to use when fitting the Telluric continuum
  3 = *
  4 = *
  5 = 11200:13300,14300:15200
  6 = 9500:10900,11700:12600
  7 = 8500:9200,9800:10700
  8 = 7500:9500

And run the gnirsPipeline with:

.. code-block:: text

 This is currently under development.

**gnirsFluxCalibrate:** To only do a flux calibration use a gnirs.cfg file like this:
*Make sure to populate scienceDirectoryList, telluricDirectoryList and calibrationDirectoryList before running!*

.. code-block:: text

  # GNIRS-Pype configuration file.

  [defaults]
  gnirs-pipeVersion = 1.0.0
  manualMode = False
  overwrite = False

  [ScienceDirectories]

  [TelluricDirectories]

  [CalibrationDirectories]

  [gnirsPipeline]
  getData = False
  sort = False
  checkData = False
  calibrationReduction = False
  scienceReduction = False
  telluricReduction = False
  combineSpectra2D = False
  extractSpectra1D = False
  telluricCorrection = False
  fluxCalibration = True
  combineOrdersXD = False
  calculateSpectrumSNR = False
  writeDataSheet = False

  [interactive]
  nsprepareInter = False
  nsflatInter = False
  nscombineInter = False
  nssdistInter = False
  nswavelengthInter = False
  nsfitcoordsInter = False
  nsextractInter = False
  hLineInter = False
  continuumInter = False
  telluricInter = False
  tempInter = False

  [getData]
  rawPath = rawData/
  program = 

  [sort]
  proprietaryCookie = 
  telluricTimeThreshold = 5400.0  ; units [seconds]

  [checkData]
  # No input parameters required

  [calibrationReduction]
  Start = 1
  Stop = 5
  cleanir_IRflats = False
  cleanir_QHflats = False
  cleanir_arcs = False
  cleanir_pinholes = False

  [scienceReduction]
  Start = 1
  Stop = 5
  cleanir = False
  radiationCorrectionMethod = fixpix
  radiationThreshold = 20
  checkPeaksMatch = True

  [telluricReduction]
  Start = 1
  Stop = 5
  cleanir = False
  radiationCorrectionMethod = fixpix
  radiationThreshold = 50

  [combineSpectra2D]
  # No input parameters required

  [extractSpectra1D]
  useApall = True
  extractionApertureRadius = 4
  checkPeaksMatch = True
  toleranceOffset = 5
  extractionFullSlit = False
  extractionStepwise = False
  extractionStepSize = 6  ; units [pixels]

  [telluricCorrection]
  Start = 1
  Stop = 5
  hLineMethod = vega

  [fluxCalibration]
  Start = 1
  Stop = 6
  telluricRA = 
  telluricDEC = 
  telluricSpectralType =  
  telluricMagnitude = 
  telluricTemperature = 

  [combineOrdersXD]
  # This section is currently under development.
  Start = 1
  Stop = 1

  [calculateSpectrumSNR]
  # This section is currently under development.
  Start = 1
  Stop = 1

  [writeDataSheet]
  # This section is currently under development.
  Start = 1
  Stop = 1

  [TelluricRegions]  ; regions to use when tweaking the Telluric line removal
  3 = *
  4 = 14200:18000
  5 = 11200:13400,14200:15200
  6 = 11000:12000
  7 = 9000:10000
  8 = 8500:9500

  [HLineRegions]  ; regions to use when tweaking the H line removal
  3 = 21400:21900
  4 = 15500:18000
  5 = 12750:12900
  6 = 10700:11200
  7 = 9000:10500
  8 = 8500:9300

  [ContinuumRegions]  ; regions to use when fitting the Telluric continuum
  3 = *
  4 = *
  5 = 11200:13300,14300:15200
  6 = 9500:10900,11700:12600
  7 = 8500:9200,9800:10700
  8 = 7500:9500

And run the gnirsPipeline with:

.. code-block:: text

 This is currently under development.

**gnirsCombineOrdersXD Orders Combining:** To only combine sifferent spectral orders use a gnirs.cfg file like this:
*Make sure to populate scienceDirectoryList, telluricDirectoryList and calibrationDirectoryList before running!*

.. code-block:: text

  # GNIRS-Pype configuration file.

  [defaults]
  gnirs-pipeVersion = 1.0.0
  manualMode = False
  overwrite = False

  [ScienceDirectories]

  [TelluricDirectories]

  [CalibrationDirectories]

  [gnirsPipeline]
  getData = False
  sort = False
  checkData = False
  calibrationReduction = False
  scienceReduction = False
  telluricReduction = False
  combineSpectra2D = False
  extractSpectra1D = False
  telluricCorrection = False
  fluxCalibration = False
  combineOrdersXD = True
  calculateSpectrumSNR = False
  writeDataSheet = False

  [interactive]
  nsprepareInter = False
  nsflatInter = False
  nscombineInter = False
  nssdistInter = False
  nswavelengthInter = False
  nsfitcoordsInter = False
  nsextractInter = False
  hLineInter = False
  continuumInter = False
  telluricInter = False
  tempInter = False

  [getData]
  rawPath = rawData/
  program = 

  [sort]
  proprietaryCookie = 
  telluricTimeThreshold = 5400.0  ; units [seconds]

  [checkData]
  # No input parameters required

  [calibrationReduction]
  Start = 1
  Stop = 5
  cleanir_IRflats = False
  cleanir_QHflats = False
  cleanir_arcs = False
  cleanir_pinholes = False

  [scienceReduction]
  Start = 1
  Stop = 5
  cleanir = False
  radiationCorrectionMethod = fixpix
  radiationThreshold = 20
  checkPeaksMatch = True

  [telluricReduction]
  Start = 1
  Stop = 5
  cleanir = False
  radiationCorrectionMethod = fixpix
  radiationThreshold = 50

  [combineSpectra2D]
  # No input parameters required

  [extractSpectra1D]
  useApall = True
  extractionApertureRadius = 4
  checkPeaksMatch = True
  toleranceOffset = 5
  extractionFullSlit = False
  extractionStepwise = False
  extractionStepSize = 6  ; units [pixels]

  [telluricCorrection]
  Start = 1
  Stop = 5
  hLineMethod = vega

  [fluxCalibration]
  Start = 1
  Stop = 6
  telluricRA = 
  telluricDEC = 
  telluricSpectralType =  
  telluricMagnitude = 
  telluricTemperature = 

  [combineOrdersXD]
  # This section is currently under development.
  Start = 1
  Stop = 1

  [calculateSpectrumSNR]
  # This section is currently under development.
  Start = 1
  Stop = 1

  [writeDataSheet]
  # This section is currently under development.
  Start = 1
  Stop = 1

  [TelluricRegions]  ; regions to use when tweaking the Telluric line removal
  3 = *
  4 = 14200:18000
  5 = 11200:13400,14200:15200
  6 = 11000:12000
  7 = 9000:10000
  8 = 8500:9500

  [HLineRegions]  ; regions to use when tweaking the H line removal
  3 = 21400:21900
  4 = 15500:18000
  5 = 12750:12900
  6 = 10700:11200
  7 = 9000:10500
  8 = 8500:9300

  [ContinuumRegions]  ; regions to use when fitting the Telluric continuum
  3 = *
  4 = *
  5 = 11200:13300,14300:15200
  6 = 9500:10900,11700:12600
  7 = 8500:9200,9800:10700
  8 = 7500:9500

And run the gnirsPipeline with:

.. code-block:: text

 This is currently under development.

Preparing the .cfg Input File
=============================

GNIRS-Pype reads data reduction parameters with the Python 2.7 built-in ConfigParser. See
https://docs.python.org/2.7/library/configparser.html for full documentation on the parser.

Interactive Input Preparation
-----------------------------

This is currently under development.

An Example Input File
---------------------

GNIRS-Pype includes a default configuration file. As of v1.0.1, it looks like this:

.. code-block:: text

  TODO(Viraja): Updated example

Data Reduction Examples
=======================

Observations of NGC 1736; GN-2011A-Q-126
----------------------------------------

This is currently under development.

Example 2
---------

This is currently under development.

Example 3
---------

This is currently under development.

Tutorials
=========

H Line Removal
--------------

This is currently under development.

Custom Telluric Corrections
---------------------------

This is currently under development.

Known Issues
============

gnirsPipeline.py
-----------------

gnirsGetData.py
-----------------

gnirsSort.py
-----------

gnirsCheckData.py
-----------------

gnirsBaselineCalibration.py
--------------------------

gnirsReduce.py
-------------

gnirsCombineSpectra2D.py
-------------

gnirsExtractSpectra1D.py
-------------

gnirsTelluric.py
-------------

gnirsFluxCalibrate.py
-------------

gnirsCombineOrdersXD.py
------------

gnirsCalculateSpectramSNR.py
------------

gnirsWriteDataSheet.py
------------

General Issues
--------------

- A longstanding bug (see `astropy <https://github.com/astropy/astropy/pull/960>`_ ) in astropy has made it
  difficult to build GNIRS-Pype as a binary executable.

Maintaining GNIRS-Pype
========================

Documentation
-------------

Right now there exists four forms of documentation.

Paper

.. Insert a paper!

README.rst


.rst Files in the docs/ directory


This file, others like it in the docs/ directory and the README are written in
reStructuredText. This markup language integrates well with Python's automatic
documentation builder (we used Sphinx) and Github as well as being human readable. You can
read more about reStructuredText `here <http://www.sphinx-doc.org/en/stable/rest.html>`_.

Comments and DocStrings in Source Code

Tests
-----

This is a TODO. Currently we do not have automated tests.

Pipeline Structure
------------------

This is currently under development.

GNIRS-Pype is built at the lowest level from Python and IRAF subroutines. It is built so that it is relatively
easy to change the implementation of the underlying tasks.

Updates
-------

To update GNIRS-Pype, do ...

This is currently under development.

Version Numbers


GNIRS-Pype uses semantic versioning(see http://semver.org/). This means version numbers come in

.. code-block:: text

  MAJOR.MINOR.PATCH

In brief, when releasing a version of Nifty that is not backward-compatible with old test recipes,
or changes break the public API, it is time to increment the MAJOR version number.

Code Conventions
----------------

For naming variables and functions, a mix of camelCase and lower_case_with_underscores was used.

Code style was influenced by the `Google Python Style Guide <https://google.github.io/styleguide/pyguide.html>`_.

GNIRS-Pype uses the Google docstring style. Examples of docstrings can be found
`here <http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html>`_.

Future Work
===========

Throughout the code, vkhatu has placed many TODO notes. These are things that
should be reviewed at some point.

Future work:

- Implement more telluric and flux calibration methods.
- Implement instrument signature removal routine.
- Implement differential atmospheric refraction correction routine.
- Implement full automatic Gemini (on server) reduction routine.
- Python 3 compatability(if possible)
- Compiling as a self-contained executable
- Full AstroConda integration

Changelog
=========
All notable changes to this project will be documented in this file.

The format is based on `Keep a Changelog <http://keepachangelog.com/en/1.0.0/>`_
and this project adheres to `Semantic Versioning <http://semver.org/spec/v2.0.0.html>`_.

Unreleased
----------
All in-development changes will be tracked here.

- Adding unit tests for each step and integration test for pipeline.

1.0.0 - 2019-08-10
------------------
This is currently under development.

API
===

Note: I didn't have time to implement this using Sphinx automodule. GNIRS-Pype has fairly good docstrings and you
can use individual steps, routines and tasks by importing them. This is a TODO.

.. TODO: implement the public API.
