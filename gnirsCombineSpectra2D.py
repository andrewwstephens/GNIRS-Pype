#!/usr/bin/env python

# MIT License

# Copyright (c) 2015, 2017 Marie Lemoine-Busserolle

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

################################################################################
#                Import some useful Python utilities/modules                   #
################################################################################

import log, os, sys, pkg_resources, glob, ConfigParser, cleanir, datetime, dateutil.parser
from astropy.io import fits
from pyraf import iraf, iraffunctions
import numpy as np

## IMPORTANT NOTE:  The command line options are not available for GNIRS as of July 2019.


def start(configfile):
    """
    This module contains all the functions needed to perform the full reduction of SCIENCE or TELLURIC data.

    Parameters are loaded from gnirs.cfg configuration file. This script will automatically detect if it is being run
    on telluric data or science data. There are 5 steps.

    INPUT FILES:
        - Configuration file
        - Science or Telluric frames
        - mdfshiftrefimage
        - masterflat
        - /database files from the appropriate calibrations directory

    OUTPUT FILES:
        - If telluric:  cleaned (optional), prepared, radiation-event corrected, reduced, spatial distortion corrected, 
          and transformed images
        - If science:  cleaned (optional), prepared, radiation-event corrected, reduced, spatial distortion corrected, 
          and transformed images

    Args:
        - kind (string): Either 'Science' or 'Telluric'
        - configfile: gnirs.cfg configuration file.
                - Paths to the Science (str), reduction truth value (boolean)
                  E.g. 'target/date/config/{Sci,Tel}_ObsID/{Calibrations,Intermediate}', True
                - Paths to the Tellurics (str), reduction truth value (boolean)
                  E.g. 'target/date/config/{Sci,Tel}_ObsID/{Calibrations,Intermediate}', True
                - manualMode (boolean): Enable optional manualModeging pauses? Default: False
                - overwrite (boolean): Overwrite old files? Default: False
                # And gnirsReduce specific settings
    """
    logger = log.getLogger('gnirsCombineSpectra2D.start')

    ###########################################################################
    ##                                                                       ##
    ##                  BEGIN - GENERAL COMBINE 2D SETUP                     ##
    ##                                                                       ##
    ###########################################################################

    # Store current working directory for later use.
    path = os.getcwd()

    logger.info('####################################################')
    logger.info('#                                                  #')
    logger.info('#        Start Combining GNIRS 2D Spectra          #')
    logger.info('#                                                  #')
    logger.info('####################################################')

    # Set up/prepare IRAF.
    iraf.gemini()
    iraf.gemtools()
    iraf.gnirs()

    # Reset to default parameters the used IRAF tasks.
    iraf.unlearn(iraf.gemini,iraf.gemtools,iraf.gnirs,iraf.imcopy)

    # From http://bishop.astro.pomona.edu/Penprase/webdocuments/iraf/beg/beg-image.html:
    # Before doing anything involving image display the environment variable stdimage must be set to the correct frame 
    # buffer size for the display servers (as described in the dev$graphcap file under the section "STDIMAGE devices") 
    # or to the correct image display device. The task GDEVICES is helpful for determining this information for the 
    # display servers.
    iraf.set(stdimage='imt1024')

    # Prepare the IRAF package for GNIRS.
    # NSHEADERS lists the header parameters used by the various tasks in the GNIRS package (excluding headers values 
    # which have values fixed by IRAF or FITS conventions).
    iraf.nsheaders("gnirs",logfile=logger.root.handlers[0].baseFilename)

    # Set clobber to 'yes' for the script. This still does not make the gemini tasks overwrite files, so: YOU WILL 
    # LIKELY HAVE TO REMOVE FILES IF YOU RE_RUN THE SCRIPT.
    us_clobber=iraf.envget("clobber")
    iraf.reset(clobber='yes')
    
    config = ConfigParser.RawConfigParser()
    config.optionxform = str  # make options case-sensitive
    config.read(configfile)
    # Read general config.
    manualMode = config.getboolean('defaults','manualMode')
    overwrite = config.getboolean('defaults','overwrite')
    # config required for combining 2D spectra
    observationSections = ['ScienceDirectories','TelluricDirectories']
    nscombineInter = config.getboolean('interactive','nscombineInter')
    calculateSpectrumSNR = config.getboolean('gnirsPipeline','calculateSpectrumSNR')

    ###########################################################################
    ##                                                                       ##
    ##                 COMPLETE - GENERAL COMBINE 2D SETUP                   ##
    ##                                                                       ##
    ###########################################################################

    # gnirsCombineSpectra2D will first check if the reduction truth value of the science and telluric directories is 
    # True -- if it is, it will then check if the required spectra to be combined are available in the directories (and
    # proceed only if it finds them there); else, it will warn the user and request to provide the files for combining.
    # If the reduction truth value of the science and telluric directories is False, the script will skip combining 2D
    # spectra in those directories.

    # Loop through all the observation (telluric and science) directories to combine 2D spectra in each one.
    for section in observationSections:
        for obspath in config.options(section):
            if config.getboolean(section,obspath):

                ###########################################################################
                ##                                                                       ##
                ##                  BEGIN - OBSERVATION SPECIFIC SETUP                   ##
                ##                                                                       ##
                ###########################################################################

                os.chdir(obspath)
                # Change the iraf directory to the current directory.
                iraffunctions.chdir(obspath)

                # Print the current directory of the 2D spectra being combined.
                logger.info("Currently working on combining 2D spectra in %s\n", obspath)           

                srclistfilename = 'src.list'
                srclist = open(srclistfilename, "r").readlines()
                srclist = [filename.strip() for filename in srclist]
                if calculateSpectrumSNR:
                    # Check if there is a list of sky images in obspath
                    skylistfilename = 'sky.list'
                    if os.path.exists(skylistfilename):
                        skylist = open(skylistfilename, "r").readlines()
                        skylist = [filename.strip() for filename in skylist]
                    else:
                        logger.warning("Parameter 'calculateSpectrumSNR' is 'True', but a list of sky images not ")
                        logger.warning("available in %s . Setting the 'calculateSpectrumSNR' parameter for ", obspath)
                        logger.warning("the current set of observations to 'False'.\n")
                        calculateSpectrumSNR = False
                nodAlistfilename = 'nodA.list'
                nodAlist = open(nodAlistfilename, "r").readlines()
                nodAlist = [filename.strip() for filename in nodAlist]
                nodBlistfilename = 'nodB.list'
                nodBlist = open(nodBlistfilename, "r").readlines()
                nodBlist = [filename.strip() for filename in nodBlist]
                
                # Check if the 2D spectra to combine avalable in the observations directory path
                logger.info("Checking if required 2D spectra available in %s", obspath)
                
                reduce_outputPrefix = 'r'
                srccomblist = sorted(glob.glob('ttf'+reduce_outputPrefix+'lnN*.fits'))
                if not srccomblist:  ## Check for spectral transformation calibration files
                    logger.warning("Required 2D source spectra not available.")
                else:
                    logger.info("Required 2D spectra available.")
                if calculateSpectrumSNR:
                    reduce_outputPrefix = 'k'
                    skycomblist = sorted(glob.glob('ttf'+reduce_outputPrefix+'lnN*.fits'))
                    if not skycomblist:  ## Check for spectral transformation calibration files
                        logger.warning("Parameter 'calculateSpectrumSNR' is 'True', but required 2D sky spectra not ")
                        logger.warning("available. Setting the 'calculateSpectrumSNR' parameter for the current set ")
                        logger.warning("of observations to 'False'.\n")
                        calculateSpectrumSNR = False
                    else:
                        logger.info("Required 2D sky spectra available.\n")
                
                logger.info("Required 2D spectra check complete.")
                
                ###########################################################################
                ##                                                                       ##
                ##                 COMPLETE - OBSERVATION SPECIFIC SETUP                 ##
                ##            BEGIN COMBINING 2D SPECTRA FOR AN OBSERVATION              ##
                ##                                                                       ##
                ###########################################################################

                if manualMode:
                    a = raw_input("About to enter combine 2D spectra.")
                
                reduce_outPrefix = 'r'
                srccombimage = 'src_comb.fits'
                crossCorrelation = 'yes'
                combineSpectra2D(srclistfilename, srclist, nscombineInter, reduce_outputPrefix, srccombimage, \
                    crossCorrelation, overwrite)
                shiftsMatch_flag = shiftsMatch(nodAlist, nodBlist, srccombimage)
                if not shiftsMatch_flag:
                    crossCorrelation = 'no'
                    combineSpectra2D(srclistfilename, srclist, nscombineInter, reduce_outputPrefix, srccombimage, \
                        crossCorrelation, overwrite)
                    shiftsMatch_flag = shiftsMatch(nodAlist, nodBlist, srccombimage)
                    if not shiftsMatch_flag:
                        logger.warning("Shifts calculated by nscombine with fl_cross='%s' not ", crossCorrelation)
                        logger.warning("within the acceptable range. The script already tried using nscombine with ")
                        logger.warning("fl_cross='yes'. This looks like the nscombine bug issue.\n")
                    else:
                        logger.info("Shifts calculated by nscombine with fl_cross='%s' within the ", crossCorrelation)
                        logger.info("acceptable range.\n")
                else:
                    logger.info("Combining 2D spectra using nscombine with fl_cross='%s' worked.\n", crossCorrelation)
                    
                # If the parameter 'calculateSpectrumSNR' is set to 'yes', the script will combine all observations 
                # reduced without sky subtraction (having reduce_outputPrefix 'k'); else, all observations reduced 
                # with sky subtraction only (having reduce_outputPrefix 'r') will be combined.
                if calculateSpectrumSNR:
                    logger.info("Combining the sky observations reduced without sky subtraction.\n")
                    reduce_outputPrefix = 'k'
                    skycombimage = 'sky_comb.fits'
                    crossCorrelation = 'no'
                    combineSpectra2D(skylistfilesname, skylist, nscombineInter, reduce_outputPrefix, skycombimage, \
                        crossCorrelation, overwrite)
            
            logger.info("##############################################################################")
            logger.info("#                                                                            #")
            logger.info("#  COMPLETE - Combining 2D spectra completed for                             #")
            logger.info("#  %s", obspath                                                                )
            logger.info("#                                                                            #")
            logger.info("##############################################################################\n")

    # Return to directory script was begun from.
    os.chdir(path)

    return

##################################################################################################################
#                                                     ROUTINES                                                   #
##################################################################################################################

def combineSpectra2D(combinlistfilename, combinlist, nscombineInter, reduce_outputPrefix, combinedimage, crossCorrelation, overwrite):
    """
    Combining the transformed science or telluric frames.
    """
    logger = log.getLogger('gnirsReduce.combineSpectra2D')

    if os.path.exists(combinedimage):
        if overwrite:
            logger.warning("Removing old %s\n", combinedimage)
            os.remove(combinedimage)
            iraf.nscombine(inimages='ttf'+reduce_outputPrefix+'ln//@'+combinlistfilename, tolerance=0.5, \
                output=combinedimage, output_suffix='', bpm="", dispaxis=1, pixscale=1., fl_cross=crossCorrelation, \
                fl_keepshift='no', fl_shiftint='yes', interptype="linear", boundary="nearest", constant=0., \
                combtype="average", rejtype="none", masktype="none", maskvalue=0., statsec="[*,*]", scale="none", \
                zero="none", weight="none", lthreshold="INDEF", hthreshold="INDEF", nlow=1, nhigh=1, nkeep=0, \
                mclip='yes', lsigma=5., hsigma=5., ron=0.0, gain=1.0, snoise="0.0", sigscale=0.1, pclip=-0.5, \
                grow=0.0, nrejfile='', fl_vardq='yes', fl_inter=nscombineInter, \
                logfile=logger.root.handlers[0].baseFilename, verbose='yes', debug='no', force='no')
        else:
            logger.warning("Old %s exists and -overwrite not set - skipping nscombine for ", combinedimage)
            logger.warning("observations in the input list %s.\n", combinlistfilename)
    else:
        iraf.nscombine(inimages='ttf'+reduce_outputPrefix+'ln//@'+combinlistfilename, tolerance=0.5, \
            output=combinedimage, output_suffix='', bpm="", dispaxis=1, pixscale=1., fl_cross=crossCorrelation, \
            fl_keepshift='no', fl_shiftint='yes', interptype="linear", boundary="nearest", constant=0., \
            combtype="average", rejtype="none", masktype="none", maskvalue=0., statsec="[*,*]", scale="none", \
            zero="none", weight="none", lthreshold="INDEF", hthreshold="INDEF", nlow=1, nhigh=1, nkeep=0, mclip='yes',\
            lsigma=5., hsigma=5., ron=0.0, gain=1.0, snoise="0.0", sigscale=0.1, pclip=-0.5, grow=0.0, nrejfile='', \
            fl_vardq='yes', fl_inter=nscombineInter, logfile=logger.root.handlers[0].baseFilename, verbose='yes', \
            debug='no', force='no')

#---------------------------------------------------------------------------------------------------------------------#

def shiftsMatch(nodAlist, nodBlist, combinedimage):
    """
    There is currently (2015) a bug in nscombine that miscalculates the shifts at certain position angles while 
    cross-correlating different frames. Unfortunately, we cannot provide a shift to nscombine directly. 
    
    If targets are bright enough, the task avoids the bug and the cross-correlation happens smoothly (set with 
    fl_cross='yes'); else, the task should be run with fl_cross='no'. XDGNIRS always uses fl_cross='yes' for telluric
    frames. For science frames, XDGNIRS checks whether the shifts used by nscombine with fl_cross='yes' were 
    reasonable; else, it runs nscombine with fl_cross = 'no'.
    """
    logger = log.getLogger('gnirsReduce.matchShifts')

    # TODO(Viraja):  Currently, the nod Q offsets are only calculated for the first set of images in the different nod
    # lists. This assumes that the Q offsets for the remaining set of images in the lists give the same nod value. A 
    # robust way of doing this will be calculating the shifts between every two consecutive image sets from the nod 
    # lists and checking if the shift is within the acceptable range; and warn the user if it is not.
    nodAheader = fits.open(nodAlist[0])[0].header
    nodBheader = fits.open(nodBlist[0])[0].header
    nodQoffset = abs(nodAheader['QOFFSET'] - nodBheader['QOFFSET'])
    if 'Long' in nodAheader['CAMERA']:
        pixelscale = 0.05
    else:
        pixelscale = 0.15
    nodQpixels = nodQoffset/pixelscale

    shifts = iraf.hselect(images=combinedimage+'[0]', fields='NSCHLX*', expr='yes', missing='INDEF', Stdout=1)
    shifts = [shift.replace('\t', ',') for shift in shifts]

    shiftsMatch_flag = True
    for shift in shifts[0].split(','):
        if shift != '0.00' and abs(float(shift) - nodQpixels) > 1.5:
            logger.warning("Shift of %s pixels in %s (NSCHLX* keyword) does not ", shift, combinedimage)
            logger.warning("match the expected value of %s pixels. Likely, the target was not bright ", nodQpixels)
            logger.warning("enough to use nscombine with fl_cross='yes'. Rerunning nscombine with fl_cross='no'.\n")
            shiftsMatch_flag = shiftsMatch_flag and False
            break
        else:
            logger.info("Shift of %s pixels in %s (NSCHLX* keyword) is 0.0 or close ", shift, combinedimage)
            logger.info("to the expected value of %s pixels.", nodQpixels)
            shiftsMatch_flag = shiftsMatch_flag and True
            continue
    
    return shiftsMatch_flag

#---------------------------------------------------------------------------------------------------------------------#

if __name__ == '__main__':
    log.configure('gnirs.log', filelevel='INFO', screenlevel='DEBUG')
    start('gnirs.cfg')
