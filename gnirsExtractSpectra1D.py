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

import log, os, sys, glob, ConfigParser, datetime, dateutil.parser
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
    logger = log.getLogger('gnirsExtractSpectra1D.start')

    ###########################################################################
    ##                                                                       ##
    ##                  BEGIN - GENERAL EXTRACT 1D SETUP                     ##
    ##                                                                       ##
    ###########################################################################

    # Store current working directory for later use.
    path = os.getcwd()

    logger.info('####################################################')
    logger.info('#                                                  #')
    logger.info('#        Start Extracting GNIRS 1D Spectra         #')
    logger.info('#                                                  #')
    logger.info('####################################################\n')

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
    # config required for extracting 1D spectra
    observationSections = ['TelluricDirectories','ScienceDirectories']
    # Above order of sections is important to later check for plausible peaks located for science targets by nsextract
    nsextractInter = config.getboolean('interactive','nsextractInter')
    calculateSpectrumSNR = config.getboolean('gnirsPipeline','calculateSpectrumSNR')
    # extract1Spectra1D specific config
    useApall = config.getboolean('extractSpectra1D','useApall')
    extractionApertureRadius = config.getfloat('extractSpectra1D','extractionApertureRadius')
    toleranceOffset = config.getfloat('extractSpectra1D','toleranceOffset')
    extractStepwise = config.getboolean('extractSpectra1D','extractStepwise')
    extractionStepsize = config.getfloat('extractSpectra1D','extractionStepsize')

    ###########################################################################
    ##                                                                       ##
    ##                 COMPLETE - GENERAL EXTRACT 1D SETUP                   ##
    ##                                                                       ##
    ###########################################################################

    # gnirsExtractSpectra1D will first check if the reduction truth value of the science and telluric directories is 
    # True -- if it is, it will then check if the required spectra to be extracted are available in the directories
    # (and proceed only if it finds them there); else, it will warn the user and request to provide the spectra for 
    # extracting. If the reduction truth value of the science and telluric directories is False, the script will skip
    # extracting 1D spectra in those directories.

    # Loop through all the observation (telluric and science) directories to extract 1D spectra in each one.
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

                # Print the current directory of observations being reduced.
                logger.info("Currently working on extracting 1D spectra in %s\n", obspath)
                
                # Check if required combined spectra available in the observations directory path
                logger.info("Checking if required combined spectra available in %s", obspath)

                srccombimage = 'src_comb.fits'
                if os.path.exists(srccombimage):
                    logger.info("Required combined source image available.")
                else:
                    logger.warning("Required combined source image not available. Please run ")
                    logger.warning("gnirsCombineSpectra2D.pyto create the combined source image or provide it ")
                    logger.warning("manually. Exiting script.\n")
                    raise SystemExit
                if calculateSpectrumSNR:
                    skycombimage = obspath+'/sky_comb.fits'
                    if os.path.exists(srccombimage):
                        logger.info("Required combined sky image available.")
                    else:
                        logger.warning("Parameter 'calculateSpectrumSNR' is 'True', but required combined sky image ")
                        logger.warning("not available. Setting the 'calculateSpectrumSNR' parameter for the current ")
                        logger.warning("set of observations to 'False'.\n")
                        calculateSpectrumSNR = False
                
                logger.info("Required combined spectra check complete.")
                
                # TODO(Viraja)?:  If 'Science' in section, check if /Tel_* directory exists in the observations path
                if 'Science' in section:
                    # Check if nsextractInter is set: if no, check if checkPeaksMatch is set: if yes, check if the
                    # required telluric extraction reference files available in the telluric /database directory; else, 
                    # warn the user that both nsextractInter and checkPeaksMatch are not set, request the user to 
                    # manually check if the science target peak identified by task nsextract might 
                    # identify a wrong peak if the science target is not bright enough. 
                    if nsextractInter:
                        pass
                    else:  ## conditions if nsextract is not run interactively
                        if checkPeaksMatch:
                            # Check if the paths to the telluric /database direcory and the telluric extraction reference 
                            # files therein exist
                            logger.info("Checking if the required telluric extraction reference files available.")

                            # TODO(Viraja):  Check with Andy if there is a better way of getting the absolute path to
                            # the telluric directory
                            telpath = (glob.glob(obspath+'/Tel_*')).pop()
                            # Print the telluric directory path
                            logger.info("The path to the telluric directory is %s\n", telpath)
                            # Check for the telluric /database directory in the current science observation directory
                            teldatabasepath = telpath+'/database'
                            if os.path.exists(teldatabasepath):  
                                logger.info("Telluric /database directory (possibly) containing telluric extraction ")
                                logger.info("reference files available.")
                                telCheck_flag = True
                            else:
                                logger.warning("Telluric /database directory (possibly) containing telluric ")
                                logger.warning("extraction reference files not available.")
                                telCheck_flag = False
                            # Check for telluric extraction aperture reference files
                            telapfiles = glob.glob(teldatabasepath+'/ap*')
                            if not telapfiles:
                                logger.warning("Reference files containing telluric extraction aperture details not ")
                                logger.warning("available in the telluric /database directory.")
                                telCheck_flag = telCheck_flag and False
                            else:
                                logger.info("Reference files containing telluric extraction aperture details ")
                                logger.info("available in the telluric /database directory.")
                                telapfileslength = len(telapfiles)
                                telCheck_flag = telCheck_flag and True

                            logger.info("Required telluric extraction reference files check complete.\n")
                            
                            if telCheck_flag:
                                logger.info("All telluric extraction reference files available in %s\n", teldatabasepath)
                            else:
                                logger.warning("Parameter 'checkPeakMatch' is set to 'True', but one or more ")
                                logger.warning("telluric extraction reference files not available in ")
                                logger.warning("%s . Setting 'checkPeaksMatch' to 'False'.\n", teldatabasepath)
                                checkPeaksMatch = False
                                # TODO(Viraja):  Can ask the user if they want to perform checkPeaksMatch and set it at 
                                # this point. This would probably need the checks to be called as a function because the 
                                # script must check for the required telluric extraction reference files once the 
                                # parameter 'checkPeaksMatch' is set to 'True'.
                        else:
                            logger.warning("Parameters 'nsextractInter' and 'checkPeaksMatch' both set to 'False'. ")
                            logger.warning("After the science spectra extraction, please check manually if nsextract ")
                            logger.warning("identified the science peaks at expected locations.\n")

                ###########################################################################
                ##                                                                       ##
                ##                 COMPLETE - OBSERVATION SPECIFIC SETUP                 ##
                ##            BEGIN EXTRACTING 1D SPECTRA FOR AN OBSERVATION             ##
                ##                                                                       ##
                ###########################################################################
                    
                if manualMode:
                    a = raw_input("About to enter extract 1D spectra.")

                if useApall:
                    # This performs a weighted extraction
                    apertureTracingColumns = 20
                    extractSpectra1D(srccombimage, nsextractInter, useApall, apertureTracingColumns, \
                        extractionApertureRadius, overwrite)
                else:
                    apertureTracingColumns = 10
                    extractSpectra1D(srccombimage, nsextractInter, useApall, apertureTracingColumns, \
                        extractionApertureRadius, overwrite)
                # If the parameter 'calculateSpectrumSNR' is set to 'yes', the script will extract spectra from 
                # the combined sky image; else, it will only extract spectra from the combined source image.
                if calculateSpectrumSNR:
                    logger.info("Extracting the combined sky spectrum reduced without sky subtraction.\n")
                    if useApall:
                        apertureTracingColumns = 20
                        extractSpectra1D(srccombimage, nsextractInter, useApall, apertureTracingColumns, \
                            extractionApertureRadius, overwrite)
                    else:
                        apertureTracingColumns = 10
                        extractSpectra1D(srccombimage, nsextractInter, useApall, apertureTracingColumns, \
                            extractionApertureRadius, overwrite)
                
                if not nsextractInter:
                    if peaksMatch:
                        logger.info("Finding palusible peaks used by nsextract.")
                        telpeaks = peaksFind(teldatabasepath, telapfileslength, telcombfilename) 
                        scipeaks = peaksFind(scidatabasepath, sciapfileslength, scicombfilename)
                        logger.info("Completed finding palusible peaks used by nsextract.")

                        peaksMatch



                        #Re-extract science target spectrum if needed
                        

                Parameter 'extractionApertureRadius' is set to 23 pixels here (-/+ 23 pixels about the aperture covers the 6.9", almost whole length of slit),
                apertureTracingColumns = 20
                extractionApertureRadius = 23

                Extraction in steps on either side of the peak
	            # Setting fl_apall='yes' this time,



            logger.info("##############################################################################")
            logger.info("#                                                                            #")
            logger.info("#  COMPLETE - Extracting 1D spectra completed for                            #")
            logger.info("#  %s", obspath                                                                )
            logger.info("#                                                                            #")
            logger.info("##############################################################################\n")

    # Return to directory script was begun from.
    os.chdir(path)

    return

##################################################################################################################
#                                                     ROUTINES                                                   #
##################################################################################################################

def extractSpectra1D(combinedimage, nsextractInter, useApall, apertureTracingColumns, extractionApertureRadius, overwrite):
    """
    Extracting 1D spectra from the combined 2D spectra using nsextract.
    """
    logger = log.getLogger('gnirsReduce.extractSpectra1D')

    if os.path.exists('v'+combinedimage):
        if overwrite:
            logger.warning("Removing old v%s", combinedimage)
            os.remove('v'+combinedimage)
            iraf.nsextract(inimages=combinedimage, outspectra='', outprefix='v', dispaxis=1, \
                database='', line=700, nsum=apertureTracingColumns, ylevel='INDEF', \
                upper=str(extractionApertureRadius), lower='-'+str(extractionApertureRadius), background='none', \
                fl_vardq='yes', fl_addvar='no', fl_skylines='yes', fl_inter=nsextractInter, fl_apall=useApall, \
                fl_trace='no', aptable='gnirs$data/apertures.fits', fl_usetabap='no', fl_flipped='yes', \
                fl_project='yes', fl_findneg='no', bgsample='*', trace='', tr_nsum=10, tr_step=10, tr_nlost=3, \
                tr_function='legendre', tr_order=5, tr_sample='*', tr_naver=1, tr_niter=0, tr_lowrej=3.0, \
                tr_highrej=3.0, tr_grow=0.0, weights='variance', logfile=logger.root.handlers[0].baseFilename, \
                verbose='yes', mode='al')
        else:
            logger.warning("Old %s exists and -overwrite not set - skipping nsextract for observations.", combinedimage)
    else:
        iraf.nsextract(inimages=combinedimage, outspectra='', outprefix='v', dispaxis=1, database='',\
            line=700, nsum=apertureTracingColumns, ylevel='INDEF', upper=str(extractionApertureRadius), \
            lower='-'+str(extractionApertureRadius), background='none', fl_vardq='yes', fl_addvar='no', \
            fl_skylines='yes', fl_inter=nsextractInter, fl_apall=useApall, fl_trace='no', \
            aptable='gnirs$data/apertures.fits', fl_usetabap='no', fl_flipped='yes', fl_project='yes', \
            fl_findneg='no', bgsample='*', trace='', tr_nsum=10, tr_step=10, tr_nlost=3, tr_function='legendre', \
            tr_order=5, tr_sample='*', tr_naver=1, tr_niter=0, tr_lowrej=3.0, tr_highrej=3.0, tr_grow=0.0, \
            weights='variance', logfile=logger.root.handlers[0].baseFilename, verbose='yes', mode='al')

#---------------------------------------------------------------------------------------------------------------------#

def peaksFind(databasepath, apfileslength, combinedimage):
    """
    Check the telluric or science extraction reference files in the telluric directory databases to find the location 
    of the respective peaks.
    """
    logger = log.getLogger('gnirsReduce.peaksFind')

    peaks = []
    for i in range(apfileslength):
        apfile = open(databasepath+'/ap'+combinedimage[:-5]+'_SCI_'+str(i)+'_', 'r') 
        for line in apfile:
            # Get the peak location, which is the number in the second column of the line beginning with 'center'
            if 'center' in line:
                peaks.append(line.split()[1])
            else:
                peaks.append('Peak not found')
    return peaks

#---------------------------------------------------------------------------------------------------------------------#

def peaksMatch(scisrccomb, telsrccomb, scipeaks, telpeaks, toleranceOffset):
    """
    Checks is NSEXTRACT located the peak of the science at a resonable location along the slit given the aperture
    center of extraction of the telluric.

    For faint targets, NSEXTRACT often finds a noise peak instead of the science peak. In such cases, it is advisable  
    to check the aperture center of extraction of the science with respect to the telluricre and re-extract at the 
    expected location. Here,
        First, check that all the extensions of the combined 2D image have been extracted.
        Third, look the science and telluric absolute Q offsets and determine if the relative location of the target 
        peak was correct.
    If not, re-extract at the expected location.
    """
    logger = log.getLogger('gnirsReduce.peaksMatch')
    
    # Find absolute Q offsets of the combined science and telluric images from their respective last acquisition images
    sciacqHeader = fits.open(sciacq)[0].header
    scisrccombHeader = fits.open(scisrccomb)[0].header
    scisrccombQoffset = abs(sciacqHeader['QOFFSET'] - scisrccombHeader['QOFFSET'])

    telacqHeader = fits.open(telacq)[0].header
    telsrccombHeader = fits.open(telsrccomb)[0].header
    telsrccombQoffset = abs(telacqHeader['QOFFSET'] - telsrccombHeader['QOFFSET'])

    pixelscale = scisrccombHeader['PIXSCALE']
    pixeldifference = (scisrccombQoffset - telsrccombQoffset)/pixelscale  ## units: [pixels]
    peaks_flag = []
    # nsextract should find the spectrum within a 'tolerance' pixels of expected location. This depends on how well the 
    # observer centred the target along the slit. Here, we use 5 pixels as a reasonable tolerance level. A more robust
    # way would be to use some measure of whether the peak found by nsextract was real, e.g. counts + FWHM. However, 
    # this information is not recorded in database.
    for i in range(len(scipeaks)):
        expectedPeak = float(telpeaks[i]) + pixeldifference
        if scipeaks[i] == 'Peak not found':
            logger.warning("In extension %d for the science, nsextract did not extract anything. Re-extracting ", i)
            logger.warning("the spectrum forcing the aperture to be at the expected peak location %.4g", expectedPeak)
            peaks_flag.append(False)
        else:    
            locatedPeak = float(scipeaks[i])
            if abs(locatedPeak - expectedPeak) < toleranceOffset:
                logging.info("In extension %d for the science, nsextract detected the spectrum close to the ", i)
                logging.info("expected peak location along slit (located position = %s vs. expected ", locatedPeak)
                logging.info("position = %s .", expectedPeak)
                peaks_flag.append(True)
            else:
                logger.warning("In extension %d for the science, nsextract extracted an unexpected peak location ", i)
                logger.warning("along the slit (located position = %s vs. expected position = ", locatedPeak)
                logger.warning("%s . It is probably extracting noise; re-extracting the spectrum ", expectedPeak)
                logger.warning("forcing the aperture to be at the expected peak location.'")
                peaks_flag.append(False)
                reExtract=True
    return peaks_flag

#---------------------------------------------------------------------------------------------------------------------#

def reExtractSpectra1D(scidatabasepath, teldatabasepath, combinedimage, peaks_flag):
    """
    """
    logger = log.getLogger('gnirsReduce.reExtractSpectra1D')

    for i in range(len(peaks_flag)):
        logger.info("Creating new aperture files in database.")   
        # There is some trouble replacing only the database files for the extracted spectra that were not well centred.
        # So, simply replacing all science database files but using the ones for which nsextract located the peaks at 
        # the right positions.
        oldapfile = scidatabasepath+'/ap'+combinedimage[:-5]+'_SCI_'+str(i)+'_'
        if os.path.exists(oldapfile):
            os.remove(oldapfile)
        apfile = open ('database/apstandard_comb_SCI_'+str(i+1)+'_', 'r')
        newapfile = open ('database/apnewstandard_comb_SCI_'+str(i+1)+'_', 'w')
        if found_spectrum[i] == 0:
            clean  = apfile.read().replace(str(std_peak[i]), str(std_peak[i] + diff)+' ').replace('standard_comb','newstandard_comb')
        else:
            clean  = apfile.read().replace(str(std_peak[i]), str(tgt_peak[i])+' ').replace('standard_comb','newstandard_comb')
        newapfile.write(clean)
        apfile.close()
        newapfile.close()
    shutil.copy(standard, 'new'+standard)
    iraf.imdelete(images='v'+target)
    #Using these settings in nsextract will force it to use the aperture size and centre in the edited "apnewstandard_comb" files
    iraf.nsextract(inimages=target,outspectra='',outprefix='v',dispaxis=1,database='',line=700,nsum=20,ylevel='INDEF',upper=aperture_plus,lower=aperture_minus,background='none',fl_vardq='yes',fl_addvar='no',fl_skylines='yes',fl_inter=nsextinter,fl_apall='yes',fl_trace='no',aptable=path_to_nsextract+'config/apertures.fits',fl_usetabap='no',fl_flipped='yes',fl_project='yes',fl_findneg='no',bgsample='*',trace='newstandard_comb',tr_nsum=10,tr_step=10,tr_nlost=3,tr_function='legendre',tr_order=5,tr_sample='*',tr_naver=1,tr_niter=0,tr_lowrej=3.0,tr_highrej=3.0,tr_grow=0.0,weights='variance',logfile='',verbose='yes',mode='al')

    #Slight complication - we occasionally find that nsextract locates the aperture too close to the end of the slit
    #Then it exits with an "Aperture too large" error and spectra aren't extracted for one or more extensions
    #But we work around that above, so when we check for errors in XDpiped.csh, we ignore this error
    #But maybe that error can happen for other reasons
    #So to be on the safe side, will check that all extensions are present in the extracted target file (should really add other files as well...)
    check_extn = iraf.gemextn(inimages=target,check='exists,mef',process='expand',index='',extname='SCI',extversion='',ikparams='',omit='',replace='',outfile='STDOUT',logfile='',Stdout=1,glogpars='',verbose='no',fail_count='0', count='20', status='0')
    if len(check_extn) != 6:
        print "ERROR: target_comb file contains only ", len(check_extn), 'extensions. Exiting script.'

#---------------------------------------------------------------------------------------------------------------------#

def stepwiseExtractSpectra1D(combinedimage, nsextractInter, useApall, apertureTracingColumns):
    """
    Extracts science spectra along (approximately) the full slit. 
    
    This method is appropriate for objects centred along length of slit (absolute Q offset for the science = 0). The 
    effect, if nsextract does not find a spectrum centred along the slit, is not known at this point.
    
    CAUTION NOTE:  From XDGNIRS, full slit and stepwise extractions have not been used or tested thoroughly. So, 
    please double check your results.
    """
    logger = log.getLogger('gnirsReduce.stepwiseExtractSpectra1D')

    iraf.nsextract(inimages=combinedimage, outspectra='', outprefix='a', dispaxis=1, database='', line=700, \
        nsum=apertureTracingColumns, ylevel='INDEF', upper=str(extractionApertureRadius), \
        lower='-'+str(extractionApertureRadius), background='none', fl_vardq='yes', fl_addvar='no', fl_skylines='yes',\
        fl_inter=nsextractInter, fl_apall=useApall, fl_trace='no', aptable='gnirs$data/apertures.fits', \
        fl_usetabap='no', fl_flipped='yes', fl_project='yes', fl_findneg='no', bgsample='*', trace='', tr_nsum=10, \
        tr_step=10, tr_nlost=3, tr_function='legendre', tr_order=5, tr_sample='*', tr_naver=1, tr_niter=0, \
        tr_lowrej=3.0, tr_highrej=3.0, tr_grow=0.0, weights='variance', logfile=logger.root.handlers[0].baseFilename, \
        verbose='yes', mode='al')

    # First trace the peak to make sure that the same part of the object is extracted in each step along the slit for  
    # all orders. This is required when there is complex, e.g., structure varying in the spectral direction in an  
    # extended object such as a galaxy; otherwise the spectra can have offsets between orders.
    
    # This first nsextract step performed outside the loop gets the trace into the science extraction database to be 
    # used during the actual stepwise extraction.
    iraf.nsextract(inimages=target, outspectra='trace_ref', outprefix='x', dispaxis=1, database='', line=700, nsum=20, ylevel='INDEF', upper=3, \
        lower=-3, background='none', fl_vardq='yes', fl_addvar='no', fl_skylines='yes', fl_inter=nsextractInter, fl_apall=useApall, fl_trace='yes', \
        aptable='gnirs$data/apertures.fits', fl_usetabap='no', fl_flipped='yes', fl_project='no', fl_findneg='no', bgsample='*', trace='', \
        tr_nsum=10, tr_step=10, tr_nlost=3, tr_function='legendre', tr_order=5, tr_sample='300:1000', tr_naver=1, tr_niter=0, tr_lowrej=3.0, \
        tr_highrej=3.0, tr_grow=0.0, weights='variance', logfile=logger.root.handlers[0].baseFilename, verbose='yes', mode='al')
    
    # This step is never done interactively, because it uses info from the previous call (and it would be very tedious)
    n = 0
    for j in range (-21,21,step):
        n = n + 1
        iraf.nsextract(inimages=target, outspectra='', outprefix='s'+str(n), dispaxis=1, database='', line=700, nsum=20, ylevel='INDEF', \
            upper=j+step, lower=j, background='none', fl_vardq='yes', fl_addvar='no', fl_skylines='yes', fl_inter=nsextractInter, fl_apall=useApall, \
            fl_trace='no', aptable='gnirs$data/apertures.fits', fl_usetabap='no', fl_flipped='yes', fl_project='yes', fl_findneg='no', \
            bgsample='*', trace='', tr_nsum=10, tr_step=10, tr_nlost=3, tr_function='legendre', tr_order=5, tr_sample='*', tr_naver=1, tr_niter=0, \
            tr_lowrej=3.0, tr_highrej=3.0, tr_grow=0.0, weights='variance', logfile=logger.root.handlers[0].baseFilename, verbose='yes', mode='al')
		
    nsteps = open ("nsteps.txt", "w")
    nsteps.write (str(n))
    nsteps.close()

#---------------------------------------------------------------------------------------------------------------------#

if __name__ == '__main__':
    log.configure('gnirs.log', filelevel='INFO', screenlevel='DEBUG')
    start('gnirs.cfg')
