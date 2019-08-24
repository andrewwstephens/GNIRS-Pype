#!/usr/bin/env python

import log
import os
import glob
import ConfigParser
import datetime
import dateutil.parser
from astropy.io import fits
from pyraf import iraf


def start(kind, configfile):
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
        - If telluric:  cleaned (optional), prepared, radiation corrected, reduced, spatial distortion corrected, and
          transformed images
        - If science:  cleaned (optional), prepared, radiation corrected, reduced, spatial distortion corrected, and
          transformed images

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
    logger = log.getLogger('gnirsReduce.start')

    # TODO(nat): Right now the pipeline will crash if you decide to skip, say, doing a bad pixel correction. This is 
    # because each step adds a prefix to the frame name, and most following steps depend on that prefix being there.
    # One way to fix this is if a step is to be skipped, iraf.copy() is called instead to copy the frame and add the 
    # needed prefix. Messy but it might work for now.

    ###########################################################################
    ##                                                                       ##
    ##                  BEGIN - GENERAL REDUCTION SETUP                      ##
    ##                                                                       ##
    ###########################################################################

    path = os.getcwd()  # Store current working directory for later use.

    logger.info('####################################################')
    logger.info('#                                                  #')
    logger.info('#  Start the GNIRS Science and Telluric Reduction  #')
    logger.info('#                                                  #')
    logger.info('####################################################\n')

    # Set up/prepare IRAF.
    iraf.gemini()
    iraf.gemtools()
    iraf.gnirs()

    # Reset to default parameters the used IRAF tasks.
    iraf.unlearn(iraf.gemini, iraf.gemtools, iraf.gnirs, iraf.imcopy)

    # Prepare the IRAF package for GNIRS.
    # NSHEADERS lists the header parameters used by the various tasks in the GNIRS package (excluding headers values 
    # which have values fixed by IRAF or FITS conventions).
    iraf.nsheaders("gnirs", logfile=logger.root.handlers[0].baseFilename)

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
    if kind == 'Science':  ## scienceReduction specific config
        observationSection = 'ScienceDirectories'
        start = config.getint('scienceReduction','Start')
        stop = config.getint('scienceReduction','Stop')
        cleanir = config.getboolean('scienceReduction','cleanir')
        radiationCorrectionMethod = config.get('scienceReduction','radiationCorrectionMethod')
        radiationThreshold = config.getfloat('scienceReduction','radiationThreshold')
        nsprepareInter = config.getboolean('interactive','nsprepareInter')
        calculateSpectrumSNR = config.getboolean('gnirsPipeline','calculateSpectrumSNR')
    elif kind == 'Telluric':  ## telluricReduction specific config
        observationSection = 'TelluricDirectories'
        start = config.getint('telluricReduction','Start')
        stop = config.getint('telluricReduction','Stop')
        cleanir = config.getboolean('telluricReduction','cleanir')
        radiationCorrectionMethod = config.get('telluricReduction','radiationCorrectionMethod')
        radiationThreshold = config.getfloat('telluricReduction','radiationThreshold')
        nsprepareInter = config.getboolean('interactive','nsprepareInter')
        calculateSpectrumSNR = config.getboolean('gnirsPipeline','calculateSpectrumSNR')
    else:
        logger.error("###########################################################################")
        logger.error("###########################################################################")
        logger.error("#                                                                         #")
        logger.error("#     ERROR in reduce: invalid kind of reduction. Please enter either     #")
        logger.error("#                      <Science> or <Telluric>. Exiting script.           #")
        logger.error("#                                                                         #")
        logger.error("###########################################################################")
        logger.error("###########################################################################\n")
        raise SystemExit  
        ## TODO(Viraja): can ask for a raw input on the command line to avoid a system exit

    ###########################################################################
    ##                                                                       ##
    ##                 COMPLETE - GENERAL REDUCTION SETUP                    ##
    ##                                                                       ##
    ###########################################################################

    # gnirsReduce will reduce observations in each science or telluric directory only if the reduction truth value for
    # that directory is True; else, it will skip the reductions (of science or telluric frames) in that dorectory. If 
    # True, the script will check if all required calibrations to reduce the science or telluric frames are available
    # in the calibrations directory -- if they are not, it will warn the user and request to provide them.

    # Loop through all the observation (telluric or science) directories to perform a reduction on each one.
    for obspath in config.options(observationSection):

        if not config.getboolean(observationSection, obspath):
            logger.debug('Skipping %s', obspath)
            continue

        ###########################################################################
        ##                                                                       ##
        ##                  BEGIN - OBSERVATION SPECIFIC SETUP                   ##
        ##                                                                       ##
        ###########################################################################

        iraf.chdir(obspath + '/Intermediate')
        logger.info("Working on reductions in %s\n", obspath)

        calpath = '../Calibrations'
        logger.info("Path to calibrations: %s\n", calpath)

        # TODO(Viraja)?:  Define a function to extract the lists by specifying their filenames
        allobsfilename = 'all.list'
        allobslist = open(allobsfilename, "r").readlines()
        allobslist = [filename.strip() for filename in allobslist]
        rawHeader = fits.open(allobslist[0])[0].header
        if calculateSpectrumSNR:
            # Check if there is a list of sky images in obspath
            skylistfilename = 'sky.list'
            if os.path.exists(skylistfilename):
                skylist = open(skylistfilename, "r").readlines()
                skylist = [filename.strip() for filename in skylist]
            else:
                logger.warning("Parameter 'calculateSpectrumSNR' is 'True', but a list of sky images not ")
                logger.warning("available in %s . Setting the 'calculateSpectrumSNR' parameter for the ", obspath)
                logger.warning("the current set of observations to 'False'.\n")
                calculateSpectrumSNR = False
        nodAlistfilename = 'nodA.list'
        nodAlist = open(nodAlistfilename, "r").readlines()
        nodAlist = [filename.strip() for filename in nodAlist]
        nodBlistfilename = 'nodB.list'
        nodBlist = open(nodBlistfilename, "r").readlines()
        nodBlist = [filename.strip() for filename in nodBlist]

        logger.debug("Checking if required calibrations available in %s\n", calpath)

        # Check for the reference image to calculate MDF
        QHflatslist = open(calpath+'/QHflats.list', "r").readlines()
        mdfshiftimage = calpath+'/n'+QHflatslist[0].strip()
        if os.path.exists(mdfshiftimage):
            logger.info("Reference image to calculate MDF information available.")
            calCheck_flag = True   # calCheck_flag defined here
        else:
            logger.warning("Reference image to calculate MDF information not available.")
            calCheck_flag = False   #  calCheck_flag defined here

        # Check for the masterflat
        masterflat = calpath+'/masterflat.fits'
        if os.path.exists(masterflat):
            logger.info("Masterflat to apply flat field correction available.")
            calCheck_flag = calCheck_flag and True
        else:
            logger.warning("Masterflat to apply flat field correction not available.")
            calCheck_flag = calCheck_flag and False

        # Check for /database directory (possibly) containing calibration files
        databasepath = calpath+'/database'
        if os.path.exists(databasepath):
            logger.info("Reference /database directory (possibly) containing calibration files available.")
            calCheck_flag = calCheck_flag and True
        else:
            logger.warning("Reference /database directory (possibly) containing calibration files not available.")
            calCheck_flag = calCheck_flag and False

        # Check for spatial distortion correction calibration files
        sdistfiles = glob.glob(databasepath+'/*_sdist')
        if not sdistfiles:
            logger.warning("Reference files to apply spatial distortion correction not available in the ")
            logger.warning("/database directory.")
            calCheck_flag = calCheck_flag and False
        else:
            logger.info("Reference files to apply spatial distortion correction available in the /database ")
            logger.info("directory.")
            sdistfileslength = len(sdistfiles)
            calCheck_flag = calCheck_flag and True

        # Check for spectral transformation calibration files
        wavecallampfiles = glob.glob(databasepath+'/*_lamp')
        if not wavecallampfiles:
            logger.warning("Reference files to apply spectral transformation not available in the /database ")
            logger.warning("directory.\n")
            calCheck_flag = calCheck_flag and False
        else:
            logger.info("Reference files to apply spectral transformation available in the /database ")
            logger.info("directory.\n")
            wavecallampfileslength = len(wavecallampfiles)
            calCheck_flag = calCheck_flag and True

        logger.info("Calibrations check complete.\n")
        if calCheck_flag:
            logger.info("All reference calibration files available in %s\n", calpath)
        else:
            logger.warning("One or more reference calibration files not available in %s ", calpath + ". Please ")
            logger.warning("run gnirsBaselineCalibration.py to create the required calibration files or update ")
            logger.warning("them manually. Exiting script.\n")
            raise SystemExit

        combinedarc = 'arc_comb.fits'  ## The combined arc filename used in the calibrations directory path

        ###########################################################################
        ##                                                                       ##
        ##                 COMPLETE - OBSERVATION SPECIFIC SETUP                 ##
        ##                BEGIN DATA REDUCTION FOR AN OBSERVATION                ##
        ##                                                                       ##
        ###########################################################################

        # Check start and stop values for reduction steps. Ask user for a correction if input is not valid.
        valindex = start
        while valindex > stop or valindex < 1 or stop > 5:
            logger.warning("#####################################################################")
            logger.warning("#####################################################################")
            logger.warning("#                                                                   #")
            logger.warning("#     WARNING in reduce: invalid start/stop values of observation   #")
            logger.warning("#                        reduction steps.                           #")
            logger.warning("#                                                                   #")
            logger.warning("#####################################################################")
            logger.warning("#####################################################################\n")

            valindex = int(raw_input("Please enter a valid start value (1 to 5, default 1): "))
            stop = int(raw_input("Please enter a valid stop value (1 to 5, default 5): "))

        while valindex <= stop:

            #############################################################################
            ##  STEP 1: Clean raw observations.                                        ##
            ##  Output: Cleaned science or telluric frames.                            ##
            #############################################################################

            if valindex == 1:
                if manualMode:
                    a = raw_input("About to enter step 1: clean raw observations.")

                if cleanir:
                    # cleanir(allobslist)
                    pass  ## this step is not modified to work with this pipeline as of July 2019
                else:
                    logger.info("######################################################################")
                    logger.info("######################################################################")
                    logger.info("#                                                                    #")
                    logger.info("#       WARNING in reduce: raw observations not cleaned.             #")
                    logger.info("#                                                                    #")
                    logger.info("######################################################################")
                    logger.info("######################################################################\n")

                logger.info("##################################################################")
                logger.info("#                                                                #")
                logger.info("#       STEP 1: Clean raw observations - COMPLETED               #")
                logger.info("#                                                                #")
                logger.info("##################################################################\n")

            ###########################################################################
            ##    STEP 2: Prepare observations (science or telluric frames)          ##
            ###########################################################################

            elif valindex == 2:
                if manualMode:
                    a = raw_input("About to enter step 2: locate the spectrum.")

                # Use the most appropriate bad pixel mask. For data taken Before the summer 2012 lens replacement,
                # use 'gnirs$data/gnirsn_2011apr07_bpm.fits'; after summer 2012, use
                # 'gnirs$data/gnirsn_2012dec05_bpm.fits'. Use keyword 'ARRAYID' in the raw file header to check
                # which camera was used for observations, and accordingly set parameter 'bpm' for the bpmfile to be
                # used from the gnirs$database.
                arrayid = rawHeader['ARRAYID'].strip()
                if arrayid == 'SN7638228.1':
                    bpmfile = 'gnirs$data/gnirsn_2011apr07_bpm.fits'
                elif arrayid == 'SN7638228.1.2':
                    bpmfile = 'gnirs$data/gnirsn_2012dec05_bpm.fits'
                else:
                    logger.error("######################################################################")
                    logger.error("######################################################################")
                    logger.error("#                                                                    #")
                    logger.error("#        ERROR in reduce: unknown array ID. Exiting script.          #")
                    logger.error("#                                                                    #")
                    logger.error("######################################################################")
                    logger.error("######################################################################\n")
                    raise SystemExit

                prepareObservations(nsprepareInter, mdfshiftimage, bpmfile, overwrite)

                logger.info("##############################################################################")
                logger.info("#                                                                            #")
                logger.info("#       STEP 2: Locate the Spectrum (prepare observations) - COMPLETED       #")
                logger.info("#                                                                            #")
                logger.info("##############################################################################\n")

            ###########################################################################
            ##  STEP 3: Prepare bad pixel masks and do radiation event correction    ##
            ###########################################################################

            elif valindex == 3:
                if manualMode:
                    a = raw_input("About to enter step 3: radiation event correction.")

                preparedHeader = fits.open('n'+allobslist[0])[0].header
                date = dateutil.parser.parse(preparedHeader['DATE-OBS'].strip() + ' ' + preparedHeader['TIME-OBS'])
                rdnoise = preparedHeader['RDNOISE']
                gain = preparedHeader['GAIN']

                min_nodA = 'min_nodA.fits'
                min_nodB = 'min_nodB.fits'
                minimagelist = [min_nodA, min_nodB]
                nodlist = [nodAlist, nodBlist]

                if radiationCorrectionMethod == 'fixpix':
                    createMinimumImage(nodAlistfilename, nodAlist, min_nodA, overwrite)
                    createMinimumImage(nodBlistfilename, nodBlist, min_nodB, overwrite)
                    radiationCorrectionFixpix(nodlist, minimagelist, radiationThreshold, rdnoise, gain, overwrite)

                elif radiationCorrectionMethod == 'dqplane':
                    createMinimumImage(nodAlistfilename, nodAlist, min_nodA, overwrite)
                    radiationCorrectionDQplane(allobslist, nodAlist, min_nodA, overwrite)
                    createMinimumImage(nodBlistfilename, nodBlist, min_nodB, overwrite)
                    radiationCorrectionDQplane(allobslist, nodBlist, min_nodB, overwrite)

                else:
                    logger.warning("########################################################################")
                    logger.warning("########################################################################")
                    logger.warning("#                                                                      #")
                    logger.warning("#   WARNING in reduce: invalid/no radiation event correction method    #")
                    logger.warning("#                      method found. Checking the observation date.    #")
                    logger.warning("#                                                                      #")
                    logger.warning("########################################################################")
                    logger.warning("########################################################################\n")

                    if date < datetime.datetime(year=2012, month=8, day=1, hour=0, minute=0, second=0):
                        logger.info("Observation date is " + str(date) + ". Performing the radiation event ")
                        logger.info("correction using the 'fixpix' method.\n")
                        radiationCorrectionMethod = 'fixpix'
                        createMinimumImage(nodAlistfilename, nodAlist, min_nodA, overwrite)
                        createMinimumImage(nodBlistfilename, nodBlist, min_nodB, overwrite)
                        radiationCorrectionFixpix(nodlist, minimagelist, radiationThreshold, rdnoise, gain, overwrite)
                    else:
                        logger.info("Observation date is " + str(date) + ". GNIRS data from 2012B not affected ")
                        logger.info("by radiation events. Skipping radiation event correction.\n")
                        pass

                logger.info("##############################################################################")
                logger.info("#                                                                            #")
                logger.info("#  STEP 3: Radiation event correction - COMPLETED                            #")
                logger.info("#                                                                            #")
                logger.info("##############################################################################\n")

            ##############################################################################
            ##  STEP 4: Flat field and subtract sky background (science or telluric)    ##
            ##############################################################################

            elif valindex == 4:
                if manualMode:
                    a = raw_input("About to enter step 4: flat fielding and sky subtraction.")

                # If the parameter 'calculateSpectrumSNR' is set to 'yes' in the configuration file, the script
                # will reduce all observations without sky subtraction (by setting the parameter 'fl_sky' in
                # nsreduce to 'no') in which case the parameter 'outprefix' in nsreduce will be 'k'; else, all
                # observations will be reduced with sky subtraction and the output prefix will be 'r'.
                skySubtraction = 'yes'
                reduce_outputPrefix = 'r'  ## output prefix assigned to the reduced observations
                reduceObservations(skySubtraction, reduce_outputPrefix, masterflat, radiationThreshold, overwrite)
                if calculateSpectrumSNR:
                    logger.info("Reducing observations without sky subtraction to generate an SNR spectrum later ")
                    logger.info("by the pipeline.\n")
                    skySubtraction = 'no'
                    reduce_outputPrefix = 'k'
                    reduceObservations(skySubtraction, reduce_outputPrefix, masterflat, radiationThreshold,
                        overwrite)

                logger.info("##############################################################################")
                logger.info("#                                                                            #")
                logger.info("#  STEP 4: Flat fielding and sky subtraction - COMPLETED                     #")
                logger.info("#                                                                            #")
                logger.info("##############################################################################\n")

            #################################################################################
            ##  STEP 5: Apply spatial distortion correction and spectral transformation    ##
            #################################################################################

            elif valindex == 5:
                if manualMode:
                    a = raw_input("About to enter step 5: spatial distortion correction and spectral transformation.")

                reduce_outputPrefix = 'r'
                SdistCorrection_SpectralTransform(databasepath, reduce_outputPrefix, allobslist, sdistfileslength,
                    wavecallampfileslength, combinedarc, overwrite)
                if calculateSpectrumSNR:
                    logger.info("Applying spatial distortion correction and spectral transformation to the ")
                    logger.info("science frames reduced without sky subtraction.\n")
                    reduce_outputPrefix = 'k'
                    SdistCorrection_SpectralTransform(databasepath, reduce_outputPrefix, allobslist,
                        sdistfileslength, wavecallampfileslength, combinedarc, overwrite)

                logger.info("##################################################################################")
                logger.info("#                                                                                #")
                logger.info("# STEP 5: spatial distortion correction and spectral transformation - COMPLETED  #")
                logger.info("#                                                                                #")
                logger.info("##################################################################################\n")

            else:
                logger.error("###########################################################################")
                logger.error("###########################################################################")
                logger.error("#                                                                         #")
                logger.error("#        ERROR in reduce: %d is not valid. Exiting script.", valindex)
                logger.error("#                                                                         #")
                logger.error("###########################################################################")
                logger.error("###########################################################################")
                raise SystemExit

            valindex += 1

        logger.info("##############################################################################")
        logger.info("#                                                                            #")
        logger.info("#  COMPLETE - Reductions completed for                                       #")
        logger.info("#  %s", obspath)
        logger.info("#                                                                            #")
        logger.info("##############################################################################\n")

    iraf.chdir(path)  # Return to directory script was begun from.

    return


##################################################################################################################
#                                                     ROUTINES                                                   #
##################################################################################################################

def prepareObservations(interactive, mdfshiftimage, bpmfile, overwrite):
    """
    Prepare raw observations (science or telluric) using nsprepare. Output prefix "n" is added to raw filenames.

    Processing with NSPREPARE will rename the data extension and add variance and data quality extensions. By default
    (see NSHEADERS) the extension names are SCI for science data, VAR for variance, and DQ for data quality (0 = good).
    NSPREPARE will add an MDF file (extension MDF) describing the GNIRS image pattern and how all images
    cross-coorelate with the first prepared QHflat provided as the "mdfshiftimage" reference with the "shiftimage"
    parameter in NSPREPARE.

    :param interactive:
    :param mdfshiftimage:
    :param bpmfile:
    :param overwrite:
    :return:
    """

    logger = log.getLogger('gnirsReduce.prepareObservations')

    # Update frames with mdf offset value and generate variance and data quality extensions. This code structure checks 
    # if iraf output files already exist. If output files exist and overwrite is specified, iraf output is overwritten.

    oldfiles = glob.glob('./nN*.fits')
    if len(oldfiles) > 0:
        if overwrite:
            logger.warning('Removing old files: %s', oldfiles)
            for f in oldfiles:
                os.remove(f)
        else:
            logger.warning("Output exists and -overwrite not set - skipping nsprepare for all observations.\n")
            return

    iraf.nsprepare(
        inimages='@all.list', rawpath='', outimages='', outprefix='n', bpm=bpmfile,
        logfile=logger.root.handlers[0].baseFilename, fl_vardq='yes', fl_cravg='no', crradius=0.0,
        fl_dark_mdf='no', fl_correct='no', fl_saturated='yes', fl_nonlinear='yes', fl_checkwcs='yes',
        fl_forcewcs='yes', arraytable='gnirs$data/array.fits', configtable='gnirs$data/config.fits',
        specsec='[*,*]', offsetsec='none', pixscale='0.15', shiftimage=mdfshiftimage, shiftx=0., shifty=0.,
        obstype='FLAT', fl_inter=interactive, verbose='yes', mode='al')

    return


# ----------------------------------------------------------------------------------------------------------------------
def createMinimumImage(inlistfilename, inlist, minimage, overwrite):
    """
    Combine images using GEMCOMBINE to create a 'minimum' image using the "minmax" combining algorithm.
    In this pipeline, the function is used to combine observations at different nod positions positions, separately,
    to create 'minimum' images at each position.

    :param inlistfilename:
    :param inlist:
    :param minimage:
    :param overwrite:
    :return:
    """

    logger = log.getLogger('gnirsReduce.createMinimumImage')

    if os.path.exists(minimage):
        if overwrite:
            logger.warning("Removing old %s", minimage)
            os.remove(minimage)
        else:
            logger.warning("Output exists and overwrite flag is not set")
            logger.warning('Skipping creating a mimimum image for the input list %s.', inlistfilename)
            return

    logger.info("Creating a relatively clean 'minimum' image of observations in the input list ")
    logger.info("%s.\n", inlistfilename)
    iraf.gemcombine(
        input='n//@'+inlistfilename, output=minimage, title="", combine="median", reject="minmax",
        offsets="none", masktype="goodvalue", maskvalue=0., scale="none", zero="none", weight="none",
        statsec="[*,*]", expname="EXPTIME", lthreshold='INDEF', hthreshold='INDEF', nlow=0,
        nhigh=len(inlist)-1, nkeep=1, mclip='yes', lsigma=3., hsigma=3., key_ron="RDNOISE", key_gain="GAIN",
        ron=0.0, gain=1.0, snoise="0.0", sigscale=0.1, pclip=-0.5, grow=0.0, bpmfile="", nrejfile="",
        sci_ext="SCI", var_ext="VAR", dq_ext="DQ", fl_vardq='yes',
        logfile=logger.root.handlers[0].baseFilename, fl_dqprop='no', verbose='yes')

    return


# ----------------------------------------------------------------------------------------------------------------------
def radiationCorrectionFixpix(nodlist, minimagelist, radiationThreshold, rdnoise, gain, overwrite):
    """
    Prepare bad pixel masks for all observations (science or telluric) using a radiation threshold value and 
    correct the bad pixels by interpolation method.

    Here, the task proto.fixpix is set to interpolate along lines to avoid bad effects when a radiation event is right
    on a skyline. Column interpolation may be better for events landing on the peak of target spectra, but effects from
    bad interpolation over skylines seem to do worse.
    """
    logger = log.getLogger('gnirsReduce.createMinimumImage')

    oldfiles_masks = glob.glob('./mask*.pl')
    oldfiles_gmasks = glob.glob('./gmask*.pl')
    oldfiles_radiationcleaned = glob.glob('./lnN*.fits')
    oldfiles_combinedmasks = glob.glob('./bgmask*.pl')
    oldfiles = oldfiles_masks + oldfiles_gmasks + oldfiles_radiationcleaned + oldfiles_combinedmasks

    if len(oldfiles) > 0:
        if overwrite:
            logger.warning('Removing old mask*, gmask*.pl, lnN*.fits, and bgmask*.pl files.')
            for f in oldfiles:
                os.remove(f)
        else:
            logger.warning("Output exists and -overwrite not set - skipping radiation correction by 'fixpix' method ")
            logger.warning("for all observations.")
            return

    logger.info("Creating initial masks for all observations.")
    for i in range(len(nodlist)):
        for image in nodlist[i]:
            iraf.imexpr("(a-b)>"+str(radiationThreshold)+"*sqrt("+str(rdnoise)+"**2+2*b/"+str(gain)+") ? 1 : 0",
                output='mask'+image[:-5]+'.pl', a='n'+image+'[SCI]', b=minimagelist[i]+'[SCI]', dims="auto",
                intype="int", outtype="auto", refim="auto", bwidth=0, btype="nearest", bpixval=0.,
                rangecheck='yes', verbose='yes', exprdb="none")

    with open('./masks.list', 'a+') as f:
        for filename in sorted(glob.glob('./mask*.pl')):
            filename = filename[filename.rfind('/')+1:]
            if filename not in f.read().split('\n'):
                f.write(filename + '\n')

    iraf.crgrow(input='@masks.list', output='g//@masks.list', radius=1.5, inval="INDEF", outval="INDEF")

    logger.info("Adding bad pixels from the [DQ] plane generated by nsprepare to the masks so that they can be ")
    logger.info("corrected as well.\n")

    for i in range(len(nodlist)):
        for image in nodlist[i]:
            iraf.copy(input='n'+image, output='ln'+image, verbose='yes')

            iraf.imexpr(
                expr="a||b", output='bgmask'+image[:-5]+'.pl', a='gmask'+image[:-5]+'.pl',
                b=minimagelist[i]+'[DQ]', dims="auto", intype="int", outtype="auto", refim="auto", bwidth=0,
                btype="nearest", bpixval=0., rangecheck='yes', verbose='yes', exprdb="none")

            iraf.proto.fixpix(images='ln'+image+'[SCI,1]', masks='bgmask'+image[:-5]+'.pl', linterp=1,
                cinterp="INDEF", verbose='yes', pixels='no')

    return


# ----------------------------------------------------------------------------------------------------------------------
def radiationCorrectionDQplane(obslist, inlist, minimage, overwrite):
    """
    Prepare bad pixel masks for all observations (science or telluric) using the data quality plane of the prepared 
    observations and correct the bad pixels by replacing the corresponding pixels from 'minimum' images.
    """
    logger = log.getLogger('gnirsReduce.radiationCorrectionDQplane')

    oldfiles_masks = glob.glob('./mask*')
    oldfiles_gmasks = glob.glob('./gmask*')
    oldfiles_radiationclean = glob.glob('./lnN*.fits')
    oldfiles = oldfiles_masks + oldfiles_gmasks + oldfiles_radiationclean

    if len(oldfiles) > 0:
        if overwrite:
            logger.warning('Removing old mask*, gmask*.pl, and lnN*.fits files.')
            for f in oldfiles:
                os.remove(f)
        else:
            logger.warning("Output files exist and -overwrite not set - skipping radiation correction by 'dqplane' ")
            logger.warning("method for all observations.\n")
            return

    logger.info("Copying the [DQ] planes generated by nsprepare as the initial masks for all observations.")
    for obs in obslist:
        iraf.imcopy(input='n'+obs+'[DQ,1]', output='mask'+obs[:-5]+'.pl', verbose='yes')

    with open('./masks.list', 'a+') as f:
        for filename in sorted(glob.glob('./mask*.pl')):
            filename = filename[filename.rfind('/')+1:]
            if filename not in f.read().split('\n'):
                f.write(filename + '\n')

    iraf.crgrow(input='@masks.list', output='g//@masks.list', radius=1.5, inval="INDEF", outval="INDEF")

    logger.info("Replacing the pixels identified as bad pixels by the final masks with pixels from the ")
    logger.info("'minimum' images.\n")
    for image in inlist:
        iraf.copy(input='n'+image, output='ln'+image, verbose='yes')
        iraf.imexpr(
            expr="c>1 ? b : a", output='ln'+image+'[SCI,overwrite]', a='n'+image+'[SCI]', b=minimage+'[SCI]',
            c='gmask'+image[:-5]+'.pl', dims="auto", intype="int", outtype="auto", refim="auto", bwidth=0,
            btype="nearest", bpixval=0., rangecheck='yes', verbose='yes', exprdb="none")

    return


# ----------------------------------------------------------------------------------------------------------------------
def reduceObservations(skySubtraction, outprefix, masterflat, radiationThreshold, overwrite):
    """
    Flat field and sky subtract observations with nsreduce. If 'calculateSpectrumSNR' is set, observations will be
    reduced without sky subtraction too and prefix 'k' will be added before the input filenames; otherwise, sky 
    subtraction will be performed and the prefix will be 'r'.

    NSREDUCE is used for basic reduction of raw data - it provides a single, unified interface to several tasks and 
    also allows the subtraction of dark frames (not done for GNIRS data) and dividing by the flat.

    :param skySubtraction:
    :param outprefix:
    :param masterflat:
    :param radiationThreshold:
    :param overwrite:
    :return:
    """
    logger = log.getLogger('gnirsReduce.reduceObservations')

    oldfiles = glob.glob('./' + outprefix + 'lnN*.fits')
    if len(oldfiles) > 0:
        if overwrite:
            logger.warning('Removing old files: %s', oldfiles)
            for f in oldfiles:
                os.remove(f)
        else:
            logger.warning("Output files exist and overwrite flag not set.")
            logger.warnign('Skipping nsreduce for all observations.')
            return

    iraf.nsreduce(
        inimages='ln//@all.list', outimages='', outprefix=outprefix, fl_cut='yes', section="",
        fl_corner='yes', fl_process_cut='yes', fl_nsappwave='no', nsappwavedb='gnirs$data/nsappwave.fits',
        crval="INDEF", cdelt="INDEF", fl_dark='no', darkimage="", fl_save_dark='no', fl_sky=skySubtraction,
        skyimages="", skysection="", combtype="median", rejtype="avsigclip", masktype="goodvalue", maskvalue=0.,
        scale="none", zero="median", weight="none", statsec="[*,*]", lthreshold="INDEF", hthreshold="INDEF",
        nlow=1, nhigh=1, nkeep=0, mclip='yes', lsigma=3., hsigma=3., snoise="0.0", sigscale=0.1, pclip=-0.5,
        grow=0.0, skyrange='INDEF', nodsize=3., fl_flat='yes', flatimage=masterflat, flatmin=0.0, fl_vardq='yes',
        logfile=logger.root.handlers[0].baseFilename, verbose='yes', debug='no', force='no')

    # Record radiation threshold value in reduced file headers
    iraf.hedit(images=outprefix+'ln//@all.list', fields='RTHRESH', value=radiationThreshold, add='yes',
        addonly='no', delete='no', verify='no', show='no', update='yes')

    return


# ----------------------------------------------------------------------------------------------------------------------
def SdistCorrection_SpectralTransform(databasepath, outprefix, obslist, sdistfileslength,
                                      wavecallampfileslength, combinedarc, overwrite):
    """
    Apply spatial distortion correction and spectral transformation. 

    In general, we would run nsfitcoords on all reduced science or telluric frames with and without sky subtraction (if 
    creating an error spectrum) to straighten the spectra and then de-tilt. nsfitcoords simply writes database info, 
    which depends only on pinholes and arcs (the *_sdist and *_lamp files that were written to the databse while 
    reducing baseline calibrations), not on on-sky data. So, here we access the nsfitcoords solution from the 
    calibrations folder and append the relevant header keywords, which will tell nstransform the database files it 
    should use, without actually calling nsfitcoords. This routine also uses the output files in the database that were
    created by the nswavelength task.
    """

    logger = log.getLogger('gnirsReduce.SdistCorrection_SpectralTransform')

    oldfiles_fitcoords = glob.glob('./f'+outprefix+'lnN*')
    oldfiles_sdisttrans = glob.glob('./tf'+outprefix+'lnN*')
    oldfiles_spectrans = glob.glob('./ttf'+outprefix+'lnN*')
    oldfiles = oldfiles_fitcoords + oldfiles_sdisttrans + oldfiles_spectrans

    if len(oldfiles) > 0:
        if overwrite:
            logger.warning("Removing old f%slnN*.fits, tf%slnN*.fits, ttf%slnN*.fits files.", outprefix,
                outprefix, outprefix)
            for f in oldfiles:
                os.remove(f)
        else:
            logger.warning("Output files exist and -overwrite not set - skipping spatial distortion correction and ")
            logger.warning("spectral transformation for all observations.")
            return

    logger.info("Editing the primary header AND science extensions.")
    for obs in obslist:

        # First, make copies of the reduced observations with new filenames and edit header to add the
        # nsfitcoords stamp.
        iraf.copy(input=outprefix+'ln'+obs, output='f'+outprefix+'ln'+obs, verbose='yes', mode="ql")

        iraf.hedit(
            images='f'+outprefix+'ln'+obs+'[0]', fields="NSFITCOO", value="copied from arc",
            add='yes', addonly='no', delete='no', verify='no', show='no', update='yes', mode='al')

        # Second, add sdist info and transform the files to make the orders vertical.
        for i in range(sdistfileslength):
            iraf.hedit(
                images='f'+outprefix+'ln'+obs+'[SCI,'+str(i+1)+']', field="FCFIT2",
                value='f'+combinedarc[:-5]+'_SCI_'+str(i+1)+'_sdist', add='yes', addonly='no', delete='no',
                verify='no', show='no', update='yes', mode='al')

        iraf.nstransform(
            inimages='f'+outprefix+'ln'+obs, outspectra='', outprefix='t', dispaxis=1,
            database=databasepath, fl_stripe='no', interptype='poly3', xlog='no', ylog='no', pixscale=1.0,
            logfile=logger.root.handlers[0].baseFilename, verbose='yes', debug='no', mode='al')

        # Third, add lamp info and transform the files to make skylines horizontal (spatial direction along
        # detector rows).
        for i in range(wavecallampfileslength):
            iraf.hedit(
                images='tf'+outprefix+'ln'+obs+'[SCI,'+str(i+1)+']', field="FCFIT1",
                value='ftf'+combinedarc[:-5]+'_SCI_'+str(i+1)+'_lamp', add='yes', addonly='no', delete='no',
                verify='no', show='no', update='yes', mode='al')

        iraf.nstransform(
            inimages='tf'+outprefix+'ln'+obs, outspectra='', outprefix='t', dispaxis=1,
            database=databasepath, fl_stripe='no', interptype='poly3', xlog='no', ylog='no', pixscale=1.0,
            logfile=logger.root.handlers[0].baseFilename, verbose='yes', debug='no', mode='al')

    return

#---------------------------------------------------------------------------------------------------------------------#

if __name__ == '__main__':
    log.configure('gnirs.log', filelevel='INFO', screenlevel='DEBUG')
    a = raw_input('Enter <Science> for science reduction or <Telluric> for telluric reduction: ')
    start(a, 'gnirs.cfg')
