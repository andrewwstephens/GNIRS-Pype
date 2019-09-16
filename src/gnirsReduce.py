#!/usr/bin/env python

import log
import os
import glob
import ConfigParser
import datetime
import dateutil.parser
from astropy.io import fits
from pyraf import iraf
import utils


def start(kind, configfile):
    """
    This module contains all the functions needed to perform the full reduction of SCIENCE or TELLURIC data.

    Parameters are loaded from gnirs-pype.cfg configuration file.

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
        - configfile: gnirs-pype.cfg configuration file.
                - Paths to the Science (str), reduction truth value (boolean)
                  E.g. 'target/date/config/{Sci,Tel}_ObsID/{Calibrations,Intermediate}', True
                - Paths to the Tellurics (str), reduction truth value (boolean)
                  E.g. 'target/date/config/{Sci,Tel}_ObsID/{Calibrations,Intermediate}', True
                - manualMode (boolean): Enable optional manualModeging pauses? Default: False
                - overwrite (boolean): Overwrite old files? Default: False
                # And gnirsReduce specific settings
    """
    logger = log.getLogger('gnirsReduce')

    path = os.getcwd()  # Store current working directory for later use.

    logger.info(' --------------------------------------------- ')
    logger.info('| Starting the Science and Telluric Reduction |')
    logger.info(' --------------------------------------------- ')

    # Set up/prepare IRAF.
    iraf.gemini()
    iraf.gemtools()
    iraf.gnirs()
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

    if kind == 'Science':  # scienceReduction specific config
        observationSection = 'ScienceDirectories'
        startstep = config.getint('scienceReduction','Start')
        stopstep = config.getint('scienceReduction','Stop')
        cleanir = config.getboolean('scienceReduction','cleanir')
        radiationCorrectionMethod = config.get('scienceReduction','radiationCorrectionMethod')
        radiationThreshold = config.getfloat('scienceReduction','radiationThreshold')

    elif kind == 'Telluric':  # telluricReduction specific config
        observationSection = 'TelluricDirectories'
        startstep = config.getint('telluricReduction','Start')
        stopstep = config.getint('telluricReduction','Stop')
        cleanir = config.getboolean('telluricReduction','cleanir')
        radiationCorrectionMethod = config.get('telluricReduction','radiationCorrectionMethod')
        radiationThreshold = config.getfloat('telluricReduction','radiationThreshold')

    else:
        logger.error('Invalid kind of reduction: %s. Please enter either Science or Telluric.', kind)
        raise SystemExit

    nsprepareInter = config.getboolean('interactive', 'nsprepareInter')
    calculateSNR = config.getboolean('gnirsPipeline', 'CalculateSNR')
    overwrite = config.getboolean('defaults', 'overwrite')
    manualMode = config.getboolean('defaults', 'manualMode')

    # gnirsReduce will reduce observations in each science or telluric directory only if the reduction truth value for
    # that directory is True; else, it will skip the reductions (of science or telluric frames) in that dorectory. If 
    # True, the script will check if all required calibrations to reduce the science or telluric frames are available
    # in the calibrations directory -- if they are not, it will warn the user and request to provide them.

    # Loop through all the observation (telluric or science) directories to perform a reduction on each one.
    for obspath in config.options(observationSection):

        if not config.getboolean(observationSection, obspath):
            logger.debug('Skipping %s', obspath)
            continue

        iraf.chdir(obspath + '/Intermediate')
        logger.info("Working in %s", obspath)

        calpath = '../Calibrations'
        logger.info("Path to calibrations: %s", calpath)

        if calculateSNR:
            # Check if there is a list of sky images in obspath
            skylistfilename = 'sky.list'
            if os.path.exists(skylistfilename):
                skylist = open(skylistfilename, "r").readlines()
                skylist = [filename.strip() for filename in skylist]
            else:
                logger.warning("Parameter 'calculateSNR' is 'True', but a list of sky images not ")
                logger.warning("available in %s . Setting the 'calculateSNR' parameter for the ", obspath)
                logger.warning("the current set of observations to 'False'.\n")
                calculateSNR = False


        logger.debug("Checking if required calibrations available in %s\n", calpath)

        # Check for the reference image to calculate MDF
        QHflatslist = open(calpath+'/QHflats.list', "r").readlines()
        mdfshiftimage = calpath + '/n' + QHflatslist[0].strip()
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

        if startstep > stopstep or startstep < 1 or stopstep > 5:
            logger.error('Invalid start/stop values')
            raise SystemExit

        for valindex in range(startstep, stopstep + 1):
            logger.debug('valindex = %d', valindex)

            if valindex == 1:
                logger.info(" -------------------------------- ")
                logger.info("| STEP 1: Clean raw observations |")
                logger.info(" -------------------------------- ")

                utils.pause(manualMode)

                if cleanir:
                    utils.clean('all.list', 'c', overwrite)
                else:
                    logger.warning("Observations not cleaned.")

            elif valindex == 2:
                logger.info(" ------------------------------- ")
                logger.info("| STEP 2: Prepare observations |")
                logger.info(" ------------------------------- ")

                utils.pause(manualMode)

                prepareObservations(nsprepareInter, mdfshiftimage, overwrite)

            elif valindex == 3:
                logger.info(" ------------------------------------------------------------- ")
                logger.info("| STEP 3: Create bad pixel masks and correct radiation events |")
                logger.info(" ------------------------------------------------------------- ")

                utils.pause(manualMode)

                if radiationCorrectionMethod not in ['fixpix', 'dqplane']:
                    logger.warning("Invalid/no radiation event correction method specified.")
                    header = fits.open('n' + utils.files_in(['all.list'])[0])[0].header
                    date = dateutil.parser.parse(header['DATE-OBS'].strip() + ' ' + header['TIME-OBS'])
                    if date < datetime.datetime(year=2012, month=8, day=1, hour=0, minute=0, second=0):
                        logger.info("Defaulting to the 'fixpix' method.")
                        radiationCorrectionMethod = 'fixpix'
                    else:
                        logger.info("Observations after August 2012 do not have radiation events; skipping correction.")
                        radiationCorrectionMethod = None

                if radiationCorrectionMethod == 'fixpix':

                    createMinimumImage('nodA.list', 'min_nodA.fits', overwrite)
                    createMinimumImage('nodB.list', 'min_nodB.fits', overwrite)
                    radiationCorrectionFixpix(radiationThreshold, overwrite)

                elif radiationCorrectionMethod == 'dqplane':

                    createMinimumImage('nodA.list', 'min_nodA.fits', overwrite)
                    createMinimumImage('nodB.list', 'min_nodB.fits', overwrite)
                    radiationCorrectionDQplane('nodA.list', 'min_nodA.fits', overwrite)
                    radiationCorrectionDQplane('nodB.list', 'min_nodB.fits', overwrite)

            elif valindex == 4:
                logger.info(" ------------------------------------- ")
                logger.info("| STEP 4: Flat field and sky subtract |")
                logger.info(" ------------------------------------- ")

                utils.pause(manualMode)

                # All observations are reduced WITH sky subtraction using an ouput file prefix of 'r'.
                # If the parameter 'calculateSNR' is set to 'yes' in the configuration file, the observations are
                # also reduced WITHOUT sky subtraction (by setting the parameter 'fl_sky' in nsreduce to 'no')
                # using the output file prefix of 'k'

                reduceObservations(masterflat, radiationThreshold, overwrite, skysub=True)

                if calculateSNR:
                    logger.info("Reducing observations without sky subtraction to generate an SNR spectrum later.")
                    reduceObservations(masterflat, radiationThreshold, overwrite, skysub=False)

            elif valindex == 5:
                logger.info(" ------------------------------------------------------------------------- ")
                logger.info("| STEP 5: Apply spatial distortion correction and spectral transformation |")
                logger.info(" ------------------------------------------------------------------------- ")

                utils.pause(manualMode)

                SdistCorrection_SpectralTransform(databasepath,  sdistfileslength, wavecallampfileslength,
                                                  combinedarc, overwrite, prefix='r')
                if calculateSNR:
                    logger.info("Processing science frames reduced without sky subtraction.")
                    SdistCorrection_SpectralTransform(databasepath, sdistfileslength, wavecallampfileslength,
                                                      combinedarc, overwrite, prefix='k')

        logger.info(" -------------------------- ")
        logger.info("| Reduction step complete. |")
        logger.info(" -------------------------- ")

    iraf.chdir(path)  # Return to directory script was begun from.

    return


# ----------------------------------------------------------------------------------------------------------------------
def prepareObservations(interactive, mdfshiftimage, overwrite):
    """
    Prepare raw observations (science or telluric) using nsprepare. Output prefix "n" is added to raw filenames.

    Processing with NSPREPARE will rename the data extension and add variance and data quality extensions. By default
    (see NSHEADERS) the extension names are SCI for science data, VAR for variance, and DQ for data quality (0 = good).
    NSPREPARE will add an MDF file (extension MDF) describing the GNIRS image pattern and how all images
    cross-coorelate with the first prepared QHflat provided as the "mdfshiftimage" reference with the "shiftimage"
    parameter in NSPREPARE.

    :param interactive:
    :param mdfshiftimage:
    :param overwrite:
    :return:
    """
    logger = log.getLogger('prepareObservations')

    infiles = utils.files_in(['all.list'])
    utils.requires(infiles)
    outfiles = ['n' + f for f in infiles]
    if utils.exists(outfiles, overwrite):
        logger.info('All files already prepared.')
        return

    bpm = utils.get_bpm(infiles[0])

    iraf.nsprepare(
        inimages='@all.list', rawpath='', outimages='', outprefix='n', bpm=bpm,
        logfile=logger.root.handlers[0].baseFilename, fl_vardq='yes', fl_cravg='no', crradius=0.0,
        fl_dark_mdf='no', fl_correct='no', fl_saturated='yes', fl_nonlinear='yes', fl_checkwcs='yes',
        fl_forcewcs='yes', arraytable='gnirs$data/array.fits', configtable='gnirs$data/config.fits',
        specsec='[*,*]', offsetsec='none', pixscale='0.15', shiftimage=mdfshiftimage, shiftx=0., shifty=0.,
        obstype='FLAT', fl_inter=interactive, verbose='yes', mode='al')

    return


# ----------------------------------------------------------------------------------------------------------------------
def createMinimumImage(filelist, minimage, overwrite):
    """
    Combine images using GEMCOMBINE to create a 'minimum' image using the "minmax" combining algorithm.
    In this pipeline, the function is used to combine observations at different nod positions positions, separately,
    to create 'minimum' images at each position.

    :param filelist:
    :param minimage:
    :param overwrite:
    :return:
    """
    logger = log.getLogger('createMinimumImage')

    infiles = ['n' + f for f in utils.files_in([filelist])]
    utils.requires(infiles)

    outfiles = [minimage]
    if utils.exists(outfiles, overwrite):
        logger.info('Minimum image already exists.')
        return

    logger.info("Creating a 'minimum' image...")
    iraf.gemcombine(
        input=','.join(infiles), output=minimage, title="", combine="median", reject="minmax",
        offsets="none", masktype="goodvalue", maskvalue=0., scale="none", zero="none", weight="none",
        statsec="[*,*]", expname="EXPTIME", lthreshold='INDEF', hthreshold='INDEF', nlow=0,
        nhigh=len(infiles)-1, nkeep=1, mclip='yes', lsigma=3., hsigma=3., key_ron="RDNOISE", key_gain="GAIN",
        ron=0.0, gain=1.0, snoise="0.0", sigscale=0.1, pclip=-0.5, grow=0.0, bpmfile="", nrejfile="",
        sci_ext="SCI", var_ext="VAR", dq_ext="DQ", fl_vardq='yes',
        logfile=logger.root.handlers[0].baseFilename, fl_dqprop='no', verbose='yes')

    return


# ----------------------------------------------------------------------------------------------------------------------
def radiationCorrectionFixpix(radiationThreshold, overwrite):
    """
    Prepare bad pixel masks for all observations (science or telluric) using a radiation threshold value and 
    correct the bad pixels by interpolation method.

    Here, the task proto.fixpix is set to interpolate along lines to avoid bad effects when a radiation event is right
    on a skyline. Column interpolation may be better for events landing on the peak of target spectra, but effects from
    bad interpolation over skylines seem to do worse.
    """
    logger = log.getLogger('radiationCorrectionFixpix')

    nodAlistfilename = 'nodA.list'
    nodAlist = utils.files_in([nodAlistfilename])
    nodBlistfilename = 'nodB.list'
    nodBlist = utils.files_in([nodBlistfilename])
    nodlist = [nodAlist, nodBlist]

    min_nodA = 'min_nodA.fits'
    min_nodB = 'min_nodB.fits'
    minimagelist = [min_nodA, min_nodB]

    files = nodAlist + nodBlist
    infiles = ['n' + f for f in files] + minimagelist
    utils.requires(infiles)

    outfiles = ['mask' + f.replace('.fits', '.pl') for f in files] + \
               ['gmask' + f.replace('.fits', '.pl') for f in files] + \
               ['bgmask' + f.replace('.fits', '.pl') for f in files] + \
               ['ln' + f for f in files]

    if utils.exists(outfiles, overwrite):
        logger.info('Radiation corrected files exist.  Skipping this step.')
        return

    header = fits.open(infiles[0])[0].header
    rdnoise = header['RDNOISE']
    gain = header['GAIN']
    logger.debug('RDNOISE: %s', rdnoise)
    logger.debug('GAIN: %s', gain)

    logger.info("Creating initial masks for all observations.")
    for i in range(len(nodlist)):
        for image in nodlist[i]:
            iraf.imexpr(
                "(a-b)>" + str(radiationThreshold) + "*sqrt(" + str(rdnoise) + "**2+2*b/" + str(gain) + ") ? 1 : 0",
                output='mask' + image[:-5] + '.pl',
                a='n' + image + '[SCI]',
                b=minimagelist[i] + '[SCI]',
                dims="auto", intype="int", outtype="auto", refim="auto", bwidth=0, btype="nearest", bpixval=0.,
                rangecheck='yes', verbose='yes', exprdb="none")

    logger.debug('Creating list of masks...')
    with open('masks.list', 'w') as masklist:
        for f in files:
            masklist.write('mask%s\n' % f.replace('.fits', '.pl'))

    logger.debug('Growing masks...')
    iraf.crgrow(input='@masks.list', output='g//@masks.list', radius=1.5, inval="INDEF", outval="INDEF")

    logger.info("Adding bad pixels from the [DQ] plane generated by nsprepare to the masks so that they can be ")
    logger.info("corrected as well.\n")

    for i in range(len(nodlist)):
        for image in nodlist[i]:
            iraf.copy(input='n' + image, output='ln' + image, verbose='yes')

            iraf.imexpr(
                expr="a||b",
                output='bgmask' + image[:-5] + '.pl',
                a='gmask'+image[:-5] + '.pl',
                b=minimagelist[i] + '[DQ]',
                dims="auto", intype="int", outtype="auto", refim="auto", bwidth=0,
                btype="nearest", bpixval=0., rangecheck='yes', verbose='yes', exprdb="none")

            iraf.proto.fixpix(
                images='ln' + image + '[SCI,1]', masks='bgmask' + image[:-5] + '.pl', linterp=1,
                cinterp="INDEF", verbose='yes', pixels='no')

    return


# ----------------------------------------------------------------------------------------------------------------------
def radiationCorrectionDQplane(inlistfile, minimage, overwrite):
    """
    Prepare bad pixel masks for all observations (science or telluric) using the data quality plane of the prepared 
    observations and correct the bad pixels by replacing the corresponding pixels from 'minimum' images.
    """
    logger = log.getLogger('radiationCorrectionDQplane')

    files = utils.files_in([inlistfile])
    infiles = ['n' + f for f in files]
    utils.requires(infiles)

    outfiles = ['mask' + f.replace('.fits', '.pl') for f in files] + \
               ['gmask' + f.replace('.fits', '.pl') for f in files] + \
               ['ln' + f for f in files]

    if utils.exists(outfiles, overwrite):
        logger.info('Radiation corrected files exist.  Skipping this step.')
        return

    logger.info("Copying the [DQ] planes generated by nsprepare as the initial masks for all observations.")
    for f in files:
        iraf.imcopy(input='n%s[DQ,1]' % f, output='mask%s' % f.replace('.fits', '.pl'), verbose='yes')

    logger.debug('Creating list of masks...')
    with open('masks.list', 'w') as masklist:
        for f in files:
            masklist.write('mask%s\n' % f.replace('.fits', '.pl'))

    logger.debug('Growing masks...')
    iraf.crgrow(input='@masks.list', output='g//@masks.list', radius=1.5, inval="INDEF", outval="INDEF")

    logger.info("Replacing the pixels identified as bad pixels by the final masks with pixels from the ")
    logger.info("'minimum' images.\n")
    for image in files:
        iraf.copy(input='n' + image, output='ln' + image, verbose='yes')
        iraf.imexpr(
            expr="c>1 ? b : a",
            output='ln'+image+'[SCI,overwrite]',
            a='n' + image + '[SCI]',
            b=minimage + '[SCI]',
            c='gmask'+image[:-5]+'.pl',
            dims="auto", intype="int", outtype="auto", refim="auto", bwidth=0,
            btype="nearest", bpixval=0., rangecheck='yes', verbose='yes', exprdb="none")

    return


# ----------------------------------------------------------------------------------------------------------------------
def reduceObservations(flat, radiationThreshold, overwrite, skysub=True):
    """
    Flat field and sky subtract observations with nsreduce. If 'calculateSNR' is set, observations will be
    reduced without sky subtraction too and prefix 'k' will be added before the input filenames; otherwise, sky 
    subtraction will be performed and the prefix will be 'r'.

    NSREDUCE is used for basic reduction of raw data - it provides a single, unified interface to several tasks and 
    also allows the subtraction of dark frames (not done for GNIRS data) and dividing by the flat.

    :param outprefix:
    :param masterflat:
    :param radiationThreshold:
    :param overwrite:
    :param skysub:
    :return:
    """
    logger = log.getLogger('reduceObservations')

    if skysub:
        outprefix = 'r'
    else:
        outprefix = 'k'

    files = utils.files_in(['all.list'])
    infiles = ['ln' + f for f in files]
    utils.requires(infiles)

    outfiles = [outprefix + f for f in infiles]
    if utils.exists(outfiles, overwrite):
        logger.info('Reduced files exist.  Skipping this step.')
        return

    iraf.nsreduce(
        inimages='ln//@all.list', outimages='', outprefix=outprefix, fl_cut='yes', section="",
        fl_corner='yes', fl_process_cut='yes', fl_nsappwave='no', nsappwavedb='gnirs$data/nsappwave.fits',
        crval="INDEF", cdelt="INDEF", fl_dark='no', darkimage="", fl_save_dark='no', fl_sky=skysub,
        skyimages="", skysection="", combtype="median", rejtype="avsigclip", masktype="goodvalue", maskvalue=0.,
        scale="none", zero="median", weight="none", statsec="[*,*]", lthreshold="INDEF", hthreshold="INDEF",
        nlow=1, nhigh=1, nkeep=0, mclip='yes', lsigma=3., hsigma=3., snoise="0.0", sigscale=0.1, pclip=-0.5,
        grow=0.0, skyrange='INDEF', nodsize=3., fl_flat='yes', flatimage=flat, flatmin=0.0, fl_vardq='yes',
        logfile=logger.root.handlers[0].baseFilename, verbose='yes', debug='no', force='no')

    logger.debug('Recording the radiation threshold value in reduced file headers...')
    iraf.hedit(
        images=outprefix + 'ln//@all.list',
        fields='RTHRESH',
        value=radiationThreshold,
        add='yes', addonly='no', delete='no', verify='no', show='no', update='yes')

    return


# ----------------------------------------------------------------------------------------------------------------------
def SdistCorrection_SpectralTransform(databasepath, sdistfileslength, wavecallampfileslength,
                                      combinedarc, overwrite, prefix='r'):
    """
    Apply spatial distortion correction and spectral transformation. 

    In general, we would run nsfitcoords on all reduced science or telluric frames with and without sky subtraction (if 
    creating an error spectrum) to straighten the spectra and then de-tilt. nsfitcoords simply writes database info, 
    which depends only on pinholes and arcs (the *_sdist and *_lamp files that were written to the database while
    reducing baseline calibrations), not on on-sky data. So here we access the nsfitcoords solution from the
    calibrations folder and append the relevant header keywords, which will tell nstransform the database files it 
    should use, without actually calling nsfitcoords. This routine also uses the output files in the database that were
    created by the nswavelength task.
    """
    logger = log.getLogger('SdistCorrection_SpectralTransform')

    # TODO: check that the database files exist in Calibrations/database/

    files = utils.files_in(['all.list'])
    infiles = [prefix + 'ln' + f for f in files]
    utils.requires(infiles)

    outfiles = ['f' + f for f in infiles] + ['tf' + f for f in infiles] + ['ttf' + f for f in infiles]
    if utils.exists(outfiles, overwrite):
        logger.info('Reduced files exist.  Skipping this step.')
        return

    logger.info("Editing the primary header AND science extensions.")
    for obs in files:

        # First, make copies of the reduced observations with new filenames and edit header to add the
        # nsfitcoords stamp.
        iraf.copy(
            input=prefix + 'ln' + obs, output='f' + prefix + 'ln' + obs, verbose='yes', mode="ql")

        iraf.hedit(
            images='f' + prefix + 'ln' + obs + '[0]', fields="NSFITCOO", value="copied from arc",
            add='yes', addonly='no', delete='no', verify='no', show='no', update='yes', mode='al')

        # Second, add sdist info and transform the files to make the orders vertical.
        for i in range(sdistfileslength):
            iraf.hedit(
                images='f' + prefix + 'ln' + obs + '[SCI,' + str(i+1) + ']',
                field="FCFIT2",
                value='f' + combinedarc[:-5] + '_SCI_' + str(i+1) + '_sdist',
                add='yes', addonly='no', delete='no', verify='no', show='no', update='yes', mode='al')

        iraf.nstransform(
            inimages='f' + prefix + 'ln' + obs, outspectra='', outprefix='t', dispaxis=1,
            database=databasepath, fl_stripe='no', interptype='poly3', xlog='no', ylog='no', pixscale=1.0,
            logfile=logger.root.handlers[0].baseFilename, verbose='yes', debug='no', mode='al')

        # Third, add lamp info and transform the files to make skylines horizontal (spatial direction along
        # detector rows).
        for i in range(wavecallampfileslength):
            iraf.hedit(
                images='tf' + prefix + 'ln' + obs + '[SCI,' + str(i+1) + ']',
                field="FCFIT1",
                value='ftf' + combinedarc[:-5] + '_SCI_' + str(i+1) + '_lamp',
                add='yes', addonly='no', delete='no', verify='no', show='no', update='yes', mode='al')

        iraf.nstransform(
            inimages='tf'+prefix+'ln'+obs, outspectra='', outprefix='t', dispaxis=1,
            database=databasepath, fl_stripe='no', interptype='poly3', xlog='no', ylog='no', pixscale=1.0,
            logfile=logger.root.handlers[0].baseFilename, verbose='yes', debug='no', mode='al')

    return


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    log.configure('gnirs-pype.log', filelevel='INFO', screenlevel='DEBUG')
    kind = raw_input('Enter <Science> for science reduction or <Telluric> for telluric reduction: ')
    start(kind, 'gnirs-pype.cfg')
