#!/usr/bin/env python

from astropy.io import fits
import ConfigParser
import glob
import log
from pyraf import iraf
import os
import utils


# ----------------------------------------------------------------------------------------------------------------------
def start(configfile):
    """
    This module combines reduced 2D science and/or telluric spectra.

    Parameters are loaded from gnirs-pype.cfg configuration file.

    Args:
        - configfile: gnirs-pype.cfg configuration file.
                - Paths to the Science (str), reduction truth value (boolean)
                  E.g. 'target/date/config/{Sci,Tel}_ObsID/{Calibrations,Intermediate}', True
                - Paths to the Tellurics (str), reduction truth value (boolean)
                  E.g. 'target/date/config/{Sci,Tel}_ObsID/{Calibrations,Intermediate}', True
                - manualMode (boolean): Enable optional manualModeging pauses? Default: False
                - overwrite (boolean): Overwrite old files? Default: False
                # And gnirsCombineSpectra2D specific settings
    """
    logger = log.getLogger('CombineSpectra2D')

    path = os.getcwd()  # Store current working directory for later use.

    # Set up/prepare IRAF.
    iraf.gemini()
    iraf.gemtools()
    iraf.gnirs()
    iraf.unlearn(iraf.gemini, iraf.gemtools, iraf.gnirs, iraf.imcopy)  # Reset parameters to the default values

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

    manualMode = config.getboolean('defaults', 'manualMode')
    overwrite = config.getboolean('defaults', 'overwrite')
    nscombineInter = config.getboolean('interactive', 'nscombineInter')
    calculateSNR = config.getboolean('gnirsPipeline', 'CalculateSNR')
    preparedPrefix = config.get('runtimeFilenames', 'preparedPrefix')
    reducedPrefix = config.get('runtimeFilenames', 'reducedPrefix')
    fitcoordsPrefix = config.get('runtimeFilenames', 'fitcoordsPrefix')
    transformPrefix = config.get('runtimeFilenames', 'transformPrefix')
    radiationCorrectedPrefix = config.get('runtimeFilenames', 'radiationCorrectedPrefix')
    noskysubReducedPrefix = config.get('runtimeFilenames', 'noskysubReducedPrefix')
    combinedsrc = config.get('runtimeFilenames', 'combinedsrc')
    combinedsky = config.get('runtimeFilenames', 'combinedsky')

    # gnirsCombineSpectra2D will first check if the reduction truth value of the science and telluric directories is 
    # True -- if it is, it will then check if the required spectra to be combined are available in the directories (and
    # proceed only if it finds them there); else, it will warn the user and request to provide the files for combining.
    # If the reduction truth value of the science and telluric directories is False, the script will skip combining 2D
    # spectra in those directories.

    # Loop through all the observation directories (telluric and science) to combine 2D spectra in each
    for section in ['ScienceDirectories', 'TelluricDirectories']:
        for obspath in config.options(section):

            if not config.getboolean(section, obspath):  # Only process directories marked True
                logger.debug('Skipping %s', obspath)
                continue

            logger.info(' ---------------------------- ')
            logger.info('| Combining GNIRS 2D Spectra |')
            logger.info(' ---------------------------- ')

            obspath += '/Intermediate'
            iraf.chdir(obspath)
            logger.info("%s", obspath)

            logger.debug('Checking for lists...')
            utils.requires(['src.list', 'nodA.list', 'nodB.list'])

            if calculateSNR:
                if not utils.exists(['sky.list'], overwrite=False):
                    logger.warning('Could not find list of sky spectra.  Setting calculateSNR = False')
                    calculateSNR = False
                if len(utils.files_in(['sky.list'])) == 0:
                    logger.warning('sky.list is empty.  Setting calculateSNR = False')
                    calculateSNR = False

            utils.pause(manualMode, 'About to start combining 2D spectra')

            prefix = transformPrefix + transformPrefix + fitcoordsPrefix + reducedPrefix + radiationCorrectedPrefix + preparedPrefix
            nscombine('src.list', nscombineInter, prefix, combinedsrc, overwrite, cross=True)

            if shiftsMatch('nodA.list', 'nodB.list', combinedsrc, preparedPrefix):
                logger.info('The nscombine shifts look good with fl_cross = yes.')

            else:
                logger.warning('The nscombine shifts look wrong; trying with fl_cross=no...')
                nscombine('src.list', nscombineInter, reducedPrefix, combinedsrc, overwrite, cross=False)

                if shiftsMatch('nodA.list', 'nodB.list', combinedsrc, preparedPrefix):
                    logger.info('The nscombine shifts look okay with fl_cross = no.')

                else:
                    logger.error('The nscombine shifts still look wrong.')
                    logger.error('I have tried with fl_cross=yes and fl_cross=no.')
                    logger.error('This could be a bug in nscombine.')
                    raise SystemExit

            # If the parameter 'calculateSNR' is set to 'yes', the script will combine all observations
            # reduced without sky subtraction (having reduce_outputPrefix 'k'); else, all observations reduced
            # with sky subtraction only (having reduce_outputPrefix 'r') will be combined.
            if calculateSNR:
                prefix = transformPrefix + transformPrefix + fitcoordsPrefix + noskysubReducedPrefix + radiationCorrectedPrefix + preparedPrefix
                logger.info('Combining the sky observations')
                nscombine('sky.list', nscombineInter, prefix, combinedsky, overwrite, cross=False)

            logger.info(' --------------------------------- ')
            logger.info('| Done combining GNIRS 2D Spectra |')
            logger.info(' --------------------------------- ')

    iraf.chdir(path)  # Return to directory script was begun from.

    return


# ----------------------------------------------------------------------------------------------------------------------
def nscombine(inlist, interact, prefix, outfile, overwrite, cross=True):
    """
    Combining the transformed science or telluric frames.
    """
    logger = log.getLogger('nscombine')
    logger.debug('inlist: %s', inlist)
    logger.debug('outfile: %s', outfile)
    logger.debug('cross: %s', cross)

    files = utils.files_in([inlist])
    infiles = [prefix + f for f in files]
    utils.requires(infiles)

    outfiles = [outfile]
    if utils.exists(outfiles, overwrite):
        logger.info('Files already combined.')
        return

    iraf.nscombine(
        inimages=prefix + '//@' + inlist, tolerance=0.5, output=outfile, output_suffix='', bpm="",
        dispaxis=1, pixscale=1., fl_cross=cross, fl_keepshift='no', fl_shiftint='yes', interptype="linear",
        boundary="nearest", constant=0., combtype="average", rejtype="none", masktype="none", maskvalue=0., 
        statsec="[*,*]", scale="none", zero="none", weight="none", lthreshold="INDEF", hthreshold="INDEF", nlow=1, 
        nhigh=1, nkeep=0, mclip='yes', lsigma=5., hsigma=5., ron=0.0, gain=1.0, snoise="0.0", sigscale=0.1, pclip=-0.5,
        grow=0.0, nrejfile='', fl_vardq='yes', fl_inter=interact, logfile=logger.root.handlers[0].baseFilename,
        verbose='yes', debug='no', force='no')

    return


# ----------------------------------------------------------------------------------------------------------------------
def shiftsMatch(nodAlist, nodBlist, combinedimage, preparedPrefix):
    """
    There is currently (2015) a bug in nscombine that miscalculates the shifts at certain position angles while 
    cross-correlating different frames. Unfortunately, we cannot provide a shift to nscombine directly. 
    
    If targets are bright enough, the task avoids the bug and the cross-correlation happens smoothly (set with 
    fl_cross='yes'); else, the task should be run with fl_cross='no'. XDGNIRS always uses fl_cross='yes' for telluric
    frames. For science frames, XDGNIRS checks whether the shifts used by nscombine with fl_cross='yes' were 
    reasonable by comparing with the shifts expected from the headers; else, it runs nscombine with fl_cross = 'no'.
    """
    logger = log.getLogger('ShiftsMatch')
    logger.debug('Starting...')

    utils.requires([combinedimage])

    # TODO:  Currently the nod Q offsets are only calculated for the first set of images in the different nod lists.
    # This assumes that the Q offsets for the remaining set of images in the lists give the same nod value. A
    # robust way of doing this will be calculating the shifts between every two consecutive image sets from the nod 
    # lists and checking if the shift is within the acceptable range; and warn the user if it is not.

    nodAs = utils.files_in([nodAlist])
    nodBs = utils.files_in([nodBlist])
    nodAheader = fits.open(preparedPrefix + nodAs[0])[0].header
    nodBheader = fits.open(preparedPrefix + nodBs[0])[0].header
    nodQoffset = abs(nodAheader['QOFFSET'] - nodBheader['QOFFSET'])
    pixelscale = nodAheader['PIXSCALE']
    nodQpixels = nodQoffset/pixelscale
    logger.debug('nodQpixels: %s', nodQpixels)

    shifts = iraf.hselect(images=combinedimage+'[0]', fields='NSCHLX*', expr='yes', missing='INDEF', Stdout=1)
    logger.debug('shifts: %s', shifts)
    shifts = [s.replace('\t', ',') for s in shifts]
    logger.debug('shifts: %s', shifts)

    match = True
    for shift in shifts[0].split(','):
        logger.debug('shift: %s', shift)
        if shift != '0.00' and abs(float(shift) - nodQpixels) > 1.5:
            logger.warning("Shift of %s pixels in %s (NSCHLX* keyword).", shift, combinedimage)
            logger.warning("This does not match the expected value of %s pixels.", nodQpixels)
            logger.warning("The target may not have been bright enough to use nscombine with fl_cross='yes'.")
            match = match and False
            break
        else:
            logger.info("Shift is %s pixels in %s (NSCHLX* keyword).", shift, combinedimage)
            logger.info("This is 0.0 or close to the expected value of %s pixels.", nodQpixels)
            match = match and True
            continue

    logger.debug('match: %s', match)
    return match


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    log.configure('gnirs-pype.log', filelevel='INFO', screenlevel='DEBUG')
    start('gnirs-pype.cfg')
