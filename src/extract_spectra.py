#!/usr/bin/env python

from astropy.io import fits
import ConfigParser
import log
import obslog
import os
from pyraf import iraf
import shutil
import utils


# ----------------------------------------------------------------------------------------------------------------------
def start(configfile):
    """
    Parameters are loaded from gnirs-pype.cfg configuration file. This script will automatically detect if it is being run
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
        - configfile: gnirs-pype.cfg configuration file.
                - Paths to the Science (str), reduction truth value (boolean)
                  E.g. 'target/date/config/{Sci,Tel}_ObsID/{Calibrations,Intermediate}', True
                - Paths to the Tellurics (str), reduction truth value (boolean)
                  E.g. 'target/date/config/{Sci,Tel}_ObsID/{Calibrations,Intermediate}', True
                - manualMode (boolean): Enable optional manualModeging pauses? Default: False
                - overwrite (boolean): Overwrite old files? Default: False
                # And gnirsReduce specific settings
    """
    logger = log.getLogger('extract_spectra')

    path = os.getcwd() # Store current working directory for later use.

    # Set up/prepare IRAF.
    iraf.gemini()
    iraf.gemtools()
    iraf.gnirs()
    iraf.unlearn(iraf.gemini, iraf.gemtools, iraf.gnirs, iraf.imcopy)  # reset parameters to default values

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

    # Order of sections is important to later check for plausible peaks located for science targets by nsextract
    nsextractInter = config.getboolean('interactive', 'nsextractInter')
    calculateSNR = config.getboolean('gnirsPipeline', 'CalculateSNR')
    
    combinedsrc = config.get('runtimeFilenames', 'combinedsrc')
    combinedsky = config.get('runtimeFilenames', 'combinedsky')
    extractRegularPrefix = config.get('runtimeFilenames', 'extractRegularPrefix')
    extractFullSlitPrefix = config.get('runtimeFilenames', 'extractFullSlitPrefix')
    extractStepwiseTracePrefix = config.get('runtimeFilenames', 'extractStepwiseTracePrefix')
    extractStepwisePrefix = config.get('runtimeFilenames', 'extractStepwisePrefix')

    # extract1Spectra1D specific config
    useApall = config.getboolean('extractSpectra1D', 'useApall')
    subtractBkg = config.get('extractSpectra1D', 'subtractBkg')
    extractApertureRadius = config.getfloat('extractSpectra1D', 'extractApertureRadius')
    checkPeaksMatch = config.getboolean('extractSpectra1D', 'checkPeaksMatch')
    toleranceOffset = config.getfloat('extractSpectra1D', 'toleranceOffset')
    extractFullSlit = config.getboolean('extractSpectra1D','extractFullSlit')
    extractStepwise = config.getboolean('extractSpectra1D','extractStepwise')
    extractionStepSize = config.getfloat('extractSpectra1D','extractStepSize')
    extractApertureWindow = config.getfloat('extractSpectra1D','extractApertureWindow')

    # gnirsExtractSpectra1D will first check if the reduction truth value of the science and telluric directories is 
    # True -- if it is, it will then check if the required spectra to be extracted are available in the directories
    # (and proceed only if it finds them there); else, it will warn the user and request to provide the spectra for 
    # extracting. If the reduction truth value of the science and telluric directories is False, the script will skip
    # extracting 1D spectra in those directories.

    # Loop through all the observation (telluric and science) directories to extract 1D spectra in each one.
    # (the Telluric standards must be done first if they are to be used as a reference)
    for section in ['TelluricDirectories', 'ScienceDirectories']:
        for obspath in config.options(section):

            if not config.getboolean(section, obspath):  # Only process directories marked True
                logger.debug('Skipping extraction of 1D spectra in %s', obspath)
                continue

            logger.info(' ----------------------------------------- ')
            logger.info('| Starting extraction of GNIRS 1D spectra |')
            logger.info(' ----------------------------------------- ')

            obspath += '/Intermediate'
            logger.info("%s\n", obspath)
            iraf.chdir(obspath) 
            utils.pause(manualMode)

            logger.debug("Checking if required combined spectra available.")
            utils.requires([combinedsrc])

            if calculateSNR:
                if not utils.exists([combinedsky], overwrite=False):
                    logger.warning('Could not find combined sky spectra.  Setting calculateSNR = False')
                    calculateSNR = False

            orders = utils.get_orders(obspath)
            extractApertureWindow = get_window(obspath)

            if nsextractInter:
                subtractBkg = 'fit'
                logger.info('Setting background subtraction method to "fit"')
            
            if useApall:
                nsum = 20
            else:
                nsum = 10

            extractSpectra(combinedsrc, extractRegularPrefix, nsextractInter, useApall, nsum, subtractBkg,
                               extractApertureRadius, overwrite)

            if calculateSNR:
                logger.info("Extracting the combined sky spectrum reduced without sky subtraction.")
                subtractBkg = 'none'
                extractSpectra(combinedsky, extractRegularPrefix, nsextractInter, useApall, nsum, subtractBkg,
                                   extractApertureRadius, overwrite)

            if 'Science' in section:
                # If the extraction was not done interactively check if checkPeaksMatch is set: if yes, check if the
                # required telluric extraction reference files available in the telluric /database directory; else,
                # warn the user that both nsextractInter and checkPeaksMatch are not set, request the user to
                # manually check if the science target peak identified by task nsextract might identify a wrong peak
                # if the science target is not bright enough.

                # Get symbolic path to the tel database directory in the sci directory
                # Relative path/link expected to be at the top level of every sci directory
                scidatabasepath = 'database'
                logger.info("Science database path: %s", scidatabasepath)
                telpath = '../Telluric/Intermediate'
                logger.info("Telluric path: %s", telpath)
                teldatabasepath = telpath + '/database'
                logger.info("Telluric database path: %s", teldatabasepath)
                sci_combinedsrc = obspath + '/' + combinedsrc
                tel_combinedsrc = telpath + '/' + combinedsrc

                if not nsextractInter:  # if nsextract was not run interactively

                    if not checkPeaksMatch:
                        logger.warning("Parameters 'nsextractInter' and 'checkPeaksMatch' are both set to False.")
                        logger.warning("Please manually verify that nsextract identified the science peaks correctly.")

                    else:
                        logger.info("Finding extraction locations for Telluric standard...")
                        telpeaks = get_peaks(teldatabasepath)
                        logger.info("Finding extration locations for Science target...")
                        scipeaks = get_peaks(scidatabasepath)

                        utils.pause(True, 'Next is the peak matching!')




                        logger.info("Comparing the science and Telluric extraction locations...")

                        stuff = compare_peaks(obspath, telpath, scipeaks, telpeaks, toleranceOffset)





                        logger.info("Completed matching the peaks found by nsextract.")

                        # Re-extract combined science 2D source spectrum if needed
                        if sciReExtract:
                            logger.info("Re-extracting one or more science spectra.")
                            newtelapfilePrefix = 're'
                            useApall = 'yes'
                            apertureTracingColumns = 20
                            reExtractSpectra1D(scidatabasepath, teldatabasepath, telpath, scitelPeaksMatched,
                                sci_combinedsrc, tel_combinedsrc, newtelapfilePrefix, orders, nsextractInter, extractRegularPrefix,
                                databaseDir, useApall, subtractBkg, apertureTracingColumns, extractApertureRadius)
                            logger.info("Completed re-extracting one or more science spectra.")
                        else:
                            logger.info("Not re-extracting any science spectra.")
                            pass







            # TODO(Viraja):  Set these up once the respective functions are ready.
            if extractFullSlit:
                logger.warning('Extract full slit should be here')

            if extractStepwise:  # Extract in steps on either side of the peak
                useApall = 'yes'
                logger.warning('Step-wise extraction should be here')

        logger.info("Extraction complete for")
        logger.info("%s", obspath)

    iraf.chdir(path)  # Return to directory script was begun from.

    return


# ----------------------------------------------------------------------------------------------------------------------
def get_window(path):
    logger = log.getLogger('get_window')
    # Determine the full-slit extraction window

    if 'LB_SXD' in path:
        # extractApertureRadius = 23 (+/-23 pixels or 6.9" covers almost the entire slit length,
        # but this is only appropriate for objects centred along length of slit (with absolute Q offset of 0).
        # [-46/2,46/2+6)   [-23.0 -17 -11 -5 1 7 13 19 23 29) warn the user if last step in extract >0.1" away
        # from the end of the slit or if extraction proceeding out of the slit
        window = 46

    elif 'LB_LXD' in path:
        window = 33    # [-33/2,33/2+6]  [-16.5 -10.5 -4.5 2.5 8.5 14.5 20.5]

    elif 'SB_SXD' in path:
        window = 46

    else:
        logger.error("Unknown GNIRS XD configuration.")
        raise SystemExit

    logger.debug('Window: %s pix', window)
    return window


# ----------------------------------------------------------------------------------------------------------------------
def extractSpectra(inimage, outprefix, interactive, apall, nsum, background, radius, overwrite):
    """
    Extracting 1D spectra from the combined 2D spectra using nsextract.
    background = Type of background to subtract (none|average|median|minimum|fit)
    """

    # This is really just a wrapper around nsextract.
    # I'm tempted to call this 'nsextract' and call this whenever I need nsextract.
    # I guess it might not have all the parameters, but they could be included as optional.

    logger = log.getLogger('extractSpectra')
    logger.debug('inimage: %s', inimage)
    logger.debug('background: %s', background)
    logger.debug('radius: %s pix', radius)
    logger.debug('nsum: %s pix', nsum)

    utils.requires([inimage])
    orders = utils.get_orders(os.getcwd())
    outfiles = [outprefix + inimage] + \
        ['database/apsrc_comb_DQ_%d_' % i for i in range(1, len(orders)+1)] + \
        ['database/apsrc_comb_SCI_%d_' % i for i in range(1, len(orders)+1)]

    if utils.exists(outfiles, overwrite):
        logger.info('Spectra already extracted.')
        return

    iraf.nsextract(
        inimages=inimage, outspectra='', outprefix=outprefix, dispaxis=1, database='', line=700,
        nsum=nsum, ylevel='INDEF', upper=radius, lower=-radius, background=background,
        fl_vardq='yes', fl_addvar='no', fl_skylines='yes', fl_inter=interactive, fl_apall=apall, fl_trace='no',
        aptable='gnirs$data/apertures.fits', fl_usetabap='no', fl_flipped='yes', fl_project='yes',  fl_findneg='no',
        bgsample='*', trace='', tr_nsum=10, tr_step=10, tr_nlost=3, tr_function='legendre', tr_order=5,
        tr_sample='*', tr_naver=1, tr_niter=0, tr_lowrej=3.0, tr_highrej=3.0, tr_grow=0.0, weights='variance',
        logfile=logger.root.handlers[0].baseFilename, verbose='yes', mode='al')

    return


# ----------------------------------------------------------------------------------------------------------------------
def get_peaks(databasepath):
    """
    Check the extraction reference files in the telluric directory databases to find the location of the peaks.
    """
    logger = log.getLogger('get_peaks')
    logger.debug('databasepath: %s', databasepath)

    orders = utils.get_orders(os.getcwd())
    infiles = ['%s/apsrc_comb_SCI_%d_' % (databasepath, i) for i in range(1, len(orders)+1)]
    utils.requires(infiles)

    peaks = []
    for i in range(1, len(orders)+1):
        apfile = '%s/apsrc_comb_SCI_%d_' % (databasepath, i)
        with open(apfile, 'r') as f:
            p = None
            for line in f:
                # The peak location is the number in the second column of the line beginning with 'center'
                if 'center' in line:
                    p = line.split()[1]
                    break
            peaks.append(p)
            if p is None:
                logger.warning('Peak not found')

    logger.debug('peaks: %s', peaks)
    return peaks


# ----------------------------------------------------------------------------------------------------------------------
def compare_peaks(scipath, telpath, scipeaks, telpeaks, tolerance):
    """
    Compare locations of the extraction locations of the science and Telluric standard.

    For faint targets, NSEXTRACT may find a noise peak instead of the science peak. In such cases, it is advisable
    to check the aperture center of extraction of the science with respect to the telluric and re-extract at the
    expected location.
    Look the science and telluric absolute Q offsets and determine if the relative location of the target
    peak was correct.  If not, we should re-extract at the expected location.
    """
    logger = log.getLogger('gnirsReduce.peaksMatch')

    # Get the absolute P,Q offsets from the obslog:
    sciinfo = obslog.readcsv(scipath + '/obslog.csv')
    telinfo = obslog.readcsv(telpath + '/obslog.csv')

    sciA = utils.files_in([scipath + '/nodA.list'])
    sciB = utils.files_in([scipath + '/nodB.list'])

    telA = utils.files_in([telpath + '/nodA.list'])
    telB = utils.files_in([telpath + '/nodB.list'])

    # Assuming that both the science and Telluric were acquired in the center of the slit, the extraction locations
    # should be the same minus any difference in the Q-offset.

    logger.debug('The telluric used offsets of: %s', telinfo[telA[0]]['Q'])
    logger.debug('The science used offsets of: %s', sciinfo[sciA[0]]['Q'])

    for f in sciA:
        logger.debug('sciA: %s', sciinfo[f]['Q'])

    for f in sciB:
        logger.debug('sciB: %s', sciinfo[f]['Q'])

    for f in telA:
        logger.debug('telA: %s', telinfo[f]['Q'])

    for f in telB:
        logger.debug('telB: %s', telinfo[f]['Q'])

    # The PIXSCALE should be in the obslog, but it's not yet, so:
    pixscale = 0.15

    # So what do we get for the science positions:

    delta = (float(sciinfo[sciA[0]]['Q']) - float(telinfo[telA[0]]['Q']) ) / pixscale
    logger.debug('delta: %s pix', delta)

    logger.debug('Predicted science extraction locations:')
    for s,t in zip(scipeaks, telpeaks):
        logger.debug('telluric: %s  science: %s  predicted sci: %s', t, s, float(t) + delta)


    utils.pause(True, '* You should probably stop here *')



    # Find absolute Q offsets of the combined science and telluric images from their respective last acquisition images
    sciacqHeader = fits.open(sciacq)[0].header
    scisrccombHeader = fits.open(scisrccomb)[0].header
    scisrccombQoffset = abs(sciacqHeader['QOFFSET'] - scisrccombHeader['QOFFSET'])

    telacqHeader = fits.open(telacq)[0].header
    telsrccombHeader = fits.open(telsrccomb)[0].header
    telsrccombQoffset = abs(telacqHeader['QOFFSET'] - telsrccombHeader['QOFFSET'])

    pixelscale = scisrccombHeader['PIXSCALE']
    pixeldifference = (scisrccombQoffset - telsrccombQoffset)/pixelscale  # units: [pixels]

    # nsextract should find the spectrum within a 'tolerance' pixels of expected location. This depends on how well the
    # observer centred the target along the slit. Here, we use 5 pixels as a reasonable tolerance level. A more robust
    # way would be to use some measure of whether the peak found by nsextract was real, e.g. counts + FWHM. However, 
    # this information is not recorded in database.

    # Is reExtract overwritten with each new order?

    matched = []
    for i in range(len(scipeaks)):
        expectedPeak = float(telpeaks[i]) + pixeldifference
        if scipeaks[i] == 'Peak not found':
            logger.warning("nsextract did not extract anything for extension %d", i)
            logger.warning("Will re-extract forcing the aperture to be at the expected location %.4g", expectedPeak)
            matched.append(False)
            reExtract = True
        else:    
            locatedPeak = float(scipeaks[i])
            if abs(locatedPeak - expectedPeak) < tolerance:
                logger.info("nsextract detected the spectrum close to the expected location for In extension %d,  ", i)
                logger.info("(located = %s vs. expected = %s", locatedPeak, expectedPeak)
                matched.append(True)
                reExtract = False
            else:
                logger.warning("nsextract extracted an unexpected peak location in extension %d", i)
                logger.warning("(located = %s vs. expected = ", locatedPeak, expectedPeak)
                logger.warning("Will re-extract forcing the aperture to be at the expected location.'")
                matched.append(False)
                reExtract = True

    # Is reExtract = not all(peaksMatched)?

    return matched

# ----------------------------------------------------------------------------------------------------------------------

def reExtractSpectra1D(scidatabasepath, teldatabasepath, telpath, scitelPeaksMatched, sci_combinedsrc,
    tel_combinedsrc, newtelapfilePrefix, orders, nsextractInter, extractRegularPrefix, databaseDir, useApall, subtractBkg,
    apertureTracingColumns, extractApertureRadius):
    """
    """
    logger = log.getLogger('gnirsReduce.reExtractSpectra1D')

    logger.info("Creating new aperture files in database.")
    for i in range(len(orders)):
        extension = i+1

        # There is some trouble replacing only the database files for the extracted spectra that were not well centred.
        # So, simply replacing all science database files but using the ones for which nsextract located the peaks at 
        # the right positions.
        oldsciapfile = scidatabasepath+'/ap'+sci_combinedsrc[:-5]+'_SCI_'+str(extension)+'_'
        if os.path.exists(oldsciapfile):
            os.remove(oldsciapfile)
        oldtelapfile = open(teldatabasepath+'ap'+tel_combinedsrc+'_SCI_'+str(extension)+'_', 'r')
        newtelapfile = open(teldatabasepath+'ap'+newtelapfilePrefix+tel_combinedsrc+'_SCI_'+str(extension)+'_', 'w')
        if not peaks_flag[i]:
            # TODO(Viraja):  Check if there is a better way to replace the peak values then how it is done below.
            replacetelapfile  = oldtelapfile.read().replace(telpeaks[i], str(float(telpeaks[i])+pixeldifference)+' ').replace(tel_combinedsrc, newtelapfilePrefix+tel_combinedsrc)
        else:
            replacetelapfile  = oldtelapfile.read().replace(telpeaks[i], telpeaks[i]+' ').replace(tel_combinedsrc, newtelapfilePrefix+tel_combinedsrc)
        newtelapfile.write(replacetelapfile)
        oldtelapfile.close()
        newtelapfile.close()
        
    shutil.copy(telpath+'/'+tel_combinedsrc, telpath+'/'+newtelapfilePrefix+tel_combinedsrc)
    os.remove('v'+sci_combinedsrc)
    
    # These settings in nsextract will force it to use the aperture size and center in the revised telluric apfiles
    iraf.nsextract(
        inimages=sci_combinedsrc, outspectra='', outprefix=extractRegularPrefix, dispaxis=1, database='', line=700,
        nsum=apertureTracingColumns, ylevel='INDEF', upper=str(extractApertureRadius),
        lower=-str(extractApertureRadius), background=subtractBkg, fl_vardq='yes', fl_addvar='no', fl_skylines='yes',
        fl_inter=nsextractInter, fl_apall=useApall, fl_trace='no', aptable='gnirs$data/apertures.fits',
        fl_usetabap='no', fl_flipped='yes', fl_project='yes', fl_findneg='no', bgsample='*',
        trace=telpath+'/'+newtelapfilePrefix+tel_combinedsrc, tr_nsum=10, tr_step=10, tr_nlost=3,
        tr_function='legendre', tr_order=5, tr_sample='*', tr_naver=1 ,tr_niter=0, tr_lowrej=3.0, tr_highrej=3.0,
        tr_grow=0.0, weights='variance', logfile=logger.root.handlers[0].baseFilename, verbose='yes', mode='al')

    # NOTE:  There is a slight complication here - we occasionally find that nsextract locates the aperture too close 
    # to the end of the slit. Due to this, it exits with an "Aperture too large" error and spectra are not extracted 
    # for one or more orders. XDGNIRS works around that in XDpiped.csh, where it ignores this error. However, the 
    # error can occur for other reasons. So, here we check if all file extensions are present for the extracted target
    # file (should really add other files as well...)  Viraja:  I believe the other files are the database files.
    extractedSpectraExtensions = iraf.gemextn(inimages=sci_combinedsrc, check='exists,mef', process='expand',
        index='', extname='SCI', extversion='', ikparams='', omit='', replace='', outfile='STDOUT',
        logfile=logger.root.handlers[0].baseFilename, glogpars='', verbose='yes', fail_count='0', count='20',
        status='0', Stdout=1)

    if len(extractedSpectraExtensions) != len(orders):
        # TODO(Viraja):  Can ask the user to change the extractApertureRadius and redo the extraction. Check with 
        # Andy if this could work. If yes, I think this is something that can also be done when the spectra are first 
        # extracted non-interactively (although I am not sure if that would make sense if we do a peak check).
        logger.error("The combined science image file contains only %d extensions.", len(extractedSpectraExtensions))
        logger.error("Please run gnirsCombineSpectra2D.py to create the combined 2D science image with the right")
        logger.error("number of extensions or provide the combined science image with the right number of extensions")
        logger.error("manually. Exiting script.")
        raise SystemExit

# ----------------------------------------------------------------------------------------------------------------------
'''
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
        nsum=apertureTracingColumns, ylevel='INDEF', upper=str(extractApertureRadius), \
        lower='-'+str(extractApertureRadius), background='none', fl_vardq='yes', fl_addvar='no', fl_skylines='yes',\
        fl_inter=nsextractInter, fl_apall=useApall, fl_trace='no', aptable='gnirs$data/apertures.fits', \
        fl_usetabap='no', fl_flipped='yes', fl_project='yes', fl_findneg='no', bgsample='*', trace='', tr_nsum=10, \
        tr_step=10, tr_nlost=3, tr_function='legendre', tr_order=5, tr_sample='*', tr_naver=1, tr_niter=0, \
        tr_lowrej=3.0, tr_highrej=3.0, tr_grow=0.0, weights='variance', logfile=logger.root.handlers[0].baseFilename, \
        verbose='yes', mode='al')

    # First trace the peak to make sure that the same part of the object is extracted in each step along the slit for  
    # all orders. This is required when the structure is complex, e.g., structure varying in the spectral direction in  
    # an extended object such as a galaxy; otherwise the spectra can have offsets between orders.
    
    # This first nsextract step performed outside the loop gets the trace into the science extraction database to be 
    # used during the actual stepwise extraction.
    extractApertureRadius = 3
    iraf.nsextract(inimages=combinedimage, outspectra='extractionStepwiseTraceReference', outprefix='x', dispaxis=1, \
        database='', line=700, nsum=apertureTracingColumns, ylevel='INDEF', upper=str(extractApertureRadius), \
        lower='-'+str(extractApertureRadius), background='none', fl_vardq='yes', fl_addvar='no', fl_skylines='yes',\
        fl_inter=nsextractInter, fl_apall=useApall, fl_trace='yes', aptable='gnirs$data/apertures.fits', \
        fl_usetabap='no', fl_flipped='yes', fl_project='no', fl_findneg='no', bgsample='*', trace='', tr_nsum=10, \
        tr_step=10, tr_nlost=3, tr_function='legendre', tr_order=5, tr_sample='300:1000', tr_naver=1, tr_niter=0, \
        tr_lowrej=3.0, tr_highrej=3.0, tr_grow=0.0, weights='variance', logfile=logger.root.handlers[0].baseFilename, \
        verbose='yes', mode='al')
    
    # This second step is never done interactively, because it uses extraction details from the previous call to 
    # nsextract
    nsextractInter = False
    for i in range(-extractionStepradius,extractionStepradius,extractionStepSize):
        iraf.nsextract(inimages=combinedimage, outspectra='', outprefix='s'+str(n), dispaxis=1, database='', line=700,\
            nsum=apertureTracingColumns, ylevel='INDEF', upper=i+extractionStepSize, lower=i, background='none', \
            fl_vardq='yes', fl_addvar='no', fl_skylines='yes', fl_inter=nsextractInter, fl_apall=useApall, \
            fl_trace='no', aptable='gnirs$data/apertures.fits', fl_usetabap='no', fl_flipped='yes', fl_project='yes', \
            fl_findneg='no', bgsample='*', trace='extractionStepwiseTraceReference', tr_nsum=10, tr_step=10, \
            tr_nlost=3, tr_function='legendre', tr_order=5, tr_sample='*', tr_naver=1, tr_niter=0, tr_lowrej=3.0, \
            tr_highrej=3.0, tr_grow=0.0, weights='variance', logfile=logger.root.handlers[0].baseFilename, \
            verbose='yes', mode='al')
'''


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    log.configure('gnirs-pype.log', filelevel='INFO', screenlevel='DEBUG')
    start('gnirs-pype.cfg')
