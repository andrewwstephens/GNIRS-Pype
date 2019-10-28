#!/usr/bin/env python

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

    path = os.getcwd()  # Store current working directory for later use.

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
    us_clobber = iraf.envget("clobber")
    iraf.reset(clobber='yes')
    
    config = ConfigParser.RawConfigParser()
    config.optionxform = str  # make options case-sensitive
    config.read(configfile)

    manualMode = config.getboolean('defaults', 'manualMode')
    overwrite = config.getboolean('defaults', 'overwrite')

    # Order of sections is important to later check for plausible peaks located for science targets by nsextract
    nsextractInter = config.getboolean('interactive', 'nsextractInter')
    combinedsrc = config.get('runtimeFilenames', 'combinedsrc')
    combinedsky = config.get('runtimeFilenames', 'combinedsky')
    extractRegularPrefix = config.get('runtimeFilenames', 'extractRegularPrefix')
    extractFullSlitPrefix = config.get('runtimeFilenames', 'extractFullSlitPrefix')
    extractStepwiseTracePrefix = config.get('runtimeFilenames', 'extractStepwiseTracePrefix')
    extractStepwisePrefix = config.get('runtimeFilenames', 'extractStepwisePrefix')
    useApall = config.getboolean('extractSpectra1D', 'useApall')
    subtractBkg = config.get('extractSpectra1D', 'subtractBkg')
    extractApertureRadius = config.getfloat('extractSpectra1D', 'extractApertureRadius')
    checkPeaksMatch = config.getboolean('extractSpectra1D', 'checkPeaksMatch')
    toleranceOffset = config.getfloat('extractSpectra1D', 'toleranceOffset')
    extractFullSlit = config.getboolean('extractSpectra1D', 'extractFullSlit')
    extractStepwise = config.getboolean('extractSpectra1D', 'extractStepwise')
    extractionStepSize = config.getfloat('extractSpectra1D', 'extractStepSize')
    extractApertureWindow = config.getfloat('extractSpectra1D', 'extractApertureWindow')

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

            logger.info(' ----------------------- ')
            logger.info('| Extracting 1D spectra |')
            logger.info(' ----------------------- ')

            obspath += '/Intermediate'
            logger.info("%s\n", obspath)
            iraf.chdir(obspath) 
            utils.pause(manualMode)

            utils.requires([combinedsrc])

            calculateSNR = config.getboolean('gnirsPipeline', 'CalculateSNR')
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

                        logger.info("Comparing the science and Telluric extraction locations...")
                        reextract, predicted = compare_peaks(obspath, telpath, scipeaks, telpeaks, toleranceOffset)

                        if any(reextract):
                            logger.warning("Re-extracting...")
                            useApall = 'yes'
                            nsum = 20
                            reExtractSpectra(reextract, scipeaks, telpeaks, predicted, obspath, telpath, nsum,
                                             extractApertureRadius, useApall, subtractBkg, nsextractInter)

                # ------------------------------------------------------------------------------------------------------
                if extractFullSlit:
                    logger.warning('Full-slit extraction is untested')
                    utils.pause(manualMode)

                    # Approx. full-slit extraction (science target only)
                    # Uses +/- 23 pix aperture (6.9", almost whole length of slit), appropriate for objects centred
                    # along length of slit (q=0).  Not sure what the effect is if nsextract finds a spectrum that's
                    # not centred along the slit.

                    iraf.nsextract(
                        inimages='src_comb', outspectra='', outprefix='a', dispaxis=1, database='', line=700,
                        nsum=20, ylevel='INDEF', upper=23, lower=-23, background='none', fl_vardq='yes', fl_addvar='no',
                        fl_skylines='yes', fl_inter=nsextractInter, fl_apall=useApall, fl_trace='no',
                        aptable='gnirs$data/apertures.fits', fl_usetabap='no', fl_flipped='yes', fl_project='yes',
                        fl_findneg='no', bgsample='*', trace='', tr_nsum=10, tr_step=10, tr_nlost=3,
                        tr_function='legendre', tr_order=5, tr_sample='*', tr_naver=1, tr_niter=0, tr_lowrej=3.0,
                        tr_highrej=3.0, tr_grow=0.0, weights='variance', logfile=logger.root.handlers[0].baseFilename,
                        verbose='yes', mode='al')

                # ------------------------------------------------------------------------------------------------------
                if extractStepwise:  # Extract in steps on either side of the peak
                    logger.warning('Step-wise extraction is untestd')
                    utils.pause(manualMode)

                    # Calling apall and tracing the peak first to make sure the same part of the object is extracted in
                    # each step along the slit for all orders (needed when there is complex, spectrally-varying
                    # structure in a galaxy, for example; otherwise the spectra can have offsets between orders)

                    # This first nsextract step, outside the loop, gets the trace into the database to be used when we
                    # do the "real" extraction

                    iraf.nsextract(
                        inimages='src_comb', outspectra='trace_ref', outprefix='x', dispaxis=1, database='', line=700,
                        nsum=20, ylevel='INDEF', upper=3, lower=-3, background='none', fl_vardq='yes', fl_addvar='no',
                        fl_skylines='yes', fl_inter=nsextractInter, fl_apall='yes', fl_trace='yes',
                        aptable='gnirs$data/apertures.fits', fl_usetabap='no', fl_flipped='yes' ,fl_project='no',
                        fl_findneg='no', bgsample='*', trace='', tr_nsum=10, tr_step=10, tr_nlost=3,
                        tr_function='legendre', tr_order=5, tr_sample='300:1000', tr_naver=1, tr_niter=0, tr_lowrej=3.0,
                        tr_highrej=3.0, tr_grow=0.0, weights='variance', logfile=logger.root.handlers[0].baseFilename,
                        verbose='yes', mode='al')

                    # This is non-interactive because it uses info from the previous call (and it would be very tedious)
                    # TODO: Make sure that the stepping range and step size results in an integer number of steps

                    step = 3
                    n = 0
                    for i in range(-21, 21, step):
                        iraf.nsextract(
                            inimages='src_comb', outspectra='', outprefix='s'+str(n), dispaxis=1, database='', line=700,
                            nsum=20, ylevel='INDEF', lower=i, upper=i+step, background='none', fl_vardq='yes',
                            fl_addvar='no', fl_skylines='yes', fl_inter='no', fl_apall='no', fl_trace='no',
                            aptable='gnirs$data/apertures.fits', fl_usetabap='no', fl_flipped='yes', fl_project='yes',
                            fl_findneg='no', bgsample='*', trace='', tr_nsum=10, tr_step=10, tr_nlost=3,
                            tr_function='legendre', tr_order=5, tr_sample='*', tr_naver=1, tr_niter=0, tr_lowrej=3.0,
                            tr_highrej=3.0, tr_grow=0.0, weights='variance',
                            logfile=logger.root.handlers[0].baseFilename, verbose='yes', mode='al')
                        n += 1

            logger.info("Extraction complete for")
            logger.info("%s", obspath)

    iraf.chdir(path)  # Return to directory script was begun from

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
    # I'm tempted to name this 'nsextract' and call this whenever I need nsextract.
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
                    p = float(line.split()[1])
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
    telA = utils.files_in([telpath + '/nodA.list'])

    # Assuming that both the science and Telluric were acquired in the center of the slit, the extraction locations
    # should be the same minus any difference in the Q-offset.
    # Here I assume that I only need to compare the "A" offset (B might be a sky):

    logger.debug('Science "A" offset: %s arcsec', sciinfo[sciA[0]]['Q'])
    logger.debug('Telluric "A" offset: %s arcsec', telinfo[telA[0]]['Q'])

    # TODO:  The PIXSCALE should be in the obslog, but it's not, so for now:
    pixscale = 0.15

    offset = (float(sciinfo[sciA[0]]['Q']) - float(telinfo[telA[0]]['Q'])) / pixscale
    logger.debug('offset: %s pix', offset)

    logger.debug('Extraction locations (pixels):')

    shifts = []
    predicted = []
    reextract = []
    for s,t in zip(scipeaks, telpeaks):
        predicted.append(t + offset)
        shifts.append(t + offset - s)
        reextract.append(abs(shifts[-1]) > tolerance)
        logger.debug('Telluric: %6.2f  Predicted sci: %6.2f  Actual sci: %6.2f', t, t + offset, s)

    logger.debug('predicted: %s', predicted)
    logger.debug('shifts: %s', shifts)
    logger.debug('reextract: %s', reextract)

    if any(reextract):
        logger.warning('Some orders are not where they were expected to be.')

    # I'm not sure if I should do the tolerance checking here, or pass back the list of shifts,
    # or a list of booleans, or do the comparison in the main program...

    # nsextract should find the spectrum within a 'tolerance' pixels of expected location. This depends on how well the
    # observer centred the target along the slit. Here, we use 5 pixels as a reasonable tolerance level. A more robust
    # way would be to use some measure of whether the peak found by nsextract was real, e.g. counts + FWHM. However, 
    # this information is not recorded in database.

    return reextract, predicted


# ----------------------------------------------------------------------------------------------------------------------
def reExtractSpectra(reextract, scipeaks, telpeaks, predicted, scipath, telpath, nsum, aperture, apall, background,
                     interactive):

    # rextract - boolean list of extensions that should be re-extracted
    # predicted - float list of new extraction locations

    logger = log.getLogger('reExtractSpectra')
    logger.debug('scipath: %s', scipath)
    logger.debug('rextract: %s', reextract)
    logger.debug('predicted: %s', predicted)

    # Rename the old extracted spectum for posterity
    # Rename the old aperture files for posterity
    # Copy in the Telluric aperture files
    # Edit the Telluric aperture files to have the predicted science spectra locations
    # Run nsextract with fl_trace=no and set the tracing reference image (trace) to the edited Telluric file
    # Test by setting the tolerance to be ~1 pix which will force some orders to be reextracted.

    logger.debug('Renaming old extracted spectra...')
    os.rename('vsrc_comb.fits', 'vsrc_comb_TRACED.fits')

    logger.debug('Generating reference files...')
    for i in range(len(reextract)):
        ext = i+1
        oldsciapfile = '%s/database/apsrc_comb_SCI_%d_' % (scipath, ext)
        os.rename(oldsciapfile, oldsciapfile + 'TRACED')
        telapfile = '%s/database/apsrc_comb_SCI_%d_' % (telpath, ext)
        refapfile = '%s/database/apref_comb_SCI_%d_' % (scipath, ext)
        shutil.copy(telapfile, refapfile)

        with open(refapfile, 'r') as f:
            data = f.read()

        # Following XDGNIRS replace the Telluric location with either the shifted Telluric or the Science location:
        with open(refapfile, 'w') as f:
            if reextract[i]:
                logger.debug('Substituting predicted position: %s', predicted[i])
                f.write(data.replace(str(telpeaks[i]), str(predicted[i])).replace('src_comb', 'ref_comb'))
            else:
                logger.debug('Substituting science position: %s', scipeaks[i])
                f.write(data.replace(str(telpeaks[i]), str(scipeaks[i])).replace('src_comb', 'ref_comb'))

    shutil.copy(telpath + '/src_comb.fits', 'ref_comb.fits')

    logger.debug('Running nsextract with the modified reference file and trace=no...')
    iraf.nsextract(
        inimages='src_comb.fits', outspectra='', outprefix='v', dispaxis=1, database='', line=700, nsum=nsum,
        ylevel='INDEF', upper=aperture, lower=-aperture, background=background, fl_vardq='yes', fl_addvar='no',
        fl_skylines='yes', fl_inter=interactive, fl_apall=apall, fl_trace='no', aptable='gnirs$data/apertures.fits',
        fl_usetabap='no', fl_flipped='yes', fl_project='yes', fl_findneg='no', bgsample='*', trace='ref_comb',
        tr_nsum=10, tr_step=10, tr_nlost=3, tr_function='legendre', tr_order=5, tr_sample='*', tr_naver=1, tr_niter=0,
        tr_lowrej=3.0, tr_highrej=3.0, tr_grow=0.0, weights='variance', logfile=logger.root.handlers[0].baseFilename,
        verbose='yes')

    # Sometimes nsextract locates the aperture too close to the end of the slit.
    # When this happens it fails with "Aperture too large" and spectra are not extracted for that order.
    # Check if all file extensions are present in the extracted target file:

    extracted_sci_extensions = iraf.gemextn(
        inimages='src_comb', check='exists,mef', process='expand', index='', extname='SCI', extversion='', ikparams='',
        omit='', replace='', outfile='STDOUT', logfile=logger.root.handlers[0].baseFilename, glogpars='', verbose='yes',
        fail_count='0', count='20', status='0', Stdout=1)

    logger.debug('extracted_sci_extensions: %s', extracted_sci_extensions)

    if len(extracted_sci_extensions) != len(reextract):
        logger.error("The combined science image file contains only %d extensions.", len(extracted_sci_extensions))
        raise SystemExit

    return


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    log.configure('gnirs-pype.log', filelevel='INFO', screenlevel='DEBUG')
    start('gnirs-pype.cfg')
