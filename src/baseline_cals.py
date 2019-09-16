#!/usr/bin/env python

import ConfigParser
import glob
import log
import os
import pkg_resources
from pyraf import iraf
import sys
import utils


# ----------------------------------------------------------------------------------------------------------------------
def start(configfile):
    """
    This module contains all the functions needed to reduce GNIRS BASELINE CALIBRATIONS.

    INPUT FILES FOR BASELINE CALIBRATION:

    INPUT FILES:
        - Configuration file                                                                                                                                                                                                                                              
        - All calibration frames
        - QH flat frames
        - IR flat frames
        - Arc frames
        - Pinhole flat frames

    OUTPUT FILES:
        - masterflat.fits
        - wavelength calibrated arc image

    Args:
        - configfile: gnirs-pype.cfg configuration file
                - Paths to the calibrations (str), reduction truth value (boolean): 
                  E.g. 'target/date/config/{Sci,Tel}_ObsID/{Calibrations,Intermediate}', True
                - manualMode (boolean): Enable optional manualModeging pauses? Default: False
                - overwrite (boolean): Overwrite old files? Default: False
                - start (int): starting step of calibration reductions. Specified at command line with -a. Default: 1
                - stop (int): stopping step of calibration reductions. Specified at command line with -z. Default: 5
                - cleanir_QHflats (boolean): Cleaning QH flat frames. Default: False
                - cleanir_IRflats (boolean): Cleaning IR flat frames. Default: False
                - cleanir_arcs (boolean): Cleaning arc frames. Default: False
                - cleanir_pinholes (boolean): Cleaning pinhole flat frames. Default: False
                - nsprepareInter (boolean): Running nsprepare interactively? Default: False
                - nscombineInter (boolean): Running nscombine interactively? Default: False
                - nssdistInter (boolean): Running nssdist interactively? Default: False
                - nswavelengthInter (boolean): Running nswavelength interactively? Default: False
                - nsfitcoordsInter (boolean): Running nsfitcoords interactively? Default: False
    """
    logger = log.getLogger('baseline_cals')

    path = os.getcwd()  # Store current working directory for later use.

    logger.info(' ------------------------------------------------')
    logger.info('| Start the GNIRS Baseline Calibration Reduction |')
    logger.info(' ------------------------------------------------')

    # Set up/prepare IRAF.
    iraf.gemini()
    iraf.gnirs()
    iraf.gemtools()

    # Reset to default parameters the used IRAF tasks.
    iraf.unlearn(iraf.gemini,iraf.gemtools,iraf.gnirs)

    # Prepare the IRAF package for GNIRS.
    # NSHEADERS lists the header parameters used by the various tasks in the GNIRS package (excluding headers values 
    # which have values fixed by IRAF or FITS conventions).
    iraf.nsheaders("gnirs", logfile=logger.root.handlers[0].baseFilename)

    # Set clobber to 'yes' for the script. This still does not make the gemini tasks overwrite files, so: YOU WILL 
    # LIKELY HAVE TO REMOVE FILES IF YOU RE_RUN THE SCRIPT.
    user_clobber=iraf.envget("clobber")
    iraf.reset(clobber='yes')

    config = ConfigParser.RawConfigParser()
    config.optionxform = str  # make options case-sensitive
    config.read(configfile)

    manualMode = config.getboolean('defaults', 'manualMode')
    overwrite = config.getboolean('defaults', 'overwrite')
    nsprepareInter = config.getboolean('interactive', 'nsprepareInter')
    nsflatInter = config.getboolean('interactive', 'nsflatInter')
    nscombineInter = config.getboolean('interactive', 'nscombineInter')
    nssdistInter = config.getboolean('interactive', 'nssdistInter')
    nswavelengthInter = config.getboolean('interactive', 'nswavelengthInter')
    nsfitcoordsInter = config.getboolean('interactive', 'nsfitcoordsInter')
    preparedPrefix = config.get('runtimeFilenames', 'preparedPrefix')
    reducedPrefix = config.get('runtimeFilenames', 'reducedPrefix')
    QHflat = config.get('runtimeFilenames', 'QHflat')
    QHflat_bpm = config.get('runtimeFilenames', 'QHflat_bpm')
    IRflat = config.get('runtimeFilenames', 'IRflat')
    IRflat_bpm = config.get('runtimeFilenames', 'IRflat_bpm')
    masterflat = config.get('runtimeFilenames', 'masterflat')
    combinedarc = config.get('runtimeFilenames', 'combinedarc')

    waveCalibPrefix = config.get('runtimeFilenames', 'waveCalibPrefix')
    fitcoordsPrefix = config.get('runtimeFilenames', 'fitcoordsPrefix')
    transformPrefix = config.get('runtimeFilenames', 'transformPrefix')
    cleanirPrefix = config.get('runtimeFilenames', 'cleanirPrefix')

    # Loop over the Calibrations directories and reduce the calibrations for the ones with caibration_direcotries True.
    for calpath in config.options('CalibrationDirectories'):
        if config.getboolean('CalibrationDirectories', calpath):

            iraf.chdir(calpath)
            logger.info("Processing calibrations in %s\n", calpath)

            orders = utils.get_orders(calpath)
            nominal_wavelengths, advertised_accuracy = utils.get_wavelengths(calpath)

            required_files = ['all.list', 'QHflats.list', 'IRflats.list', 'pinholes.list', 'arcs.list']
            for f in required_files:
                if os.path.exists(f):
                    logger.info('Found %s', f)
                else:
                    logger.error('Could not find %s', f)
                    logger.error("Please run make_lists.py or provide it manually in %s", calpath)
                    raise SystemExit

            sdistrefimage = utils.files_in(['pinholes.list'])[0]  # use the first pinhole image
            logger.info('S-Distortion reference image: %s', sdistrefimage)

            startstep = config.getint('calibrationReduction', 'Start')
            stopstep = config.getint('calibrationReduction', 'Stop')
            if startstep > stopstep or startstep < 1 or stopstep > 5:
                logger.error('Invalid start/stop values')
                raise SystemExit

            for valindex in range(startstep, stopstep + 1):
                logger.debug('valindex = %d', valindex)

                if valindex == 1:
                    logger.info(" ----------------------------------")
                    logger.info("| STEP 1: Clean calibration frames |")
                    logger.info(" ----------------------------------")

                    utils.pause(manualMode)

                    if config.getboolean('calibrationReduction','cleanir_QHflats'):
                        utils.clean('QHflats.list', cleanirPrefix, overwrite)
                    else:
                        logger.warning('QH Flats not cleaned')

                    if config.getboolean('calibrationReduction','cleanir_IRflats'):
                        utils.clean('IRflats.list', cleanirPrefix, overwrite)
                    else:
                        logger.warning('IR Flats not cleaned')

                    if config.getboolean('calibrationReduction','cleanir_arcs'):
                        utils.clean('arcs.list', cleanirPrefix, overwrite)
                        pass
                    else:
                        logger.warning('Arcs not cleaned')

                    if config.getboolean('calibrationReduction','cleanir_pinholes'):
                        utils.clean('pinholes.list', cleanirPrefix, overwrite)
                        pass
                    else:
                        logger.warning('Pinholes not cleaned')

                    # NOTE:  At this point in XDGNIRS, StatsAndPrepartXD.py checks for any deviations in statistical
                    # parameters between the raw and the cleaned images. Since cleaning is not incorporated yet in 
                    # this pipeline, this check is skipped here.

                    # NOTE:  StatsAndPrepartXD.py also records any deviations in the mean values of cleaned QHflats
                    # and IRflats before preparing calibrations for reduction (here, step 2 prepareCalibrations). If 
                    # cleaning is skipped, statistics is done on the raw flat frames. This check is skipped in this 
                    # pipeline for now, but it could be useful to do.

                elif valindex == 2:

                    logger.info(" ------------------------------------")
                    logger.info("| STEP 2: Prepare calibration frames |")
                    logger.info(" ------------------------------------")
                    utils.pause(manualMode)

                    prepareCalibrations('QHflats.list', 'IRflats.list', 'arcs.list', 'pinholes.list',
                                        nsprepareInter, preparedPrefix, reducedPrefix, overwrite)

                elif valindex == 3:

                    logger.info(" ---------------------------")
                    logger.info("| STEP 3: Create flat field |")
                    logger.info(" ---------------------------")
                    utils.pause(manualMode)

                    makeFlat(preparedPrefix, reducedPrefix, 'QHflats.list', 'IRflats.list', QHflat, QHflat_bpm,
                             IRflat, IRflat_bpm, masterflat, nsflatInter, overwrite)

                elif valindex == 4:

                    logger.info(" ---------------------------------------------------------")
                    logger.info("| STEP 4: Trace spatial curvature and spectral distortion |")
                    logger.info(" ---------------------------------------------------------")
                    utils.pause(manualMode)
                    # Output: Pinhole flat for which the spatial curvature and S-distortion has been determined.

                    if 'LB_SXD' in calpath or 'LB_LXD' in calpath:
                        pinhole_coordlist = 'gnirs$data/pinholes-long-dense-north.lis'
                        pinholenumber = 9
                    elif 'SB_SXD' in calpath:
                        pinhole_coordlist = 'gnirs$data/pinholes-short-dense-north.lis'
                        pinholenumber = 6
                    else:
                        logger.error('Unknown camera.')
                        raise SystemExit

                    makeSdistortion(nssdistInter, sdistrefimage, pinhole_coordlist, pinholenumber,
                                    reducedPrefix, preparedPrefix, overwrite)

                elif valindex == 5:

                    logger.info(" ----------------------------------------------------------------------- ")
                    logger.info("| STEP 5: Combine arc frames, determine the wavelength solution, and    |")
                    logger.info("| create the wavelength referenced arc corrected for spatial distortion |")
                    logger.info(" ----------------------------------------------------------------------- ")
                    # OUTPUT: Combined and wavelength calibrated arc image.

                    utils.pause(manualMode, 'Starting combine_arcs')
                    combine_arcs(preparedPrefix, reducedPrefix, combinedarc, nscombineInter, overwrite)

                    if '10' in calpath and '32' in calpath:
                        coordlist = 'gnirs$data/lowresargon.dat'
                    elif '110' in calpath:
                        coordlist = 'gnirs$data/argon.dat'
                    else:
                        logger.error("Uknown grating.")
                        raise SystemExit

                    utils.pause(manualMode, 'Starting makeWaveCal')
                    makeWaveCal(preparedPrefix, reducedPrefix, combinedarc, fitcoordsPrefix, transformPrefix,
                                waveCalibPrefix, nsfitcoordsInter, nswavelengthInter, sdistrefimage,
                                coordlist, orders, overwrite)

                    utils.pause(manualMode, 'Starting check_wavelengths.')
                    check_wavelengths(combinedarc, fitcoordsPrefix, transformPrefix, orders, nominal_wavelengths,
                                      advertised_accuracy, overwrite)

            logger.info(" -------------------------------------- ")
            logger.info("| Baseline Calibrations step complete. |")
            logger.info(" -------------------------------------- ")

    iraf.chdir(path)  # Return to where we began
        
    return


# ----------------------------------------------------------------------------------------------------------------------
def prepareCalibrations(QHflats, IRflats, arcs, pinholes, interactive, preparedPrefix,
                        reducedPrefix, overwrite):
    """
    Prepare calibration frames for further processing.

    Use NSPREPARE on all calibration frames (QHflats, IRflats, arcs, pinholes, and darks to update the raw data headers
    and attach the mask definition file (MDF) as a binary table on all files. Note that dark frames will not have an
    MDF attached by default. Instead, the appropriate MDF is added in NSREDUCE or NSFLAT to match the data being
    reduced.

    Use NSREDUCE to cut the calibration frame spectra to the size specified by the MDF, placing different spectral
    orders in separate image extensions. Darks will not be used for further reductions after this point.

    :param QHflats: filename of QH Flats list
    :param IRflats:
    :param arcs:
    :param pinholes:
    :param interactive:
    :param preparedPrefix:
    :param reducedPrefix:
    :param overwrite:
    :return:
    """
    logger = log.getLogger('prepareCalibrations')

    infiles = utils.files_in([QHflats, IRflats, arcs, pinholes])
    utils.requires(infiles)
    outfiles = [preparedPrefix + f for f in infiles] + [reducedPrefix + preparedPrefix + f for f in infiles]
    if utils.exists(outfiles, overwrite):
        logger.info('All calibrations already prepared.')
        return

    bpm = utils.get_bpm(infiles[0])

    # Update all calibration frames with mdf offset value and generate variance and data quality extensions.
    iraf.nsprepare(
        inimages='@' + QHflats + ',@' + IRflats + ',@' + arcs + ',@' + pinholes, rawpath='',
        outimages='', outprefix=preparedPrefix, bpm=bpm, logfile=logger.root.handlers[0].baseFilename,
        fl_vardq='yes', fl_cravg='no', crradius=0.0, fl_dark_mdf='no', fl_correct='no', fl_saturated='yes', 
        fl_nonlinear='yes', fl_checkwcs='yes', fl_forcewcs='yes', arraytable='gnirs$data/array.fits',
        configtable='gnirs$data/config.fits', specsec='[*,*]', offsetsec='none', pixscale='0.15', shiftimage='',
        shiftx='INDEF', shifty='INDEF', obstype='FLAT', fl_inter=interactive, verbose='yes', mode='al')

    # Cut all calibration frames according to the size specified by the MDFs.
    iraf.nsreduce(
        inimages=preparedPrefix + '//@' + QHflats + ',' + preparedPrefix + '//@' + IRflats + ',' +
                 preparedPrefix + '//@' + arcs + ',' + preparedPrefix + '//@' + pinholes,
        outimages='', outprefix=reducedPrefix, fl_cut='yes', section='', fl_corner='yes', fl_process_cut='yes',
        fl_nsappwave='no', nsappwavedb='gnirs$data/nsappwave.fits', crval='INDEF', cdelt='INDEF', fl_dark='no',
        darkimage='', fl_save_dark='no', fl_sky='no', skyimages='', skysection='', combtype='median',
        rejtype='avsigclip', masktype='goodvalue', maskvalue=0.0, scale='none', zero='median', weight='none',
        statsec='[*,*]', lthreshold='INDEF', hthreshold='INDEF', nlow=1, nhigh=1, nkeep=0, mclip='yes', lsigma=3.0,
        hsigma=3.0, snoise='0.0', sigscale=0.1, pclip=-0.5, grow=0.0, skyrange='INDEF', nodsize=3.0, fl_flat='no',
        flatimage='', flatmin=0.0, fl_vardq='yes', logfile=logger.root.handlers[0].baseFilename, verbose='yes',
        debug='no', force='no', mode='al')

    return


# ----------------------------------------------------------------------------------------------------------------------
def makeFlat(preparedPrefix, reducedPrefix, QHflatsfilename, IRflatsfilename, QHflat, QHflat_bpm, IRflat, IRflat_bpm,
    masterflat, nsflatInter, overwrite):
    """
    Make flat field and bad pixel masks.

    Use NSFLAT to generate a normalized flat field for the QHflats and IRflats (for each cross-dispersed order).
    A bad-pixel mask (BPM) will also be generated by thresholding - this can be used to flag bad pixels in other data.

    Use FXCOPY followed by FXINSERT to produce the final flat by grouping order 3 of IRflats and orders 4-8 of
    QHflats. The output from this task is used as the flat field for further reductions.

    TODO(Viraja):  Sometimes nsflat crashes with a fixpix or other, unknown error. In such a situation, tweaking the 
    lthresh parameter sometimes helps. XDGNIRS loops through a fixed list of lthresh values until it (hopefully) runs
    without an error in three trials after which it fails and exits. The nsflat function in this pipeline is currently
    not set to run through a loop; it simply uses the specified lthresh. However, it could be made to run through a loop
    of random lthresh values in a specified range until it runs without an error to avoid a system exit.
    """
    logger = log.getLogger('makeFlat')

    infiles = ['n%s' % f for f in utils.files_in([QHflatsfilename, IRflatsfilename])]
    utils.requires(infiles)
    outfiles = [QHflat, QHflat_bpm, IRflat, IRflat_bpm, masterflat]
    if utils.exists(outfiles, overwrite):
        logger.info('Flat fields and BPMs exist.  Skipping this step.')
        return

    logger.info('Processing QH flats...')
    iraf.nsflat(
        lampson=reducedPrefix + preparedPrefix + '//@' + QHflatsfilename, darks='', flatfile=QHflat, darkfile='',
        fl_corner='yes', fl_save_darks='no', flattitle='default', bpmtitle='default', bpmfile=QHflat_bpm, process="fit",
        statsec='MDF', fitsec='MDF', thr_flo=0.35, thr_fup=4.0, thr_dlo=-20, thr_dup=100, fl_inter=nsflatInter,
        fl_range='no', fl_fixbad='yes', fixvalue=1.0, function='spline3', order=5, normstat='midpt',
        combtype='default', rejtype='ccdclip', masktype='goodvalue', maskvalue=0.0, scale='median', zero='none',
        weight='none', lthreshold=50.0, hthreshold='INDEF', nlow=1, nhigh=1, nkeep=0, mclip='yes', lsigma=3.0,
        hsigma=3.0, snoise='0.0', sigscale=0.1, pclip=-0.5, grow=0.0, box_width=20, box_length=1, trace='',
        traceproc='none', threshold=100.0, aptable='gnirs$data/apertures.fits', database='', apsum=10, tr_step=10,
        tr_nlost=3, tr_function='legendre', tr_order=5, tr_naver=1, tr_niter=0, tr_lowrej=3.0, tr_highrej=3.0,
        tr_grow=0.0, ap_lower=-30, ap_upper=30, fl_vardq='yes', logfile=logger.root.handlers[0].baseFilename,
        verbose='yes', mode='al')

    logger.info('Processing IR flats...')
    iraf.nsflat(
        lampson=reducedPrefix + preparedPrefix + '//@' + IRflatsfilename, darks='', flatfile=IRflat, darkfile='',
        fl_corner='yes', fl_save_darks='no', flattitle='default', bpmtitle='default', bpmfile=IRflat_bpm, process="fit",
        statsec='MDF', fitsec='MDF', thr_flo=0.35, thr_fup=1.5, thr_dlo=-20, thr_dup=100, fl_inter=nsflatInter,
        fl_range='no', fl_fixbad='yes', fixvalue=1.0, function='spline3', order=10, normstat='midpt',
        combtype='default', rejtype='ccdclip', masktype='goodvalue', maskvalue=0.0, scale='none', zero='none',
        weight='none', lthreshold=50.0, hthreshold='INDEF', nlow=1, nhigh=1, nkeep=0, mclip='yes', lsigma=3.0,
        hsigma=3.0, snoise='0.0', sigscale=0.1, pclip=-0.5, grow=0.0, box_width=20, box_length=1, trace='',
        traceproc='none', threshold=100.0, aptable='gnirs$data/apertures.fits', database='', apsum=10, tr_step=10,
        tr_nlost=3, tr_function='legendre', tr_order=5, tr_naver=1, tr_niter=0, tr_lowrej=3.0, tr_highrej=3.0,
        tr_grow=0.0, ap_lower=-30, ap_upper=30, fl_vardq='yes', logfile=logger.root.handlers[0].baseFilename,
        verbose='yes', mode='al')

    # Group order 3 of IRflats and orders 4-8 of QHflats to create the final flat field image.
    # The output flat field image will be masterflat.fits.
    logger.info('Combining flats...')
    iraf.fxcopy(input=IRflat, output=masterflat, group="0-3", new_file='yes', verbose='no', mode='ql')
    iraf.fxinsert(input=QHflat, output=masterflat+'[3]', groups="4-18", verbose='no', mode='ql')

    return

# ----------------------------------------------------------------------------------------------------------------------
def makeSdistortion(interactive, sdistrefimage, pinhole_coordlist, pinholenumber, reducedPrefix,
    preparedPrefix, overwrite):
    """
    Establish Spatial distortion calibration.

    NSSDIST uses the information in the pinhole calibration images to calibrate the spatial distorion of the GNIRS 
    field. The pinhole frame is a dispersed flat field image with a slit-mask in the field so that the illumination on 
    GNIRS is in a pattern of 6 different holes that are stacked in the y-dimension on the field. Proper alignment of
    the slit across the pattern can be used for spatial rectification of the on-sky science data. The spatial solution
    determined by NSSSDIST is linked to the science data in NSFITCOORDS.
    """
    logger = log.getLogger('makeSdistortion')

    infiles = ['rn' + sdistrefimage]
    utils.requires(infiles)

    # FIXME: The number of files should not be hard-coded:
    outfiles = ['database/id' + reducedPrefix + preparedPrefix + utils.nofits(sdistrefimage) + '_SCI_%d_' % i for i in range(1,7)]
    if utils.exists(outfiles, overwrite):
        logger.info('S-Distortion files exist.  Skipping this step.')
        return

    logger.info('Running nssdist...')
    iraf.nssdist(
        inimages='rn' + sdistrefimage, outsuffix='_sdist', pixscale=1.0, dispaxis=1, database='',
        firstcoord=0.0, coordlist=pinhole_coordlist, aptable='gnirs$data/apertures.fits', fl_inter=interactive,
        fl_dbwrite='yes', section='default', nsum=30, ftype='emission', fwidth=10.0, cradius=10.0,
        threshold=1000.0, minsep=5.0, match=-6.0, function="legendre", order=5, sample='', niterate=3,
        low_reject=5.0, high_reject=5.0, grow=0.0, refit='yes', step=10, trace='no', nlost=0, aiddebug='',
        logfile=logger.root.handlers[0].baseFilename, verbose='yes', debug='no', force='no', mode='al')

    logger.info("Checking the number of pinholes identified...")
    for idrnfile in outfiles:
        with open(idrnfile, 'r') as f:
            for line in f:  # Read up to the first occurrence of "features" in the file
                if "features" in line:
                    npinholes = int(line.split()[1])
                    logger.debug('Number of pinholes: %d', npinholes)
                    if npinholes != pinholenumber:
                        logger.warning("Expected %d pinholes to be detected by nssdist, but found ", pinholenumber)
                        logger.warning("%s pinholes in extension %d. This can cause problems. ", line.split()[1], i)
                        logger.warning("Please check the transformed data files (ttf*.fits) later and look for ")
                        logger.warning("inter-order offsets in 'orders.pdf' file created by gnirsCombineOrdersXD.py.")
                        break
                    else:
                        logger.info('Correct number of pinholes detected in %s', idrnfile)
                        break
    logger.info("Number of pinholes check complete.\n")

    return


# ----------------------------------------------------------------------------------------------------------------------
def combine_arcs(preparedPrefix, reducedPrefix, combinedarc, interactive, overwrite):
    # Combine arc frames.
    logger = log.getLogger('combine_arcs')

    infiles = ['rn%s' % f for f in utils.files_in(['arcs.list'])]
    utils.requires(infiles)

    if utils.exists([combinedarc], overwrite):
        logger.info('Combined arc exists.  Skipping combining arcs.')
        return

    iraf.nscombine(
        inimages=reducedPrefix + preparedPrefix + '//@arcs.list', tolerance=0.5, output=combinedarc,
        output_suffix='_comb', bpm='', dispaxis=1, pixscale=1.0, fl_cross='no', fl_keepshift='no',
        fl_shiftint='yes', interptype='linear', boundary='nearest', constant=0.0, combtype='average',
        rejtype='sigclip', masktype='goodvalue', maskvalue=0.0, statsec='[*,*]', scale='none', zero='none',
        weight='none', lthreshold='INDEF', hthreshold='INDEF', nlow=1, nhigh=1, nkeep=0, mclip='yes',
        lsigma=5.0, hsigma=5.0, ron=0.0, gain=1.0, snoise='0.0', sigscale=0.1, pclip=-0.5, grow=0.0,
        nrejfile='', fl_vardq='yes', fl_inter=interactive, logfile=logger.root.handlers[0].baseFilename,
        verbose='yes', debug='no', force='no', mode='al')

    return


# ----------------------------------------------------------------------------------------------------------------------
def makeWaveCal(preparedPrefix, reducedPrefix, combinedarc, fitcoordsPrefix, transformPrefix,
    waveCalibPrefix, nsfitcoordsInter, nswavelengthInter, sdistrefimage, coordlist, orders, overwrite):
    """
    Determine the wavelength solution of the combined arc.

    Uses NSWAVELENGTH to calibrate arcs (after cutting and optionally applying a flatfield with NSREDUCE previously).

    Note that better RMS fits can be obtained by running the wavelength calibration interactively and identifying all 
    of the lines manually. Tedious, but will give more accurate results than the automatic mode (i.e., fl_inter-).
    """
    logger = log.getLogger('makeWaveCal')

    infiles = [combinedarc]
    utils.requires(infiles)

    # FIXME: The number of files should not be hard-coded:
    outfiles = ['farc_comb.fits', 'tfarc_comb.fits', 'wtfarc_comb.fits', 'ftfarc_comb.fits', 'tftfarc_comb.fits'] + \
               ['database/fcfarc_comb_SCI_%d_sdist' % i for i in range(1, 7)] + \
               ['database/fcftfarc_comb_SCI_%d_lamp' % i for i in range(1, 7)] + \
               ['database/idwtfarc_comb_SCI_%d_' % i for i in range(1, 7)]
    if utils.exists(outfiles, overwrite):
        logger.info('All the wavelength calibration files exist.  Skipping this step.')
        return

    # WHY run nsfitcoords twice?

    logger.info("Running nsfitcoords and nstransform on the combined arc")
    logger.info("to spatially rectify it before running nswavelength.")
    iraf.nsfitcoords(
        inimages=combinedarc, outspectra='', outprefix=fitcoordsPrefix, lamptransf='',
        sdisttransf='rn' + sdistrefimage, dispaxis=1, database='', fl_inter=nsfitcoordsInter, fl_align='no',
        function='chebyshev', lxorder=2, lyorder=4, sxorder=4, syorder=4, pixscale=1.0, 
        logfile=logger.root.handlers[0].baseFilename, verbose='yes', debug='no', force='no', mode='al')

    iraf.nstransform(
        inimages='f'+combinedarc, outspectra='', outprefix='t', dispaxis=1, database='',
        fl_stripe='no', interptype='poly3', xlog='no', ylog='no', pixscale=1.0, 
        logfile=logger.root.handlers[0].baseFilename, verbose='yes', debug='no', mode='al')
        
    logger.info("The arc lines might still be a bit tilted because of the detector rotation.")
    logger.info("Running nswavelength on the transformed arc with fl_median-.")
    iraf.nswavelength(
        lampspectra=transformPrefix+fitcoordsPrefix+combinedarc, outspectra='',
        outprefix=waveCalibPrefix, crval='INDEF', cdelt='INDEF', crpix='INDEF', dispaxis=2, database='',
        coordlist=coordlist, fl_inter=nswavelengthInter, nsappwavedb='gnirs$data/nsappwave.fits',
        fl_median='no', sdist='', sdorder=4, xorder=2, yorder=2, aptable='gnirs$data/apertures.fits', 
        section='auto', nsum=10, ftype='emission', fwidth=5.0, cradius=5.0, threshold=300.0, minsep=2.0, 
        match=-6.0, function='chebyshev', order=4, sample='*', niterate=10, low_reject=3.0, high_reject=3.0, 
        grow=0.0, refit='yes', step=2, trace='no', nlost=10, fl_overwrite='yes', aiddebug='', fmatch=0.2, 
        nfound=6, sigma=0.05, rms=0.1, logfile=logger.root.handlers[0].baseFilename, verbose='yes', 
        debug='no', mode='al')
        
    logger.info("Running nsfitcoords to produce the combined arc with horizontal lines")
    logger.info("and then nstransform on the transformed arc using the output from nswavelength")
    logger.info("for the 'lamptrans' parameter.")
    iraf.nsfitcoords(
        inimages=transformPrefix+fitcoordsPrefix+combinedarc, outspectra='', outprefix=fitcoordsPrefix,
        lamptrans=waveCalibPrefix+transformPrefix+fitcoordsPrefix+combinedarc, sdisttransf='', dispaxis=1, 
        database='', fl_inter=nsfitcoordsInter, fl_align='no', function='chebyshev', lxorder=2, lyorder=4,
        sxorder=4, syorder=4, pixscale=1.0, logfile=logger.root.handlers[0].baseFilename, verbose='yes', 
        debug='no', force='no', mode='al')
    iraf.nstransform(
        inimages=fitcoordsPrefix+transformPrefix+fitcoordsPrefix+combinedarc, outspectra='',
        outprefix=transformPrefix, dispaxis=1, database='', fl_stripe='no', interptype='poly3', xlog='no', ylog='no',
        pixscale=1.0, logfile=logger.root.handlers[0].baseFilename, verbose='no', debug='no', mode='al')

    return


# ----------------------------------------------------------------------------------------------------------------------
def check_wavelengths(combinedarc, fitcoordsPrefix, transformPrefix, orders, nominal_wavelengths, advertised_accuracy, overwrite):

    # Check if wavelength calibration looks reasonable by extracting a single column approximately down the middle
    # of each order, then combining orders

    logger = log.getLogger('check_wavelengths')

    logger.info("Creating a wavelength calibrated arc spectrum.")
    # TODO(Viraja):  Check with Andy if it is a good idea to randomize the column selection, within a fixed range

    prefix = transformPrefix + fitcoordsPrefix + transformPrefix + fitcoordsPrefix

    infiles = [prefix + combinedarc]
    utils.requires(infiles)

    outfiles = [prefix + utils.nofits(combinedarc) + '_order%d.fits' % i for i in range(1, 7)] + \
               ['arc_wavelength_test.fits']
    utils.exists(outfiles, overwrite=True)  # Force over-writing here?  I don't want to skip the wavelength check later.

    columns = [88, 77, 65, 54, 53, 92]
    for i in range(len(orders)):
        extension = i + 1
        try:
            iraf.imcopy(
                input=prefix + utils.nofits(combinedarc) + '[SCI,' + str(extension) + '][' + str(columns[i]) + ',*]',
                output=prefix + utils.nofits(combinedarc) + '_order' + str(extension),
                verbose='yes', mode='ql')

        except(TypeError, iraf.IrafError) as e:
            logger.error(e)
            logger.error("A problem was encountered in making the wavelength solution. First, check that the raw ")
            logger.error("and the reduced arcs look sensible. Second, check that the spectral orders are in the ")
            logger.error("same place in the array (in the x direction) in all the raw files to within a few ")
            logger.error("pixels. Third, try running nssdist, nsfitcoords, and nswavelength interactively by ")
            logger.error("setting nssdistInter, nsfitcoordsInter, and nswavelengthInter, respectively. Exiting ")
            logger.error("script.\n")
            sys.exit(0)

    iraf.odcombine(
        input=','.join([prefix + utils.nofits(combinedarc) + '_order%d' % i for i in range(1, 7)]),
        output='arc_wavelength_test.fits',
        headers='', bpmasks='', rejmask='', nrejmasks='', expmasks='', sigma='',
        logfile=logger.root.handlers[0].baseFilename, apertures='', group='all', first='no', w1='INDEF', w2='INDEF',
        dw='INDEF', nw='INDEF', log='no', combine='average', reject='none', outtype='real', outlimits='',
        smaskformat='bpmspectrum', smasktype='none', smaskvalue=0.0, blank=0.0, scale='none', zero='none',
        weight='none', statsec='', expname='', lthreshold='INDEF', hthreshold='INDEF', nlow=1, nhigh=1, nkeep=1,
        mclip='yes', lsigma=3.0, hsigma=3.0, rdnoise='0.0', gain='1.0', snoise='0.0', sigscale=0.1, pclip=-0.5,
        grow=0.0, offsets='physical', masktype='none', maskvalue=0.0, mode='al')

    # ------------------------------------------------------------------------------------------------------------------
    # Check if the starting and ending wavelengths of the orders are reasonable. As of January 2016, we advertised, 
    # that "wavelength settings ... are accurate to better than 5 percent of the wavelength coverage." (see, note a 
    # at the bottom of the table here: http://www.gemini.edu/sciops/instruments/gnirs/spectroscopy). 
    # Take the nominal starting and ending wavelengths from the OT, and check if the ones from the wavelength 
    # calibration are the same to within the advertised tolerance.

    logger.info("Checking that the wavelength solutions are reasonable.")
    logger.debug('Advertised accuracy: %s', advertised_accuracy)
    for i in range(len(orders)):
        extension = i+1

        waveStart = float(iraf.hselect(images=prefix + combinedarc[:combinedarc.rfind('.')]+'[SCI,'+str(extension)+']',
            fields='CRVAL2', expr='yes', missing='INDEF', Stdout=1)[0])
        logger.debug('wavestart: %s', waveStart)

        waveDelt = float(iraf.hselect(images=prefix + combinedarc[:combinedarc.rfind('.')]+'[SCI,'+str(extension)+']',
            fields='CDELT2', expr='yes', missing='INDEF', Stdout=1)[0])
        logger.debug('waveDelt: %s', waveDelt)

        waveEnd = waveStart + (1022 * waveDelt)
        logger.debug('waveEnd: %s', waveEnd)

        logger.debug('Nominal wavelengths for order %d:  %.3f - %.3f',
                     extension, nominal_wavelengths[extension][0], nominal_wavelengths[extension][1])

        allowed_error = (nominal_wavelengths[extension][1] - nominal_wavelengths[extension][0]) * advertised_accuracy/100.0
        logger.debug('Allowed error: %s nm', allowed_error)

        if abs(waveStart - nominal_wavelengths[extension][0]) < allowed_error and \
                abs(waveEnd - nominal_wavelengths[extension][1]) < allowed_error:
            logger.info('%d:  Expected: %.2f - %.2f   Actual: %.2f - %.2f nm', extension, waveStart, waveEnd,
                        nominal_wavelengths[extension][0], nominal_wavelengths[extension][1])
        else:
            logger.warning('%d:  Expected: %.2f - %.2f   Actual: %.2f - %.2f nm', extension, waveStart, waveEnd,
                           nominal_wavelengths[extension][0], nominal_wavelengths[extension][1])

    logger.info("Wavelength solutions check complete.")
    return


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    log.configure('gnirs-pype.log', filelevel='INFO', screenlevel='DEBUG')
    start('gnirs-pype.cfg')
