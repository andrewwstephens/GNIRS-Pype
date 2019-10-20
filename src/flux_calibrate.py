#!/usr/bin/env python

from astropy.io import fits
import ConfigParser
import glob
import header
import log
from pyraf import iraf
import obslog
import os
import utils


# ----------------------------------------------------------------------------------------------------------------------
def start(configfile):
    """
    Do a flux calibration.
    """
    logger = log.getLogger('flux_calibrate')

    path = os.getcwd()  # Store current working directory for later use.

    # Set up/prepare IRAF.
    iraf.gemini()
    iraf.gemtools()
    iraf.gnirs()
    iraf.onedspec()
    iraf.imutil()

    # Reset to default parameters the used IRAF tasks.
    iraf.unlearn(iraf.gemini, iraf.gemtools, iraf.gnirs, iraf.imutil)

    logger.info("Parameters read from %s", configfile)

    config = ConfigParser.RawConfigParser()
    config.optionxform = str  # make options case-sensitive
    config.read(configfile)

    manualMode = config.getboolean('defaults', 'manualMode')
    overwrite = config.getboolean('defaults', 'overwrite')
    extractFullSlit = config.getboolean('extractSpectra1D', 'extractFullSlit')
    extractStepwise = config.getboolean('extractSpectra1D', 'extractStepwise')
    extractStepSize = config.getfloat('extractSpectra1D', 'extractStepSize')
    combinedsrc = config.get('runtimeFilenames', 'combinedsrc')
    combinedsky = config.get('runtimeFilenames', 'combinedsky')  # Viraja:  Not sure if this is required; surprising though
    extractRegularPrefix = config.get('runtimeFilenames', 'extractRegularPrefix')
    extractFullSlitPrefix = config.get('runtimeFilenames', 'extractFullSlitPrefix')
    extractStepwisePrefix = config.get('runtimeFilenames', 'extractStepwisePrefix')
    hLinePrefix = config.get('runtimeFilenames', 'hLinePrefix')
    dividedTelContinuumPrefix = config.get('runtimeFilenames', 'dividedTelContinuumPrefix')
    telluricPrefix = config.get('runtimeFilenames', 'telluricPrefix')
    bb_unscaled = config.get('runtimeFilenames', 'bb_unscaled')
    bb_scaled = config.get('runtimeFilenames', 'bb_scaled')
    fluxCalibPrefix = config.get('runtimeFilenames', 'fluxCalibPrefix')
    fluxCalibrationMethod = config.get('fluxCalibration', 'fluxCalibrationMethod')
    zeroMagFlux = utils.dictify(config.items('zeroMagnitudeFluxes'))

    for scipath in config.options("ScienceDirectories"):

        if not config.getboolean("ScienceDirectories", scipath):  # only process directories that are marked True
            logger.info('Skipping flux calibaration in %s', scipath)
            continue

        logger.info(' ------------------ ')
        logger.info('| Flux Calibration |')
        logger.info(' ------------------ ')

        scipath += '/Intermediate'
        iraf.chdir(scipath)
        logger.info('%s', scipath)

        if fluxCalibrationMethod == 'fluxcalibrator':
            stdpath = '../Standard/Intermediate'
            if not os.path.exists(stdpath):
                logger.error('fluxCalibrationMethod is fluxcalibrator but could not find Standard star directory')
                raise SystemExit
        elif fluxCalibrationMethod == 'telluricapproximate':
            stdpath = '../Telluric/Intermediate'
            if not os.path.exists(stdpath):
                logger.error('fluxCalibrationMethod is telluricapproximate but could not find Telluric directory')
                raise SystemExit
        else:
            logger.error('Uknown fluxCalibrationMethod')
            raise SystemExit
        logger.debug('stdpath: %s', stdpath)

        orders = utils.get_orders(scipath)
        sciroot = scipath + '/' + dividedTelContinuumPrefix + telluricPrefix + extractRegularPrefix + utils.nofits(combinedsrc)
        stdroot = stdpath + '/' + dividedTelContinuumPrefix + hLinePrefix + extractRegularPrefix + utils.nofits(combinedsrc)
        scifiles = ['%s_order%d.fits' % (sciroot, o) for o in orders]
        stdfiles = ['%s_order%d.fits' % (stdroot, o) for o in orders]
        logger.debug('scifiles: %s', scifiles)
        logger.debug('stdfiles: %s', stdfiles)
        utils.requires(scifiles + stdfiles)

        std_obslog = obslog.readcsv(stdpath + '/obslog.csv')
        firstfile = std_obslog.keys()[0]
        logger.debug('firstfile: %s', firstfile)
        std_name = std_obslog[firstfile]['OBJECT']
        logger.debug('Standard: %s', std_name)

        std_pars = utils.dictify(config.items(std_name))

        temperature = std_pars['Temperature']
        logger.debug('Std Temp: %s', temperature)

        sci_exptime = fits.getheader(scifiles[0])['EXPTIME']
        std_exptime = fits.getheader(stdfiles[0])['EXPTIME']
        logger.debug('Sci Exptime: %s s', sci_exptime)
        logger.debug('Std Exptime: %s s', std_exptime)

        # TODO:  if extractionStepwise:

        if fluxCalibrationMethod == 'fluxcalibrator':
            logger.error('Method not implemented')
            raise SystemExit

        elif fluxCalibrationMethod == 'telluricapproximate':

            # Convert magnitude to flux density for the telluric.
            # Derived spectrum (FLambda) for the telluric for each order.
            # Make a blackbody at the Telluric temperature for each order.
            # Get the scale for the blackbody to the telluric derived spectrum.
            # Blackbody scale (float) at each order.
            # Scale the blackbody to the telluric derived spectrum.
            # Scaled blackbody at each order.

            infile = [os.path.basename(f) for f in scifiles]
            bbody = ['bb_unscaled_order%d.fits' % o for o in orders]
            bbody_scaled = ['bb_scaled_order%d.fits' % o for o in orders]
            outfile = [fluxCalibPrefix + os.path.basename(f) for f in scifiles]

            outfiles = bbody + bbody_scaled + outfile
            if utils.exists(outfiles, overwrite):
                logger.info('All files already flux calibrated.')
                continue  # to the next science directory

            for i in range(len(orders)):
                logger.debug('Order: %s', orders[i])
                logger.debug("Input: %s", infile[i])
                logger.debug("Output: %s", outfile[i])
                logger.debug('Blackbody: %s', bbody[i])
                band = utils.band(orders[i])
                logger.debug('Band: %s', band)
                stdmag = std_pars[band]
                logger.debug('Standard star mag: %s', stdmag)

                logger.info('Flux calibrating %s', infile)

                # FIXME: We should never do abs flux cal in some orders and relative flux cal in other orders

                if stdmag is not None:  # if the magnitude value in the configuration file is not empty
                    # Convert magnitude to flux density (erg/cm2/s/A) for a rough flux scaling
                    logger.info("Performing absolute flux calibration")
                    flambda = 10**(-stdmag/2.5) * zeroMagFlux[utils.band(orders[i])] * std_exptime / sci_exptime
                    absolute_fluxcalib = True

                else:   # if the magnitude value in the config file is empty, no absolute flux calibration is performed
                    logger.info("Performing relative flux calibration based on the ratio of exposure times.")
                    flambda = std_exptime / sci_exptime
                    absolute_fluxcalib = False

                logger.debug('Flambda: %s', flambda)

                wstart, wend, ndim = get_wave_range(stdfiles[i])

                logger.info('Making a %dK blackbody...', temperature)
                iraf.mk1dspec(
                    bbody[i], ap=1, rv=0.0, z='no', title='', ncols=ndim, naps=1, header='', wstart=wstart, wend=wend,
                    continuum=1000, slope=0.0, temperature=temperature, fnu='no', lines='', nlines=0,
                    profile='gaussian', peak=-0.5, gfwhm=20.0, lfwhm=20.0, seed=1, comments='yes', mode='ql')

                if orders[i] == 3:  # Scale the blackbody to the Telluric K-band magnitude

                    bbmean = float(iraf.imstat(
                        images=bbody[i], fields="mean", lower='INDEF', upper='INDEF', nclip=0, lsigma=3.0, usigma=3.0,
                        binwidth=0.1, format='yes', cache='no', mode='al', Stdout=1)[1])
                    logger.debug('bbody mean: %s', bbmean)
                    bb_scale = flambda / bbmean

                else:  # for the other orders match the region of overlap with the previous order

                    logger.debug('Calculating the region of overlap with the previous order...')

                    # The list of orders [3, 4, 5, ...] goes from longer wavelengths to shorter wavelengths, so the
                    # region of overlap between order N and N-1 is the region betwen wstart(N-1) and wend(N).

                    wstart_previous, wend_previous, ndim_previous = get_wave_range(stdfiles[i-1])

                    logger.debug('Overlap: %s - %s', wstart_previous, wend)

                    # wstart_previous should be < wend. If it's not, then there's something weird with the wavelength
                    # range of the data this isn't necessarily a big problem, scientifically, but this script can't
                    # handle it so check for this and exit if there is no overlap between any two orders:

                    if wstart_previous > wend:
                        logger.error("Orders %d and %d do not overlap in wavelength.", orders[i-1], orders[i])
                        logger.error("This is unusual and suggests that the grating was not at the expected position.")
                        logger.error("This may not be a problem for the scientific use of the data, but this script")
                        logger.error("cannot handle this and is not able to flux calibrate the spectral orders.")
                        logger.error("Please plot INTERMEDIATE/calibrated_arc.fits to see if the data cover the")
                        logger.error("wavelength range you need.")
                        raise SystemExit

                    # Find the mean in the overlapping wavelength region

                    [os.remove(filename) for filename in ['overlap_bbody.fits', 'overlap_bbscaled.fits']]

                    iraf.scopy(input=bbody[i], output='overlap_bbody', w1=wstart_previous, w2=wend,
                               apertures='', beams='', bands='', apmodulus=0, format='multispec', renumber='no',
                               offset=0, clobber='no', merge='no', rebin='yes', verbose='yes', mode='ql')

                    iraf.scopy(input=bbody_scaled[i-1], output='overlap_bbscaled', w1=wstart_previous, w2=wend,
                               apertures='', beams='', bands='', apmodulus=0, format='multispec', renumber='no',
                               offset=0, clobber='no', merge='no', rebin='yes', verbose='yes', mode='ql')

                    bbody_mean = float(iraf.imstat(
                        images='overlap_bbody', fields='mean', lower='INDEF', upper='INDEF', nclip=0,
                        lsigma=3.0, usigma=3.0, binwidth=0.1, format='yes', cache='no', mode='al', Stdout=1)[1])
                    logger.debug('bbody_mean: %s', bbody_mean)

                    bbody_scaled_mean = float(iraf.imstat(
                        images='overlap_bbscaled', fields='mean', lower='INDEF', upper='INDEF', nclip=0,
                        lsigma=3.0, usigma=3.0, binwidth=0.1, format='yes', cache='no', mode='al', Stdout=1)[1])
                    logger.debug('bbody_scaled_mean: %s', bbody_scaled_mean)

                    bb_scale = bbody_scaled_mean / bbody_mean

                logger.info("Blackbody scale factor for order %d: %s", orders[i], bb_scale)

                iraf.imarith(
                    operand1=bbody[i], op="*", operand2=bb_scale, result=bbody_scaled[i], title='',
                    divzero=0.0, hparams='', pixtype='', calctype='', verbose='yes', noact='no', mode='al')

                # TODO: Is this correct?  This is not what I expected...

                iraf.imarith(
                    operand1=infile[i], op="*", operand2=bbody_scaled[i], result=outfile[i], title='',
                    divzero=0.0, hparams='', pixtype='', calctype='', verbose='yes', noact='no', mode='al')

                if absolute_fluxcalib:
                    add2header(outfile[i], 'FUNITS', 'erg/cm^2/s/A')
                else:
                    add2header(outfile[i], 'FUNITS', 'Flambda, relative')

                #iraf.wmef(input=outfile[i], output=fluxCalibrationOutput_MEF, extnames=', phu=infile[i], verbose='yes', mode='al')

                if extractFullSlit:
                    logger.error('Unsupported mode')
                    raise SystemExit

                if extractStepwise:
                    logger.error('Unsupported mode')
                    raise SystemExit

        logger.info(' --------------------------- ')
        logger.info('| Flux calibration complete |')
        logger.info(' --------------------------- ')

    os.chdir(path)  # Return to directory script was begun from.

    return


# ----------------------------------------------------------------------------------------------------------------------
def add2header(image, field, value):
    iraf.hedit(images=image, fields=field, value=value, add='yes', addonly='no',
               delete='no', verify='no', show='no', update='yes')
    return


# ----------------------------------------------------------------------------------------------------------------------
def get_wave_range(image):
    # Get the wavelength range of a spectrum
    logger = log.getLogger('get_wave_range')
    refpix = float(iraf.hselect(images=image, fields='CRPIX1', expr='yes', missing='INDEF', mode='al', Stdout=1)[0])
    logger.debug('refpix: %s', refpix)
    refval = float(iraf.hselect(images=image, fields='CRVAL1', expr='yes', missing='INDEF', mode='al', Stdout=1)[0])
    logger.debug('refval: %s', refval)
    dwave = float(iraf.hselect(images=image, fields='CD1_1', expr='yes', missing='INDEF', mode='al', Stdout=1)[0])
    logger.debug('dwave: %s', dwave)
    ndim = int(iraf.hselect(images=image, fields='NAXIS1', expr='yes', missing='INDEF', mode='al', Stdout=1)[0])
    logger.debug('ndim: %s', ndim)
    wstart = refval - (refpix-1) * dwave
    wend = wstart + float(ndim) * dwave
    logger.debug('Wavelength Range: %s - %s', wstart, wend)
    return wstart, wend, ndim


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    log.configure('gnirs-pype.log', filelevel='INFO', screenlevel='DEBUG')
    start('gnirs-pype.cfg')
