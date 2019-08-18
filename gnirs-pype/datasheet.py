#!/usr/bin/env python

from astropy.io import fits
import ConfigParser
import log
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
from matplotlib.pyplot import table
from pyraf import iraf
import numpy
import os


# ----------------------------------------------------------------------------------------------------------------------
def write(configfile):

    logger = log.getLogger('datasheet.write')

    # iraf.gemini(_doprint=0, motd="no")
    # iraf.gnirs(_doprint=0)
    # iraf.imutil(_doprint=0)
    iraf.onedspec()  # (_doprint=0)
    # iraf.nsheaders('gnirs')  # , Stdout='/dev/null')

    config = ConfigParser.RawConfigParser()
    config.optionxform = str  # make options case-sensitive
    config.read(configfile)

    for path, process in config.items("ScienceDirectories"):  # Returns list of (variable, value) pairs
        logger.debug('%s = %s', path, process)

        if not process:
            logger.debug('Skipping %s', path)
            continue

        sci_data = imexam(path)
        tel_data = imexam(path + '/Telluric')

        # sci_data['SNR'] = estimate_snr('?')                         # XDGNIRS -> flam1.fits
        tel_data['SNR'] = estimate_snr(path + '/Telluric/Intermediate/hvsrc_comb_order1_SEF.fits')  # ftell_nolines1

        sci_data['PARANGLE'] = parallactic(dec=float(sci_data['DEC']),
                                           ha=hms2deg(sci_data['HA']),
                                           lat=location(sci_data['OBSERVAT'])['latitude'],
                                           az=float(sci_data['AZIMUTH']), units='degrees')

        tel_data['PARANGLE'] = parallactic(dec=float(tel_data['DEC']),
                                           ha=hms2deg(tel_data['HA']),
                                           lat=location(tel_data['OBSERVAT'])['latitude'],
                                           az=float(tel_data['AZIMUTH']), units='degrees')

        print 'SCI:', sci_data
        print 'TEL:', tel_data

    return


# ----------------------------------------------------------------------------------------------------------------------
def imexam(path, ypos=340):
    """
    Measure the spectrum peak and FWHM
    :param path: target path of image to measure
    :param ypos: Y-position to perform measurements [340 pix]
    :return: dictionary of measurements {? overkill for only 2 values ?}
    """
    logger = log.getLogger('datasheet.imexam')

    fits_file = 'Intermediate/src_comb.fits'
    ap_file = 'Intermediate/database/apsrc_comb_SCI_1_'

    original_path = os.getcwd()
    iraf.chdir(path)  # This shouldn't be necessary, but IRAF has path length limits

    with open(ap_file, 'r') as f:
        for line in f.readlines():
            if 'center' in line:
                xpos = float(line.split()[1])
                break
    logger.debug('Spectrum X-position: %.2f pix', xpos)

    cursor = 'tmp.cur'  # Write a cursor file for imexam
    with open(cursor, 'w') as f:
        f.write('%.3f %.3f\n' % (xpos, ypos))

    logger.info('Running IRAF imexam to measure the spectrum peak and FHWM...')
    iraf.unlearn(iraf.imexam)
    iraf.unlearn(iraf.jimexam)    # jimexam = 1-dimensional gaussian line fit
    # iraf.jimexam.naverage = 50  # Number of lines, columns, or width perpendicular to a vector to average
    # iraf.jimexam.width = 100    # Width of background region for background subtraction (pix)
    # iraf.jimexam.rplot = 100    # Radius to which the radial profile or 1D profile fits are plotted (pix)
    # iraf.jimexam.sigma =        # Initial sigma (pix)
    # iraf.imexam.graphics = 'stgkern' #  Force the use of the standard IRAF graphics kernel
    logger.debug('iraf.jimexam.sigma = %0.3f', iraf.jimexam.sigma)
    logger.debug('iraf.jimexam.naverage = %0.3f', iraf.jimexam.naverage)
    logfile = 'tmp.log'

    iraf.imexam(
        input=fits_file + '[SCI,1]', frame=1, output='', logfile=logfile, keeplog='yes', defkey='j',
        ncstat=5, nlstat=5, imagecur=cursor, use_display='no', Stdout=1)

    logger.debug('Parsing imexam results from the log file...')
    peak = None
    fwhm = None
    with open(logfile) as f:
        for line in f:
            if '#' in line:
                continue
            logger.debug('%s', line.strip())
            vals = line.replace('=', ' ').split()
            if vals[0] == 'Lines':      # record measure of x center
                center = float(vals[3])
                peak = float(vals[5])
                fwhm = float(vals[9])
                break
    logger.debug('center = %s  peak = %s  fwhm = %s', center, peak, fwhm)
    data = {'peak': peak, 'fwhm': fwhm}

    logger.debug('Cleaning up...')
    for f in [cursor, logfile]:
        os.remove(f)

    logger.debug('Reading some FITS header keywords...')
    header = fits.open(fits_file)[0].header
    for key in ['OBJECT', 'GEMPRGID', 'AIRMASS', 'RA', 'DEC', 'HA', 'AZIMUTH', 'PA',
                'OBSERVAT', 'RAWIQ', 'RAWCC', 'RAWWV', 'RAWBG', 'DATE-OBS']:
        try:
            data[key] = header[key].strip() if isinstance(header[key], str) else header[key]
        except:
            logger.warning('%s[%s] is undefined', f, key)
            data[key] = None

    iraf.chdir(original_path)

    logger.debug('data: %s', data)
    return data


# ----------------------------------------------------------------------------------------------------------------------
def estimate_snr(onedspectrum, wav1=21000, wav2=22000, interactive=False):
    """
    Estimate Signal-to-Noise ratio

    :param onedspectrum: input one-dimensional (extracted) spectrum
    :param wav1: starting wavelength rof ange to fit and measure
    :param wav2: ending wavelength of range to fit and measure
    :param interactive:
    :return: signal-to-noise ratio (float)
    """

    logger = log.getLogger('datasheet.snr')
    logger.info('Estimating S/N...')

    output = 'tmp.fits'
    stdout = 'tmp.out'
    cursor = 'tmp.cur'
    logfile = 'tmp.log'

    with open(cursor, 'w') as f:  # Generate a cursor file for bplot
        f.write('%d 0 1 m\n' % wav1)
        f.write('%d 0 1 m\n' % wav2)
        f.write('q')

    logger.debug('continuum input: %s', onedspectrum)
    logger.debug('sample: %d:%d', wav1, wav2)

    iraf.sfit.logfile = logfile
    iraf.continuum(
        input=onedspectrum, output=output, lines='*', bands='1', type='ratio', replace=False, wavescale=True,
        logscale=False, override=False, logfile=logfile, interactive=interactive, sample='%d:%d' % (wav1, wav2),
        naverage=1, function='spline3', order=3, low_rej=2, high_rej=3, niterate=5, grow=1)

    iraf.splot.save_file = logfile
    iraf.bplot(
        images=output, apertures="", band=1, cursor=cursor, next_image="",
        new_image="", overwrite="no", spec2="", constant=0.0, wavelength=0.0, linelist="",
        wstart=0.0, wend=0.0, dw=0.0, boxsize=2, Stdout=stdout)  # graphics="stgkern", StdoutG="dev$null")

    logger.debug('Parsing output...')
    snr = None
    with open(stdout, 'r') as f:
        for line in f.readlines():
            if 'snr' in line:
                snr = float(line.split()[-1])
    logger.debug('SNR: %s', snr)

    for f in [cursor, logfile, output, stdout]:
        os.remove(f)

    return snr


# ----------------------------------------------------------------------------------------------------------------------
def parallactic(dec, ha, lat, az, units='degrees'):
    """
    Compute the parallactic angle
    :param dec:  target declination
    :param ha:   hour angle
    :param lat:  observatory latitude
    :param az:   target azimuth
    :param units:  degrees or radians for input and ouput quantities
    :return:     parallactic angle (float)
    """
    logger = log.getLogger('parallactic')

    if units == 'degrees':
        dec *= numpy.pi / 180.
        ha *= numpy.pi / 180.
        lat *= numpy.pi / 180.
        az *= numpy.pi / 180.

    if numpy.cos(dec) != 0.0:
        sinp = -1.0*numpy.sin(az)*numpy.cos(lat)/numpy.cos(dec)
        cosp = -1.0*numpy.cos(az)*numpy.cos(ha)-numpy.sin(az)*numpy.sin(ha)*numpy.sin(lat)
        pa = numpy.arctan2(sinp, cosp)
    else:
        if lat > 0.0:
            pa = numpy.pi
        else:
            pa = 0.0

    if units == 'degrees':
        pa *= 180. / numpy.pi

    logger.debug('Parallactic Angle: %.3f %s', pa, units)
    return pa


# ----------------------------------------------------------------------------------------------------------------------
def hms2deg(angle):
    """Convert sexagesimal HH:MM:SS.sss to decimal degrees"""
    h, m, s = angle.split(':')
    hours = float(h) + float(m)/60. + float(s)/3600.
    return hours / 24. * 360.


# ----------------------------------------------------------------------------------------------------------------------
def location(observatory):
    """Return the observatory location as a dictionary"""
    if observatory == 'Gemini-North':
        latitude = 297.35709        # 19:49:25.7016
        longitude = -155.46906      # -155:28:08.616
        elevation = 4213            # meters
    elif observatory == 'Gemini-South':
        latitude = 453.61125        # -30:14:26.700
        longitude = -70.7366933333  # -70:44:12.096
        elevation = 2722            # meters
    else:
        raise SystemExit('Unknown observatory')

    return {'latitude': latitude, 'longitude': longitude, 'elevation': elevation}


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    log.configure('gnirs.log', filelevel='INFO', screenlevel='DEBUG')
    write('gnirs.cfg')
