#!/usr/bin/env python

# Incorporated Daniel Ruschel Dutra's code into XDGNIRS, July 2014 -REM

################################################################################
# CHANGE LOG                                                                   #
#                                                                              #
# 23 Jul 2014 - * Moved the nc term to the end of the SNR equation to increase #
#               readability.                                                   #
#               * Removed the 1/2 factor multiplying the nc term, because the  #
#               number of images in NSCHL0* cards is already the amount of     #
#               source exposures.                                              #
#               * Removed the import statement for spectools. Since the        #
#               the functions in this script no longer write ASCII files,      #
#               there is no need to build an array of wavelengths.             #
#                                                                              #
# KNOWN ISSUES                                                                 #
#                                                                              #
# - Limit to 99 combined images: The nc term is the result of a search in the  #
# header for NSCHL0* cards, which end in a two digit counter, e.g. NSCHL001.   # 
# Therefore, there is an artificial limitation to 99 combined source images    #
# until a better way of determining the number of combined exposures is        #
# implemented. Nevertheless, this should no be a problem since the number of   #
# combined exposures is usually around 10.                                     #
#                                                                              #
################################################################################

from astropy.io import fits
import ConfigParser
from copy import deepcopy
import log
import numpy
import re
import os.path
from pyraf import iraf
from scipy.interpolate import interp1d


# ----------------------------------------------------------------------------------------------------------------------
def get_wl(image, dimension=0, hdrext=0, dataext=0, dwlkey='CD1_1', wl0key='CRVAL1', pix0key='CRPIX1'):
    """
    Obtains the wavelength coordinates from the header keywords of the
    FITS file image. The default keywords are CD1_1 for the delta lambda,
    CRVAL for the value of the first pixel and CRPIX1 for the number of
    the first pixel. These keywords are the standard for GEMINI images.

    The function is prepared to work with Multi-Extesion FITS (MEF) files.

    Parameters:
    -----------
    image : string
        Name of the FITS file containing the spectrum
    dimension : integer
        Dimension of the dispersion direction
    hdrext : number
        Extension that contains the header
    dataext : number
        Extension that contains the actual spectral data
    dwlkey : string
        Header keyword for the interval between two data points
    wl0key : string
        Header keyword for the first pixel value
    pix0key : string
        Header keyword for the first pixel coordinate

    Returns:
    --------
    wl : numpy.array
        Wavelength coordinates for each data point
    """

    h = fits.getheader(image, ext=hdrext)
    dwl, wl1, pix0 = [float(h[i]) for i in [dwlkey, wl0key, pix0key]]
    npoints = numpy.shape(fits.getdata(image, dataext))[dimension]
    wl = wl1 + (numpy.arange(1, npoints + 1) - pix0) * dwl

    return wl


# ----------------------------------------------------------------------------------------------------------------------
def build_unc(inimage, skyspec, outimage=None, clobber=True, writefits=True, aperture_width=12, skyglow=True,
              readnoise=None, dark=0.15, verbose=False):
    """
  Builds a Poisson noise estimate for the all the extensions in the image ''inimage''.

  Parameters:
  -----------
  inimage : string
    Name of the input image.
  outimage : string
    Name of the output image.
  clobber : boolean
    Overwrite existing images.
  writefits : boolean
    Writes the results in a FITS image.
  aperture_width : float
    Default aperture width to use when there is no
    APNUM1 card in the header.
  skyglow : boolean
    Includes the emission lines from the atmosphere
    in the noise estimate.
  readnoise : float
    Reading noise in electrons per pixel. If None
    its value will be read from the 'rdnoise' card in
    the image header.
  dark : float
    Detector's dark current in electrons per pixel
    per second. The standard value comes from the manual
    of GNIRS, and is 0.15.
  skyspec : string
    Atmospheric emission spectrum in units of ADU per second
    per pixel.
  
  Returns:
  --------
  snrlst : list
    A list containing all the SNR spectra, beginning with the
    one with longer wavelength (K).
  Actually I changed this to return the name of the output fits spectrum for use in the order-joining function -REM

  Description:
  ------------
  This function takes a FITS file with 18 extensions and
  writes a copy of it with every 2+3n extension replaced
  by a S/N estimate.

  The estimate relies on the 1+3n extensions being the
  spectrum in electrons s^-1 units. According to Steve B. Howell,
  in Handbook of CCD Astronomy (2nd ed. page 73), and also
  Mortara & Fowler 1981 (SPIE Conference Proceedings, Vol. 290.,
  page 28), the signal to noise ratio is given by

  S/N = (Nx)/sqrt(Nx + npix*(Nd*t + Ns*t + Nr**2))

  where Nx is the number of electrons from the astronomical
  source, t is the exposure time, npix is the number of pixels
  in the aperture, Nd is the dark current in electrons per pixel
  per second, Ns is the number of electrons from the
  atmospheric emission per pixel and per second and Nr is the
  read noise times the square root of the number of exposures.

  By comparison with data from the GNIRS Integration Time Calculator
  it has been found that the data at the 'vtarget_comb.fits' stage
  is already in units of electrons, therefore no gain multiplication
  is being performed by this function.
  """

    im = fits.open(inimage)
    h = im[0].header
    nc = len(h['NSCHL0*'])  # number of combined frames
    # sky = pf.getdata(skyspec)
    sky = fits.open(skyspec)
    # sky[sky < 1.e-3] = 1.e-3  # to avoid divisions by zero
    snrlst = []

    for i in numpy.arange(1, 17, 3):
        hext = im[i].header
        t = float(hext['exptime'])
        s = deepcopy(im[i].data)
        sky2 = deepcopy(sky[i].data)  # this
        if readnoise is None:
            readnoise = 7.  # assuming the best possible read noise from GNIRS' manual
        rn = readnoise
        s[s < 1.e-3] = 1.e-3  # to avoid divisions by zero
        sky2[sky2 < 1.e-3] = 1.e-3  # this

        # getting the aperture width
        try:
            npix = numpy.diff([float(k) for k in hext['apnum1'].split()[2:]])
        except KeyError:
            if verbose:
                print 'WARNING! APNUM1 card not found in {:s}.'.format(inimage + '[' + str(i) + ']')
                print 'Proceeding with default value of {:.1f} pixels.'.format(aperture_width)
            npix = aperture_width

        if skyglow:
            # for sky spectrum extracted by XDGNIRS from science data (same aperture and exposure time)
            snr = s / numpy.sqrt(s + sky2 + npix * (+dark * t + rn ** 2)) * numpy.sqrt(nc)
            # canned sky spectrum in electrons/second/pix
            # snr = s/sqrt(s+npix*(sky[(i-1)/3]*t+dark*t+rn**2)) * sqrt(nc)
        else:
            snr = s / numpy.sqrt(s + npix * (dark * t + rn ** 2)) * numpy.sqrt(nc)

        # Temporarily writing out K-band input spectrum, sky, and snr for plotting
        if i == 1:
            writeout = open('noise_info.txt', 'w')
            for k in range(0, len(snr)):
                writeout.write(str(s[k]) + "  " + str(sky2[k]) + "  " + str(snr[k]) + '\n')
            writeout.close()

        im[i + 1].data = snr
        snrlst.append(snr)
        if verbose:
            print 'npix = {:.2f}; nc = {:d}; Combined read noise = {:.2f}; Texp = {:.2f}'.format(float(npix), nc, rn, t)

    if writefits:
        if outimage is None:
            outimage = inimage[:inimage.find('.fits')] + '_unc.fits'
        im.writeto(outimage, output_verify='silentfix', clobber=clobber)

    # Outimage contains the input source spectrum in every 3n+1 extension, and its SNR spectrum in every 3n+2 extension -REM
    return outimage


# ----------------------------------------------------------------------------------------------------------------------
def join_extensions(infile, reference=None, outjoined=None, outname=None):
    """
    Joins in to a single spectrum the different orders in a GNIRS XD
    spectrum, with one order from each FITS extension. This functions
    is specifically aimed at bringing together the SNR spectrum with
    the flux calibrated one.

    Parameters:
    -----------
    infile : string
        Input FITS file with SNR spectra in the 2+3n extensions.
    reference : string
        Name of the FITS file containing the combined, flux calibrated and (if available) redshift-corrected spectrum.
    outjoined : string
        Name of the "iraf.noao.onedspec.odcombine" output.
    output : string
        Name of the final output file, which consists of an ASCII with
        three columns: wavelength, flux density, error in flux density.

    Returns:
    --------
    Name of the output spectrum -REM
    """
    logger = log.getLogger('join_extensions')

    s = infile[:infile.find('_unc.fits')]

    if reference is None:
        try:
            srest = '../PRODUCTS/' + s + '_rest.fits'  # FIXME
        except IOError:
            srest = srest.replace('NGC', 'ngc')
    else:
        srest = reference

    h = fits.getheader(srest)['IMC*']
    hfull = fits.getheader(srest)

    if fits.getheader(infile)['object'] == fits.getheader(srest)['object']:
        pass
    else:
        logger.error('Objects in input and reference files do not match.')
        print 'Input object: {:s}'.format(fits.getheader(infile)['object'])
        print 'Reference object: {:s}'.format(fits.getheader(srest)['object'])
        print 'Exiting function now.'
        return

    gb = [str(h[i])[str(h[i]).find('[') + 1:str(h[i]).find(']')] for i in range(6)]

    with open('listod.lst', 'w') as f:
        for i in numpy.arange(6):
            f.write(infile + '[' + str(2 + 3 * i) + '][' + gb[i] + ']\n')

    iraf.noao()
    iraf.onedspec()
    if outjoined is None:
        outjoined = srest.replace('.fits', '') + '_unc_joined.fits'
    if os.path.isfile(srest.replace('.fits', '') + '_unc_joined.fits'):
        iraf.delete(files=srest.replace('.fits', '') + '_unc_joined.fits', verify="no", Stdout="/dev/null")
        iraf.onedspec.interp = 'linear'

    iraf.odcombine(
        input='@listod.lst', output=outjoined, headers='', bpmasks='', rejmasks='', nrejmasks='',
        expmasks='', sigmas='', logfile='STDOUT', apertures='', group='all', first='no', w1='INDEF',
        w2='INDEF', dw='INDEF', nw='INDEF', log='no', combine='average', reject='none', outtype='real',
        outlimits='', smaskformat='bpmspectrum', smasktype='none', smaskvalue=0.0, blank=0.0, scale='none',
        zero='none', weight='none', statsec='', expname='', lthreshold='INDEF', hthreshold='INDEF', nlow=1,
        nhigh=1, nkeep=1, mclip='yes', lsigma=3.0, hsigma=3.0, rdnoise='0.0', gain='1.0', snoise='0.0',
        sigscale=0.1, pclip=-0.5, grow=0.0, offsets='physical', masktype='none', maskvalue=0.0, mode='al')

    srestdata = fits.getdata(srest)
    suncdata = fits.getdata(outjoined)

    if len(srestdata) != len(suncdata):
        print '*****************'
        print 'ERROR! Length of flux calibrated spectrum does not match the'
        print 'length of combined uncertainty spectrum.'
        print '{:s} : {:d}'.format(srest, len(srestdata))
        print '{:s} : {:d}'.format(outjoined, len(suncdata))
        print '*****************'
        return

    if outname is None:
        outname = srest.replace('.fits', '') + '_unc.fits'

    # savetxt(output,column_stack([wl,srestdata,srestdata/suncdata]),fmt='%.6e')
    #outh = deepcopy(hfull)
    #outh['CTYPE2'] = 'LINEAR'
    #outh['CRVAL2'] = 1
    #outh['CRPIX2'] = 1
    #outh['CDELT2'] = 1
    #outh['CD2_2'] = 1
    #outh['WCSDIM'] = 2
    #outh['LTM2_2'] = 1
    #outh['APNUM2'] = '2 2'
    #outh['BANDID1'] = 'Spectrum'
    #outh['BANDID2'] = 'Flux uncertainty'
    #pf.writeto(outname, data=row_stack([srestdata, srestdata / suncdata]), header=outh)

    # Replaced nice python code with nasty iraf code as original was crashing with "keyword not found" errors -REM

    if os.path.isfile(outname):
        iraf.delete(files=outname, verify="no", Stdout="/dev/null")

    # Divide the signal by the SNR to get the noise spectrum in erg/cm2/s/A
    iraf.imarith(operand1=srest, op='/', operand2=outjoined, result='snoise', title='Noise', divzero=0.0, hparams='',
                 pixtype='', calctype='', verbose='no', noact='no')

    # Put spectrum, noise, and snr into three extensions of a FITS file
    iraf.fxcopy(input=srest, output=outname, group="", new_file='yes', verbose='yes', mode='ql')
    iraf.fxinsert(input='snoise', output=outname + '[1]', groups="", verbose='yes', mode='ql')
    iraf.fxinsert(input=outjoined, output=outname + '[2]', groups="", verbose='yes', mode='ql')

    # Write text files
    iraf.wspectext(input=outname + '[0]', output=outname.replace('_unc.fits', '_spec.txt'), header='no', wformat='')
    iraf.wspectext(input=outname + '[1]', output=outname.replace('_unc.fits', '_noise.txt'), header='no', wformat='')
    iraf.wspectext(input=outname + '[2]', output=outname.replace('_unc.fits', '_snr.txt'), header='no', wformat='')

    return outname


# ----------------------------------------------------------------------------------------------------------------------
def gen_sn_spec(target, reference, sky):
    """
    Generate single S/N spectrum, all extensions combined using same pixels as in combined flux spectrum.
    :param target: the combined science target data file
    :param reference: the reference spectrum file (to determine which pixels went into the combined science spectrum)
    :param sky: the sky spectrum file
    :return:
    """
    logger = log.getLogger('gen_sn_spec')
    logger.debug('target: %s', target)
    logger.debug('reference: %s', reference)
    logger.debug('sky: %s', sky)
    sn_spec = build_unc(target, sky)
    final_spec = join_extensions(sn_spec, reference=reference)
    logger.info('%s contains the flux-calibrated source spectrum in ext[0],'
                'its noise spectrum in ext[1], and the SNR spectrum in ext[2]', final_spec)
    return


# ----------------------------------------------------------------------------------------------------------------------
def start(configfile):

    logger = log.getLogger('noise_spectrum')

    config = ConfigParser.RawConfigParser()
    config.optionxform = str  # make config file options case-sensitive
    config.read(configfile)

    for sdir in config.options('ScienceDirectories'):
        if config.getboolean('ScienceDirectories', sdir):  # Only use directories marked True

            # Here we nned to figure out if there are nods to sky, i.e. if the nod size is > slit size; see inslit()
            nod_to_sky = True  # FIXME

            if nod_to_sky:

                obj = re.sub('[^a-zA-Z0-9]', '', fits.getheader(sdir + '/Intermediate/vsrc_comb.fits')['OBJECT'])

                if config.getboolean('combineOrdersXD', 'shiftToRestframe'):
                    refspec = sdir + '/Final/' + obj + '_rest.fits'
                else:
                    refspec = sdir + '/Final/' + obj + '.fits'

                gen_sn_spec(
                    target=sdir + '/Intermediate/vsrc_comb.fits',
                    reference=refspec,
                    sky=sdir + '/Intermediate/vsky_comb.fits')

            else:
                logger.warning('Calculation of the noise spectrum is not possible without nods to sky.')

    return


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    log.configure('gnirs-pype.log', filelevel='INFO', screenlevel='DEBUG')
    start('gnirs-pype.cfg')
