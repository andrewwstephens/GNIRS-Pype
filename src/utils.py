#!/usr/bin/env python

from astropy.io import fits
import log
from matplotlib import pyplot
import os
import re


# ----------------------------------------------------------------------------------------------------------------------
def band(order):
    # Give the order to band mapping.  band[4] = H
    return {1: 'M', 2: 'L', 3: 'K', 4: 'H', 5: 'J', 6: 'J', 7: 'J', 8: 'J'}[order]


# ----------------------------------------------------------------------------------------------------------------------
def boxit(text, character='-', center=True):
    lines = text.split('\n')
    max_length = max([len(line) for line in lines])
    border = character * (max_length + 4)
    print("%s" % border)

    for line in lines:
        dif = (max_length-len(line))
        if dif == 0:
            print(character + " %s " % line + character)
        else:
            if center:
                print(character + ' ' * (dif/2) + ' %s ' % line + ' ' * (dif/2) + character)
            else:
                print(character + " %s " % line + ' ' * dif + character)

    print("%s" % border)
    return


# ----------------------------------------------------------------------------------------------------------------------
def clean(filelist, outputPrefix, overwrite):
    """
    Run cleanir to remove any pattern noise.
    :param filelist: name of input file list
    :param outputPrefix:
    :param overwrite:
    :return:
    """
    logger = log.getLogger('clean')

    infiles = files_in(filelist)
    requires(infiles)

    outfiles = ['c' + f for f in infiles]
    if exists(outfiles, overwrite):
        logger.info('All files already cleaned')
        return

    logger.error('CLEANIR NOT IMPLEMENTED')  # FIXME
    # cleanir.py inlist

    return


# ----------------------------------------------------------------------------------------------------------------------
def dictify(itemlist, fmt='str'):
    # Convert a ConfigParser items list [(item1,value1),(item2,value2),...] to a dictionary {item1:value1,...}
    # fmt specifies format of the KEY (not the value)
    logger = log.getLogger('dictify')
    logger.debug('itemlist: %s', itemlist)
    d = {}
    for key, value in itemlist:
        try:
            v = float(value)
        except:
            v = value
        if fmt == 'int':
            d[int(key)] = v
        else:
            d[key] = v
    logger.debug('dict: %s', d)
    return d


# ----------------------------------------------------------------------------------------------------------------------
def exists(inlist, overwrite=False):
    """
    Check for the existence of the files in the input list.
    If overwrite == True then delete any existing files and return False
    if overwrite == False:
    - if no files exist then return False
    - if all the files exist then return True, indicating that this reduction step has been done and can be skipped
    - if some (but not all) of the files exist then throw an error because something failed somewhere
    The goal is to be able to say something like:  if exists() then skip this step, or if not exists() do stuff.
    :param inlist: list of files
    :param overwrite: boolean
    :return: boolean
    """
    logger = log.getLogger('exists')
    logger.debug('inlist: %s', inlist)
    logger.debug('overwrite: %s', overwrite)

    exist = [os.path.exists(f) for f in inlist]
    logger.debug('exists: %s', exist)

    if overwrite:
        for i in range(len(inlist)):
            if exist[i]:
                logger.warning('Removing the old %s', inlist[i])
                os.remove(inlist[i])
        return False

    else:
        if all(exist):
            logger.debug("All files exist")
            return True

        elif not any(exist):
            logger.debug('No files exist')
            return False

        else:
            logger.error('Some (but not all) the files exist and overwrite = False')
            raise SystemExit


# ----------------------------------------------------------------------------------------------------------------------
def requires(filelist):
    """
    Check that all the required files exist, and throw an error if any are missing.
    :param filelist: list of file names
    :return:
    """
    logger = log.getLogger('requires')
    logger.debug('Files: %s', filelist)

    exist = [os.path.exists(f) for f in filelist]
    logger.debug('Exist: %s', exist)

    if not all(exist):
        logger.error('Some required files are missing:')
        for f, e in zip(filelist, exist):
            if not e:
                logger.error('Cannot find %s', f)
        raise SystemExit
    return


# ----------------------------------------------------------------------------------------------------------------------
def get_bpm(filename):
    """
    Find the most appropriate bad pixel mask. For data taken before the summer 2012 lens replacement use
    gnirs$data/gnirsn_2011apr07_bpm.fits.  After summer 2012 use gnirs$data/gnirsn_2012dec05_bpm.fits'.
    Use keyword 'ARRAYID' in the raw file header to check which camera was used for observations:
    """
    logger = log.getLogger('get_bpm')
    arrayid = fits.getheader(filename)['ARRAYID']
    logger.debug('ARRAYID: %s', arrayid)
    if arrayid == 'SN7638228.1':
        bpm = 'gnirs$data/gnirsn_2011apr07_bpm.fits'
    elif arrayid == 'SN7638228.1.2':
        bpm = 'gnirs$data/gnirsn_2012dec05_bpm.fits'
    else:
        logger.error('Unknown array ID.')
        raise SystemExit
    logger.info('BPM: %s', bpm)
    return bpm


# ----------------------------------------------------------------------------------------------------------------------
def get_orders(path):
    logger = log.getLogger('get_orders')

    if 'LB_SXD' in path:
        orders = [3, 4, 5]
    elif 'LB_LXD' in path:
        orders = [3, 4, 5, 6, 7, 8]
    elif 'SB_SXD' in path:
        orders = [3, 4, 5, 6, 7, 8]
    else:
        logger.error('Unknown GNIRS configuration.  Cannot determine the XD orders.')
        raise SystemExit

    logger.debug('orders: %s', orders)
    return orders


# ----------------------------------------------------------------------------------------------------------------------
def get_wavelengths(path):
    logger = log.getLogger('get_wavelengths')

    if 'LB_SXD' in path:
        wavelengths = None
        accuracy = None
    elif 'LB_LXD' in path:
        wavelengths = None
        accuracy = None
    elif 'SB_SXD' in path:
        wavelengths = ([0,0], [18690,25310], [14020,18980], [11220,15180], [9350,12650], [8020,10840], [7020,9480])
        accuracy = 5  # percent
    else:
        logger.error("Unknown GNIRS configuration.  Cannot determine the nominal wavelenghts.")
        raise SystemExit

    logger.debug('wavelengths: %s', wavelengths)
    logger.debug('accuracy: %s', accuracy)
    return wavelengths, accuracy


# ----------------------------------------------------------------------------------------------------------------------
def files_in(filelist):
    """
    Return a list of all the files in the list of supplied file list names.
    :param filelist: list of file-lists, e.g. ['pinholes.list', 'arcs.list']
    :return: list of all the filts, e.g. ['N20190102S0123', 'N20190102S0124', ...]
    """
    logger = log.getLogger('files_in')
    logger.debug('Filelist: %s', filelist)
    allfiles = []
    for fl in filelist:
        with open(fl, 'r') as f:
            allfiles.extend(f.read().splitlines())
    logger.debug('Files: %s', allfiles)
    return allfiles


# ----------------------------------------------------------------------------------------------------------------------
def make_list(prefix, orders=None, suffix='_MEF[1]', regions=None):
    """
    Take a filename and return a list of all the orders (filename_orderN)
    # If the regions dictionary is specified then append the appropriate order, e.g. filename_orderN[xxx:yyy]
    # The orders are only needed if the regions are not specified in which case we don't know the orders.
    :param prefix: string filename prefix
    :param orders: list of spectral orders, e.g. [3,4,5,6,7,8]
    :param suffix: file name suffix
    :param regions: dictionary of all orders and corresponding regions, e.g. {3:"123:456"}
    :return: list of files for each order with optional regions, e.g. [..., prefix_orderN_suffix[region], ...]
    """
    logger = log.getLogger('make_list')
    logger.debug('prefix: %s', prefix)
    logger.debug('orders: %s', orders)
    logger.debug('suffix: %s', suffix)
    logger.debug('regions: %s', regions)

    if not orders and not regions:
        logger.error('Must specify orders or regions')
        raise SystemExit

    output = []
    if regions:
        for o in sorted(regions.keys()):
            output.append('%s_order%d%s[%s]' % (prefix, o, suffix, regions[o]))
    else:
        for o in orders:
            output.append('%s_order%d%s' % (prefix, o, suffix))

    logger.debug('list: %s', output)
    return output


# ----------------------------------------------------------------------------------------------------------------------
def nofits(filename):
    return filename.replace('.fits', '')


# ----------------------------------------------------------------------------------------------------------------------
def pause(active, message=None):
    # This is a replacement for the original "Manual Mode."
    # If active == True then print the optional message and wait for the user to hit <enter>, else do nothing.
    if active:
        if message:
            print message
        answer = raw_input('Continue? [y]/n ')
        if 'n' in answer.lower():
            raise SystemExit

    return


# ----------------------------------------------------------------------------------------------------------------------
def plot(points1=None, points2=None, line1=None, line2=None, label1=None, label2=None, title=None):
    # Plot up to two spectra given the filenames with the option of specifying the extension in IRAF format
    logger = log.getLogger('plot')
    pyplot.figure()

    def read(filename):
        if '[' in filename:  # filename.fits[SCI,N]
            i1 = filename.find('[')
            i2 = filename.find(',')
            extname = filename[i1+1:i2]
            extver = int(filename[i2+1:filename.rfind(']')])
            return fits.getdata(filename[:i1], extname=extname, extver=extver)
        else:
            return fits.getdata(filename)

    if points1:
        logger.debug('Plotting %s', points1)
        pyplot.plot(read(points1), linestyle='', marker='o', label=label1 if label1 else points1)

    if points2:
        logger.debug('Plotting %s', points2)
        pyplot.plot(read(points2), linestyle='', marker='o', label=label2 if label2 else points2)

    if line1:
        logger.debug('Plotting %s', line1)
        pyplot.plot(read(line1), linestyle='-', marker='', label=label1 if label1 else line1)

    if line2:
        logger.debug('Plotting %s', line2)
        pyplot.plot(read(line2), linestyle='-', marker='', label=label2 if label2 else line2)

    pyplot.legend(loc='best', numpoints=1, fancybox=True)
    pyplot.grid(linewidth=0.25)
    if title:
        pyplot.title(title)
    pyplot.show()
    return


# ----------------------------------------------------------------------------------------------------------------------
def get_target(fitsfile):
    logger = log.getLogger('gettarget')
    obj = fits.getheader(fitsfile)['OBJECT']
    logger.debug('Object: %s', obj)
    target = re.sub('[^a-zA-Z0-9]', '', obj)  # replace non-alphanumeric characters
    logger.debug('Target: %s', target)
    return target


# ----------------------------------------------------------------------------------------------------------------------
def joinorders_noresampling(inlist, merged_spec, outspec, writefits, snrlist=None):

    """
    Join the different spectral orders from the cross-dispersed exposure
    into a single spectrum, having the wavelength coordinates as a
    lookup table in the header.

    Author: Daniel Ruschel Dutra, 16 Apr 2015

    Parameters:
    -----------
    inlist : list
    List of the names of FITS files containing the 1D spectra of each order.
    In the standard nomenclature of the XDGNIRS pipeline these would be the
    final?.fits files.
    merged_spec : string
    Name of FITS file containing the final merged spectrum, combined with
ODCOMBINE. This is required to get the wavelength limits of each order.
outspec : string
Name of the final output spectrum. The code also writes an ASCII file with
the same root name, but with '.fits' replaced with '.txt'. Therefore
outspec needs to have '.fits' at the end.
writefits : boolean
Writes the resulting stitched spectrum in a FITS file.
snrlist : list or string
If snrlist is a list, it is expected to have six elements, each one
being the name of a FITS 1D signal-to-noise ratio spectrum. The order of
this list has equal to the inlist, from longer to shorter wavelengths.
snrlist can also be the name of a MEF FITS file, where the signal-to-noise
ratio is written in the 3n+2 extension, with n from 1 to 6.


Returns:
--------
jspec : numpy.ndarray
jspec[:,0] - Wavelength coordinates
jspec[:,1] - Flux coordinates

Description:
------------
This function was developed to produce a stitched spectrum spanning the full
wavelength range available in the cross dispersed exposure of GNIRS without
resampling the different orders with the same wavelength step. Instead each
order carries the dispersion it had when it was flux calibrated, and the
final stitched spectrum has the wavelength coordinates written in a lookup
table at the image header, with 6 differente dispersions.

-- Motivation --
The ODCOMBINE routine that is employed by the pipeline can only produce a
spectrum with a linear dispersion function, defined by a fixed step and
starting coordinate. This of course leads to oversampling of the regions
with lower dispersion, or longer wavelengths. Originally the oversampling
was done with the highest dispersion between the input spectra, which led
to a sampling frequency of roughly 3 at the 6th order, and still smaller
frequencies for every other order. Such a low frequency caused the flux
values to be 'smoothed', specially when using high order polynomials for the
interpolation.

-- Methodology --
The main problem with stitching the different difraction orders was the
appropriate treatment of the overlaping regions. The method employed here
consist in storing each spectrum as a scipy.interpolate.interp1d function
and averaging both spectra at the wavelength coordinates of the one with the
higher dispersion.
"""

    if type(snrlist) == list:
        snrmef = False
    elif type(snrlist) == str:
        if len(pf.open(snrlist)) == 19:
            snrmef = True
        else:
            print 'ERROR! Expecting 19 extension MEF as snrlist.'
            return

    #
    #  Reading general info from the input files
    #
    for i in inlist:
        print i
    cdelts = array([pf.getheader(i)['CD1_1'] for i in inlist])
    # cdelts is not actually needed, but who knows?
    wls = array([get_wl(i) for i in inlist])

    #
    # Storing each spectrum as a linear interpolation function
    #
    f = array([interp1d(wls[i], pf.getdata(inlist[i])) for i in range(6)])
    if snrlist != None:
        if snrmef:
            n = array([interp1d(wls[i], pf.getdata(snrlist, ext=(3 * i + 2))) for i in range(6)])
        else:
            n = array([interp1d(wls[i], pf.getdata(snrlist[i])) for i in range(6)])

    #
    # Getting the limits of each order
    #
    h = pf.getheader(merged_spec)['IMC*']

    # gb stands for stringed good bits
    # If it takes longer than one line it is wrong.

    gb = array([h[i][h[i].find('[') + 1:h[i].find(']')].split(':') for i in range(6)], dtype='int')
    wllims = array([[wls[i][gb[i, 0]], wls[i][gb[i, 1]]] for i in range(6)], dtype='int')

    # js stands for joined spectrum
    js = array([])
    jwl = array([])

    if snrlist != None:
        jsnr = array([])

    for i in arange(6)[::-1]:

        # shortest wavelength
        if i == 5:
            singlecondition = (wls[i] > wllims[i, 0]) & (wls[i] < wllims[i - 1, 0])
            overlapcondition = (wls[i] > wllims[i - 1, 0]) & (wls[i] < wllims[i, 1])

        # longest wavelength
        elif i == 0:
            # only the single spectrum part
            singlecondition = (wls[i] > wllims[i + 1, 1]) & (wls[i] < wllims[i, 1])
            js = append(js, f[i](wls[i][singlecondition]))
            jwl = append(jwl, wls[i][singlecondition])
            if snrlist != None:
                jsnr = append(jsnr, n[i](wls[i][singlecondition]))
            break

        # everything else
        else:
            singlecondition = (wls[i] > wllims[i + 1, 1]) & (wls[i] < wllims[i - 1, 0])
            overlapcondition = (wls[i] > wllims[i - 1, 0]) & (wls[i] < wllims[i, 1])

        js = append(js, f[i](wls[i][singlecondition]))
        js = append(js, average([f[i](wls[i][overlapcondition]), f[i - 1](wls[i][overlapcondition])], 0))

        jwl = append(jwl, wls[i][singlecondition])
        jwl = append(jwl, wls[i][overlapcondition])

        if snrlist != None:
            jsnr = append(jsnr, n[i](wls[i][singlecondition]))
            # this should be revised as soon as the code starts working
            jsnr = append(jsnr,
                          average([n[i](wls[i][overlapcondition]) + n[i - 1](wls[i][overlapcondition])], 0) / sqrt(2.))

    #
    #  Writing the 1D spectrum to a FITS file with the aid of noao.onedspec.rspectext
    #
    outtext = outspec.replace('.fits', '.txt')

    if snrlist == None:
        jspec = column_stack([jwl, js])
        if writefits:
            savetxt(outtext, jspec, fmt='%.6e\t%.6e')
    else:
        jspec = column_stack([jwl, js, jsnr])
        if writefits:
            savetxt(outtext, jspec, fmt='%.6e\t%.6e\t%.6e')

    if writefits:
        iraf.noao()
        iraf.onedspec()
        iraf.rspectext(outtext, outspec, dtype='nonlinear', flux='yes', title=pf.getheader(merged_spec)['OBJECT'])

    return jspec


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    log.configure('gnirs-pype.log', filelevel='INFO', screenlevel='DEBUG')
    boxit('Some text that needs a star box around it', '*')
    boxit('Some text that needs a line box around it', '-')
    boxit('A long centered multi-line message\nthat needs to wrap.', '#')

