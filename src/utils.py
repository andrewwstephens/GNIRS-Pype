#!/usr/bin/env python
import log
import os


# ----------------------------------------------------------------------------------------------------------------------
def boxit(text, character, center=True):
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
def exists(inlist, overwrite):
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

    exist = [os.path.exists(f) for f in inlist]
    logger.debug('Exist: %s', exist)

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
            # We're stuck.  We don't have all the files required to skip this step,
            # but we can't erase files and try again.  Give up and let the user figure it out.
            logger.error('Some (but not all) the files exist and overwrite = False')
            raise SystemExit


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
        raise SystemExit("Unknown GNIRS configuration.  Cannot determine the XD orders.")
    logger.debug('orders: %s', orders)
    return orders


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

