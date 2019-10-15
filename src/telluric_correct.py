#!/usr/bin/env python

from astropy.io import fits
import ConfigParser
import log
import os
from pyraf import iraf
import utils

# TODO: Should this module be split into two (like it is for XDGNIRS):
# part 1 iterates through all the Telluric paths and does steps 1,2,3
# part 2 iterates through the science paths and does steps 4,5


# ----------------------------------------------------------------------------------------------------------------------
def start(configfile):
    """
    Do a Telluric correction using the IRAF TELLURIC task.
    """
    logger = log.getLogger('gnirsTelluric')

    path = os.getcwd()  # Store current working directory for later use

    # Set up/prepare IRAF.
    iraf.gemini()
    iraf.gemtools()
    iraf.gnirs()

    # Reset to default parameters the used IRAF tasks.
    iraf.unlearn(iraf.gemini,iraf.gemtools,iraf.gnirs)

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
    runtimedata = config.get('defaults', 'runtimeData')
    hLineInter = config.getboolean('interactive', 'hLineInter')
    continuumInter = config.getboolean('interactive', 'continuumInter')
    telluricInter = config.getboolean('interactive', 'telluricInter')
    extractFullSlit = config.getboolean('extractSpectra1D', 'extractFullSlit')
    extractStepwise = config.getboolean('extractSpectra1D', 'extractStepwise')
    extractStepSize = config.getfloat('extractSpectra1D', 'extractStepSize')

    combinedsrc = config.get('runtimeFilenames', 'combinedsrc')
    extractRegularPrefix = config.get('runtimeFilenames', 'extractRegularPrefix')
    extractFullSlitPrefix = config.get('runtimeFilenames', 'extractFullSlitPrefix')
    extractStepwiseTracePrefix = config.get('runtimeFilenames', 'extractStepwiseTracePrefix')
    extractStepwisePrefix = config.get('runtimeFilenames', 'extractStepwisePrefix')
    hLinePrefix = config.get('runtimeFilenames', 'hLinePrefix')
    fitTelContinuumPrefix = config.get('runtimeFilenames', 'fitTelContinuumPrefix')
    dividedTelContinuumPrefix = config.get('runtimeFilenames', 'dividedTelContinuumPrefix')
    telluricPrefix = config.get('runtimeFilenames', 'telluricPrefix')

    startstep = config.getint('telluricCorrection', 'Start')
    stopstep = config.getint('telluricCorrection', 'Stop')
    hLineMethod = config.get('telluricCorrection', 'hLineMethod')

    hLineRegions = utils.dictify(config.items('hLineRegions'), fmt='int')
    continuumRegions = utils.dictify(config.items('continuumRegions'), fmt='int')
    telluricRegions = utils.dictify(config.items('telluricRegions'), fmt='int')
    telluricFitOrders = utils.dictify(config.items('TelluricFitOrders'), fmt='int')

    logger.debug('hLineRegions: %s', hLineRegions)
    logger.debug('continuumRegions: %s', continuumRegions)
    logger.debug('telluricRegions: %s', telluricRegions)
    logger.debug('telluricFitOrders: %s', telluricFitOrders)

    for scipath in config.options("ScienceDirectories"):

        if not config.getboolean("ScienceDirectories", scipath):  # only process directories that are marked True
            logger.debug('Skipping %s', scipath)
            continue

        logger.info(' --------------------- ')
        logger.info('| Telluric Correction |')
        logger.info(' --------------------- ')

        scipath += '/Intermediate'
        logger.info("%s", scipath)
        iraf.chdir(scipath)

        telpath = '../Telluric/Intermediate'
        logger.info("Telluric directory: %s", telpath)
        logger.info("Runtime data path: %s", runtimedata)

        sci_spec = scipath + '/' + extractRegularPrefix + combinedsrc
        tel_spec = extractRegularPrefix + combinedsrc
        vega_spec = runtimedata + 'vega_ext.fits'
        utils.requires([sci_spec, tel_spec, vega_spec])

        # TODO: support step-wise extraction

        orders = utils.get_orders(scipath)

        tel_h_corrected = utils.nofits(telpath + '/' + hLinePrefix + os.path.basename(tel_spec))
        telluric_fitContinuum = utils.nofits(telpath + '/' + fitTelContinuumPrefix + os.path.basename(tel_h_corrected))
        telluric_dividedContinuum = utils.nofits(telpath + '/' + dividedTelContinuumPrefix + os.path.basename(tel_h_corrected))
        science_dividedTelluricLines = utils.nofits(telluricPrefix + os.path.basename(sci_spec))
        science_correctedTelluric = utils.nofits(dividedTelContinuumPrefix + os.path.basename(science_dividedTelluricLines))

        for valindex in range(startstep, stopstep + 1):
            logger.debug('valindex = %d', valindex)

            if valindex == 1:
                logger.info(" ------------------------------- ")
                logger.info("| STEP 1: Remove hydrogen lines |")
                logger.info(" ------------------------------- ")
                # Output: H line corrected Telluric 1D spectrum
                utils.pause(manualMode)

                if utils.exists([tel_h_corrected + '_order%d.fits' % o for o in orders], overwrite):
                    logger.info('Telluric standard already had the hydrogen lines removed.  Skipping this step.')
                    continue

                if hLineMethod == 'None':
                    logger.warning('Skipping H line removal')
                    logger.warning('Copying files to have the right names for later use in the pipeline.')
                    for i in range(1, len(orders)+1):
                        iraf.imcopy(input=tel_spec + '[SCI,' + str(i) + ']',
                                    output=tel_h_corrected + '_order%d.fits' % orders[i], verbose='yes')

                elif 'vega' in hLineMethod.lower():
                    hcor_vega(tel_spec, vega_spec, tel_h_corrected, hLineRegions, orders, hLineInter)

                elif 'lorentz' in hLineMethod.lower():
                    hcor_lorentz(tel_spec, tel_h_corrected, orders, runtimedata)

                elif 'manual' in hLineMethod.lower():
                    hcor_manual(tel_spec, tel_h_corrected, orders)

                else:
                    logger.error("Unrecognized H line removal method: %s", hLineMethod)
                    raise SystemExit

                if 'tweak' in hLineMethod.lower():  # Vega-Tweak or Lorentz-Tweak
                    hcor_manual(tel_h_corrected, tel_h_corrected, orders, tweak=True)

            elif valindex == 2:
                logger.info(" ------------------------------------------ ")
                logger.info("| STEP 2: Normalize the Telluric continuum |")
                logger.info(" ------------------------------------------ ")
                # Output: 1D Fit to the telluric continuum
                # Output: Continuum-divided, H line removed telluric 1D source spectra
                utils.pause(manualMode)

                normalize(tel_h_corrected, telluric_fitContinuum, telluric_dividedContinuum,
                          continuumInter, orders, continuumRegions, telluricFitOrders, overwrite)

            elif valindex == 4:
                logger.info(" -------------------------------------------- ")
                logger.info("| STEP 4: Divide the science by the Telluric |")
                logger.info(" -------------------------------------------- ")
                # Output: Telluric line removed science 1D source spectra
                utils.pause(manualMode)

                telluricCorrect(sci_spec, telluric_dividedContinuum, science_dividedTelluricLines, telluricInter,
                                orders, telluricRegions, 'science_telluricInfo.txt', overwrite)

            elif valindex == 5:
                logger.info(" -------------------------------------------- ")
                logger.info("| STEP 5: reintroduceTelluricContinuum       |")
                logger.info(" -------------------------------------------- ")
                # Output: Continuum-divided, telluric corrected science 1D source spectra.
                utils.pause(manualMode)

                reintroduceContinuum(science_dividedTelluricLines, telluric_fitContinuum,
                                     science_correctedTelluric, orders, overwrite)

        logger.info(" --------------------------------------------------- ")
        logger.info("| Telluric correction complete for this observation |")
        logger.info(" --------------------------------------------------- ")

    iraf.chdir(path)  # Return to the original directory

    return


# ----------------------------------------------------------------------------------------------------------------------
def hcor_vega(infile, calfile, outfile, regions, order, interactive):
    # Use IRAF TELLURIC to remove H absorption lines by scaling and shifting a Vega spectrum.

    logger = log.getLogger('hcor_vega')
    logger.debug('--------------------------------------------------')
    logger.debug('infile: %s', infile)
    logger.debug('calfile: %s', calfile)
    logger.debug('outfile: %s', outfile)
    logger.debug('regions: %s', regions)
    logger.debug('interactive: %s', interactive)
    utils.requires([infile, calfile])
    airmass = fits.getheader(infile)['AIRMASS']
    info = open('../Telluric/Intermediate/telluric_hlines.txt', 'w')  # What is this for ?

    for i in range(len(order)):
        extension = i+1
        inspec = '%s[SCI,%d]' % (infile, extension)
        calspec = '%s[%d]' % (calfile, extension)
        outspec = '%s_order%d.fits' % (outfile, order[i])
        sample = regions[order[i]]
        logger.debug('-------------------------')
        logger.debug('inspec: %s', inspec)
        logger.debug('calspec: %s', calspec)
        logger.debug('outspec: %s', outspec)

        iraf.hedit(images=inspec, fields='AIRMASS', value=airmass, add='yes', addonly='no', delete='no',
                   verify='no', show='no', update='yes')

        results = iraf.telluric(
            input=inspec, output=outspec, cal=calspec, sample=sample, ignoreaps='yes', xcorr='yes', tweakrms='yes',
            interactive=interactive, threshold=0.1, lag=3, shift=0., dshift=0.05, scale=1., dscale=0.05, offset=0,
            smooth=1, cursor='', answer='yes', mode='al', Stdout=1)

        logger.debug('iraf.telluric screen output: %s', results)
        info.write(str(results) + '\n')  # Record shift and scale results for future reference

        if 'limits' in results[-1].split()[-1]:      # if there is a warning about pixels ouside limits
            normalization = results[-2].split()[-1]  # the normalization value is second to last
        else:
            normalization = results[-1].split()[-1]  # otherwise the normalization value is last
        logger.debug('Normalization: %s', normalization)

        logger.debug('Undoing the normalization introduced by iraf.telluric...')
        iraf.imarith(
            operand1=outspec, operand2=normalization, op='/', result=outspec, title='',
            divzero=0.0, hparams='', pixtype='', calctype='', verbose='yes', noact='no', mode='al')

        # utils.plot(line1=inspec, line2=outspec, label1='in', label2='out',
        #           title='Telluric before and after H-line removal')

    info.close()

    # Comment from nifsTelluric.py: There are subtle bugs in iraf mean imarith does not work. So we use an
    # astropy/numpy solution.  See NIFTY code for the workaround if necessary.

    return


# ----------------------------------------------------------------------------------------------------------------------
def hcor_lorentz(infile, outfile, order, cursordir):
    # Automatially fit Lorentz profiles to the hydrogen lines pre-defined in lorentz_orderN.cur files.
    logger = log.getLogger('hcor_lorentz')
    logger.debug('--------------------------------------------------')
    utils.requires([infile])

    for i in range(len(order)):
        ext = i+1
        inspec = infile + '[SCI,' + str(ext) + ']'
        outspec = outfile + '_order%d.fits' % order[i]
        cursor1 = cursordir + 'lorentz_order%d.cur' % order[i]
        cursor2 = '../Telluric/Intermediate/lorentz_fitting_input_order%d' % order[i]
        stdout = '../Telluric/Intermediate/lorentz_fitting_output_order%d' % order[i]
        logger.debug('-------------------------')
        logger.debug('inspec: %s', inspec)
        logger.debug('outspec: %s', outspec)
        logger.debug('cursor1: %s', cursor1)

        utils.requires([cursor1])
        utils.exists([stdout, cursor2], overwrite=True)

        # Input the hydrogen line wavelengths and output pixel positions:  (is there another way to do this?)
        specpos = iraf.bplot(images=inspec, cursor=cursor1, Stdout=1, StdoutG='/dev/null')
        specpos = str(specpos).split("'x,y,z(x):")
        logger.debug('specpos: %s', specpos)

        # Write a new cursor file with line x,y info and bplot commands to fit Lorentz profiles:
        write_line_positions(cursor2, specpos)

        # Fit and subtract Lorentz profiles and write output to file for posterity:
        iraf.bplot(images=inspec, cursor=cursor2, new_image=outspec, overwrite="yes",
                   StdoutG='/dev/null', Stdout=stdout)

        #utils.plot(line1=inspec, line2=outspec, label1='in', label2='out',
        #           title='Telluric before and after H-line removal')

    return


def write_line_positions(filename, var):
    with open(filename, "w") as curfile:
        for i in range(len(var)):
            if i != 0:
                var[i] = var[i].split()
                var[i][2] = var[i][2].replace("',", "").replace("']", "")
            if not i % 2 and i != 0:
                # even number, means RHS of H line
                # write x and y position to file, also "k" key
                curfile.write(var[i][0]+" "+var[i][2]+" 1 k \n")
                # LHS of line, write info + "l" key to file
                curfile.write(var[i-1][0]+" "+var[i-1][2]+" 1 l \n")
                # now repeat but writing the "-" key to subtract the fit
                curfile.write(var[i][0]+" "+var[i][2]+" 1 - \n")
                curfile.write(var[i-1][0]+" "+var[i-1][2]+" 1 - \n")
        curfile.write("0 0 1 i\n")
        curfile.write("0 0 1 q\n")
    return


# ----------------------------------------------------------------------------------------------------------------------
def hcor_manual(infile, outfile, order, tweak=False):
    # Run IRAF splot so the user can manually remove hydrogen lines by fitting and subtracting thier favorite profiles.
    logger = log.getLogger('hcor_manual')
    logger.debug('--------------------------------------------------')

    if not tweak:
        utils.requires([infile])

    logger.info('Running IRAF splot so that you can fit and remove hydrogen lines.')
    logger.info('Type "?" to get help.')
    logger.info('Type "i" to save your line-free spectrum')

    for i in range(len(order)):
        ext = i+1
        if tweak:  # the input files will be simple FITS
            inspec = infile + '_order%d.fits' % order[i]
        else:      # the input files will be MEFs
            inspec = infile + '[SCI,' + str(ext) + ']'
        outspec = outfile + '_order%d.fits' % order[i]
        stdout = '../Telluric/Intermediate/manual_fitting_output_order%d' % order[i]
        logger.debug('-------------------------')
        logger.debug('inspec: %s', inspec)
        logger.debug('outspec: %s', outspec)
        utils.exists([stdout], overwrite=True)
        iraf.splot(images=inspec, new_image=outspec, save_file=stdout, overwrite='yes')

        # It's easy to forget to use the 'i' key to save the line-free spectrum, so check that it exists:
        # With the 'tweak' options, the line-free spectrum will already exist, so the user can simply type 'q'
        # and move on without editing.  Too bad if they edit and forget to hit 'i'  :-(

        while True:
            if os.path.exists(outspec):
                break
            else:
                logger.error('Cannot find the output file.')
                logger.error('You may have fogotten to type "i" to save the spectrum')
                logger.error('Re-running splot...')
                iraf.splot(images=inspec, new_image=outspec, save_file=stdout, overwrite='yes')


# ----------------------------------------------------------------------------------------------------------------------
def normalize(inroot, fitroot, outroot, interactive, orders, samples, fitorders, overwrite):
    # Normalize the H line corrected Telluric 1D spectra using IRAF CONTINUUM.
    # inroot = file name root of spectra to be fit
    # fitroot = file name root of continuum fits
    # outroot = file name root of normalized spectra
    logger = log.getLogger('fitContinuum')

    infiles = ['%s_order%d.fits' % (inroot, o) for o in orders]
    fitfiles = ['%s_order%d.fits' % (fitroot, o) for o in orders]
    outfiles = ['%s_order%d.fits' % (outroot, o) for o in orders]

    logger.debug('infiles: %s', infiles)
    logger.debug('fitfiles: %s', fitfiles)
    logger.debug('outfiles: %s', outfiles)
    logger.debug('samples: %s', samples)
    logger.debug('fit orders: %s', fitorders)

    utils.requires(infiles)
    if utils.exists(fitfiles + outfiles, overwrite):
        logger.info('Telluric standard already normalized.')
        return

    for o in orders:
        inspec = inroot + '_order%d.fits' % o
        fitspec = fitroot + '_order%d.fits' % o
        outspec = outroot + '_order%d.fits' % o
        sample = samples[o]
        fitorder = fitorders[o]

        # Here, iraf.continuum is used with type='fit' to write out the continuum fit rather than type='ratio'.

        # NOTE:  As of August 2019, the help documentation for iraf.continuum does not show 'ask' as a
        # non-optional input parameter; however, <lpar continuum> lists ask as the third non-optional paranter
        # to be given with the task. If the input is a list of spectra as against an individual image and 'ask'
        # is set to 'yes', the user will be asked, while proceeding to the next spectrum in the input list,
        # whether they would like to run the task interactively. A "YES" would run the task interactively for
        # all spectra in the input list.

        logger.info('Fitting %s in the region [%s] with an order %s spline...', inspec, sample, fitorder)

        iraf.continuum(
            input=inspec, output=fitspec, ask='yes', lines='*', bands='1', type='fit', replace='no', wavescale='yes',
            logscale='no', override='no', listonly='no', logfiles=logger.root.handlers[0].baseFilename,
            inter=interactive, sample=sample, naverage=1, func='spline3', order=fitorder, low_reject=1.0,
            high_reject=3.0, niterate=2, grow=1.0, markrej='yes', graphics='stdgraph', cursor='', mode='ql')

        #utils.plot(line1=inspec, line2=fitspec, label1='in', label2='fit', title='Telluric continuum and fit')

        iraf.imarith(
            operand1=inspec, operand2=fitspec, op='/', result=outspec, title='', divzero=0.0, hparams='',
            pixtype='', calctype='', verbose='yes', noact='no', mode='al')

        #utils.plot(line1=outspec, title='Normalized Telluric')

    return


# ----------------------------------------------------------------------------------------------------------------------
def telluricCorrect(sci_input, telluric, sci_output, interactive, orders, regions, infofile, overwrite):
    """
    Run IRAF telluric task on each order
    :param sci_input: input multi-extension fits file name
    :param telluric: input calibration spectrum file name root
    :param sci_output: output spectrum file name root
    :param interactive: interactive parameter for iraf.telluric
    :param orders: list of orders
    :param regions: list of regions to fit in each order [(order3,region1), (order4,region2), ...]
    :param infofile: filename for output results that are normally printed to the terminal
    :param overwrite: overwrite previous results
    """
    logger = log.getLogger('telluricCorrect')

    infiles = [sci_input] + ['%s_order%d.fits' % (telluric, o) for o in orders]
    logger.debug('infiles: %s', infiles)

    outfiles = ['%s_order%d.fits' % (sci_output, o) for o in orders]
    logger.debug('outfiles: %s', outfiles)

    utils.requires(infiles)

    if utils.exists(outfiles, overwrite):
        logger.info('Science spectra already Telluric corrected.')
        return

    airmass = fits.getheader(sci_input)['AIRMASS']
    info = open(infofile, 'w')

    for i in range(len(orders)):
        extension = i+1
        inspec = '%s[SCI,%d]' % (sci_input, extension)
        calspec = '%s_order%d.fits' % (telluric, orders[i])
        outspec = '%s_order%d.fits' % (sci_output, orders[i])
        logger.debug('-------------------------')
        logger.debug('inspec: %s', inspec)
        logger.debug('calspec: %s', calspec)
        logger.debug('outspec: %s', outspec)

        results = iraf.telluric(
            input=inspec, output=outspec, cal=calspec, airmass=airmass, answer='yes', ignoreaps='yes', xcorr='yes',
            tweakrms='yes', interactive=interactive, sample=regions[orders[i]], threshold=0.1, lag=10, shift=0.0,
            dshift=1.0, scale=1.0, dscale=0.2, offset=1, smooth=1, cursor='', Stdout=1)

        logger.debug('iraf.telluric screen output: %s', results)
        info.write(str(results) + '\n')  # Record shift and scale results for future reference

        if "limits" in results[-1].split()[-1]:
            index = -2
        else:
            index = -1
        normalization = results[index].split()[-1]
        scale = float(results[index].split()[-4].replace(',', ''))
        shift = float(results[index].split()[-7].replace(',', ''))
        logger.info('Normalization: %s', normalization)

        logger.info('Shift: %.2f pixels', shift)
        if abs(shift) > 1.0:
            logger.warning("Shift shift > 1 pixel!")

        logger.info("Scale: %.3f", scale)
        if abs(scale - 1.0) > 0.1:
            logger.warning("Scaling is > 10%")

        logger.debug('Undoing the normalization introduced by iraf.telluric...')
        iraf.imarith(
            operand1=outspec, operand2=normalization, op='/', result=outspec, title='', divzero=0.0,
            hparams='', pixtype='', calctype='', verbose='yes', noact='no', mode='al')

        #utils.plot(line1=inspec, line2=outspec, label1='in', label2='out',
        # title='Order %d with Telluric Correction' % orders[i])

    info.close()
    return


# ----------------------------------------------------------------------------------------------------------------------
def reintroduceContinuum(inroot, fitroot, outroot, orders, overwrite):
    # Re-introduce telluric continuum shape into Telluric-corrected science spectra.
    logger = log.getLogger('reintroduceContinuum')

    infiles = ['%s_order%d.fits' % (inroot, o) for o in orders]
    fitfiles = ['%s_order%d.fits' % (fitroot, o) for o in orders]
    outfiles = ['%s_order%d.fits' % (outroot, o) for o in orders]
    logger.debug('infiles: %s', infiles)
    logger.debug('fitfiles: %s', fitfiles)
    logger.debug('outfiles: %s', outfiles)

    utils.requires(infiles)
    if utils.exists(outfiles, overwrite):
        logger.info('Step already performed.')
        return

    for o in orders:
        inspec = '%s_order%d.fits' % (inroot, o)
        fitspec = '%s_order%d.fits' % (fitroot, o)
        outspec = '%s_order%d.fits' % (outroot, o)

        logger.info('Re-introducing Telluric continuum back into order %d', o)
        iraf.imarith(
            operand1=inspec, operand2=fitspec, op="/", result=outspec, title='', divzero=0.0,
            hparams='', pixtype='', calctype='', verbose='yes', noact='no', mode='al')

        utils.plot(line1=outspec, title='Order %d with Telluric continuum' % o)

    return


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    log.configure('gnirs-pype.log', filelevel='INFO', screenlevel='DEBUG')
    start('gnirs-pype.cfg')
