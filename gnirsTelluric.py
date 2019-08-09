#!/usr/bin/env python

import ConfigParser
import log
import os
import glob
from   astropy.io import fits
import matplotlib.pyplot as plt
from   pyraf import iraf


# ----------------------------------------------------------------------------------------------------------------------
def start(configfile):
    """
    Do a Telluric correction using the IRAF TELLURIC task.
    """
    logger = log.getLogger('gnirsTelluric.start')

    path = os.getcwd()  # Store current working directory for later use

    logger.info('#################################################')
    logger.info('#                                               #')
    logger.info('#       Start GNIRS Telluric Correction         #')
    logger.info('#                                               #')
    logger.info('#################################################')

    # Set up/prepare IRAF.
    iraf.gemini()
    iraf.gemtools()
    iraf.gnirs()

    # Reset to default parameters the used IRAF tasks.
    iraf.unlearn(iraf.gemini,iraf.gemtools,iraf.gnirs)

    # Prepare the IRAF package for GNIRS.
    # NSHEADERS lists the header parameters used by the various tasks in the GNIRS package (excluding headers values 
    # which have values fixed by IRAF or FITS conventions).
    iraf.nsheaders("gnirs",logfile=logger.root.handlers[0].baseFilename)

    # Set clobber to 'yes' for the script. This still does not make the gemini tasks overwrite files, so: YOU WILL 
    # LIKELY HAVE TO REMOVE FILES IF YOU RE_RUN THE SCRIPT.
    us_clobber=iraf.envget("clobber")
    iraf.reset(clobber='yes')
    
    config = ConfigParser.RawConfigParser()
    config.optionxform = str  # make options case-sensitive
    config.read(configfile)
    manualMode = config.getboolean('defaults', 'manualMode')
    overwrite = config.getboolean('defaults', 'overwrite')
    hLineInter = config.getboolean('interactive', 'hLineInter')
    continuumInter = config.getboolean('interactive', 'continuumInter')
    telluricInter = config.getboolean('interactive', 'telluricInter')
    tempInter = config.getboolean('interactive', 'tempInter')
    extractionStepwise = config.getboolean('extractSpectra1D', 'extractionStepwise')
    calculateSpectrumSNR = config.getboolean('gnirsPipeline', 'calculateSpectrumSNR')
    # telluricCorrection specific config
    start = config.getint('telluricCorrection', 'Start')
    stop = config.getint('telluricCorrection', 'Stop')
    hLineMethod = config.get('telluricCorrection', 'hLineMethod')


    for scipath in config.options("ScienceDirectories"):
        if config.getboolean("ScienceDirectories",scipath):

            ###########################################################################
            ##                                                                       ##
            ##                  BEGIN - OBSERVATION SPECIFIC SETUP                   ##
            ##                                                                       ##
            ###########################################################################

            #os.chdir(scipath)   # Why?
            iraf.chdir(scipath)  # Change the iraf directory to the current directory.

            # Print the current directory of observations being reduced.
            logger.info("Currently working on telluric correction in %s\n", scipath)
            tempScipath = scipath.split(os.sep)
            # Get symbolic path to the telluric directory within the science directory
            # TODO(Viraja):  Check with Andy if there is a better way of getting the absolute path to the telluric
            # directory
            telpath = (glob.glob("/".join(tempScipath[:-1])+'/Tel_*')).pop()      # This may not work as expected
            logger.info("Telluric directory: %s\n", telpath)
            # Get the runtime data directory path
            runtimedatapath = "/".join(tempScipath[:-5])+'/runtimeData'
            # Print the runtime data directory path
            logger.info("Runtime data path: %s\n", runtimedatapath)

            # Check if required runtime data files and extracted spectra available in their respective paths
            logger.info("Checking if required runtime data files and extracted spectra available in %s", runtimedatapath)
            logger.info("and %s, respectively.\n", scipath)
            
            telluric_hLineRegions_filename = 'telluric_hLineRegions.dat'
            if os.path.exists(runtimedatapath+'/'+telluric_hLineRegions_filename):
                logger.info("Required telluric H line regions reference file available.")
                telluric_hLineRegions = open(runtimedatapath+'/'+telluric_hLineRegions_filename, "r").readlines()
            else:
                logger.warning("Required telluric H line regions reference file not available. Please provide the ")
                logger.warning("file in %s. Exiting script.\n", runtimedatapath)
                raise SystemExit

            telluric_continuumRegions_filename = 'telluric_continuumRegions.dat'
            if os.path.exists(runtimedatapath+'/'+telluric_continuumRegions_filename):
                logger.info("Required telluric continuum regions reference file available.")
                telluric_continuumRegions = open(runtimedatapath+'/'+telluric_continuumRegions_filename, "r").readlines()
            else:
                logger.warning("Required telluric continuum regions reference file not available. Please provide the ")
                logger.warning("file in %s. Exiting script.\n", runtimedatapath)
                raise SystemExit

            if os.path.exists(runtimedatapath + '/' + 'telluric_telluricRegions.dat'):
                logger.info("Required telluric regions reference file available.")
                telluricRegions = {}
                with open(runtimedatapath+'/' + 'telluric_telluricRegions.dat', "r") as f:
                    for line in f.readlines():
                        chunk = line.split()
                        telluricRegions[chunk[0]] = chunk[1]
                logger.debug('TelluricRegions: %s', telluricRegions)
            else:
                logger.warning("Required telluric regions reference file not available. Please provide the file in")
                logger.warning("%s. Exiting script.\n", runtimedatapath)
                raise SystemExit

            vega_spectrum = runtimedatapath+'/vega_ext.fits'
            if os.path.exists(vega_spectrum):
                logger.info("Required vega spectrum available.")
            else:
                logger.warning("Required vega spectrum not available. Please provide the vega spectrum in")
                logger.warning("%s. Exiting script.\n", runtimedatapath)
                raise SystemExit
            
            # Assigning different variable names to "src_comb.fits" and "sky_comb.fits" is not necessary with the
            # current file/directory arrangement for this pipeline (where this works as the telluric and science images
            # are not in the same directory), but it can be useful if the user decides to change the filenames of the
            # combined images.
            scisrccombimage = 'src_comb.fits'  
            scisrcextractedspectrum = scipath+'/v'+scisrccombimage
            sciHeader = fits.open(scisrcextractedspectrum)[0].header
            sciAirmass = sciHeader['AIRMASS']
            if os.path.exists(scisrcextractedspectrum):
                logger.info("Required extracted science source spectrum available.")
            else:
                logger.warning("Required extracted science source spectrum not available. Please run")
                logger.warning("gnirsExtarctSpectra1D.py to create the extracted spectrum or provide it manually in")
                logger.warning("%s. Exiting script.\n", scipath)
                raise SystemExit
            
            telsrccombimage = 'src_comb.fits'
            telsrcextractedspectrum = telpath+'/v'+telsrccombimage
            telHeader = fits.open(telsrcextractedspectrum)[0].header
            telAirmass = telHeader['AIRMASS']
            if os.path.exists(telsrcextractedspectrum):
                logger.info("Required extracted telluric 1D source spectrum available.")
            else:
                logger.warning("Required extracted telluric 1D source spectrum not available. Please run")
                logger.warning("gnirsExtarctSpectra1D.py to create the extracted spectrum or provide it manually in")
                logger.warning("%s. Exiting script.\n", telpath)
                raise SystemExit

            '''
            if extractionStepwise:
            '''

            if calculateSpectrumSNR:
                sciskycombimage = 'sky_comb.fits'
                sciskyextractedspectrum = scipath+'/v'+sciskycombimage
                if os.path.exists(sciskyextractedspectrum):
                    logger.info("Required extracted science sky spectrum available.")
                else:
                    logger.warning("Parameter 'calculateSpectrumSNR' is 'True', but required extracted science sky")
                    logger.warning("spectrum not available. Setting 'calculateSpectrumSNR' for the current set of")
                    logger.warning("observations to 'False'.\n")
                    calculateSpectrumSNR = False
                
                telskycombimage = 'sky_comb.fits'
                telskyextractedspectrum = telpath+'/v'+telskycombimage
                if os.path.exists(telskyextractedspectrum):
                    logger.info("Required extracted telluric sky spectrum available.")
                else:
                    logger.warning("Parameter 'calculateSpectrumSNR' is 'True', but required extracted telluric sky")
                    logger.warning("spectrum not available. Setting 'calculateSpectrumSNR' for the current set of")
                    logger.warning("observations to 'False'.\n")
                    calculateSpectrumSNR = False
            
            logger.info("Required extracted spectra check complete.")

            # Record the right number of order expected according to the GNIRS XD configuration.
            if 'Long' in scipath and 'SXD' in scipath:
                orders = [3, 4, 5]
            elif 'Long' in scipath and 'LXD' in scipath:
                orders = [3, 4, 5, 6, 7, 8]
            elif 'Short' in scipath and 'SXD' in scipath:
                orders = [3, 4, 5, 6, 7, 8]
            else:
                logger.error("#############################################################################")
                logger.error("#############################################################################")
                logger.error("#                                                                           #")
                logger.error("#     ERROR in telluric: unknown GNIRS XD configuration. Exiting script.    #")
                logger.error("#                                                                           #")
                logger.error("#############################################################################")
                logger.error("#############################################################################\n")
                raise SystemExit
            
            # Output filenames with different prefixes added at different stages of script but without the extension 
            # '.fits'
            telluric_hLineCorrectedSpectrum = telpath+'/h'+telsrcextractedspectrum[telsrcextractedspectrum.rfind('/')+1:]
            telluric_hLineCorrectedSpectrum = telluric_hLineCorrectedSpectrum[:telluric_hLineCorrectedSpectrum.rfind('.')]
            
            telluric_fitContinuum = telpath+'/fit'+telluric_hLineCorrectedSpectrum[telluric_hLineCorrectedSpectrum.rfind('/')+1:]

            telluric_dividedContinuum = telpath+'/d'+telluric_hLineCorrectedSpectrum[telluric_hLineCorrectedSpectrum.rfind('/')+1:]

            science_dividedTelluricLines = scipath+'/u'+scisrcextractedspectrum[scisrcextractedspectrum.rfind('/')+1:]
            science_dividedTelluricLines = science_dividedTelluricLines[:science_dividedTelluricLines.rfind('.')]

            science_correctedTelluric = scipath+'/d'+science_dividedTelluricLines[science_dividedTelluricLines.rfind('/')+1:]

            ###########################################################################
            ##                                                                       ##
            ##                 COMPLETE - OBSERVATION SPECIFIC SETUP                 ##
            ##             BEGIN TELLURIC CORRECTION FOR AN OBSERVATION              ##
            ##                                                                       ##
            ###########################################################################

            valindex = start
            while valindex > stop or valindex < 1 or stop > 6:
                logger.warning("#####################################################################")
                logger.warning("#####################################################################")
                logger.warning("#                                                                   #")
                logger.warning("#   WARNING in telluric: invalid start/stop values of telluric      #")
                logger.warning("#                        correction steps.                          #")
                logger.warning("#                                                                   #")
                logger.warning("#####################################################################")
                logger.warning("#####################################################################")

                valindex = int(raw_input("Please enter a valid start value (1 to 5, default 1): "))
                stop = int(raw_input("Please enter a valid stop value (1 to 5, default 5): "))

            while valindex <= stop:

                #############################################################################
                ##  STEP 1: H line removal.                                                ##
                ##  Output: H line corrected telluric 1D source spectra.                   ##
                #############################################################################

                if valindex == 1:
                    
                    telhlineinfofile = telpath+'/telluric_hLineInfo.txt'
                    hLineRemoval(telsrcextractedspectrum, telluric_hLineCorrectedSpectrum, hLineInter, orders,
                        hLineMethod, telluric_hLineRegions, telAirmass, vega_spectrum, tempInter, telhlineinfofile,
                        overwrite)

                    logger.info("##################################################################")
                    logger.info("#                                                                #")
                    logger.info("#       STEP 1: H line removal - COMPLETED                       #")
                    logger.info("#                                                                #")
                    logger.info("##################################################################")

                #############################################################################
                ##  STEP 2: Fit telluric continuum.                                        ##
                ##  Output: 1D Fit to the telluric continuum.                              ##
                #############################################################################

                elif valindex == 2:
                    
                    fitTelluricContinuum(telluric_hLineCorrectedSpectrum, telluric_fitContinuum, continuumInter,
                        orders, telluric_continuumRegions, tempInter, overwrite)
                    
                    logger.info("##################################################################")
                    logger.info("#                                                                #")
                    logger.info("#       STEP 2: Fit telluric continuum - COMPLETED               #")
                    logger.info("#                                                                #")
                    logger.info("##################################################################")

                #############################################################################################
                ##  STEP 3: Division of telluric by the telluric continuum.                                ##
                ##  Output: Continuum-divided, H line removed telluric 1D source spectra.                  ##
                #############################################################################################

                elif valindex == 3:
                    
                    divideTelluricContinuum(telluric_hLineCorrectedSpectrum, telluric_fitContinuum,
                        telluric_dividedContinuum, orders, overwrite)
                    
                    logger.info("##############################################################################")
                    logger.info("#                                                                            #")
                    logger.info("#    STEP 3: Division of telluric by the telluric continuum  - COMPLETED     #")
                    logger.info("#                                                                            #")
                    logger.info("##############################################################################")

                #############################################################################
                ##  STEP 4: Telluric line removal.                                         ##
                ##  Output: Telluric line removed science 1D source spectra.               ##
                #############################################################################

                elif valindex == 4:

                    telluric(infile=scisrcextractedspectrum, calib=telluric_dividedContinuum,
                             outfile=science_dividedTelluricLines, interactive=telluricInter,
                             regions=config.items('TelluricRegions'), orders=orders,
                             airmass=sciAirmass, results_file='science_telluricInfo.txt', overwrite=overwrite)
                    
                    logger.info("##################################################################")
                    logger.info("#                                                                #")
                    logger.info("#       STEP 4: Telluric line removal - COMPLETED                #")
                    logger.info("#                                                                #")
                    logger.info("##################################################################")
                
                #################################################################################
                ##  STEP 5: Division of telluric by the telluric continuum.                    ##
                ##  Output: Continuum-divided, telluric corrected science 1D source spectra.   ##
                #################################################################################
                
                elif valindex == 5:
                    
                    reintroduceTelluricContinuum(science_dividedTelluricLines, science_correctedTelluric, orders,
                        overwrite)
                    
                    logger.info("##################################################################################")
                    logger.info("#                                                                                #")
                    logger.info("#       STEP 5: Division of telluric by the telluric continuum - COMPLETED       #")
                    logger.info("#                                                                                #")
                    logger.info("##################################################################################")

                else:
                    logger.error("###########################################################################")
                    logger.error("###########################################################################")
                    logger.error("#                                                                         #")
                    logger.error("#      ERROR in telluric: %d is not valid. Exiting script.", valindex       )
                    logger.error("#                                                                         #")
                    logger.error("###########################################################################")
                    logger.error("###########################################################################")
                    raise SystemExit
            
                valindex += 1

            logger.info("##############################################################################")
            logger.info("#                                                                            #")
            logger.info("#  COMPLETE - Telluric correction completed for                              #")
            logger.info("#  %s", scipath)
            logger.info("#                                                                            #")
            logger.info("##############################################################################")

    # Return to directory script was begun from.
    os.chdir(path)
    return

##################################################################################################################
#                                                     ROUTINES                                                   #
##################################################################################################################

def hLineRemoval(telsrcextractedspectrum, telluric_hLineCorrectedSpectrum, hLineInter, orders, hLineMethod,
    telluric_hLineRegions, telAirmass, vega_spectrum, tempInter, telhlineinfofile, overwrite):
    """
    Remove hydrogen (H) absorption lines from the extracted 1D telluric source spectrum. Output filename prefix is 'h'.

    Reads:
        Extracted 1D telluric source spectrum.

    Writes:
        H line corrected 1D telluric source spectrum.
    """
    logger = log.getLogger('gnirsTelluric.hLineRemoval')

    for i in range(len(orders)):
        # TODO(Viraja):  Can provide regions to be used by iraf.telluric in the code rather than reading from a file.
        for line in telluric_hLineRegions:
            if 'order'+str(orders[i]) in line:
                telluric_hLineRegions_sample = line.split()[1]
                break
            else:
                pass  
        
        extension = i+1
        telluric_hLineRemovalInput = telsrcextractedspectrum+'[SCI,'+str(extension)+']'
        calibrationInput = vega_spectrum+'['+str(extension)+']'
        telluric_hLineRemovalOutput_SEF = telluric_hLineCorrectedSpectrum+'_order'+str(extension)+'_SEF.fits'
        telluric_hLineRemovalOutput_MEF = telluric_hLineCorrectedSpectrum+'_order'+str(extension)+'_MEF.fits'

        oldfiles = glob.glob(telluric_hLineCorrectedSpectrum+'_order'+str(extension)+'*.fits')
        if not oldfiles:
            iraf.hedit(images=telluric_hLineRemovalInput, fields='AIRMASS', value=telAirmass, add='yes',
                addonly='no', delete='no', verify='no', show='no', update='yes', mode='al')
            
            if hLineMethod == '':
                # Need to copy files so have right names for later use
                logger.warning("Parameter 'hLineMethod' in the configuration file is empty. Skipping H line")
                logger.warning("removal for the telluric 1D source spectrum, extension %d.", extension)
                logger.info("Copying files to have the right names for later use in the pipeline.")
                iraf.imcopy(input=telluric_hLineRemovalInput, output=telluric_hLineRemovalOutput_SEF, verbose='yes')

            elif hLineMethod:
                logger.info("Removing H lines from the telluric 1D source spectrum, extension %d.", extension)
                if hLineMethod == 'vega':
                    vega(telluric_hLineRemovalInput, calibrationInput, telluric_hLineRemovalOutput_SEF, hLineInter,
                        telAirmass, telluric_hLineRegions_sample, telhlineinfofile, overwrite)

            else:
                logger.error("Unrecognized H line removal method encountered in the configuration file. Exiting")
                logger.error("script.")
                raise SystemExit
            
            # TODO(Viraja):  Check while testing this script if this is the right position to use iraf.wmef to add the
            # primary header back to the spectra files. Might change depending on whether the other H line removal
            # methods retain the primary headers (they probably do not).
            iraf.wmef(input=telluric_hLineRemovalOutput_SEF, output=telluric_hLineRemovalOutput_MEF, extnames='',
                phu=telsrcextractedspectrum, verbose='yes', mode='al')
            
            # TODO(Viraja):  Temporarily, pause to display and take a look at the telluric 1D source spectra before and 
            # after H line correction. I am not sure if this can be more useful than just looking at the spectra. If 
            # this has to be done, it would require adding another science extension to the MEF created in the previous 
            # step using iraf.wmef task BEFORE starting the H line correction in this loop.
            if tempInter:
                # Plot the telluric 1D source spectrum with and without H line correction
                telluric_hLineUncorrected = fits.open()  ## not completed
                telluric_hLineCorrected = fits.open()  ## not completed
                plt.title('Telluric 1D Source Spectrum Before and After H Line Removal')
                plt.plot(telluric_hLineUncorrected)
                plt.plot(telluric_hLineCorrected)
                plt.show()
        else:
            if overwrite:
                logger.warning("Removing old %s"+"_order"+str(extension)+"*.fits files.", telluric_hLineCorrectedSpectrum)
                [os.remove(filename) for filename in oldfiles]

                iraf.hedit(images=telluric_hLineRemovalInput, fields='AIRMASS', value=telAirmass, add='yes',
                    addonly='no', delete='no', verify='no', show='no', update='yes', mode='al') 

                if hLineMethod == '':
                    # Need to copy files so have right names for later use
                    logger.warning("Parameter 'hLineMethod' in the configuration file is empty. Skipping H line")
                    logger.warning("removal for the telluric 1D source spectrum, extension %d.", extension)
                    logger.info("Copying files to have the right names for later use in the pipeline.")
                    iraf.imcopy(input=telluric_hLineRemovalInput, output=telluric_hLineRemovalOutput_SEF,
                        verbose='yes')

                elif hLineMethod:
                    logger.info("Removing H lines from the telluric 1D source spectrum, extension %d.", extension)
                
                    if hLineMethod == 'vega':
                        vega(telluric_hLineRemovalInput, calibrationInput, telluric_hLineRemovalOutput_SEF,
                            hLineInter, telAirmass, telluric_hLineRegions_sample, telhlineinfofile, overwrite)
                    
                    # TODO(Viraja):  Update the functions for the rest of the H line removal options (commented below).
                    # NOTE: Untested because interactive scripted iraf tasks are broken... Once ready, add them to the 
                    # part of the script where the output did not already exists and the task had to be run to generate
                    # files for the first time. 
                    '''
                    if hLineMethod == "lineFitAuto" and not no_hLine:
                        lineFitAuto(combined_extracted_1d_spectra, grating)

                    if hLineMethod == "lineFitManual" and not no_hLine:
                        lineFitManual(combined_extracted_1d_spectra+'[sci,1]', grating)

                    if hLineMethod == "vega_tweak" and not no_hLine:
                        # First, do H line removal using the 'vega' method automatically, and then give user a chance  
                        # to interact with the spectrum
                        vega(combined_extracted_1d_spectra, path, hLineInter, telluric_shift_scale_record, overwrite)
                        lineFitManual("final_tel_no_hLines_no_norm", grating)

                    if hLineMethod == "lineFit_tweak" and not no_hLine:
                        # First, do H line removal using the 'lineFitAuto' method automatically, then give user a 
                        # chance to interact with the spectrum
                        lineFitAuto(combined_extracted_1d_spectra,grating)
                        lineFitManual("final_tel_no_hLines_no_norm", grating)
                    '''
                else:
                    logger.error("Unrecognized H line removal method encountered in the configuration file. Exiting")
                    logger.error("script.")
                    raise SystemExit

                # TODO(Viraja):  Check while testing this script if this is the right position to use iraf.wmef to add
                # the primary header back to the spectra files. Might change depending on whether the other H line
                # remova methods retain the primary headers (they probably do not).
                iraf.wmef(input=telluric_hLineRemovalOutput_SEF, output=telluric_hLineRemovalOutput_MEF, extnames='',
                    phu=telsrcextractedspectrum, verbose='yes', mode='al')

                # TODO(Viraja):  Temporarily, pause to display and take a look at the telluric 1D source spectra
                # before and after H line correction. I am not sure if this can be more useful than just looking at the
                # spectra. If this has to be done, it would require adding another science extension to the MEF created
                # in the previous step using iraf.wmef task BEFORE starting the H line correction in this loop.
                if tempInter:
                    # Plot the telluric 1D source spectrum with and without H line correction
                    telluric_hLineUncorrected = fits.open()  ## not completed
                    telluric_hLineCorrected = fits.open()  ## not completed
                    plt.title('Telluric 1D Source Spectrum Before and After H Line Removal')
                    plt.plot(telluric_hLineUncorrected)
                    plt.plot(telluric_hLineCorrected)
                    plt.show()
            else:
                logger.warning("Output exists and -overwrite not set - skipping H line removal for the telluric 1D")
                logger.warning("source spectrum, extension %d.", extension)


# ----------------------------------------------------------------------------------------------------------------------
def vega(inputtelspectrum, inputcalspectrum, outputtelspectrum, hLineInter, tel_airmass, tel_hline_regions,
    tel_hline_infofile, overwrite):
    """
    Use IRAF TELLURIC task to remove H absorption lines from the telluric 1D source spectrum using the 'vega' method, 
    then divide normalization introduced by TELLURIC with IMARITH.
    """
    logger = log.getLogger('gnirsTelluric.vega')

    # For each order, this will be done interactively if parameter 'hLineInter' in the configuration file is 'yes'
    telluric_hLineInfo = iraf.telluric(input=inputtelspectrum, output=outputtelspectrum, cal=inputcalspectrum,
        ignoreaps='yes', xcorr='yes', tweakrms='yes', interactive=hLineInter, sample=tel_hline_regions, threshold=0.1,
        lag=3, shift=0., dshift=0.05, scale=1., dscale=0.05, offset=0, smooth=1, cursor='', airmass=tel_airmass,
        answer='yes', mode='al', Stdout=1)

    # Record shift and scale info for future reference in the script
    # NOTE: If function "hLineremoval" is ran more than once, the function call to "vega" will automatically overwrite 
    # the existing tel_hline_infofile.
    tel_hline_info = open(tel_hline_infofile,'w')
    tel_hline_info.write(str(telluric_hLineInfo)+'\n')
    tel_hline_info.close()
    
    # This loop identifies telluric output containing warning about pixels outside calibration limits (different 
    # formatting)  Viraja:  This comment is from HlinesXD.py in XDGNIRS. I am not sure what the "different formatting"
    # means.
    if 'limits' in telluric_hLineInfo[-1].split()[-1]:
        normalization = telluric_hLineInfo[-2].split()[-1]
    else:
        normalization = telluric_hLineInfo[-1].split()[-1]
    
    # Divide normalization introduced by iraf.telluric
    iraf.imarith(operand1=outputtelspectrum, operand2=normalization, op='/', result=outputtelspectrum, title='',
        divzero=0.0, hparams='', pixtype='', calctype='', verbose='yes', noact='no', mode='al')

    '''
    # Comment from nifsTelluric.py -- There are subtle bugs in iraf mean imarith does not work. So we use an 
    # astropy/numpy solution. Open the image and the scalar we will be dividing it by.
    # Viraja:  I am not sure what "... subtle bugs in iraf mean imarith does not work" means. The comment is not clear.
    # TODO(Viraja):  Incorporate this commented section later if the imarith way from XDGNIRS is found to not work 
    # properly.
    operand1 = fits.open(outputtelspectrum)[0].data
    operand2 = float(normalization)
    # Create a new data array
    if operand2 != 0:
        operand1 = operand1 / operand2
    else:
        operand1 = 1
        # Viraja:  Why is operand1 set to 1? I would imagine it is set to itself as the divisor is 0.
    '''


# ----------------------------------------------------------------------------------------------------------------------
def fitTelluricContinuum(telluric_hLineCorrectedSpectrum, telluric_fitContinuum, continuumInter, orders,
    telluric_continuumRegions, tempInter, overwrite):
    """
    Fit the continua of the H line corrected telluric 1D source spectra using IRAF CONTINUUM task to normalize them.
    """
    logger = log.getLogger('gnirsTelluric.fitTelluricContinuum')

    # Continuum fitting order depends on the spectral order.
    for i in range(len(orders)):
        for line in telluric_continuumRegions:
            if 'order'+str(orders[i]) in line:
                telluric_continuumRegions_sample = line.split()[1]
                break
            else:
                pass
        
        extension = i+1
        telluric_fitContinuumInput = telluric_hLineCorrectedSpectrum+'_order'+str(extension)+'_MEF.fits[1]'
        telluric_fitContinuumOutput = telluric_fitContinuum+'_order'+str(extension)+'.fits'

        if orders[i] == 3:
            fitcontinuumorder = 5
        elif orders[i] == 4:
            fitcontinuumorder = 2
        elif orders[i] == 5:
            fitcontinuumorder = 3
        elif orders[i] == 6:
            fitcontinuumorder = 5
        elif orders[i] == 7:
            fitcontinuumorder = 5
        elif orders[i] == 8:
            fitcontinuumorder = 5
        else:
            logger.error("Unrecognized spectral order for GNIRS encountered during telluric continuum fitting.")
            logger.error("Exiting script.")
            raise SystemExit
    
        if os.path.exists(telluric_fitContinuumOutput):
            if overwrite:
                logger.warning("Removing old %s", telluric_fitContinuumOutput)
                os.remove(telluric_fitContinuumOutput)

                # Here, iraf.continuum is used with type='fit' to write out the continuum fit rather than type='ratio'.
                logger.info("Fitting telluric continuum for the H line corrected spectrum, extension %d.", extension)
                
                # NOTE:  As of August 2019, the help documentation for iraf.continuum does not show 'ask' as a 
                # non-optional input parameter; however, <lpar continuum> lists ask as the third non-optional paranter
                # to be given with the task. If the input is a list of spectra as against an individual image and 'ask' 
                # is set to 'yes', the user will be asked, while proceeding to the next spectrum in the input list, 
                # whether they would like to run the task interactively. A "YES" would run the task interactively for
                # all spectra in the input list.

                # NOTE: In TelluricXD.py in XDGNIRS, parameter 'logfiles', where to write the power series coefficients, 
                # was NULL (""). In nifsTelluric.py in NIFTY, parameter 'logfiles' was set to the main logging file for 
                # the data reduction. 
            
                iraf.continuum(input=telluric_fitContinuumInput, output=telluric_fitContinuumOutput, ask='yes',
                    lines='*', bands='1', type="fit", replace='no', wavescale='yes', logscale='no', override='no',
                    listonly='no', logfiles=logger.root.handlers[0].baseFilename, inter=continuumInter,
                    sample=telluric_continuumRegions_sample, naverage=1, func='spline3', order=fitcontinuumorder,
                    low_reject=1.0, high_reject=3.0, niterate=2, grow=1.0, markrej='yes', graphics='stdgraph',
                    cursor='', mode='ql')
                
                if tempInter:
                    # Plot H line corrected telluric 1D source spectrum without continuum fitting and the continuum fit
                    telluric_withContinuumFitting = astropy.io.fits.open()
                    telluric_ContinuumFit = astropy.io.fits.open()
                    plt.title('H Line Corrected Telluric 1D Source Spectrum and Continuum Fit for Normalization')
                    plt.plot(telluric_withContinuumFitting)
                    plt.plot(telluric_ContinuumFit)
                    plt.show()
            else:
                logger.warning("Output exists and -overwrite not set - skipping telluric continuum fitting for the H")
                logger.warning("line corrected spectrum, extension %d.", extension)
        else:
            logger.info("Fitting telluric continuum for the H line corrected spectrum, extension %d.", extension)
            iraf.continuum(input=telluric_fitContinuumInput, output=telluric_fitContinuumOutput, ask='yes', lines='*',
                bands='1', type="fit", replace='no', wavescale='yes', logscale='no', override='no', listonly='no',
                logfiles=logger.root.handlers[0].baseFilename, inter=continuumInter,
                sample=telluric_continuumRegions_sample, naverage=1, func='spline3', order=fitcontinuumorder,
                low_reject=1.0, high_reject=3.0, niterate=2, grow=1.0, markrej='yes', graphics='stdgraph', cursor='',
                mode='ql')
            
            if tempInter:
                # Plot H line corrected telluric 1D source spectrum without continuum fitting and the continuum fit
                telluric_withContinuumFitting = fits.open()  ## not completed
                telluric_ContinuumFit = fits.open()  ## not completed
                plt.title('H Line Corrected Telluric 1D Source Spectrum and Continuum Fit for Normalization')
                plt.plot(telluric_withContinuumFitting)
                plt.plot(telluric_ContinuumFit)
                plt.show()


# ----------------------------------------------------------------------------------------------------------------------
def divideTelluricContinuum(telluric_hLineCorrectedSpectrum, telluric_fitContinuum, telluric_dividedContinuum, orders,
    overwrite):
    """
    Divide the H line corrected telluric 1D source spectra by its continuum fit to normalize it.
    """
    logger = log.getLogger('gnirsTelluric.divideTelluricContinuum')

    for i in range(len(orders)):
        extension = i+1
        telluric_divideContinuumInput = telluric_hLineCorrectedSpectrum+'_order'+str(extension)+'_MEF.fits[1]'
        telluricContinuumFit = telluric_fitContinuum+'_order'+str(extension)+'.fits'
        telluric_divideContinuumOutput_SEF = telluric_dividedContinuum+'_order'+str(extension)+'_SEF.fits'
        telluric_divideContinuumOutput_MEF = telluric_dividedContinuum+'_order'+str(extension)+'_MEF.fits'

        oldfiles = glob.glob(telluric_dividedContinuum+'_order'+str(extension)+'*.fits')
        if not oldfiles:
            logger.info("Dividing the H line corrected telluric 1D source spectrum, extension %d, by the", extension)
            logger.info("telluric continuum fit.")
            
            iraf.imarith(operand1=telluric_divideContinuumInput, operand2=telluricContinuumFit, op='/',
                result=telluric_divideContinuumOutput_SEF, title='', divzero=0.0, hparams='', pixtype='', calctype='',
                verbose='yes', noact='no', mode='al')

            '''
            # Comment from nifsTelluric.py -- There are subtle bugs in iraf mean imarith does not work. So we use an
            # astropy/numpy solution. Open the image and the scalar we will be dividing it by.
            # TODO(Viraja):  Incorporate this commented section later if the imarith way from XDGNIRS is found to not
            # work properly.
            operand1 = fits.open(telluric_divideContinuumInput)[1].data
            operand2 = fits.open(telluricContinuumFit)[0].data
            # Create a new data array
            multiplied = np.array(operand1, copy=True)
            for i in range(len(telluric_divideContinuumOutput_SEF)):
                if operand2[1] != 0:
                    multiplied[1] = operand1[1] / operand2[1]
                else:
                    multiplied[i] = 0.0
                    # Viraja:  Why is operand3 set to 0.0?
            table = Table(multiplied)
            table.write(telluric_divideContinuumOutput_SEF, format='fits')
            '''
            # TODO(Viraja):  Check if a new MEF with the same filename as the previous one-extension image can be 
            # created sucessfully.
            iraf.wmef(input=telluric_divideContinuumOutput_SEF, output=telluric_divideContinuumOutput_MEF,
                extnames='', phu=telluric_hLineCorrectedSpectrum+'_order'+str(extension)+'_MEF', verbose='yes',
                mode='al')
        else:
            if overwrite:
                logger.warning("Removing old %s"+"_order"+str(extension)+"*.fits,", telluric_dividedContinuum)
                [os.remove(filename) for filename in oldfiles]
                
                # Divide the H line corrected telluric 1D source spectra by the continuum fit.
                logger.info("Dividing the H line corrected telluric 1D source spectrum, extension %d, by", extension)
                logger.info("the telluric continuum fit.")
    
                iraf.imarith(operand1=telluric_divideContinuumInput, operand2=telluricContinuumFit, op='/',
                    result=telluric_divideContinuumOutput_SEF, title='', divzero=0.0, hparams='', pixtype='',
                    calctype='', verbose='yes', noact='no', mode='al')
                '''
                # Comment from nifsTelluric.py -- There are subtle bugs in iraf mean imarith does not work. So we use 
                # an astropy/numpy solution. Open the image and the scalar we will be dividing it by.
                # TODO(Viraja):  Incorporate this commented section later if the imarith way from XDGNIRS is found to 
                # not work properly.
                operand1 = fits.open(telluric_divideContinuumInput)[1].data
                operand2 = fits.open(telluricContinuumFit)[0].data
                # Create a new data array
                telluric_divideContinuumOutput_SEF = np.array(operand1, copy=True)
                for i in range(len(telluric_divideContinuumOutput_SEF)):
                    if operand2[1] != 0:
                        telluric_divideContinuumOutput_SEF[1] = operand1[1] / operand2[1]
                    else:
                        telluric_divideContinuumOutput_SEF[i] = 0.0
                        # Viraja:  Why is operand3 set to 0.0?
                '''
                # TODO(Viraja):  Check if a new MEF with the same filename as the previous one-extension image can be 
                # created sucessfully.
                iraf.wmef(input=telluric_divideContinuumOutput_SEF, output=telluric_divideContinuumOutput_MEF,
                    extnames='', phu=telluric_hLineCorrectedSpectrum+'_order'+str(extension)+'_MEF', verbose='yes',
                    mode='al')
            else:
                logger.warning("Output exists and -overwrite not set - skipping division of the H line corrected")
                logger.warning("telluric 1D source spectrum, extension %d, by the telluric continuum fit.", extension)  


# ----------------------------------------------------------------------------------------------------------------------
def telluric(infile, calib, outfile, interactive, orders, airmass, regions, results_file, overwrite):
    """
    Run IRAF telluric task on each order
    :param infile: input multi-extension fits file name
    :param calib:
    :param outfile:
    :param interactive:
    :param orders: list of orders
    :param airmass: float airmass
    :param regions: list of regions to fit in each order [(order1,region1), (order2,region2), ...]
    :param results_file: filename for output results that are normally printed to the terminal
    :param overwrite:
    """
    logger = log.getLogger('gnirsTelluric.telluric')

    sample = {}
    for order, samp in regions:
        sample[int(order)] = samp
    logger.debug('Sample: %s', sample)

    for i in range(len(orders)):

        extension = i+1
        output = outfile + '_order' + str(extension)+'_SEF.fits'
        if os.path.exists(output):
            if overwrite:
                logger.debug('Removing %s', output)
                os.remove(output)
            else:
                logger.warning("Output exists and overwrite flag not set.")
                logger.warning("Skipping division of the science 1D source")
                logger.warning("spectrum, extension %d, by the H line corrected, continuum-divided", extension)
                logger.warning("telluric spectrum, extension %d.", extension)
                continue

        logger.debug('\n')
        logger.info("Running IRAF telluric...")
        logger.debug('Order: %s', orders[i])
        logger.debug("Input: %s", infile + '[SCI,' + str(extension) + ']')
        logger.debug("Calib: %s", calib + '_order' + str(extension) + '_MEF.fits[1]')
        logger.debug("Output: %s", output)
        logger.debug("Sample: %s", sample[orders[i]])
        logger.debug("Airmass: %s", airmass)

        results = iraf.telluric(
            input=infile + '[SCI,' + str(extension) + ']',
            output=output,
            cal=calib + '_order' + str(extension) + '_MEF.fits[1]',
            airmass=airmass,
            answer='yes',
            ignoreaps='yes',
            xcorr='yes',
            tweakrms='yes',
            interactive=interactive,
            sample=sample[orders[i]],
            threshold=0.1,
            lag=10,
            shift=0.0, dshift=1.0,
            scale=1.0, dscale=0.2,
            offset=1, smooth=1, cursor='', Stdout=1)
        logger.debug('Results: %s', results)

        with open(results_file,'w') as f:
            f.write(str(results)+'\n')

        if "limits" in results[-1].split()[-1]:
            index = -2
        else:
            index = -1
        normalization = results[index].split()[-1]
        scale = float(results[index].split()[-4].replace(',', ''))
        shift = float(results[index].split()[-7].replace(',', ''))
        logger.info('Shift: %.2f pixels', shift)
        if abs(float(shift) > 1.0):
            logger.warning("Shift shift > 1 pixel!")

        logger.info("Scale: %.3f", scale)
        if abs(scale - 1.0) > 0.1:
            logger.warning("Scaling is > 10%")

        # Undo normalization introduced by iraf.telluric:
        # iraf.imarith(operand1=science_divideTelluricLinesOutput_SEF, operand2=normalization, op='/',
        #    result=science_divideTelluricLinesOutput_SEF, title='', divzero=0.0, hparams='', pixtype='',
        #    calctype='', verbose='yes', noact='no', mode='al')

        # TODO(Viraja):  Check if a new MEF with the same filename as the previous one-extension image can be
        # created sucessfully.
        # iraf.wmef(input=science_divideTelluricLinesOutput_SEF, output=science_divideTelluricLinesOutput_MEF,
        #    extnames='', phu=input, verbose='yes', mode='al')

    return


# ----------------------------------------------------------------------------------------------------------------------
def reintroduceTelluricContinuum(science_dividedTelluricLines, science_correctedTelluric, orders, overwrite):
    """
    Re-introduce telluric continuum into telluric lines removed science 1D source spectra.
    """
    logger = log.getLogger('gnirsTelluric.reintroduceTelluricContinuum')

    # Start with regularly extracted science spectra (not using full slit or stepwise extraction methods)
    for i in range(len(orders)):
        extension = i+1
        science_correctedTelluricInput = scisrcextractedspectrum[scisrcextractedspectrum.rfind('.')+1:]+'[SCI,'+\
            str(extension)+']'
        telluric_fitContinuumOutput = telluric_fitContinuumOutput+'_order'+str(extension)+'.fits'
        science_correctedTelluricOutput = science_correctedTelluric+'_order'+str(extension)+'.fits'
        
        if os.path.exists(science_correctedTelluricOutput):
            if overwrite:
                logger.warning("Removing old %s", science_correctedTelluricOutput)
                os.remove(science_correctedTelluricOutput)
                
                logger.info("Re-introducing telluric continuum into th telluric lines removed science 1D source")
                logger.info("spectrum, extension %d.", extension)
                iraf.imarith(operand1=science_correctedTelluricInput, operand2=telluric_fitContinuumOutput, op="/",
                    result=science_correctedTelluricOutput, title='', divzero=0.0, hparams='', pixtype='',
                    calctype='', verbose='yes', noact='no', mode='al')
                '''
                # Comment from nifsTelluric.py -- There are subtle bugs in iraf mean imarith does not work. So we use  
                # an astropy/numpy solution. Open the image and the scalar we will be dividing it by.
                # TODO(Viraja):  Incorporate this commented section later if the imarith way from XDGNIRS is found to 
                # not work properly.
                operand1 = fits.open(science_correctedTelluricInput)[0].data
                operand2 = float(normalization)
                # Create a new data array
                if operand2 != 0:
                    operand3 = operand1 / operand2
                else:
                    operand3 = 1
                    # Viraja:  I am not sure what operand3 would be in the above line??
                '''
                # TODO(Viraja):  Check if a new MEF with the same filename as the previous one-extension image can be 
                # created sucessfully.
                iraf.wmef(input=science_correctedTelluricOutput, output=science_correctedTelluricOutput, extnames='',
                    phu=science_correctedTelluricInput, verbose='yes', mode='al')
            else:
                logger.warning("Output exists and -overwrite not set - skipping re-introduction of telluric continuum")
                logger.warning("into the telluric lines removed science 1D source spectrum, extension %d.", extension)
        else:
            logger.info("Re-introducing telluric continuum into th telluric lines removed science 1D source")
            logger.info("spectrum, extension %d.", extension)
            iraf.imarith(operand1=science_correctedTelluricInput, operand2=telluric_fitContinuumOutput, op="/",
                result=science_correctedTelluricOutput, title='', divzero=0.0, hparams='', pixtype='',
                calctype='', verbose='yes', noact='no', mode='al')
            '''
            # Comment from nifsTelluric.py -- There are subtle bugs in iraf mean imarith does not work. So we use  
            # an astropy/numpy solution. Open the image and the scalar we will be dividing it by.
            # TODO(Viraja):  Incorporate this commented section later if the imarith way from XDGNIRS is found to 
            # not work properly.
            operand1 = fits.open(science_correctedTelluricInput)[0].data
            operand2 = float(normalization)
            # Create a new data array
            if operand2 != 0:
                operand3 = operand1 / operand2
            else:
                operand3 = 1
                # Viraja:  I am not sure what operand3 would be in the above line??
            '''
            # TODO(Viraja):  Check if a new MEF with the same filename as the previous one-extension image can be 
            # created sucessfully.
            iraf.wmef(input=science_correctedTelluricOutput, output=science_correctedTelluricOutput, extnames='',
                phu=science_correctedTelluricInput, verbose='yes', mode='al')
'''
#---------------------------------------------------------------------------------------------------------------------#

def getShiftScale(rawFrame, telluricInter, log, over):
    """
    Use iraf.telluric() to get the best shift and scale of a telluric correction spectrum.

    Use IRAF TELLURIC to get the best shift and scale of the H line corrected, continuum-divided telluric 1D source 
    spectrum.

    Writes:
        "6_shiftScale"+rawFrame+".txt" :
    """
    if os.path.exists('5_oneDCorrected'+rawFrame+'.fits') and os.path.exists("6_shiftScale"+rawFrame+".txt"):
        if over:
            os.remove('5_oneDCorrected'+rawFrame+'.fits')
            # TODO(nat): implement logging for this
            iraf.chdir(os.getcwd())
            tell_info = iraf.telluric(input='4_cubeslice'+rawFrame+'.fits[0]',output='5_oneDCorrected'+rawFrame+'.fits',cal="3_chtel"+rawFrame+'.fits[0]',airmass=1.0,answer='yes',ignoreaps='yes',xcorr='yes',tweakrms='yes',inter=telluricInter,sample="*",threshold=0.1,lag=3,shift=0.,dshift=0.1,scale=1.0,dscale=0.1, offset=1,smooth=1,cursor='',mode='al',Stdout=1)
        else:
            logging.info("\nOutput exists and -over not set - skipping get shift scale of telluric correction and fit")
            return
    else:
        iraf.chdir(os.getcwd())
        tell_info = iraf.telluric(input='4_cubeslice'+rawFrame+'.fits[0]',output='5_oneDCorrected'+rawFrame+'.fits',cal="3_chtel"+rawFrame+'.fits[0]',airmass=1.0,answer='yes',ignoreaps='yes',xcorr='yes',tweakrms='yes',inter=telluricInter,sample="*",threshold=0.1,lag=3,shift=0.,dshift=0.1,scale=1.0,dscale=0.1, offset=1,smooth=1,cursor='',mode='al',Stdout=1)
    # Get shift and scale from the list of values iraf.telluric() returns.
    # Sample tell_info:
    # ['cubeslice.fits[0]: norm.fits[1]: cubeslice.fits[0]: dshift 5.', 'window:again:window:window:again:window:window:again:window:TELLURIC:',
    # '  Output: vtella - HE1353-1917', '  Input: cubeslice.fits[0] - HE1353-1917', '
    # Calibration: norm.fits[1] - Hip70765', '  Tweak: shift = 59.12, scale = 1.323,
    # normalization = 0.9041', '  WARNING: 3 pixels outside of calibration limits']
    tellshift = 0.
    scale = 1.0
    for i in range(len(tell_info)):
        # Now string looks like '  Tweak: shift = 59.12, scale = 1.323, normalization = 0.9041'
        if "Tweak" in tell_info[i]:
            # Remove the first 9 characters,
            temp = tell_info[i][9:]
            # Split into a list; now it looks like '['shift', '=', '59.12,', 'scale', '=', '1.323,', 'normalization', '=', '0.9041']'
            temp = temp.split()
            # Index two is the shift value with a trailing comma, index 5 is the scale value with a trailing comma.
            # Remove trailing comma.
            tellshift = temp[2].replace(',', '')
            # Turn it into a float.
            tellshift = float(tellshift) # Convert to a clean float
            # Do the same for the scale.
            scale = temp[5].replace(',', '')
            scale = float(scale)
    with open("6_shiftScale"+rawFrame+".txt", "w") as text_file:
        text_file.write("Shift: {} Scale: {} \n".format(tellshift, scale))

#---------------------------------------------------------------------------------------------------------------------#

def shiftScaleSpec(rawFrame, inPrefix, outPrefix, log, over):
    """
    Shifts and scales a spectrum using scipy.
    Replaces overflow with 1.
    """
    spectrum = astropy.io.fits.open(inPrefix+rawFrame+'.fits')
    spectrumData = spectrum[0].data
    try:
        with open("6_shiftScale"+rawFrame+".txt", "r") as f:
            line = f.readlines()
        line = line[0].strip().split()
        tellshift = line[1]
        scale = line[3]
        logger.info("\nRead a shift for "+inPrefix+" spectrum for " + str(tellshift))
        logger.info("\nRead a scale of "+inPrefix+" spectrum for  " + str(scale))
    except IOError:
        logger.info("\nNo shiftScale file found for " + rawFrame + " in " + str(os.getcwd() + ". Skipping."))
        return

    if os.path.exists(outPrefix+rawFrame+'.fits'):
        if over:
            os.remove(outPrefix+rawFrame+'.fits')
            # Shift using SciPy, substituting 1 where data overflows.
            # TODO(nat): doesn't look like interpolation is happening but could be tested more.
            # Works but it's gross. The int(round(float())) is a funny way to turn "-0.02" into 0
            spectrumData = scipy.ndimage.interpolation.shift(spectrumData, -1*int(round(float(tellshift))), cval=1.)
            # Scale by simple multiplication; 1D spectrum times a scalar.
            spectrumData = spectrumData * float(scale)
            spectrum[0].data = spectrumData
            spectrum.writeto(outPrefix+rawFrame+'.fits')
        else:
            logger.info("\nOutput exists and -over not set - skipping shift and scale of " + inPrefix)
    else:
        # Shift using SciPy, substituting 1 where data overflows.
        # TODO(nat): doesn't look like interpolation is happening but could be tested more.
        # Works but it's gross. The int(round(float())) is a funny way to turn "-0.02" into 0
        spectrumData = scipy.ndimage.interpolation.shift(spectrumData, -1*int(round(float(tellshift))), cval=1.)
        # Scale by simple multiplication; 1D spectrum times a scalar.
        spectrumData = spectrumData * float(scale)
        spectrum[0].data = spectrumData
        spectrum.writeto(outPrefix+rawFrame+'.fits')

#---------------------------------------------------------------------------------------------------------------------#

# TODO(nat): linefitAuto and linefitManual could be useful at some point.
def lineFitAuto(spectrum, grating):
    """
    Automatically fit Lorentz profiles to lines defined in existing cur* files. Go to x position in cursor file and use 
    space bar to find spectrum at each of those points.
    """
    logger = log.getLogger('gnirsTelluric.lineFitAuto')

    specpos = iraf.bplot(images=spectrum+'[SCI,1]', cursor='cur'+grating, Stdout=1, StdoutG='/dev/null')
    specpose = str(specpos).split("'x,y,z(x):")
    nextcur = 'nextcur'+grating+'.txt'
    # Write line x,y info to file containing Lorentz fitting commands for bplot
    write_line_positions(nextcur, specpos)
    iraf.delete('final_tel_no_hLines_no_norm.fits',ver="no",go_ahead='yes',Stderr='/dev/null')
    # Fit and subtract Lorentz profiles. Might as well write output to file.
    iraf.bplot(images=spectrum+'[sci,1]',cursor='nextcur'+grating+'.txt', new_image='final_tel_no_hLines_no_norm', overwrite="yes",StdoutG='/dev/null',Stdout='Lorentz'+grating)

#---------------------------------------------------------------------------------------------------------------------#

def lineFitManual(spectrum, grating):
    """ 
    Enter splot so the user can fit and subtract lorentz (or, rather any) profiles.
    """
    logger = log.getLogger('gnirsTelluric.lineFitManual')

    iraf.splot(images=spectrum, new_image='final_tel_no_hLines_no_norm', save_file='../PRODUCTS/lorentz_hLines.txt', overwrite='yes')
    # it's easy to forget to use the 'i' key to actually write out the line-free spectrum, so check that it exists:
    # with the 'tweak' options, the line-free spectrum will already exists, so this lets the user simply 'q' and move on w/o editing (too bad if they edit and forget to hit 'i'...)
    while True:
        try:
            with open("final_tel_no_hLines_no_norm.fits") as f: pass
            break
        except IOError as e:
            logger.info("It looks as if you didn't use the i key to write out the lineless spectrum. We'll have to try again. --> Re-entering splot")
            iraf.splot(images=spectrum, new_image='final_tel_no_hLines_no_norm', save_file='../PRODUCTS/lorentz_hLines.txt', overwrite='yes')
'''

# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    log.configure('gnirs.log', filelevel='INFO', screenlevel='DEBUG')

    config = ConfigParser.RawConfigParser()
    config.optionxform = str  # make config file options case-sensitive
    config.read('gnirs.cfg')

    start('gnirs.cfg')

    #telluric(infile='/tmp/vsrc_comb', calib='/tmp/dhvsrc_comb', outfile='testoutput.fits', interactive='no',
    #         regions=config.items('TelluricRegions'), orders=[3, 4, 5, 6, 7, 8], airmass=1.25,
    #         results_file='results.txt', overwrite=True)
