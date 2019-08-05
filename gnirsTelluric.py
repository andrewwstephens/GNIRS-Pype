# MIT License

# Copyright (c) 2015, 2017 Marie Lemoine-Busserolle

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

################################################################################
#                Import some useful Python utilities/modules                   #
################################################################################

import log, glob, shutil, os, glob, re, traceback, ConfigParser, scipy.ndimage.interpolation
import astropy.io import fits
import astropy.coordinates as coord
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astroquery.simbad import Simbad
from pyraf import iraf, iraffunctions


def start():
    """
    Do a telluric correction using the IRAF TELLURIC task.
    """
    logger = log.getLogger('gnirsTelluric.start')

    # Store current working directory for later use.
    path = os.getcwd()

    logging.info('#################################################')
    logging.info('#                                               #')
    logging.info('#       Start GNIRS Telluric Correction         #')
    logging.info('#                                               #')
    logging.info('#################################################\n')

    # Set up/prepare IRAF.
    iraf.gemini()
    iraf.gemtools()
    iraf.gnirs()

    # Reset to default parameters the used IRAF tasks.
    iraf.unlearn(iraf.gemini,iraf.gemtools,iraf.gnirs,)

    # From http://bishop.astro.pomona.edu/Penprase/webdocuments/iraf/beg/beg-image.html:
    # Before doing anything involving image display the environment variable stdimage must be set to the correct frame 
    # buffer size for the display servers (as described in the dev$graphcap file under the section "STDIMAGE devices") 
    # or to the correct image display device. The task GDEVICES is helpful for determining this information for the 
    # display servers.
    iraf.set(stdimage='imt1024')

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
    # Read general config.
    manualMode = config.getboolean('defaults','manualMode')
    overwrite = config.getboolean('defaults','overwrite')
    # config required for extracting 1D spectra
#    observationSections = ['ScienceDirectories','TelluricDirectories']
    hLineInter = config.getboolean('interactive','hLineInter')
    continuumInter = config.getboolean('interactive','continuumInter')
    telluricInter = config.getboolean('interactive','telluricInter')
    tempInter = config.getboolean('interactive','tempInter')
    calculateSpectrumSNR = config.getboolean('gnirsPipeline','calculateSpectrumSNR')
    # telluricCorrection specific config
    start = config.getint('telluricCorrection','Start')
    stop = config.getint('telluricCorrection','Stop')
    hLineMethod = config.get('telluricCorrection','hLineMethod')
    telluricRA = config.get('telluricCorrection','telluricRA')
    telluricDEC = config.get('telluricCorrection','telluricDEC')
    telluricSpectralType = config.get('telluricCorrection','telluricSpectralType')
    telluricMagnitude = config.get('telluricCorrection','telluricMagnitude')
    telluricTemperature = config.get('telluricCorrection','telluricTemperature')
#    telluricBand = config.get('telluricCorrection','telluricBand')

    for scipath in config.options("ScienceDirectories"):
        if config.getboolean("ScienceDirectories",scipath):

            ###########################################################################
            ##                                                                       ##
            ##                  BEGIN - OBSERVATION SPECIFIC SETUP                   ##
            ##                                                                       ##
            ###########################################################################

            os.chdir(scipath)
            # Change the iraf directory to the current directory.
            iraffunctions.chdir(scipath)

            # Print the current directory of observations being reduced.
            logger.info("Currently working on telluric correction in %s\n", scipath)
            # Get symbolic path to the telluric directory within the science directory
            # TODO(Viraja):  Check with Andy if there is a better way of getting the absolute path to the telluric
            # directory
            telpath = (glob.glob(obspath+'/Tel_*')).pop()
            # Print the telluric directory path
            logger.info("The symbolic path to the telluric directory is %s\n", telpath)
            
            reference_starstable_filename = 'stars_spectraltypes_temperature.txt'
            reference_starstable = open(reference_starstable_filename, "r").readlines()
            reference_starstable.close()

            # Check if required combined spectra available in the observations directory path
            logger.info("Checking if required extracted spectra available in %s", scipath)

            scisrccombimage = scipath+'/src_comb.fits'
            scisrcextractedspectrum = scipath+'/v'+scisrccombimage[scisrccombimage.rfind('/')+1:]
            if os.path.exists(scisrcextractedspectrum):
                logger.info("Required extracted science source spectrum available.")
            else:
                logger.warning("Required extracted science source spectrum not available. Please run ")
                logger.warning("gnirsExtarctSpectra1D.py to create the extracted spectrum or provide it manually.")
                logger.warning("Exiting script.\n")
                raise SystemExit
            
            telsrccombimage = telpath+'/src_comb.fits'
            telsrcextractedspectrum = telpath+'/v'+telsrccombimage[telsrccombimage.rfind('/')+1:]
            if os.path.exists(telsrcextractedspectrum):
                logger.info("Required extracted telluric source spectrum available.")
            else:
                logger.warning("Required extracted telluric source spectrum not available. Please run ")
                logger.warning("gnirsExtarctSpectra1D.py to create the extracted spectrum or provide it manually.")
                logger.warning("Exiting script.\n")
                raise SystemExit
            
            if calculateSpectrumSNR:
                sciskycombimage = scipath+'/sky_comb.fits'
                sciskyextractedspectrum = scipath+'/v'+sciskycombimage[sciskycombimage.rfind('/')+1:]
                if os.path.exists(sciskyextractedspectrum):
                    logger.info("Required extracted science sky spectrum available.")
                else:
                    logger.warning("Parameter 'calculateSpectrumSNR' is 'True', but required extracted science sky ")
                    logger.warning("spectrum not available. Setting 'calculateSpectrumSNR' for the current set of")
                    logger.warning("observations to 'False'.\n")
                    calculateSpectrumSNR = False
                
                telskycombimage = telpath+'/sky_comb.fits'
                telskyextractedspectrum = telpath+'/v'+telskycombimage[telskycombimage.rfind('/')+1:]
                if os.path.exists(telskyextractedspectrum):
                    logger.info("Required extracted telluric sky spectrum available.")
                else:
                    logger.warning("Parameter 'calculateSpectrumSNR' is 'True', but required extracted telluric sky ")
                    logger.warning("spectrum not available. Setting 'calculateSpectrumSNR' for the current set of ")
                    logger.warning("observations to 'False'.\n")
                    calculateSpectrumSNR = False
            
            logger.info("Required extracted spectra check complete.")

            # Record the right number of order expected according to the GNIRS XD configuration.
            if 'Long' in calpath and 'SXD' in scipath:
                orders = [3, 4, 5]
            elif 'Long' in calapth and 'LXD' in scipath:
                orders = [3, 4, 5, 6, 7, 8]
            elif 'Short' in calpath and 'SXD' in scipath:
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

            ###########################################################################
            ##                                                                       ##
            ##                 COMPLETE - OBSERVATION SPECIFIC SETUP                 ##
            ##             BEGIN TELLURIC CORRECTION FOR AN OBSERVATION              ##
            ##                                                                       ##
            ###########################################################################

            valindex = start
            while valindex > stop or valindex < 1 or stop > 5:
                logger.warning("#####################################################################")
                logger.warning("#####################################################################")
                logger.warning("#                                                                   #")
                logger.warning("#   WARNING in telluric: invalid start/stop values of telluric      #")
                logger.warning("#                        correction steps.                          #")
                logger.warning("#                                                                   #")
                logger.warning("#####################################################################")
                logger.warning("#####################################################################\n")

                valindex = int(raw_input("Please enter a valid start value (1 to 7, default 1): "))
                stop = int(raw_input("Please enter a valid stop value (1 to 7, default 7): "))

            while valindex <= stop:

                #############################################################################
                ##  STEP 1: Get telluric information by querrying SIMBAD if not obtained   ## 
                ##          from the configuration file.                                   ##
                #############################################################################

                if valindex == 1:
                    
                    telluricinfofilename = 'telluric_info.txt'
                    getTelluricInfo(telsrcextractedspectrum, telluricrRA, telluricDEC, telluricSpectralType, \
                        telluricMagnitude, telluricTemperature, telluricinfofilename, overwrite)

                    logger.info("##################################################################")
                    logger.info("#                                                                #")
                    logger.info("#       STEP 1: Get telluric information - COMPLETED             #")
                    logger.info("#                                                                #")
                    logger.info("##################################################################\n")

                #############################################################################
                ##  STEP 1: Clean raw observations.                                        ##
                ##  Output: Cleaned science or telluric frames.                            ##
                #############################################################################

                elif valindex == 2:
                    
                    hLineCorrection(rawFrame, grating, hLineInter, hLineMethod, tempInter, log, over)

                    logger.info("##################################################################")
                    logger.info("#                                                                #")
                    logger.info("#       STEP 1: Clean raw observations - COMPLETED               #")
                    logger.info("#                                                                #")
                    logger.info("##################################################################\n")

                #############################################################################
                ##  STEP 1: Clean raw observations.                                        ##
                ##  Output: Cleaned science or telluric frames.                            ##
                #############################################################################

                elif valindex == 3:
                    
                    fitContinuum(rawFrame, grating, continuumInter, tempInter, log, over)
                    
                    logger.info("##################################################################")
                    logger.info("#                                                                #")
                    logger.info("#       STEP 1: Clean raw observations - COMPLETED               #")
                    logger.info("#                                                                #")
                    logger.info("##################################################################\n")

                #############################################################################
                ##  STEP 1: Clean raw observations.                                        ##
                ##  Output: Cleaned science or telluric frames.                            ##
                #############################################################################

                elif valindex == 4:
                    
                    divideByContinuum(rawFrame, log, over)
                    
                    logger.info("##################################################################")
                    logger.info("#                                                                #")
                    logger.info("#       STEP 1: Clean raw observations - COMPLETED               #")
                    logger.info("#                                                                #")
                    logger.info("##################################################################\n")

                #############################################################################
                ##  STEP 1: Clean raw observations.                                        ##
                ##  Output: Cleaned science or telluric frames.                            ##
                #############################################################################

                elif valindex == 5:
                    
                    getShiftScale(rawFrame, telluricInter, log, over)
                    
                    logger.info("##################################################################")
                    logger.info("#                                                                #")
                    logger.info("#       STEP 1: Clean raw observations - COMPLETED               #")
                    logger.info("#                                                                #")
                    logger.info("##################################################################\n")
                
                #############################################################################
                ##  STEP 1: Clean raw observations.                                        ##
                ##  Output: Cleaned science or telluric frames.                            ##
                #############################################################################
                # Shift and scale the telluric correction spectrum and continuum fit to the telluric correction spectrum.
                elif valindex == 6:
                    
                    shiftScaleSpec(rawFrame, "2_fit", "6_shiftedFit", log, over)
                    shiftScaleSpec(rawFrame, "3_chtel", "7_schtel", log, over)
                    
                    logger.info("##################################################################")
                    logger.info("#                                                                #")
                    logger.info("#       STEP 1: Clean raw observations - COMPLETED               #")
                    logger.info("#                                                                #")
                    logger.info("##################################################################\n")

                #############################################################################
                ##  STEP 1: Clean raw observations.                                        ##
                ##  Output: Cleaned science or telluric frames.                            ##
                #############################################################################

                elif valindex == 7:
                    
                    divideCubebyTel(rawFrame, log, over)
                    
                    logger.info("##################################################################")
                    logger.info("#                                                                #")
                    logger.info("#       STEP 1: Clean raw observations - COMPLETED               #")
                    logger.info("#                                                                #")
                    logger.info("##################################################################\n")
                
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
            logger.info("#  %s", scipath                                                                )
            logger.info("#                                                                            #")
            logger.info("##############################################################################\n")

    # Return to directory script was begun from.
    os.chdir(path)
    
    return

##################################################################################################################
#                                                     ROUTINES                                                   #
##################################################################################################################

def getTelluricInfo(telsrcextractedspectrum, telluricRA, telluricDEC, telluricSpectralType, telluricMagnitude, telluricTemperature, telluricinfofilename, overwrite):
    """
    Find telluric star spectral type, temperature and/or magnitude, and exposure time. Based on XDGNIRS code. Modified 
    from nifsTelluric code.

    Executes a SIMBAD query and parses the resulting html to find spectal type, temperature and/or magnitude.

    Reads:
        Extracted 1D telluric source spectrum.
    """
    logger = log.getLogger('gnirsTelluric.getTelluricInfo')

    # If user did not specify a telluricMagnitude, telluricTemperature, telluricRA, or telluricDEC, get telluricRA and 
    # telluricDEC from the extracted telluric spectrum header. Use SIMBAD to look up telluricSpectralType, 
    # telluricMagnitude, and telluricTemperature.
    if (not telluricMagnitude or not telluricTemperature) and (not telluricRA or not telluricDEC):
        telHeader = fits.open(extractedtelluricspectrum)[0].header
        if not telluricRA:
            telluricRA = telHeader['RA']
        if not telluricDEC:
            telluricDEC = telHeader['DEC']
        '''
        # Format RA and DEC to pass to SIMBAD.
        if '-' in str(telluricDEC):
            telluricCoordinates = str(telluricRA)+'d'+str(telluricDEC)+'d'
        else:
            telluricCoordinates = str(telluricRA)+'d+'+str(telluricDEC)+'d'
        '''
        telExptime = str(telHeader['EXPTIME'])
    else:
        # Get the telluric exposure time anyways
        telHeader = fits.open(extractedtelluricspectrum)[0].header
        telExptime = str(telHeader['EXPTIME'])

    # Check to see if a spectral type or temperature has been given
    if telluricTemperature:
        spectraltypeFind = False
        temperatureFind = False
    else:
        spectraltypeFind = True
        temperatureFind = True

    if telluricMagnitude:
        magnitudeFind = False
    else:
        magnitudeFind = True

    if spectraltypeFind or temperatureFind or magnitudeFind:
        # Construct URL based on telluric coordinates and execute SIMBAD query to find the spectral type
        Simbad.add_votable_fields('flux(K)', 'sp')  ## Viraja:  Why does it query only the K magnitude?
        simbad_startable = Simbad.query_region(coord.SkyCoord(ra=telluricRA, dec=telluricDEC, unit=(u.deg, u.deg), \
            frame='fk5'), radius=0.1 * u.deg)  ## Viraja:  How are the RA and DEC formatted in XDpiped.csh??

        if spectraltypeFind:
            # Get spectral type -- only the first 3 characters (strip off end of types like AIVn as they are not in the 
            # 'reference_starstable.txt')
            telluricSpectralType = simbad_startable['SP_TYPE'][0][0:3]
        else:
            logger.error("Cannot locate the spectral type of the telluric in the table generated by the SIMBAD query.")
            logger.error("Please update the parameter 'telluricSpectralType' in the configuration file.")
            raise SystemExit
        
        if magnitudeFind:
            telluricMagnitude = str(simbad_startable['FLUX_K'][0])
        else:
            logger.error("Cannot find a the K magnitude for the telluric in the table generated by the SIMBAD query.")
            logger.error("Please update the parameter 'telluricMagnitude' in the configuration file. Exiting script.")
            raise SystemExit

        if temperatureFind:
            # Find temperature for the spectral type in 'reference_starstable.txt'
            count = 0
            for line in reference_starstable:
                if '#' in line:
                    continue
                else:
                    if telluricSpectralType in line.split()[0]:
                        telluricTemperature = line.split()[1]
                        count = 0
                        break
                    else:
                        count += 1
            if count > 0:  ## Viraja:  I am wondering why this condition is given and why is this an error??
                logger.error("Cannot find a temperature for spectral type %s of the telluric", telluricSpectralType)
                logger.error("Please update the parameter 'telluricTemperature' in the configuration file.")
                raise SystemExit

    telluricInfoTable = open(telluricinfofilename,'w')
    telluricInfoTable.write('k K '+telluricMagnitude+' '+telluricTemperature+'\n')
    #not sure if the rest of these lines are necessary now we're not using "sentelluricInfoTableunc" etc. for flux calibration
    telluricInfoTable.write('h H '+telluricMagnitude+' '+telluricTemperature+'\n')  ## Viraja: Why the same as the K for the rest of them? --> comment in mag2mass.py in XDGNIRS 
    telluricInfoTable.write('j J '+telluricMagnitude+' '+telluricTemperature+'\n')
    telluricInfoTable.write('j J '+telluricMagnitude+' '+telluricTemperature+'\n')
    telluricInfoTable.write('i J '+telluricMagnitude+' '+telluricTemperature+'\n')
    telluricInfoTable.write('i J '+telluricMagnitude+' '+telluricTemperature+'\n')
    telluricInfoTable.close()

    with open(telluricinfofilename,'r'):
        telluricInfoTableLines = telluricInfoTable.readlines()
    logger.info("Contents of %s:", telluricinfofilename)
    for line in telluricInfoTableLines:
        logger.info(line)

#---------------------------------------------------------------------------------------------------------------------#

def hLineCorrection(rawFrame, grating, hLineInter, hLineMethod, tempInter, log, over):
    """
    Remove hydrogen lines from the spectrum of a telluric standard,
    using a model of vega's atmosphere.

    Args:

    Returns:

    Reads:
        0_tel + rawFrame + .fits  : Combined, extracted one d standard star spectrum

    Writes:
        1_htel + rawFrame + .fits : Hline corrected standard star spectrum
    """

    # Remove H lines from standard star correction spectrum
    if os.path.exists("1_htel" + rawFrame + ".fits"):
        if over:
            os.remove("1_htel" + rawFrame + ".fits")
            if hLineMethod == "vega":
                vega(rawFrame, grating, hLineInter, log, over)
        else:
            logging.info("Output file exists and -over- not set - skipping H line removal")
    else:
        if hLineMethod == "vega":
            vega(rawFrame, grating, hLineInter, log, over)

    #if hLineMethod == "linefitAuto" and not no_hLine:
    #    linefitAuto(combined_extracted_1d_spectra, grating)

    # Disabled and untested because interactive scripted iraf tasks are broken...
    #if hLineMethod == "linefitManual" and not no_hLine:
    #    linefitManual(combined_extracted_1d_spectra+'[sci,1]', grating)

    #if hLineMethod == "vega_tweak" and not no_hLine:
        #run vega removal automatically first, then give user chance to interact with spectrum as well
    #    vega(combined_extracted_1d_spectra,grating, path, hLineInter, telluric_shift_scale_record, log, over)
    #    linefitManual("final_tel_no_hLines_no_norm", grating)

    #if hLineMethod == "linefit_tweak" and not no_hLine:
        #run Lorentz removal automatically first, then give user chance to interact with spectrum as well
    #    linefitAuto(combined_extracted_1d_spectra,grating)
    #    linefitManual("final_tel_no_hLines_no_norm", grating)

    if tempInter:
        # Plot the non-hLine corrected spectrum and the h-line corrected spectrum.
        uncorrected = astropy.io.fits.open("1_htel" + rawFrame + ".fits")[1].data
        corrected = astropy.io.fits.open("1_htel" + rawFrame + ".fits")[0].data
        plt.title('Before and After HLine Correction')
        plt.plot(uncorrected)
        plt.plot(corrected)
        plt.show()

#---------------------------------------------------------------------------------------------------------------------#

def vega(extractedspectrum, orders, hLineInter, overwrite):
    """
    Use IRAF TELLURIC task to remove H lines from the telluric star, then remove normalization added by telluric with 
    IMARITH.
    """
    logger = log.getLogger('gnirsTelluric.vega')

    # This task will be done interactively if parameter 'hLineInter' in the configuration file is 'yes'.
    tel_hLineInfo = iraf.telluric(input='v'+spectrum+'[sci,'+str(ext)+']',output='tell_nolines'+str(ext),cal='vega_ext.fits['+str(ext)+']',airmass=airmass_standard,answer='yes',ignoreaps='yes',xcorr='yes',tweakrms='yes',inter=telluricinter,sample=tell_sample1,threshold=0.1,lag=3,shift=0.,dshift=0.0,scale=1.,dscale=0.0,offset=0,smooth=1,cursor='',mode='al',Stdout=1)
    #record shift & scale info for future reference
    t1.write(str(tell_info) + '\n')
    #need this loop to identify telluric output containing warning about pix outside calibration limits (different formatting)
    if "limits" in tell_info[-1].split()[-1]:
        norm = tell_info[-2].split()[-1]
    else:
        norm = tell_info[-1].split()[-1]
    iraf.imarith (operand1='tell_nolines'+str(ext),op="/",operand2=norm,result='ftell_nolines'+str(ext),title='',divzero=0.0,hparams='',pixtype='',calctype='',verbose='no',noact='no',mode='al')




    if grating=='K':
        ext = '1'
        sample = "21537:21778"
        scale = 0.8
    elif grating=='H':
        ext = '2'
        sample = "16537:17259"
        scale = 0.7
    elif grating=='J':
        ext = '3'
        sample = "11508:13492"
        scale = 0.885
    elif grating=='Z':
        ext = '4'
        sample = "*"
        scale = 0.8
    else:
        logging.info("\nWARNING: invalid standard star band. Exiting this correction.")
        return
    if os.path.exists("1_htel" + rawFrame + ".fits"):
            if over:
                os.remove("1_htel" + rawFrame + ".fits")
                iraf.chdir(os.getcwd())
                tell_info = iraf.telluric(input="0_tel" + rawFrame + ".fits[1]", output="1_htel"+rawFrame, cal= RUNTIME_DATA_PATH+'vega_ext.fits['+ext+']', xcorr='yes', tweakrms='yes', airmass=1.0, inter=hLineInter, sample=sample, threshold=0.1, lag=3, shift=0., dshift=0.05, scale=scale, dscale=0.05, offset=0., smooth=1, cursor='', mode='al', Stdout=1)
            else:
                logging.info("Output file exists and -over not set - skipping H line correction")
                return
    else:
        iraf.chdir(os.getcwd())
        tell_info = iraf.telluric(input="0_tel" + rawFrame + ".fits[1]", output="1_htel"+rawFrame, cal= RUNTIME_DATA_PATH+'vega_ext.fits['+ext+']', xcorr='yes', tweakrms='yes', airmass=1.0, inter=hLineInter, sample=sample, threshold=0.1, lag=3, shift=0., dshift=0.05, scale=scale, dscale=0.05, offset=0., smooth=1, cursor='', mode='al', Stdout=1)

    # need this loop to identify telluric output containing warning about pix outside calibration limits (different formatting)
    if "limits" in tell_info[-1].split()[-1]:
        norm=tell_info[-2].split()[-1]
    else:
        norm=tell_info[-1].split()[-1]

    if os.path.exists("final_tel_no_hLines_no_norm.fits"):
        if over:
            # Subtle bugs in iraf mean imarith doesn't work. So we use an astropy/numpy solution.
            # Open the image and the scalar we will be dividing it by.
            operand1 = astropy.io.fits.open("1_htel" + rawFrame+'.fits')[0].data
            operand2 = float(norm)
            # Create a new data array
            multiplied = np.array(operand1, copy=True)
            # Don't forget to include the original header! If you don't later IRAF tasks get confused.
            header = astropy.io.fits.open("1_htel" + rawFrame+'.fits')[0].header
            for i in range(len(multiplied)):
                if operand2 != 0:
                    multiplied[i] = operand1[i] / operand2
                else:
                    multiplied[i] = 1
            # Set the data and header of the in-memory image
            hdu = astropy.io.fits.PrimaryHDU(multiplied)
            hdu.header = header
            # Finally, write the new image to a new .fits file. It only has one extension; zero, with a header and data.
            hdu.writeto('final_tel_no_hLines_no_norm.fits')
            #iraf.imarith(operand1="1_htel" + rawFrame, op='/', operand2=norm, result='final_tel_no_hLines_no_norm', title='', divzero=0.0, hparams='', pixtype='', calctype='', verbose='yes', noact='no', mode='al')
        else:
            logging.info("Output file exists and -over not set - skipping H line normalization correction")
    else:
        #iraf.imarith(operand1="1_htel" + rawFrame, op='/', operand2=norm, result='final_tel_no_hLines_no_norm', title='', divzero=0.0, hparams='', pixtype='', calctype='', verbose='yes', noact='no', mode='al')
        operand1 = astropy.io.fits.open("1_htel" + rawFrame+'.fits')[0].data
        operand2 = float(norm)
        multiplied = np.array(operand1, copy=True)
        header = astropy.io.fits.open("1_htel" + rawFrame+'.fits')[0].header
        for i in range(len(multiplied)):
            if operand2 != 0:
                multiplied[i] = operand1[i] / operand2
            else:
                multiplied[i] = 1
        hdu = astropy.io.fits.PrimaryHDU(multiplied)
        hdu.header = header
        hdu.writeto('final_tel_no_hLines_no_norm.fits')

    if os.path.exists('final_tel_no_hLines_no_norm.fits'):
        os.remove("1_htel" + rawFrame + ".fits")
        shutil.move('final_tel_no_hLines_no_norm.fits', "1_htel" + rawFrame + ".fits")

#---------------------------------------------------------------------------------------------------------------------#

def fitContinuum(rawFrame, grating, continuumInter, tempInter, log, over):
    """
    Fit a continuum to the telluric correction spectrum to normalize it. The continuum
    fitting regions were derived by eye and can be improved.

    Results are in fit<Grating>.fits
    """
    # These were found to fit the curves well by hand. You can probably improve them; feel free to fiddle around!
    if grating == "K":
        order = 5
        sample = "20279:20395,20953:24283"
    elif grating == "J":
        order = 5
        sample = "11561:12627,12745:12792,12893:13566"
    elif grating == "H":
        order = 5
        sample = "*"
    elif grating == "Z":
        order = 5
        sample = "9453:10015,10106:10893,10993:11553"
    if os.path.exists('2_fit'+rawFrame+'.fits'):
        if over:
            os.remove('2_fit'+rawFrame+'.fits')
            iraf.continuum(input='1_htel'+rawFrame,output='2_fit'+rawFrame,ask='yes',lines='*',bands='1',type="fit",replace='no',wavescale='yes',logscale='no',override='no',listonly='no',logfiles=log,inter=continuumInter,sample=sample,naverage=1,func='spline3',order=order,low_rej=1.0,high_rej=3.0,niterate=2,grow=1.0,markrej='yes',graphics='stdgraph',cursor='',mode='ql')
        else:
            logging.info("\nOutput exists and -over not set - skipping continuum fit to telluric correction")
    else:
        iraf.continuum(input='1_htel'+rawFrame,output='2_fit'+rawFrame,ask='yes',lines='*',bands='1',type="fit",replace='no',wavescale='yes',logscale='no',override='no',listonly='no',logfiles=log,inter=continuumInter,sample=sample,naverage=1,func='spline3',order=order,low_rej=1.0,high_rej=3.0,niterate=2,grow=1.0,markrej='yes',graphics='stdgraph',cursor='',mode='ql')
    if os.path.exists('../products_fluxcal_AND_telluric_corrected/0_fit'+rawFrame+'.fits'):
        if over:
            os.remove('../products_fluxcal_AND_telluric_corrected/0_fit'+rawFrame+'.fits')
            shutil.copy('2_fit'+rawFrame+'.fits', '../products_fluxcal_AND_telluric_corrected/0_fit'+rawFrame+'.fits')
        else:
            logging.info("\nOutput exists and -over not set - skipping copy of fit to products_fluxcal_AND_telluric_corrected")
    else:
        shutil.copy('2_fit'+rawFrame+'.fits', '../products_fluxcal_AND_telluric_corrected/0_fit'+rawFrame+'.fits')

    if tempInter:
        # Plot the telluric correction spectrum with the continuum fit.
        final_tel_no_hLines_no_norm = astropy.io.fits.open('1_htel'+rawFrame+'.fits')[0].data
        fit = astropy.io.fits.open('2_fit'+rawFrame+'fit.fits')[0].data
        plt.title('Unnormalized Telluric Correction and Continuum fit Used to Normalize')
        plt.plot(final_tel_no_hLines_no_norm)
        plt.plot(fit)
        plt.show()

#---------------------------------------------------------------------------------------------------------------------#

def divideByContinuum(rawFrame, log, over):
    """
    Divide the standard star spectrum by the continuum to normalize it.
    """
    if os.path.exists("3_chtel"+rawFrame+'.fits'):
        if over:
            os.remove("3_chtel"+rawFrame+'.fits')
            # This is related to issue #3
            #iraf.imarith("1_htel"+rawFrame+'.fits', "/", "2_fit"+rawFrame+'.fits', result="3_chtel"+rawFrame+'.fits',title='',divzero=0.0,hparams='',pixtype='',calctype='',verbose='no',noact='no',mode='al')
            operand1 = astropy.io.fits.open("1_htel"+rawFrame+'.fits')[0].data
            operand2 = astropy.io.fits.open("2_fit"+rawFrame+'.fits')[0].data
            header = astropy.io.fits.open("1_htel"+rawFrame+'.fits')[0].header
            multiplied = np.array(operand1, copy=True)
            for i in range(len(multiplied)):
                if operand2[i] != 0:
                    multiplied[i] = operand1[i] / operand2[i]
                else:
                    multiplied[i] = 0.0
            hdu = astropy.io.fits.PrimaryHDU(multiplied)
            hdu.header = header
            hdu.writeto("3_chtel"+rawFrame+".fits")
            logging.info("\nDivided telluric correction by continuum")
        else:
            logging.info("\nOutput exists and -over not set - skipping division by continuum")
    else:
        # This is related to issue #3
        #iraf.imarith('1_htel'+rawFrame+'.fits', "/", '2_fit'+rawFrame+'.fits', result='3_chtel'+rawFrame+'.fits',title='',divzero=0.0,hparams='',pixtype='',calctype='',verbose='no',noact='no',mode='al')
        operand1 = astropy.io.fits.open("1_htel"+rawFrame+'.fits')[0].data
        operand2 = astropy.io.fits.open("2_fit"+rawFrame+'.fits')[0].data
        header = astropy.io.fits.open("1_htel"+rawFrame+'.fits')[0].header
        multiplied = np.array(operand1, copy=True)
        for i in range(len(multiplied)):
            if operand2[i] != 0:
                multiplied[i] = operand1[i] / operand2[i]
            else:
                multiplied[i] = 0.0
        hdu = astropy.io.fits.PrimaryHDU(multiplied)
        hdu.header = header
        hdu.writeto("3_chtel"+rawFrame+".fits")
        logging.info("\nDivided telluric correction by continuum")

#---------------------------------------------------------------------------------------------------------------------#

def get1dSpecFromCube(rawFrame, log, over):
    """
    Turn a cube into a 1D spec, used to find shift and scale values of telluric spectrum.
    Currently: Extracts 1D spectra from center of cube.
    """
    cube = astropy.io.fits.open('ctfbrsn'+rawFrame+'.fits')
    cubeheader = cube[1].header
    cubeslice = cube[1].data[:,30,30]
    # Create a PrimaryHDU object to encapsulate the data and header.
    hdu = astropy.io.fits.PrimaryHDU(cubeslice)
    # Modify the cd1_1 and CRVAL1 values; this adds the wavelength calibration to the correct cube dimension.
    hdu.header = cubeheader
    hdu.header['CRVAL1'] = cubeheader['CRVAL3']
    hdu.header['CD1_1'] = cubeheader['CD3_3']
    hdu.header['CRPIX1'] = 1.
    if os.path.exists('4_cubeslice'+rawFrame+'.fits'):
        if over:
            os.remove('4_cubeslice'+rawFrame+'.fits')
            # Write the spectrum and header to a new .fits file.
            hdu.writeto('4_cubeslice'+rawFrame+'.fits', output_verify="ignore")
        else:
            logging.info("\nOutput exists and -over not set - skipping extraction of single cube slice")
    else:
        # Write the spectrum and header to a new .fits file.
        hdu.writeto('4_cubeslice'+rawFrame+'.fits', output_verify="ignore")

#---------------------------------------------------------------------------------------------------------------------#

def getShiftScale(rawFrame, telluricInter, log, over):
    """
    Use iraf.telluric() to get the best shift and scale of a telluric correction spectrum.

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
        logging.info("\nRead a shift for "+inPrefix+" spectrum for " + str(tellshift))
        logging.info("\nRead a scale of "+inPrefix+" spectrum for  " + str(scale))
    except IOError:
        logging.info("\nNo shiftScale file found for " + rawFrame + " in " + str(os.getcwd() + ". Skipping."))
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
            logging.info("\nOutput exists and -over not set - skipping shift and scale of " + inPrefix)
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

def divideCubebyTel(rawFrame, log, over):
    """
    Divide every element of a data cube by the derived telluric correction spectrum.
    """
    # Open the uncorrected data cube.
    cube = astropy.io.fits.open('ctfbrsn'+rawFrame+'.fits')
    # Open the shifted, scaled telluric correction spectrum.
    telluricSpec = astropy.io.fits.open('7_schtel'+rawFrame+'.fits')
    if os.path.exists("actfbrsn"+rawFrame+'.fits'):
        if over:
            os.remove("actfbrsn"+rawFrame+'.fits')
            # Divide each slice of cube by telluric correction spectrum.
            for i in range(cube[1].header['NAXIS2']):         # NAXIS2 is the y axis of the final cube.
                for j in range(cube[1].header['NAXIS1']):     # NAXIS1 is the x axis of the final cube.
                    cube[1].data[:,i,j] /= (telluricSpec[0].data)
            # Write the telluric corrected cube to a new file.
            cube.writeto("actfbrsn"+rawFrame+'.fits', output_verify='ignore')
        else:
            logging.info("\nOutput exists and -over not set - skipping application of telluric correction to cube")
    else:
        for i in range(cube[1].header['NAXIS2']):         # NAXIS2 is the y axis of the final cube.
            for j in range(cube[1].header['NAXIS1']):     # NAXIS1 is the x axis of the final cube.
                cube[1].data[:,i,j] /= (telluricSpec[0].data)
        cube.writeto("actfbrsn"+rawFrame+'.fits', output_verify='ignore')

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
            logging.info("It looks as if you didn't use the i key to write out the lineless spectrum. We'll have to try again. --> Re-entering splot")
            iraf.splot(images=spectrum, new_image='final_tel_no_hLines_no_norm', save_file='../PRODUCTS/lorentz_hLines.txt', overwrite='yes')
