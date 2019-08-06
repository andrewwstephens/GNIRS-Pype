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

    logger.info('#################################################')
    logger.info('#                                               #')
    logger.info('#       Start GNIRS Telluric Correction         #')
    logger.info('#                                               #')
    logger.info('#################################################\n')

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
            # Print the runtime data directory path
            tempScipath = scipath.split(os.sep)        
            runtimedatapath = "/".join(tempScipath[:-5])+'/runtimeData'
            logger.info("The runtime data path is %s\n", runtimedatapath)

            # Check if required runtime data files and extracted spectra available in their respective paths
            logger.info("Checking if required runtime data files and extracted spectra available in %s", runtimedatapath)
            logger.info("and %s, respectively.", scipath)

            reference_startable_filename = 'stars_spectraltypes_temperatures.txt'
            if os.path.exists(runtimedatapath+'/'+reference_startable_filename):
                logger.info("Required reference star table of spectral types and temperatures available.")
                reference_startable = open(runtimedatapath+'/'+reference_startable_filename, "r").readlines()
                reference_startable.close()
            else:
                logger.warning("Required reference star table of spectral types and temperatures not available.")
                logger.warning("Please provide the star table in %s. Exiting script.\n", runtimedatapath)
                raise SystemExit
            
            telluric_hLineRegions_filename = runtimedatapath+'/telluric_hLineRegions.dat'
            if os.path.exists(runtimedatapath+'/'+telluric_hLineRegions_filename):
                logger.info("Required telluric H line regions reference file available.")
                telluric_hLineRegions = open(runtimedatapath+'/'+telluric_hLineRegions_filename, "r").readlines()
                telluric_hLineRegions.close()
            else:
                logger.warning("Required telluric H line regions reference file not available. Please provide the ")
                logger.warning("file in %s. Exiting script.\n", runtimedatapath)
                raise SystemExit

            telluric_continuumRegions_filename = runtimedatapath+'/telluric_continuumRegions.dat'
            if os.path.exists(runtimedatapath+'/'+telluric_continuumRegions_filename):
                logger.info("Required telluric continuum regions reference file available.")
                telluric_continuumRegions = open(runtimedatapath+'/'+telluric_continuumRegions_filename, "r").readlines()
                telluric_continuumRegions.close()
            else:
                logger.warning("Required telluric continuum regions reference file not available. Please provide the ")
                logger.warning("file in %s. Exiting script.\n", runtimedatapath)
                raise SystemExit

            telluric_telluricRegions_filename = runtimedatapath+'/telluric_telluricRegions.dat'
            if os.path.exists(runtimedatapath+'/'+telluric_telluricRegions_filename):
                logger.info("Required telluric regions reference file available.")
                telluric_telluricRegions = open(runtimedatapath+'/'+telluric_telluricRegions_filename, "r").readlines()
                telluric_telluricRegions.close()
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
            
            # Output filenames without the extension '.fits'
            telluric_hLineCorrectedSpectrum = telpath+'/h'+telsrcextractedspectrum[telsrcextractedspectrum.rfind('/')+1:]
            telluric_hLineCorrectedSpectrum = telluric_hLineCorrectedSpectrum[telluric_hLineCorrectedSpectrum.rfind('.')+1:]
            
            telluric_fitContinuum = telpath+'/fith'+telsrcextractedspectrum[telsrcextractedspectrum.rfind('/')+1:]
            telluric_fitContinuum = telluric_fitContinuum[telluric_fitContinuum.rfind('.')+1:]

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
                ##  Output: An ascii file containing telluric information.                 ##
                #############################################################################

                if valindex == 1:
                    
                    telinfofilename = 'telluric_info.txt'
                    getTelluricInfo(telHeader, telluricrRA, telluricDEC, telluricSpectralType, telluricMagnitude, \
                        telluricTemperature, telluricinfofilename, overwrite)

                    logger.info("##################################################################")
                    logger.info("#                                                                #")
                    logger.info("#       STEP 1: Get telluric information - COMPLETED             #")
                    logger.info("#                                                                #")
                    logger.info("##################################################################\n")

                #############################################################################
                ##  STEP 2: H line removal.                                                ##
                ##  Output: Telluric source spectrum corrected for H line absorption.      ##
                #############################################################################

                elif valindex == 2:
                    
                    hLineRemoval(telsrcextractedspectrum, telluric_hLineCorrectedSpectrum, hLineInter, orders, \
                        hLineMethod, telluric_hLineRegions, telAirmass, vega_spectrum, tempInter, overwrite)

                    logger.info("##################################################################")
                    logger.info("#                                                                #")
                    logger.info("#       STEP 2: H line removal - COMPLETED                       #")
                    logger.info("#                                                                #")
                    logger.info("##################################################################\n")

                #############################################################################
                ##  STEP 3: Fit telluric continuum.                                        ##
                ##  Output: Cleaned science or telluric frames.                            ##
                #############################################################################

                elif valindex == 3:
                    
                    fitTelluricContinuum(telluric_hLineCorrectedSpectrum, telluric_fitContinuum, continuumInter, \
                        orders, tempInter, overwrite)
                    
                    logger.info("##################################################################")
                    logger.info("#                                                                #")
                    logger.info("#       STEP 3: Fit telluric continuum - COMPLETED               #")
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

def getTelluricInfo(telHeader, telluricRA, telluricDEC, telluricSpectralType, telluricMagnitude, telluricTemperature, telinfofilename, overwrite):
    """
    Find telluric star spectral type, temperature, magnitude, and exposure time. Based on XDGNIRS code. Modified from
    nifsTelluric code.

    Executes a SIMBAD query and parses the resulting html to find spectal type, temperature, and magnitude.

    Reads:
        Extracted 1D telluric source spectrum.
    
    Writes:
        Telluric information file in the science directory path.
    """
    logger = log.getLogger('gnirsTelluric.getTelluricInfo')

    # This code structure checks if output files already exist. If output files exist and overwrite is specified, all 
    # output is overwritten.
    if os.path.exists(telinfofilename):
        if overwrite:
            logger.info("Removing old %s", telinfofilename)
            os.remove(telinfofilename)

            # If user did not specify a telluricMagnitude, telluricTemperature, telluricRA, or telluricDEC, get telluricRA and 
            # telluricDEC from the extracted telluric spectrum header. Use SIMBAD to look up telluricSpectralType, 
            # telluricMagnitude, and telluricTemperature.
            if (not telluricMagnitude or not telluricTemperature) and (not telluricRA or not telluricDEC):
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
                    # 'reference_startable.txt')
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
                    # Find temperature for the spectral type in 'reference_startable.txt'
                    count = 0
                    for line in reference_startable:
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

            telinfotable = open(telinfofilename,'w')
            telinfotable.write('k K '+telluricMagnitude+' '+telluricTemperature+'\n')
            # Viraja: Why the same as the K for the rest of them?? --> comment in mag2mass.py in XDGNIRS -- "Not sure if the 
            # rest of these lines are necessary now we're not using "sensfunc" etc. for flux calibration
            telinfotable.write('h H '+telluricMagnitude+' '+telluricTemperature+'\n')
            telinfotable.write('j J '+telluricMagnitude+' '+telluricTemperature+'\n')
            telinfotable.write('j J '+telluricMagnitude+' '+telluricTemperature+'\n')
            telinfotable.write('i J '+telluricMagnitude+' '+telluricTemperature+'\n')
            telinfotable.write('i J '+telluricMagnitude+' '+telluricTemperature+'\n')
            telinfotable.close()

            with open(telinfofilename,'r'):
                telinfotable = telinfotable.readlines()
            logger.info("Contents of %s:", telinfofilename)
            for line in telinfotable:
                logger.info("%s", line)  ## Viraja:  Check if this works (only if an entire line is read as a string?)
        else:
            logger.info("Output file exists and -overwrite not set - skipping getting telluric information.\n")
    else:
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
            telHeader = fits.open(telsrcextractedspectrum)[0].header
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
                # 'reference_startable.txt')
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
                # Find temperature for the spectral type in 'reference_startable.txt'
                count = 0
                for line in reference_startable:
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

        telinfotable = open(telinfofilename,'w')
        telinfotable.write('k K '+telluricMagnitude+' '+telluricTemperature+'\n')
        # Viraja: Why the same as the K for the rest of them?? --> comment in mag2mass.py in XDGNIRS -- "Not sure if the 
        # rest of these lines are necessary now we're not using "sensfunc" etc. for flux calibration
        telinfotable.write('h H '+telluricMagnitude+' '+telluricTemperature+'\n')
        telinfotable.write('j J '+telluricMagnitude+' '+telluricTemperature+'\n')
        telinfotable.write('j J '+telluricMagnitude+' '+telluricTemperature+'\n')
        telinfotable.write('i J '+telluricMagnitude+' '+telluricTemperature+'\n')
        telinfotable.write('i J '+telluricMagnitude+' '+telluricTemperature+'\n')
        telinfotable.close()

        with open(telinfofilename,'r'):
            telinfotable = telinfotable.readlines()
        logger.info("Contents of %s:", telinfofilename)
        for line in telinfotable:
            logger.info("%s", line)  ## Viraja:  Check if this works (only if an entire line is read as a string?)

#---------------------------------------------------------------------------------------------------------------------#

def hLineRemoval(telsrcextractedspectrum, telluric_hLineCorrectedSpectrum, hLineInter, orders, hLineMethod, telluric_hLineRegions, telAirmass, vega_spectrum, tempInter, overwrite):
    """
    Remove hydrogen (H) absorption lines from the extracted 1D telluric source spectrum. Output filename prefix is 'h'.

    Reads:
        Extracted 1D telluric source spectrum.

    Writes:
        H line corrected 1D telluric source spectrum.
    """
    logger = log.getLogger('gnirsTelluric.hLineRemoval')

    for i in range(len(orders)):
        for line in telluric_hLineRegions:
            if 'order'+str(orders[i]) in line:
                telluric_hLineRegions_sample = line.split()[1]
            else:
                pass
        
        extension = i+1
        telluric_hLineCorrectedInput = telsrcextractedspectrum[telsrcextractedspectrum.rfind('.')+1:]+'[SCI,'+\
            str(extension)+']'
        calibrationInput = vega_spectrum+'['+str(extension)+']'
        telluric_hLineCorrectedOutput = telluric_hLineCorrectedSpectrum+'_order'+str(extension)+'.fits'

        if os.path.exists(telluric_hLineCorrectedOutput):
            if overwrite:
                os.remove(telluric_hLineCorrectedOutput)
                logger.warning("Removing old %s", telluric_hLineCorrectedOutput)

                iraf.hedit(images=telluric_hLineCorrectedInput, fields='AIRMASS', value=telAirmass, add='yes', \
                    addonly='no', delete='no', verify='no', show='no', update='yes', mode='al') 

                logger.info("Removing H lines from the telluric 1D source spectrum, extension %s.", str(extension))
                if hLineMethod == "vega":
                    vega(telluric_hLineCorrectedInput, calibrationInput, telluric_hLineCorrectedOutput, hLineInter, \
                        vega_spectrum, telAirmass, telluric_hLineRegions_sample, overwrite)
                
                # TODO(Viraja):  Check with Marie about the no_hLine option and update the functions for the rest of
                # the H line removal options (commented out below).
                ''' 
                # Untested because interactive scripted iraf tasks are broken...
                if hLineMethod == "lineFitAuto" and not no_hLine:
                    lineFitAuto(combined_extracted_1d_spectra, grating)

                if hLineMethod == "lineFitManual" and not no_hLine:
                    lineFitManual(combined_extracted_1d_spectra+'[sci,1]', grating)

                if hLineMethod == "vega_tweak" and not no_hLine:
                    # First, do H line removal using the 'vega' method automatically, and then give user a chance to  
                    # interact with the spectrum
                    vega(combined_extracted_1d_spectra, path, hLineInter, telluric_shift_scale_record, overwrite)
                    lineFitManual("final_tel_no_hLines_no_norm", grating)

                if hLineMethod == "lineFit_tweak" and not no_hLine:
                    # First, do H line removal using the 'lineFitAuto' method automatically, then give user a chance to 
                    # interact with the spectrum
                    lineFitAuto(combined_extracted_1d_spectra,grating)
                    lineFitManual("final_tel_no_hLines_no_norm", grating)
                '''

                # TODO(Viraja):  Check while testing this script if this is the right position to use iraf.wmef to add
                # the primary header back to the spectra files. Might change depending on whether the other H line
                # remova methods retain the primary headers (they probably do not).
                iraf.wmef(input=telluric_hLineCorrectedOutput+'[SCI,'+str(extension)+']', \
                    output=telluric_hLineCorrectedOutput, extnames='', phu=telsrcextractedspectrum, verbose='yes', \
                    mode='al')

                # TODO(Viraja):  Incorporate a temporary pause to display and take a look at the telluric 1D source 
                # spectra before and after correction. I am not sure if this can be more useful than just looking at  
                # the spectrum. If this has to be done, it would require adding another science extension to the MEFs  
                # created in the previous step using iraf.wmef task BEFORE starting the H line correction in this loop.
                if tempInter:
                    # Plot the telluric 1D source spectrum with and without H line correction
                    telluric_hLineUncorrected = fits.open()  ## not completed
                    telluric_hLineCorrected = fits.open()  ## not completed
                    plt.title('Telluric 1D Source Spectrum Before and After H Line Removal')
                    plt.plot(telluric_hLineUncorrected)
                    plt.plot(telluric_hLineCorrected)
                    plt.show()
            else:
                logger.info("Output file exists and -overwrite not set - skipping H line removal.\n")
        else:
            iraf.hedit(images=telluric_hLineCorrectedInput, fields='AIRMASS', value=telAirmass, add='yes', \
                addonly='no', delete='no', verify='no', show='no', update='yes', mode='al')

            if hLineMethod == "vega":
                vega(telluric_hLineCorrectedInput, calibrationInput, telluric_hLineCorrectedOutput, hLineInter, \
                    vega_spectrum, telAirmass, telluric_hLineRegions_sample, overwrite)

            iraf.wmef(input=telluric_hLineCorrectedOutput+'[SCI,'+str(extension)+']', \
                output=telluric_hLineCorrectedOutput, extnames='', phu=telsrcextractedspectrum, verbose='yes', \
                mode='al')
            
            # TODO(Viraja):  Incorporate a temporary pause to display and take a look at the telluric 1D source spectra
            # before and after correction. I am not sure if this can be more useful than just looking at the spectrum. 
            # If this has to be done, it would require adding another science extension to the MEFs created in the 
            # previous step using iraf.wmef task BEFORE starting the H line correction in this loop.
            if tempInter:
                # Plot the telluric 1D source spectrum with and without H line correction
                telluric_hLineUncorrected = fits.open()  ## not completed
                telluric_hLineCorrected = fits.open()  ## not completed
                plt.title('Telluric 1D Source Spectrum Before and After H Line Removal')
                plt.plot(telluric_hLineUncorrected)
                plt.plot(telluric_hLineCorrected)
                plt.show()

#---------------------------------------------------------------------------------------------------------------------#

def vega(telluric_hLineCorrectedInput, calibrationInput, telluric_hLineCorrectedOutput, hLineInter, telAirmass, telluric_hLineRegions_sample, overwrite):
    """
    Use IRAF TELLURIC task to remove H absorption lines from the telluric 1D source spectrum using the 'vega' method, 
    then remove normalization added by TELLURIC with IMARITH.
    """
    logger = log.getLogger('gnirsTelluric.vega')

    # For each order, this will be done interactively if parameter 'hLineInter' in the configuration file is 'yes'
    telluric_hLineInfo = iraf.telluric(input=telluric_hLineCorrectedInput, \
        output=telluric_hLineCorrectedOutput, cal=calibrationInput, ignoreaps='yes', xcorr='yes', tweakrms='yes', \
        interactive=hLineInter, sample=telluric_hLineRegions_sample, threshold=0.1, lag=3, shift=0., dshift=1.0, \
        scale=1., dscale=0.2, offset=0, smooth=1, cursor='', airmass=telAirmass, answer='yes', interp=poly5, \
        mode='al', Stdout=1)
    
    # This loop identifies telluric output containing warning about pixels outside calibration limits (different 
    # formatting)  Viraja:  This comment is from HlinesXD.py in XDGNIRS. I am not sure what the "different formatting"
    # means.
    if 'limits' in telluric_hLineInfo[-1].split()[-1]:
        normalization = telluric_hLineInfo[-2].split()[-1]
    else:
        normalization = telluric_hLineInfo[-1].split()[-1]
    # Remove normalization added by iraf.telluric
    iraf.imarith(operand1=telluric_hLineCorrectedOutput, operand2=normalization, op='/', \
        result=telluric_hLineCorrectedOutput, title='', divzero=0.0, hparams='', pixtype='', calctype='', \
        verbose='yes', noact='no', mode='al')
    '''
    # Comment from nifsTelluric -- There are subtle bugs in iraf mean imarith does not work. So we use an astropy/numpy
    # solution. Open the image and the scalar we will be dividing it by.
    # Viraja:  I am not sure what "... subtle bugs in iraf mean imarith does not work" means. The comment is not clear.
    This commented section can be implemented if the imarith way from XDGNIRS is later found to not work properly.
    operand1 = fits.open(telluric_hLineCorrectedOutput)[0].data
    operand2 = float(normalization)
    # Create a new data array
    if operand2 != 0:
        operand1 = operand1 / operand2
    else:
        operand1 = 1  ## Viraja:  Why is operand1 set to 1? I would imagine it is set to itself as the divisor is 0.
    # Set the data and header of the in-memory image
    '''

#---------------------------------------------------------------------------------------------------------------------#

def fitTelluricContinuum(telluric_hLineCorrectedSpectrum, telluric_fitContinuum, continuumInter, orders, tempInter, overwrite):
    """
    Fit the continua of the H line corrected telluric 1D source spectra using IRAF CONTINUUM task to normalize them.
    """
    logger = log.getLogger('gnirsTelluric.fitTelluricContinuum')

    # Continuum fitting order depends on the spectral order. Here, iraf.continuum is used with type='fit' to write out 
    # the continuum fit.
    for i in range(len(orders)):
        for line in telluric_continuumRegions:
            if 'order'+str(orders[i]) in line:
                telluric_continuumRegions_sample = line.split()[1]
            else:
                pass
        
        extension = i+1
        telluric_fitContinuumInput = telluric_hLineCorrectedSpectrum+'_order'+str(extension)+'.fits'
        telluric_fitContinuumOutput = telluric_fitContinuumOutput+'_order'+str(extension)+'.fits'

        if os.path.exists(telluric_fitContinuumOutput):
            if overwrite:
                os.remove(telluric_fitContinuumOutput)
                logger.warning("Removing old %s", telluric_fitContinuumOutput)

                logger.info("Fitting telluric continuum for the H line corrected spectra, extension %s.", str(extension))
                iraf.continuum(input=telluric_fitContinuumInput, output=telluric_fitContinuumOutput, lines='*',bands='1',type="fit",replace='no',wavescale='yes',logscale='no',override='no',listonly='no',logfiles='',inter=continuuminter,sample=cont_sample,naverage=1,ask='yes',func='spline3',order=2,low_rej=1.0,high_rej=3.0,niterate=2,grow=1.0,markrej='yes',graphics='stdgraph',cursor='',mode='ql')

    if os.path.exists('2_fit'+rawFrame+'.fits'):
        if over:
            os.remove('2_fit'+rawFrame+'.fits')
            iraf.continuum(input='1_htel'+rawFrame,output='2_fit'+rawFrame,ask='yes',lines='*',bands='1',type="fit",replace='no',wavescale='yes',logscale='no',override='no',listonly='no',logfiles=log,inter=continuumInter,sample=sample,naverage=1,func='spline3',order=order,low_rej=1.0,high_rej=3.0,niterate=2,grow=1.0,markrej='yes',graphics='stdgraph',cursor='',mode='ql')
        else:
            logger.info("\nOutput exists and -over not set - skipping continuum fit to telluric correction")
    else:
        iraf.continuum(input='1_htel'+rawFrame,output='2_fit'+rawFrame,ask='yes',lines='*',bands='1',type="fit",replace='no',wavescale='yes',logscale='no',override='no',listonly='no',logfiles=log,inter=continuumInter,sample=sample,naverage=1,func='spline3',order=order,low_rej=1.0,high_rej=3.0,niterate=2,grow=1.0,markrej='yes',graphics='stdgraph',cursor='',mode='ql')
    if os.path.exists('../products_fluxcal_AND_telluric_corrected/0_fit'+rawFrame+'.fits'):
        if over:
            os.remove('../products_fluxcal_AND_telluric_corrected/0_fit'+rawFrame+'.fits')
            shutil.copy('2_fit'+rawFrame+'.fits', '../products_fluxcal_AND_telluric_corrected/0_fit'+rawFrame+'.fits')
        else:
            logger.info("\nOutput exists and -over not set - skipping copy of fit to products_fluxcal_AND_telluric_corrected")
    else:
        shutil.copy('2_fit'+rawFrame+'.fits', '../products_fluxcal_AND_telluric_corrected/0_fit'+rawFrame+'.fits')

    if tempInter:
        # Plot the telluric correction spectrum with the continuum fit.
        final_tel_no_hLines_no_norm = astropy.io.fits.open()
        fit = astropy.io.fits.open()
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
            logger.info("\nDivided telluric correction by continuum")
        else:
            logger.info("\nOutput exists and -over not set - skipping division by continuum")
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
        logger.info("\nDivided telluric correction by continuum")

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
            logger.info("\nOutput exists and -over not set - skipping extraction of single cube slice")
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
            # TODO(nat): implement logger for this
            iraf.chdir(os.getcwd())
            tell_info = iraf.telluric(input='4_cubeslice'+rawFrame+'.fits[0]',output='5_oneDCorrected'+rawFrame+'.fits',cal="3_chtel"+rawFrame+'.fits[0]',airmass=1.0,answer='yes',ignoreaps='yes',xcorr='yes',tweakrms='yes',inter=telluricInter,sample="*",threshold=0.1,lag=3,shift=0.,dshift=0.1,scale=1.0,dscale=0.1, offset=1,smooth=1,cursor='',mode='al',Stdout=1)
        else:
            logger.info("\nOutput exists and -over not set - skipping get shift scale of telluric correction and fit")
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
            logger.info("\nOutput exists and -over not set - skipping application of telluric correction to cube")
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
            logger.info("It looks as if you didn't use the i key to write out the lineless spectrum. We'll have to try again. --> Re-entering splot")
            iraf.splot(images=spectrum, new_image='final_tel_no_hLines_no_norm', save_file='../PRODUCTS/lorentz_hLines.txt', overwrite='yes')
