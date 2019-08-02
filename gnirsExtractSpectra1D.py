#!/usr/bin/env python

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

import log, os, sys, glob, ConfigParser, datetime, dateutil.parser
from astropy.io import fits
from pyraf import iraf, iraffunctions
import numpy as np

## IMPORTANT NOTE:  The command line options are not available for GNIRS as of July 2019.


def start(configfile):
    """
    This module contains all the functions needed to perform the full reduction of SCIENCE or TELLURIC data.

    Parameters are loaded from gnirs.cfg configuration file. This script will automatically detect if it is being run
    on telluric data or science data. There are 5 steps.

    INPUT FILES:
        - Configuration file
        - Science or Telluric frames
        - mdfshiftrefimage
        - masterflat
        - /database files from the appropriate calibrations directory

    OUTPUT FILES:
        - If telluric:  cleaned (optional), prepared, radiation-event corrected, reduced, spatial distortion corrected, 
          and transformed images
        - If science:  cleaned (optional), prepared, radiation-event corrected, reduced, spatial distortion corrected, 
          and transformed images

    Args:
        - kind (string): Either 'Science' or 'Telluric'
        - configfile: gnirs.cfg configuration file.
                - Paths to the Science (str), reduction truth value (boolean)
                  E.g. 'target/date/config/{Sci,Tel}_ObsID/{Calibrations,Intermediate}', True
                - Paths to the Tellurics (str), reduction truth value (boolean)
                  E.g. 'target/date/config/{Sci,Tel}_ObsID/{Calibrations,Intermediate}', True
                - manualMode (boolean): Enable optional manualModeging pauses? Default: False
                - overwrite (boolean): Overwrite old files? Default: False
                # And gnirsReduce specific settings
    """
    logger = log.getLogger('gnirsExtractSpectra1D.start')

    ###########################################################################
    ##                                                                       ##
    ##                  BEGIN - GENERAL EXTRACT 1D SETUP                     ##
    ##                                                                       ##
    ###########################################################################

    # Store current working directory for later use.
    path = os.getcwd()

    logger.info('####################################################')
    logger.info('#                                                  #')
    logger.info('#        Start Extracting GNIRS 1D Spectra         #')
    logger.info('#                                                  #')
    logger.info('####################################################\n')

    # Set up/prepare IRAF.
    iraf.gemini()
    iraf.gemtools()
    iraf.gnirs()

    # Reset to default parameters the used IRAF tasks.
    iraf.unlearn(iraf.gemini,iraf.gemtools,iraf.gnirs,iraf.imcopy)

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
    observationSections = ['ScienceDirectories','TelluricDirectories']
    nsextractInter = config.getboolean('interactive','nsextractInter')
    calculateSpectrumSNR = config.getboolean('gnirsPipeline','calculateSpectrumSNR')
    # extract1Spectra1D specific config
    useApall = config.getboolean('extractSpectra1D','useApall')
    extractionApertureRadius = config.getfloat('extractSpectra1D','extractionApertureRadius')
    toleranceOffset = config.getfloat('extractSpectra1D','toleranceOffset')
    extractStepwise = config.getboolean('extractSpectra1D','extractStepwise')
    extractionStepsize = config.getfloat('extractSpectra1D','extractionStepsize')

    ###########################################################################
    ##                                                                       ##
    ##                 COMPLETE - GENERAL EXTRACT 1D SETUP                   ##
    ##                                                                       ##
    ###########################################################################

    # gnirsExtractSpectra1D will first check if the reduction truth value of the science and telluric directories is 
    # True -- if it is, it will then check if the required spectra to be extracted are available in the directories
    # (and proceed only if it finds them there); else, it will warn the user and request to provide the spectra for 
    # extracting. If the reduction truth value of the science and telluric directories is False, the script will skip
    # extracting 1D spectra in those directories.

    # Loop through all the observation (telluric and science) directories to extract 1D spectra in each one.
    for section in observationSections:
        for obspath in config.options(section):
            if config.getboolean(section,obspath):

                ###########################################################################
                ##                                                                       ##
                ##                  BEGIN - OBSERVATION SPECIFIC SETUP                   ##
                ##                                                                       ##
                ###########################################################################

                os.chdir(obspath)
                # Change the iraf directory to the current directory.
                iraffunctions.chdir(obspath)

                # Print the current directory of observations being reduced.
                logger.info("Currently working on extracting 1D spectra in %s\n", obspath)
                
                # Check if required combined spectra available in the observations directory path
                logger.info("Checking if required combined spectra available in %s", obspath)

                srccombimage = 'src_comb.fits'
                if os.path.exists(srccombimage):
                    logger.info("Required combined source image available.")
                else:
                    logger.warning("Required combined source image not available. Please run ")
                    logger.warning("gnirsCombineSpectra2D.pyto create the combined source image or provide it ")
                    logger.warning("manually. Exiting script.\n")
                    raise SystemExit
                if calculateSpectrumSNR:
                    skycombimage = obspath+'/sky_comb.fits'
                    if os.path.exists(srccombimage):
                        logger.info("Required combined sky image available.")
                    else:
                        logger.warning("Parameter 'calculateSpectrumSNR' is 'True', but required combined sky image ")
                        logger.warning("not available. Setting the 'calculateSpectrumSNR' parameter for the current ")
                        logger.warning("set of observations to 'False'.\n")
                        calculateSpectrumSNR = False
                
                logger.info("Required combined spectra check complete.")
                

                if 'Telluric' in section:
                    calpath = obspath+'/Calibrations'
                    logger.debug("obspath: %s", obspath)
                    logger.debug("calpath: %s", calpath)
                


                '''
                # Make an extraction database directory in the observations directory path
                extractdatabasepath = os.mkdir('/database')
                # Print the extraction database directory path
                logger.info("The path to the extraction database is %s\n", extractiondatabasepath)
                
                ###########################################################################
                ##                                                                       ##
                ##                 COMPLETE - OBSERVATION SPECIFIC SETUP                 ##
                ##            BEGIN EXTRACTING 1D SPECTRA FOR AN OBSERVATION             ##
                ##                                                                       ##
                ###########################################################################
                    
                if manualMode:
                    a = raw_input("About to enter extract 1D spectra.")

                if useApall:
                    # This performs a weighted extraction
                    apertureTracingColumns = 20
                    extractSpectra1D(srccombimage, nsextractInter, useApall, apertureTracingColumns, \
                        extractionApertureRadius, extractdatabasepath, overwrite)
                else:
                    apertureTracingColumns = 10
                    extractSpectra1D(srccombimage, nsextractInter, useApall, apertureTracingColumns, \
                        extractionApertureRadius, extractdatabasepath, overwrite)
                # If the parameter 'calculateSpectrumSNR' is set to 'yes', the script will extract spectra from 
                # the combined sky image; else, it will only extract spectra from the combined source image.
                if calculateSpectrumSNR:
                    logger.info("Extracting the combined sky spectrum reduced without sky subtraction.\n")
                    if useApall:
                        apertureTracingColumns = 20
                        extractSpectra1D(srccombimage, nsextractInter, useApall, apertureTracingColumns, \
                            extractionApertureRadius, extractdatabasepath, overwrite)
                    else:
                        apertureTracingColumns = 10
                        extractSpectra1D(srccombimage, nsextractInter, useApall, apertureTracingColumns, \
                            extractionApertureRadius, extractdatabasepath, overwrite)
                '''
            logger.info("##############################################################################")
            logger.info("#                                                                            #")
            logger.info("#  COMPLETE - Extracting 1D spectra completed for                            #")
            logger.info("#  %s", obspath                                                                )
            logger.info("#                                                                            #")
            logger.info("##############################################################################\n")

    # Return to directory script was begun from.
    os.chdir(path)

    return

##################################################################################################################
#                                                     ROUTINES                                                   #
##################################################################################################################

def extractSpectra1D(combinedimage, nsextractInter, useApall, apertureTracingColumns, extractionApertureRadius, extractdatabasepath, overwrite):
    """
    Extracting 1D spectra from the combined 2D spectra using nsextract.
    """
    logger = log.getLogger('gnirsReduce.extractSpectra1D')

    if os.path.exists('v'+combinedimage):
        if overwrite:
            logger.warning("Removing old v%s", combinedimage)
            os.remove('v'+combinedimage)
            iraf.nsextract(inimages=combinedimage, outspectra='', outprefix='v', dispaxis=1, \
                database=extractdatabasepath, line=700, nsum=apertureTracingColumns, ylevel='INDEF', \
                upper=str(extractionApertureRadius), lower='-'+str(extractionApertureRadius), background='none', \
                fl_vardq='yes', fl_addvar='no', fl_skylines='yes', fl_inter=nsextractInter, fl_apall=useApall, \
                fl_trace='no', aptable='gnirs$data/apertures.fits', fl_usetabap='no', fl_flipped='yes', \
                fl_project='yes', fl_findneg='no', bgsample='*', trace='', tr_nsum=10, tr_step=10, tr_nlost=3, \
                tr_function='legendre', tr_order=5, tr_sample='*', tr_naver=1, tr_niter=0, tr_lowrej=3.0, \
                tr_highrej=3.0, tr_grow=0.0, weights='variance', logfile=logger.root.handlers[0].baseFilename, \
                verbose='yes', mode='al')
        else:
            logger.warning("Old %s exists and -overwrite not set - skipping nsextract for observations.", combinedimage)
    else:
        iraf.nsextract(inimages=combinedimage, outspectra='', outprefix='v', dispaxis=1, database=extractdatabasepath,\
            line=700, nsum=apertureTracingColumns, ylevel='INDEF', upper=str(extractionApertureRadius), \
            lower='-'+str(extractionApertureRadius), background='none', fl_vardq='yes', fl_addvar='no', \
            fl_skylines='yes', fl_inter=nsextractInter, fl_apall=useApall, fl_trace='no', \
            aptable='gnirs$data/apertures.fits', fl_usetabap='no', fl_flipped='yes', fl_project='yes', \
            fl_findneg='no', bgsample='*', trace='', tr_nsum=10, tr_step=10, tr_nlost=3, tr_function='legendre', \
            tr_order=5, tr_sample='*', tr_naver=1, tr_niter=0, tr_lowrej=3.0, tr_highrej=3.0, tr_grow=0.0, \
            weights='variance', logfile=logger.root.handlers[0].baseFilename, verbose='yes', mode='al')

#---------------------------------------------------------------------------------------------------------------------#
'''
def func():
    if nsextinter== 'no':
        std_peak = column_lookup('database/apstandard_comb_SCI_') 
        tgt_peak = column_lookup('database/aptarget_comb_SCI_')
        
        posfile = open ('../LISTS/qoffsets.txt')
        tgtq = float(posfile.readline())
        stdq = float(posfile.readline())
        obs = pyfits.open(target)
        pixscale = obs[0].header["PIXSCALE"]
        diff = (tgtq - stdq)/pixscale 
        reExtract = False
        found_spectrum = []
        #nsextract should find spectrum within 'tolerance' pixels of expected location
        #depends on how well the observer centred the target along the slit. 5 pix is an initial guess at what's reasonable.
        #(rather than an offset, would be better to use some measure of whether the peak found by nsextract was real, e.g. counts + FWHM. Not recorded in database, though.)
        tolerance = 5 
        for i in range (0,6):
            expected = '%s' % float('%.4g' % (std_peak[i] + diff))
            if tgt_peak[i] == 'not found':
                print 'In extension '+str(i+1)+' nscombine did not extract anything. Re-extracting with the aperture forced to be at '+expected
                found_spectrum.append(0)
                reExtract=True
            else:    
                found = '%s' % float('%.4g' % tgt_peak[i])
                if abs((tgt_peak[i] - diff) - std_peak[i]) < tolerance:
                    print '***CHECK: In extension '+str(i+1)+' nscombine detected the target spectrum close to the expected location along slit (x = '+found+' vs expected x = '+expected+')'
                    found_spectrum.append(1)
                else:
                    print '***WARNING: In extension '+str(i+1)+' nscombine extracted an unexpected location along the slit (x = '+found+' vs expected x = '+expected+'). It is probably extracting noise; re-extracting with the aperture forced to be at the expected location.'
                    found_spectrum.append(0)
                    reExtract=True
            
        if reExtract:        
            #Re-extract science target spectrum if needed
            for i in range (0,6):
                #create new aperture files in database   
                #ran into trouble when only replacing ones that weren't well centred, so just replacing them all 
                #(but using nsextract peak for target when it seemed to be in the right place)
                if os.path.isfile('database/aptarget_comb_SCI_'+str(i+1)+'_'):
                    os.remove ('database/aptarget_comb_SCI_'+str(i+1)+'_')
                apfile = open ('database/apstandard_comb_SCI_'+str(i+1)+'_', 'r')
                newapfile = open ('database/apnewstandard_comb_SCI_'+str(i+1)+'_', 'w')
                if found_spectrum[i] == 0:
                    clean  = apfile.read().replace(str(std_peak[i]), str(std_peak[i] + diff)+' ').replace('standard_comb','newstandard_comb')
                else:
                    clean  = apfile.read().replace(str(std_peak[i]), str(tgt_peak[i])+' ').replace('standard_comb','newstandard_comb')
                newapfile.write(clean)
                apfile.close()
                newapfile.close()
            shutil.copy(standard, 'new'+standard)
            iraf.imdelete(images='v'+target)
            #Using these settings in nsextract will force it to use the aperture size and centre in the edited "apnewstandard_comb" files
            iraf.nsextract(inimages=target,outspectra='',outprefix='v',dispaxis=1,database='',line=700,nsum=20,ylevel='INDEF',upper=aperture_plus,lower=aperture_minus,background='none',fl_vardq='yes',fl_addvar='no',fl_skylines='yes',fl_inter=nsextinter,fl_apall='yes',fl_trace='no',aptable=path_to_nsextract+'config/apertures.fits',fl_usetabap='no',fl_flipped='yes',fl_project='yes',fl_findneg='no',bgsample='*',trace='newstandard_comb',tr_nsum=10,tr_step=10,tr_nlost=3,tr_function='legendre',tr_order=5,tr_sample='*',tr_naver=1,tr_niter=0,tr_lowrej=3.0,tr_highrej=3.0,tr_grow=0.0,weights='variance',logfile='',verbose='yes',mode='al')

        #Slight complication - we occasionally find that nsextract locates the aperture too close to the end of the slit
        #Then it exits with an "Aperture too large" error and spectra aren't extracted for one or more extensions
        #But we work around that above, so when we check for errors in XDpiped.csh, we ignore this error
        #But maybe that error can happen for other reasons
        #So to be on the safe side, will check that all extensions are present in the extracted target file (should really add other files as well...)
        check_extn = iraf.gemextn(inimages=target,check='exists,mef',process='expand',index='',extname='SCI',extversion='',ikparams='',omit='',replace='',outfile='STDOUT',logfile='',Stdout=1,glogpars='',verbose='no',fail_count='0', count='20', status='0')
        if len(check_extn) != 6:
            print "ERROR: target_comb file contains only ", len(check_extn), 'extensions. Exiting script.'
'''
#---------------------------------------------------------------------------------------------------------------------#

if __name__ == '__main__':
    log.configure('gnirs.log', filelevel='INFO', screenlevel='DEBUG')
    start('gnirs.cfg')
