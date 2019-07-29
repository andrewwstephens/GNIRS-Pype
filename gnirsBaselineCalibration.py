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

import log, os, pkg_resources, glob, shutil, ConfigParser, cleanir
from astropy.io import fits
from pyraf import iraf, iraffunctions


def start(configfile):
    """
    This module contains all the functions needed to reduce GNIRS GENERAL BASELINE CALIBRATIONS.

    INPUT FILES FOR BASELINE CALIBRATION:

    INPUT FILES:
        - Configuration file                                                                                                                                                                                                                                              
        - All calibration frames
        - QH flat frames
        - IR flat frames
        - Arc frames
        - Pinhole flat frames

    OUTPUT FILES:
        - masterflat.fits
        - Wavelength-calibrated arc image

    Args:
        # Loaded from gnirs.cfg
        Paths to the calibrations (str): E.g. ['Target/Date/Configuration/Calibrations'].
        manualMode (boolean): Enable optional manualModeging pauses. Default: False.
        overwrite (boolean): Overwrite old files. Default: False.
        start (int): starting step of calibration reductions. Specified at command line with -a. Default: 1.  ## this 
                     option not updated for GNIRS as of July 2019
        stop (int): stopping step of calibration reductions. Specified at command line with -z. Default: 5.  ## this 
                     option not updated for GNIRS as of July 2019
        cleanir_QHflats (boolean): Cleaning QH flat frames. Default: False.
        cleanir_IRflats (boolean): Cleaning IR flat frames. Default: False.
        cleanir_arcs (boolean): Cleaning arc frames. Default: False.
        cleanir_pinholes (boolean): Cleaning pinhole flat frames. Default: False.
    """
    logger = log.getLogger('gnirsBaselineCalibration.start')

    # Store current working directory for later use.
    path = os.getcwd()

    logger.info('##################################################')
    logger.info('#                                                #')
    logger.info('# Start the GNIRS Baseline Calibration Reduction #')
    logger.info('#                                                #')
    logger.info('##################################################\n')

    logger.info("Parameters read from %s", configfile)
    
    # Set up/prepare IRAF.
    iraf.gemini()
    iraf.gnirs()
    iraf.gemtools()

    # Reset to default parameters the used IRAF tasks.
    iraf.unlearn(iraf.gemini,iraf.gemtools,iraf.gnirs)

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
    user_clobber=iraf.envget("clobber")
    iraf.reset(clobber='yes')

    config = ConfigParser.RawConfigParser()
    config.optionxform = str  # make options case-sensitive
    config.read(configfile)
    # Read general config.
    manualMode = config.getboolean('defaults','manualMode')
    overwrite = config.getboolean('defaults','overwrite')
    # Read calibrationReduction specfic config.
    start = config.getint('calibrationReduction','Start')
    stop = config.getint('calibrationReduction','Stop')
    cleanir_QHflats = config.getboolean('calibrationReduction','cleanir_QHflats')
    cleanir_IRflats = config.getboolean('calibrationReduction','cleanir_IRflats')
    cleanir_arcs = config.getboolean('calibrationReduction','cleanir_arcs')
    cleanir_pinholes = config.getboolean('calibrationReduction','cleanir_pinholes')

    ################################################################################
    # Define Variables, Reduction Lists AND identify/run number of reduction steps #
    ################################################################################

    # Loop over the Calibrations directories and reduce the calibrations for the ones with caibration_direcotries True.
    for calpath in config.options('CalibrationDirectories'):
        if config.getboolean('CalibrationDirectories', calpath):
            
            os.chdir(calpath)

            firstfilename = sorted(glob.glob(calpath+'/N*.fits'))[0]
            header = fits.open(firstfilename)[0].header

            iraffunctions.chdir(calpath)

            # However, don't do the reduction for a Calibrations_"configuration" directory without associated telluric
            # or science data. Check that a Calibrations_"configuration" directory exists at the same level as the 
            # Calibrations_"configuration" directory. If not, skip the reduction of calibrations in that 
            # Calibrations_"configuration" directory. "configuration" should be the last directory of calpath.
            ''' TODO(Viraja):  If required, check if the science and telluric directories present for that specific
                               configuration.
            configuration = os.path.split(calpath)[-1]
            print(configuration)
            if not os.path.exists("../"+configuration):

                logger.error("#######################################################################")
                logger.error("#######################################################################")
                logger.error("                                                                       ")
                logger.error("    No configuration directory (including science or telluric data)    ")
                logger.error("    found for  %s.", calpath)
                logger.error("    Skipping reduction of calibrations in that directory.              ")
                logger.error("                                                                       ")
                logger.error("#######################################################################")
                logger.error("#######################################################################")

                continue
            '''
            # Create lists of each type of calibration from textfiles in Calibrations directory.
            allcallist = open("all.list", "r").readlines()
            QHflatlist = open("QHflats.list", "r").readlines()
            IRflatlist = open("IRflats.list", "r").readlines()
            pinholelist = open("pinholes.list", "r").readlines()
            arclist = open("arcs.list", "r").readlines()
            # Write the name of the first flat in pinholelist for which spatial distortion will be calculated.
            sdistrefimage = 'rn'+pinholelist[0].strip()
        
            # Check start and stop values for reduction steps. Ask user for a correction if input is not valid.
            valindex = start
            while valindex > stop  or valindex < 1 or stop > 5:
                logger.warning("######################################################################")
                logger.warning("######################################################################")
                logger.warning("#                                                                    #")
                logger.warning("#   WARNING in calibrate: invalid start/stop values of calibration   #")
                logger.warning("#                         reduction steps.                           #")
                logger.warning("#                                                                    #")
                logger.warning("######################################################################")
                logger.warning("######################################################################\n")

                valindex = int(raw_input("\nPlease enter a valid start value (1 to 5, default 1): "))
                stop = int(raw_input("\nPlease enter a valid stop value (1 to 5, default 5): "))

            # Print the current directory of calibrations being processed.
            logger.info("Currently working on calibrations in %s\n", calpath)

            while valindex <= stop:
            
                #############################################################################
                ##  STEP 1: Clean calibration frames.                                      ##
                ##  Output: Cleaned QHflats, IRflats, arcs, and pinholes.                  ##
                #############################################################################

                if valindex == 1:
                    if manualMode:
                        a = raw_input("About to enter step 1: clean calibration frames.")
                    '
                    if cleanir_arcs:
#                        cleanir(arclist)
                        pass  ## this step is not modifoed to work with this pipeline as of July 2019
                    else:
                        logger.info("######################################################################")
                        logger.info("######################################################################")
                        logger.info("#                                                                    #")
                        logger.info("#               WARNING in calibrate: arcs not cleaned.              #")
                        logger.info("#                                                                    #")
                        logger.info("######################################################################")
                        logger.info("######################################################################\n")
                
                    if cleanir_QHflats:
#                        cleanir(QHflatlist)
                        pass  ## this step is not modifoed to work with this pipeline as of July 2019
                    else:
                        logger.info("######################################################################")
                        logger.info("######################################################################")
                        logger.info("#                                                                    #")
                        logger.info("#              WARNING in calibrate: QHflats not cleaned.            #")
                        logger.info("#                                                                    #")
                        logger.info("######################################################################")
                        logger.info("######################################################################\n")
                
                    if cleanir_IRflats:
#                        cleanir(IRflatlist)
                        pass  ## this step is not modifoed to work with this pipeline as of July 2019
                    else:
                        logger.info("######################################################################")
                        logger.info("######################################################################")
                        logger.info("#                                                                    #")
                        logger.info("#               WARNING in calibrate: IRflats not cleaned.           #")
                        logger.info("#                                                                    #")
                        logger.info("######################################################################")
                        logger.info("######################################################################\n")

                    if cleanir_pinholes:
#                        cleanir(pinholelist)
                        pass  ## this step is not modifoed to work with this pipeline as of July 2019
                    else:
                        logger.info("######################################################################")
                        logger.info("######################################################################")
                        logger.info("#                                                                    #")
                        logger.info("#              WARNING in calibrate: pinholes not cleaned.           #")
                        logger.info("#                                                                    #")
                        logger.info("######################################################################")
                        logger.info("######################################################################\n")
                    
                    logger.info("##################################################################")
                    logger.info("#                                                                #")
                    logger.info("#          STEP 1: Clean calibration frames - COMPLETED          #")
                    logger.info("#                                                                #")
                    logger.info("##################################################################\n")
                    
                #############################################################################
                ##  STEP 2: Preparing all calibration frames for further processing.       ##
                ##  Output: All nsprepare'd and nsreduce'd calibration images.             ##
                #############################################################################

                elif valindex == 2:
                    if manualMode:
                        a = raw_input("About to enter step 2: preparing all calibration frames.")
                    
                    arrayid = header['ARRAYID'].strip()
                    if arrayid == 'SN7638228.1':
                        bpmfile = 'gnirs$data/gnirsn_2011apr07_bpm.fits'
                    elif arrayid == 'SN7638228.1.2':
                        bpmfile = 'gnirs$data/gnirsn_2012dec05_bpm.fits'
                    else:
                        logger.error("######################################################################")
                        logger.error("######################################################################")
                        logger.error("#                                                                    #")
                        logger.error("#       ERROR in calibrate: unknown array ID. Exiting script.        #")
                        logger.error("#                                                                    #")
                        logger.error("######################################################################")
                        logger.error("######################################################################\n")
                        raise SystemExit
                    
                    prepareCalibrations(bpmfile, overwrite)

                    logger.info("####################################################################")
                    logger.info("#                                                                  #")
                    logger.info("#    STEP 2: Preparing all calibration frames (QHflats, IRflats    #")
                    logger.info("#            arcs, and pinholes) - COMPLETED                       #")
                    logger.info("#                                                                  #")
                    logger.info("####################################################################\n")

                #############################################################################
                ##  STEP 3: Create flat field images for QHflats and IR flats.             ##
                ##          Create the final flat by grouping order 3 of IRflats and       ## 
                ##          orders 4-18 of QHflats.                                        ##
                ##  Output: Flat Field image.                                              ##
                #############################################################################

                elif valindex == 3:
                    if manualMode:
                        a = raw_input("About to enter step 3: create flat field.")

                    makeFlat(overwrite)

                    logger.info("###################################################################")
                    logger.info("#                                                                 #")
                    logger.info("#      STEP 3: Create flat field (QHflat, IRflat, and the         #")
                    logger.info("#              final flat) - COMPLETED                            #")
                    logger.info("#                                                                 #")
                    logger.info("###################################################################\n")

                ########################################################################################
                ##  Step 4: Trace the spatial curvature and spectral distortion in the pinhole flat.  ##
                ##  Output: Pinhole flat for which the spatial curvature and S-distortion has been    ##
                ##          determined.                                                               ##
                ########################################################################################

                elif valindex == 4:
                    if manualMode:
                        a = raw_input("About to enter step 4: spatial curvature and spectral distortion.")
                
                    if 'Long' in calpath:
                        pinhole_coordlist = 'gnirs$data/pinholes-long-dense-north.lis'
                    elif 'Short' in calpath:
                        pinhole_coordlist = 'gnirs$data/pinholes-short-dense-north.lis'
                    else:
                        logger.error("#######################################################################")
                        logger.error("#######################################################################")
                        logger.error("#                                                                     #")
                        logger.error("#         ERROR in calibrate: unknown camera. Exiting script.         #")
                        logger.error("#                                                                     #")
                        logger.error("#######################################################################")
                        logger.error("#######################################################################\n")
                        raise SystemExit

                    makeSdistortion(pinholelist, pinhole_coordlist, overwrite)
                
                    logger.info("####################################################################")
                    logger.info("#                                                                  #")
                    logger.info("#     Step 4: Spatial curvature and spectral distortion (trace     #")
                    logger.info("#             the spatial curvature and spectral distortion in     #") 
                    logger.info("#             the pinhole flat) - COMPLETED                        #")
                    logger.info("#                                                                  #")
                    logger.info("####################################################################\n")

                ############################################################################
                ##  STEP 5: Combine arc frames.                                           ##
                ##          Determine the wavelength solution and create the wavelength-  ##
                ##          calibrated arc image corrected for spatial distortion.        ##
                ##  OUTPUT: Combined and wavelength-calibrated arc image.                 ##
                ############################################################################

                elif valindex == 5:
                    if manualMode:
                        a = raw_input("About to enter step 5: wavelength solution.")
                    
                    if '10' in calpath and '32' in calpath:
                        waveCal_coordlist = 'gnirs$data/lowresargon.dat'
                    elif '110' in calpath:
                        waveCal_coordlist = 'gnirs$data/argon.dat'
                    else:
                        logger.error("########################################################################")
                        logger.error("########################################################################")
                        logger.error("#                                                                      #")
                        logger.error("#         ERROR in calibrate: unknown grating. Exiting script.         #")
                        logger.error("#                                                                      #")
                        logger.error("########################################################################")
                        logger.error("########################################################################\n")
                        raise SystemExit

                    makeWaveCal(sdistrefimage, waveCal_coordlist, overwrite)
                    
                    logger.info("####################################################################")
                    logger.info("#                                                                  #")
                    logger.info("#   STEP 5: Wavelength solution (Combine arc frames, determine     #")
                    logger.info("#           the wavelength solution, and create the wavelength     #")
                    logger.info("#           referenced arc image) - COMPLETED                      #")
                    logger.info("#                                                                  #")
                    logger.info("####################################################################\n")

                else:
                    logger.error("########################################################################")
                    logger.error("########################################################################")
                    logger.error("#                                                                      #")
                    logger.error("#         ERROR in calibrate: ", valindex, " is not valid.             #")
                    logger.error("#                                                                      #")
                    logger.error("########################################################################")
                    logger.error("########################################################################\n")
                    raise SystemExit

                valindex += 1

            logger.info("###########################################################################")
            logger.info("#                                                                         #")
            logger.info("#             COMPLETE - Calibration reductions completed for             #")
            logger.info("#  %s", calpath                                                             )
            logger.info("#                                                                         #")
            logger.info("###########################################################################\n")
        
        else:
            logger.warning("########################################################################")
            logger.warning("########################################################################")
            logger.warning("#                                                                      #")
            logger.warning("#    GNIRS baseline calibration turned off. Not reducing baseline      #")
            logger.warning("#      calibrations in %s", calpath                                      )
            logger.warning("#                                                                      #")
            logger.warning("########################################################################")
            logger.warning("########################################################################\n")

    # Return to directory script was begun from.
    os.chdir(path)
        
    return

#####################################################################################
#                                        FUNCTIONS                                  #
#####################################################################################

def prepareCalibrations(bpmfile, overwrite):
    """
    Prepare calibration frames for further processing.

    Use NSPREPARE on all calibration frames (QHflats, IRflats, arcs, pinholes, and darks to update the raw data headers 
    and attach the mask definition file (MDF) as a binary table on all files. Note that dark frames will not have an 
    MDF attached by default. Instead, the appropriate MDF is added in NSREDUCE or NSFLAT to match the data being 
    reduced.

    Use NSREDUCE to cut the calibration frame spectra to the size specified by the MDF, placing different spectral
    orders in separate image extensions. Darks will not be used for further reductions after this point.
    """
    logger = log.getLogger('gnirsBaselineCalibration.prepareCalibrations')
    
    # Update all calibration frames with mdf offset value and generate variance and data quality extensions. This code
    # structure checks if iraf output files already exist. If output files exist and overwrite is specified, iraf
    # output is overwritten.
    oldfiles = glob.glob('./nN*.fits')
    if not oldfiles:
        iraf.nsprepare(inimages='@QHflats.list,@IRflats.list,@arcs.list,@pinholes.list', rawpath='', outimages='', \
            outprefix='n', bpm=bpmfile, logfile=logger.root.handlers[0].baseFilename, fl_vardq='yes', fl_cravg='no', \
            crradius=0.0, fl_dark_mdf='no', fl_correct='no', fl_saturated='yes', fl_nonlinear='yes', \
            fl_checkwcs='yes', fl_forcewcs='yes', arraytable='gnirs$data/array.fits', \
            configtable='gnirs$data/config.fits', specsec='[*,*]', offsetsec='none', pixscale='0.15', shiftimage='', \
            shiftx='INDEF', shifty='INDEF', obstype='FLAT', fl_inter='no', verbose='yes', mode='al')
    else:
        if overwrite:
            logger.warning('Removing old nN*.fits files.')
            [os.remove(filename) for filename in oldfiles]
            iraf.nsprepare(inimages='@QHflats.list,@IRflats.list,@arcs.list,@pinholes.list', rawpath='', outimages='',\
                outprefix='n', bpm=bpmfile, logfile=logger.root.handlers[0].baseFilename, fl_vardq='yes', \
                fl_cravg='no', crradius=0.0, fl_dark_mdf='no', fl_correct='no', fl_saturated='yes', \
                fl_nonlinear='yes', fl_checkwcs='yes', fl_forcewcs='yes', arraytable='gnirs$data/array.fits', \
                configtable='gnirs$data/config.fits', specsec='[*,*]', offsetsec='none', pixscale='0.15', \
                shiftimage='', shiftx='INDEF', shifty='INDEF', obstype='FLAT', fl_inter='no', verbose='yes', mode='al')
        else:
            logger.warning("Output files exist and -overwrite not set - skipping nsprepare for all calibration frames.")
    
    # Cut all calibration frames according to the size specified by the MDFs.
    oldfiles = glob.glob('./rnN*.fits')
    if not oldfiles:
        iraf.nsreduce(inimages='n//@QHflats.list,n//@IRflats.list,n//@arcs.list,n//@pinholes.list', outimages='', \
            outprefix='r', fl_cut='yes', section='', fl_corner='yes', fl_process_cut='yes', fl_nsappwave='no', \
            nsappwavedb='gnirs$data/nsappwave.fits', crval='INDEF', cdelt='INDEF', fl_dark='no', darkimage='', \
            fl_save_dark='no', fl_sky='no', skyimages='', skysection='', combtype='median', rejtype='avsigclip', \
            masktype='goodvalue', maskvalue=0.0, scale='none', zero='median', weight='none', statsec='[*,*]', \
            lthreshold='INDEF', hthreshold='INDEF', nlow=1, nhigh=1, nkeep=0, mclip='yes', lsigma=3.0, hsigma=3.0, \
            snoise='0.0', sigscale=0.1, pclip=-0.5, grow=0.0, skyrange='INDEF', nodsize=3.0, fl_flat='no', \
            flatimage='', flatmin=0.0, fl_vardq='yes', logfile=logger.root.handlers[0].baseFilename, verbose='yes', \
            debug='no', force='no', mode='al')
    else:
        if overwrite:
            logger.warning('Removing old rnN*.fits files.')
            [os.remove(filename) for filename in oldfiles]
            iraf.nsreduce(inimages='n//@QHflats.list,n//@IRflats.list,n//@arcs.list,n//@pinholes.list', outimages='', \
                outprefix='r', fl_cut='yes', section='', fl_corner='yes', fl_process_cut='yes', fl_nsappwave='no', \
                nsappwavedb='gnirs$data/nsappwave.fits', crval='INDEF', cdelt='INDEF', fl_dark='no', darkimage='', \
                fl_save_dark='no', fl_sky='no', skyimages='', skysection='', combtype='median', rejtype='avsigclip', \
                masktype='goodvalue', maskvalue=0.0, scale='none', zero='median', weight='none', statsec='[*,*]', \
                lthreshold='INDEF', hthreshold='INDEF', nlow=1, nhigh=1, nkeep=0, mclip='yes', lsigma=3.0, hsigma=3.0,\
                snoise='0.0', sigscale=0.1, pclip=-0.5, grow=0.0, skyrange='INDEF', nodsize=3.0, fl_flat='no', \
                flatimage='', flatmin=0.0, fl_vardq='yes', logfile=logger.root.handlers[0].baseFilename, \
                verbose='yes', debug='no', force='no', mode='al')
        else:
            logger.info("Output files exist and -overwrite not set - skipping nsreduce for all calibration frames.")

#---------------------------------------------------------------------------------------------------------------------#

def makeFlat(overwrite):
    """
    Make flat field and bad pixel masks.

    Use NSFLAT to generate a normalized flat field for the QHflats and IRflats (for each cross-dispersed order). A mask 
    (BPM) will also be generated by thresholding - this can be used to flag bad pixels in other data.

    Use FXCOPY followed by FXINSERT to produce the final flat by grouping order 3 of IRflats and orders 4-18 of 
    QHflats. The output from this task is used as the flatfield image for further reduction.

    NOTE:  Sometimes nsflat crashes with a fixpix or other, unknown error. In such a situation, tweaking the lthresh 
    parameter sometimes helps. XDGNIRS loops through a fixed list of lthresh values until it (hopefully) runs without
    an error in three trials after which it fails and exits. The nsflat function in this pipeline is currently not set 
    to run through a loop; it simply uses the specified lthresh. However, it can be made to run through a loop of
    random lthresh values in a specified range until it runs without an error to avoid the system exit.
    """
    logger = log.getLogger('gnirsBaselineCalibration.makeFlat')

    # Generate normalized flat field images for QHflats and IRflats.
    QHflat = 'QHflat.fits'
    QHflat_bpm = 'QHflat_bpm.pl'
    if os.path.exists(QHflat):
        if overwrite:
            logger.warning('Removing old %s', QHflat)
            os.remove(QHflat)
            QHflat_flag = False
        else:
            logger.warning("Old %s exists and -overwrite not set - using existing output file for further ")
            logger.warning("reduction.", QHflat)
            QHflat_flag = True
    if os.path.exists(QHflat_bpm):
        if overwrite:
            logger.warning('Removing old %s', QHflat_bpm)
            os.remove(QHflat_bpm)
            QHflat_bpm_flag = False
        else:
            logger.warning("Old %s exists and -overwrite not set - using existing output for further ")
            logger.warning("reduction.", QHflat_bpm)
            QHflat_bpm_flag = True
    if QHflat_flag and QHflat_bpm_flag:
        pass
    else:
        iraf.nsflat(lampson='rn//@QHflats.list', darks='', flatfile=QHflat, darkfile='', fl_corner='yes', \
            fl_save_darks='no', flattitle='default', bpmtitle='default', bpmfile=QHflat_bpm, process="fit", \
            statsec='MDF', fitsec='MDF', thr_flo=0.35, thr_fup=4.0, thr_dlo=-20, thr_dup=100, fl_inter='no', \
            fl_range='no', fl_fixbad='yes', fixvalue=1.0, function='spline3', order=5, normstat='midpt', \
            combtype='default', rejtype='ccdclip', masktype='goodvalue', maskvalue=0.0, scale='median', zero='none', \
            weight='none', lthreshold=50.0, hthreshold='INDEF', nlow=1, nhigh=1, nkeep=0, mclip='yes', lsigma=3.0, \
            hsigma=3.0, snoise='0.0', sigscale=0.1, pclip=-0.5, grow=0.0, box_width=20, box_length=1, trace='', \
            traceproc='none', threshold=100.0, aptable='gnirs$data/apertures.fits', database='', apsum=10, tr_step=10,\
            tr_nlost=3, tr_function='legendre', tr_order=5, tr_naver=1, tr_niter=0, tr_lowrej=3.0, tr_highrej=3.0, \
            tr_grow=0.0, ap_lower=-30, ap_upper=30, fl_vardq='yes', logfile=logger.root.handlers[0].baseFilename, \
            verbose='yes', mode='al')

    IRflat = 'IRflat.fits'
    IRflat_bpm = 'IRflat_bpm.pl'
    if os.path.exists(IRflat):
        if overwrite:
            logger.warning('Removing old %s', IRflat)
            os.remove(IRflat)
            IRflat_flag = False
        else:
            logger.warning("Old %s exists and -overwrite not set - using existing output for further ")
            logger.warning("reduction.", IRflat)
            IRflat_flag = True
    if os.path.exists(IRflat_bpm):
        if overwrite:
            logger.warning('Removing old %s', IRflat_bpm)
            os.remove(IRflat_bpm)
            IRflat_bpm_flag = False
        else:
            logger.warning("Old %s exists and -overwrite not set - using existing output for further ")
            logger.warning("reduction.", IRflat_bpm)
            IRflat_bpm_flag = True
    if IRflat_flag and IRflat_bpm_flag:
        pass
    else:
        iraf.nsflat(lampson='rn//@IRflats.list', darks='', flatfile=IRflat, darkfile='', fl_corner='yes', \
            fl_save_darks='no', flattitle='default', bpmtitle=IRflat_bpm, bpmfile='default', process="fit", \
            statsec='MDF', fitsec='MDF', thr_flo=0.35, thr_fup=1.5, thr_dlo=-20, thr_dup=100, fl_inter='no', \
            fl_range='no', fl_fixbad='yes', fixvalue=1.0, function='spline3', order=10, normstat='midpt', \
            combtype='default', rejtype='ccdclip', masktype='goodvalue', maskvalue=0.0, scale='none', zero='none', \
            weight='none', lthreshold=50.0, hthreshold='INDEF', nlow=1, nhigh=1, nkeep=0, mclip='yes', lsigma=3.0, \
            hsigma=3.0, snoise='0.0', sigscale=0.1, pclip=-0.5, grow=0.0, box_width=20, box_length=1, trace='', \
            traceproc='none', threshold=100.0, aptable='gnirs$data/apertures.fits', database='', apsum=10, tr_step=10,\
            tr_nlost=3, tr_function='legendre', tr_order=5, tr_naver=1, tr_niter=0, tr_lowrej=3.0, tr_highrej=3.0, \
            tr_grow=0.0, ap_lower=-30, ap_upper=30, fl_vardq='yes', logfile=logger.root.handlers[0].baseFilename, \
            verbose='yes', mode='al')

    # Group order 3 of IRflats and orders 4-18 of QHflats to create the final flat field image. Output flat field file 
    # will have name "masterflat.fits".
    masterflat = "masterflat.fits"
    if os.path.exists(masterflat):
        if overwrite:
            logger.warning('Removing old %s', masterflat)
            os.remove(masterflat)
            iraf.fxcopy(input=IRflat, output=masterflat, group="0-3", new_file='yes', verbose='no', mode='ql')
            iraf.fxinsert(input=QHflat, output=masterflat+'[3]', groups="4-18", verbose='no', mode='ql')
        else:
            logger.warning("Output exists and -overwrite not set - skipping fxcopy and fxinsert to create the ")
            logger.warning("masterflat.")
    else:
        iraf.fxcopy(input=IRflat, output=masterflat, group="0-3", new_file='yes', verbose='no', mode='ql')
        iraf.fxinsert(input=QHflat, output=masterflat+'[3]', groups="4-18", verbose='no', mode='ql')

#---------------------------------------------------------------------------------------------------------------------#

def makeSdistortion(pinholelist, pinhole_coordlist, overwrite):
    """
    Establish Spatial-distortion calibration.

    NSSDIST uses the information in the pinhole calibration images to calibrate the spatial distorion of the GNIRS 
    field. The pinhole frame is a dispersed flat field image with a slit-mask in the field so that the illumination on 
    GNIRS is in a pattern of 6 different holes that are stacked in the y-dimension on the field. Proper alignment of
    the slit across the pattern can be used for spatial rectification of the on-sky science data. The spatial solution
    determined by NSSSDIST is linked to the science data in NSFITCOORDS.
    """
    logger = log.getLogger('gnirsBaselineCalibration.makeSdistortion')
   
    oldfiles = glob.glob('./database/idrn*')
    if not oldfiles:
        iraf.nssdist(inimages='rn'+pinholelist[0].strip(), outsuffix='_sdist', pixscale=1.0, dispaxis=1, database='', \
            firstcoord=0.0, coordlist=pinhole_coordlist, aptable='gnirs$data/apertures.fits', fl_inter='no', \
            fl_dbwrite='yes', section='default', nsum=30, ftype='emission', fwidth=10.0, cradius=10.0, \
            threshold=1000.0, minsep=5.0, match=-6.0, function="legendre", order=5, sample='', niterate=3, \
            low_reject=5.0, high_reject=5.0, grow=0.0, refit='yes', step=10, trace='no', nlost=0, aiddebug='', \
            logfile=logger.root.handlers[0].baseFilename, verbose='no', debug='no', force='no', mode='al')
    else:
        if overwrite:
            logger.warning('Removing old /database/idrn* files.')
            [os.remove(filename) for filename in oldfiles]
            iraf.nssdist(inimages='rn'+pinholelist[0].strip(), outsuffix='_sdist', pixscale=1.0, dispaxis=1, \
                database='', firstcoord=0.0, coordlist=pinhole_coordlist, aptable='gnirs$data/apertures.fits', \
                fl_inter='no', fl_dbwrite='yes', section='default', nsum=30, ftype='emission', fwidth=10.0, \
                cradius=10.0, threshold=1000.0, minsep=5.0, match=-6.0, function="legendre", order=5, sample='', \
                niterate=3, low_reject=5.0, high_reject=5.0, grow=0.0, refit='yes', step=10, trace='no', nlost=0, \
                aiddebug='', logfile=logger.root.handlers[0].baseFilename, verbose='no', debug='no', force='no', \
                mode='al')
        else:
            logger.warning("Output exists and -overwrite not set - skipping the spatial-distortion correction ")
            logger.warning("calculation for pinholes.")

#---------------------------------------------------------------------------------------------------------------------#

def makeWaveCal(sdistrefimage, waveCal_coordlist, overwrite):
    """
    Determine the wavelength solution of the combined arc.

    Uses NSWAVELENGTH to calibrate arc data (after cutting and optionally applying a flatfield with NSREDUCE 
    previously).

    DATA REDUCTION HINT -
    For the nswavelength call, the different wavelength settings use different values for some of the parameters. For 
    optimal auto results, use:

    ????

    Note that better RMS fits can be obtained by running the wavelength calibration interactively and identifying all 
    of the lines manually. Tedious, but will give more accurate results than the automatic mode (i.e., fl_inter-). Use 
    fl_inter+ for manual mode.
    """
    logger = log.getLogger('gnirsBaselineCalibration.makeWaveCal')
    
    # Combine arc frames. Output combined file will have the name "arc_comb.fits".
    combinedarc = 'arc_comb.fits'
    if os.path.exists(combinedarc):
        if overwrite:
            logger.warning('Removing old %s', combinedarc)
            os.remove(combinedarc)
            iraf.nscombine(inimages='rn//@arcs.list', tolerance=0.5, output=combinedarc, output_suffix='_comb', \
                bpm='', dispaxis=1, pixscale=1.0, fl_cross='no', fl_keepshift='no', fl_shiftint='yes', \
                interptype='linear', boundary='nearest', constant=0.0, combtype='average', rejtype='sigclip', \
                masktype='goodvalue', maskvalue=0.0, statsec='[*,*]', scale='none', zero='none', weight='none', \
                lthreshold='INDEF', hthreshold='INDEF', nlow=1, nhigh=1, nkeep=0, mclip='yes', lsigma=5.0, hsigma=5.0,\
                ron=0.0, gain=1.0, snoise='0.0', sigscale=0.1, pclip=-0.5, grow=0.0, nrejfile='', fl_vardq='yes', \
                fl_inter='no', logfile=logger.root.handlers[0].baseFilename, verbose='yes', debug='no', force='no', \
                mode='al')
        else:
            logger.warning("Output file exists and -overwrite not set - skipping combining arcs.")
    else:
        iraf.nscombine(inimages='rn//@arcs.list', tolerance=0.5, output=combinedarc, output_suffix='_comb', bpm='', \
            dispaxis=1, pixscale=1.0, fl_cross='no', fl_keepshift='no', fl_shiftint='yes', interptype='linear', \
            boundary='nearest', constant=0.0, combtype='average', rejtype='sigclip', masktype='goodvalue', \
            maskvalue=0.0, statsec='[*,*]', scale='none', zero='none', weight='none', lthreshold='INDEF', \
            hthreshold='INDEF', nlow=1, nhigh=1, nkeep=0, mclip='yes', lsigma=5.0, hsigma=5.0, ron=0.0, gain=1.0, \
            snoise='0.0', sigscale=0.1, pclip=-0.5, grow=0.0, nrejfile='', fl_vardq='yes', fl_inter='no', \
            logfile=logger.root.handlers[0].baseFilename, verbose='yes', debug='no', force='no', mode='al')

    oldfiles_combinedarc = glob.glob('./*f'+combinedarc)
    oldfiles_waveCal = glob.glob('./database/idwtf*')
    oldfiles_sdist = glob.glob('./database/*_sdist')
    oldfiles_lamp = glob.glob('./database/*_lamp')
    oldfiles = oldfiles_combinedarc + oldfiles_waveCal + oldfiles_sdist + oldfiles_lamp
    if not oldfiles:
        logger.info("Running nsfitcoords and nstransform on the combined arc to spatially rectify it before running ") 
        logger.info("nswavelength.")
        iraf.nsfitcoords(inimages=combinedarc, outspectra='', outprefix='f', lamptransf='', \
            sdisttransf=sdistrefimage, dispaxis=1, database='', fl_inter='no', fl_align='no', function='chebyshev', \
            lxorder=2, lyorder=4, sxorder=4, syorder=4, pixscale=1.0, logfile=logger.root.handlers[0].baseFilename, \
            verbose='yes', debug='no', force='no', mode='al')
        iraf.nstransform(inimages='f'+combinedarc, outspectra='', outprefix='t', dispaxis=1, database='', \
            fl_stripe='no', interptype='poly3', xlog='no', ylog='no', pixscale=1.0, \
            logfile=logger.root.handlers[0].baseFilename, verbose='yes', debug='no', mode='al')
            
        logger.info("The arc lines might still be a bit tilted because of the detector rotation. Running ")
        logger.info("nswavelength on the transformed arc with fl_median-.")
        iraf.nswavelength(lampspectra='tf'+combinedarc, outspectra='', outprefix='w', crval='INDEF', cdelt='INDEF', \
            crpix='INDEF', dispaxis=2, database='', coordlist=waveCal_coordlist, fl_inter='no', \
            nsappwavedb='gnirs$data/nsappwave.fits', fl_median='no', sdist='', sdorder=4, xorder=2, yorder=2, \
            aptable='gnirs$data/apertures.fits', section='auto', nsum=10, ftype='emission', fwidth=5.0, cradius=5.0, \
            threshold=300.0, minsep=2.0, match=-6.0, function='chebyshev', order=4, sample='*', niterate=10, \
            low_reject=3.0, high_reject=3.0, grow=0.0, refit='yes', step=2, trace='no', nlost=10, fl_overwrite='yes', \
            aiddebug='', fmatch=0.2, nfound=6, sigma=0.05, rms=0.1, logfile=logger.root.handlers[0].baseFilename, \
            verbose='yes', debug='no', mode='al')
           
        logger.info("Calling nsfitcoords to produce the combined arc with horizontal lines and then nstransform on ")
        logger.info("the transformed arc using the output from nswavelength for the 'lamptrans' parameter.")
        iraf.nsfitcoords(inimages='tf'+combinedarc, outspectra='', outprefix='f', lamptrans='wtf'+combinedarc, \
            sdisttransf='', dispaxis=1, database='', fl_inter='no', fl_align='no', function='chebyshev', lxorder=2, \
            lyorder=4, sxorder=4, syorder=4, pixscale=1.0, logfile=logger.root.handlers[0].baseFilename, \
            verbose='yes', debug='no', force='no', mode='al')
        iraf.nstransform(inimages='ftf'+combinedarc, outspectra='', outprefix='t', dispaxis=1, database='', \
            fl_stripe='no', interptype='poly3', xlog='no', ylog='no', pixscale=1.0, \
            logfile=logger.root.handlers[0].baseFilename, verbose='no', debug='no', mode='al')
    else:
        if overwrite:
            logger.warning('Removing old files.')
            [os.remove(filename) for filename in oldfiles]
            logger.info('Removing old files.')

            logger.info("Running nsfitcoords and nstransform on the combined arc to spatially rectify it before ") 
            logger.info("running nswavelength.")
            iraf.nsfitcoords(inimages=combinedarc, outspectra='', outprefix='f', lamptransf='', \
                sdisttransf=sdistrefimage, dispaxis=1, database='', fl_inter='no', fl_align='no', \
                function='chebyshev', lxorder=2, lyorder=4, sxorder=4, syorder=4, pixscale=1.0, \
                logfile=logger.root.handlers[0].baseFilename, verbose='yes', debug='no', force='no', mode='al')
            iraf.nstransform(inimages='f'+combinedarc, outspectra='', outprefix='t', dispaxis=1, database='', \
                fl_stripe='no', interptype='poly3', xlog='no', ylog='no', pixscale=1.0, \
                logfile=logger.root.handlers[0].baseFilename, verbose='yes', debug='no', mode='al')
            
            logger.info("The arc lines might still be a bit tilted because of the detector rotation. Running ")
            logger.info("nswavelength on the transformed arc with fl_median-.")
            iraf.nswavelength(lampspectra='tf'+combinedarc, outspectra='', outprefix='w', crval='INDEF', \
                cdelt='INDEF', crpix='INDEF', dispaxis=2, database='', coordlist=waveCal_coordlist, fl_inter='no', \
                nsappwavedb='gnirs$data/nsappwave.fits', fl_median='no', sdist='', sdorder=4, xorder=2, yorder=2, \
                aptable='gnirs$data/apertures.fits', section='auto', nsum=10, ftype='emission', fwidth=5.0, \
                cradius=5.0, threshold=300.0, minsep=2.0, match=-6.0, function='chebyshev', order=4, sample='*', \
                niterate=10, low_reject=3.0, high_reject=3.0, grow=0.0, refit='yes', step=2, trace='no', nlost=10, \
                fl_overwrite='yes', aiddebug='', fmatch=0.2, nfound=6, sigma=0.05, rms=0.1, \
                logfile=logger.root.handlers[0].baseFilename, verbose='yes', debug='no', mode='al')
            
            logger.info("Calling nsfitcoords to produce the combined arc with horizontal lines and then nstransform ")
            logger.info("on the transformed arc using the output from nswavelength for the 'lamptrans' parameter.")
            iraf.nsfitcoords(inimages='tf'+combinedarc, outspectra='', outprefix='f', lamptrans='wtf'+combinedarc, \
                sdisttransf='', dispaxis=1, database='', fl_inter='no', fl_align='no', \
                function='chebyshev', lxorder=2, lyorder=4, sxorder=4, syorder=4, pixscale=1.0, \
                logfile=logger.root.handlers[0].baseFilename, verbose='yes', debug='no', force='no', mode='al')
            iraf.nstransform(inimages='ftf'+combinedarc, outspectra='', outprefix='t', dispaxis=1, database='', \
                fl_stripe='no', interptype='poly3', xlog='no', ylog='no', pixscale=1.0, \
                logfile=logger.root.handlers[0].baseFilename, verbose='no', debug='no', mode='al')
        else:
            logger.warning("Output files exist and -overwrite not set - skipping wavelength calibration and ")
            logger.warning("spatial-distortion correction of arcs.")

#---------------------------------------------------------------------------------------------------------------------#

if __name__ == '__main__':
    log.configure('gnirs.log', filelevel='INFO', screenlevel='DEBUG')
    start('gnirs.cfg')
