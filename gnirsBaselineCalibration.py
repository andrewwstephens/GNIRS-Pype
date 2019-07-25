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
    This module contains all the functions needed to reduce GNIRS GENERAL BASELINE CALIBRATIONS

    INPUT FILES FOR EACH BASELINE CALIBRATION:

    Raw files:
        - QH flat frames
        - IR flat frames
        - Arc frames
        - Pinhole flat frames

    OUTPUT FILES:
        - Shift file. Eg: sCALFLAT.fits
        - Bad Pixel Mask. Eg: rgnCALFLAT_sflat_bmp.pl
        - Flat field. Eg: rgnCALFLAT_flat.fits
        - Reduced arc frame. Eg: wrgnARC.fits
        - Reduced ronchi mask. Eg: rgnRONCHI.fits
        - Reduced dark frame. Eg: rgnARCDARK.fits

    Args:
        # Loaded from gnirs.cfg
        calibrationDirectoryList: list of paths to calibrations. ['path/obj/date/Calibrations_grating']
        overwrite (boolean): overwrite old files. Default: False.
        start (int): starting step of daycal reduction. Specified at command line with -a. Default: 1.
        stop (int): stopping step of daycal reduction. Specified at command line with -z. Default: 6.
        manualMode (boolean): enable optional manualModeging pauses. Default: False.
    """
    logger = log.getLogger('BaselineCalibrations')

    # Store current working directory for later use.
    path = os.getcwd()

#    log = os.getcwd()+'/gnirs.log'

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
    iraf.nsheaders("gnirs",logfile=log)

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
    start = config.getint('calibrationReductionConfig','Start')
    stop = config.getint('calibrationReductionConfig','Stop')
    cleanir_arcs = config.getboolean('calibrationReductionConfig','cleanir_arcs')
    cleanir_QHflats = config.getboolean('calibrationReductionConfig','cleanir_QHflats')
    cleanir_IRflats = config.getboolean('calibrationReductionConfig','cleanir_IRflats')
    cleanir_pinholes = config.getboolean('calibrationReductionConfig','cleanir_pinholes')


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

            # However, don't do the reduction for a Calibrations_"configuration" directory without associated telluric or
            # science data. Check that a Calibrations_"configuration" directory exists at the same level as the 
            # Calibrations_"configuration" directory. If not, skip the reduction of calibrations in that 
            # Calibrations_"configuration" directory. "configuration" should be the last directory of calpath.
            configuration = os.path.split(calpath)[-1]
            if not os.path.exists("../"+configuration):

                logger.info("#######################################################################")
                logger.info("                                                                       ")
                logger.info("    No configuration directory (including science or telluric data)    ")
                logger.info("    found for  %s", calpath)
                logger.info("    Skipping reduction of calibrations in that directory.              ")
                logger.info("                                                                       ")
                logger.info("#######################################################################")

                continue

            # Create lists of each type of calibration from textfiles in Calibrations directory.
            allcallist = open('all.list', "r").readlines()
            QHflatlist = open('QHflats.list', "r").readlines()
            IRflatlist = open("IRflats.list", "r").readlines()
            pinholelist = open("pinholes.list", "r").readlines()
            arclist = open("arcs.list", "r").readlines()

            # Write the name of the first nsprepare'd QHflat in allcallist (to be used as the reference image to shift  
            # the MDFs for the rest of the images taken on that same day) into a text file called  
            # "mdfshiftrefimagefile.txt" to be used by the pipeline later.
            open("mdfshiftrefimagefile.txt", "w").write('n'+QHflatlist[0].strip())
            # Write the name of the first flat in pinholelist on which S-distortion will be applied and will be used 
            # as the reference image for S-distortion correction by the pipeline later.
            sdistrefimage = 'rn'+pinholelist[0].strip()
            open("sdistrefimagefile.txt", "w").write('rn'+pinholelist[0].strip())

            masterflat = "masterflat.fits"  ## name of the final flat having order 3 of IRflats and orders 4-18 of 
                                            ## QHflats grouped together
            combinedarc = "arc_comb"  ## name of the arc image combined using gscombine
        
            # Check start and stop values for reduction steps. Ask user for a correction if input is not valid.
            valindex = start
            while valindex > stop  or valindex < 1 or stop > 5:
                logger.info("######################################################################")
                logger.info("######################################################################")
                logger.info("                                                                      ")
                logger.info("    WARNING in calibrate: invalid start/stop values of calibration    ")
                logger.info("                           reduction steps.                           ")
                logger.info("                                                                      ")
                logger.info("######################################################################")
                logger.info("#####################################################################\n")

                valindex = int(raw_input("\nPlease enter a valid start value (1 to 5, default 1): "))
                stop = int(raw_input("\nPlease enter a valid stop value (1 to 5, default 5): "))

            # Print the current directory of calibrations being processed.
            logger.info("##########################################################################")
            logger.info("                                                                          ")
            logger.info("                  Currently working on calibrations in                    ")
            logger.info("  %s", calpath)
            logger.info("                                                                          ")
            logger.info("##########################################################################\n")

            while valindex <= stop:
            
                #############################################################################
                ##  STEP 1: Clean calibration frames.                                      ##
                ##  Output: Cleaned QHflats, IRflats, arcs, and pinholes.                  ##
                #############################################################################

                if valindex == 1:
                    if manualMode:
                        a = raw_input("About to enter step 1: clean calibration frames.")
                    '''
                    if cleanir_arcs:
                        cleanir(arclist)
                    else:
                        logger.info("######################################################################")
                        logger.info("######################################################################")
                        logger.info("                                                                      ")
                        logger.info("               WARNING in calibrate: arcs not cleaned.                ")
                        logger.info("                                                                      ")
                        logger.info("######################################################################")
                        logger.info("###############################################@######################\n")
                
                    if cleanir_QHflats:
                        cleanir(QHflatlist)
                    else:
                        logger.info("######################################################################")
                        logger.info("######################################################################")
                        logger.info("                                                                      ")
                        logger.info("              WARNING in calibrate: QHflats not cleaned.              ")
                        logger.info("                                                                      ")
                        logger.info("######################################################################")
                        logger.info("######################################################################\n")
                
                    if cleanir_IRflats:
                        cleanir(IRflatlist)
                    else:
                        logger.info("######################################################################")
                        logger.info("######################################################################")
                        logger.info("                                                                      ")
                        logger.info("               WARNING in calibrate: IRflats not cleaned.             ")
                        logger.info("                                                                      ")
                        logger.info("######################################################################")
                        logger.info("######################################################################\n")

                    if cleanir_pinholes:
                        cleanir(pinholelist)
                    else:
                        logger.info("######################################################################")
                        logger.info("######################################################################")
                        logger.info("                                                                      ")
                        logger.info("              WARNING in calibrate: pinholes not cleaned.             ")
                        logger.info("                                                                      ")
                        logger.info("######################################################################")
                        logger.info("######################################################################\n")
                    
                    logger.info("##################################################################")
                    logger.info("                                                                  ")
                    logger.info("           STEP 1: Clean calibration frames - COMPLETED           ")
                    logger.info("                                                                  ")
                    logger.info("##################################################################\n")
                    '''
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
                        logger.error("                                                                      ")
                        logger.error("       ERROR in calibrate: unknown array ID. Exiting script.          ")
                        logger.error("                                                                      ")
                        logger.error("######################################################################")
                        logger.error("#####################################################################\n")
                    
                    prepareCalibrations(allcallist, bpmfile, overwrite)

                    logger.info("###################################################################")
                    logger.info("                                                                  ")
                    logger.info("    STEP 2: Preparing all calibration frames (QHflats, IRflats    ")
                    logger.info("            arcs, and pinholes) - COMPLETED                       ")
                    logger.info("                                                                  ")
                    logger.info("##################################################################\n")


                #############################################################################
                ##  STEP 3: Create flat field images for QHflats and IR flats.             ##
                ##          Create the final flat by grouping order 3 of IRflats and       ## 
                ##          orders 4-18 of QHflats.                                        ##
                ##  Output: Flat Field image.                                              ##
                #############################################################################

                elif valindex == 3:
                    if manualMode:
                        a = raw_input("About to enter step 3: create flat field.")

                    makeFlat(QHflatlist, IRflatlist, masterflat, overwrite)

                    logger.info("###################################################################")
                    logger.info("                                                                  ")
                    logger.info("    STEP 3: Create flat field (QHflat, IRflat, and the            ")
                    logger.info("            final flat) - COMPLETED                               ")
                    logger.info("                                                                  ")
                    logger.info("##################################################################\n")

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
                        logger.error("######################################################################")
                        logger.error("######################################################################")
                        logger.error("                                                                      ")
                        logger.error("         ERROR in calibrate: unknown camera. Exiting script.          ")
                        logger.error("                                                                      ")
                        logger.error("######################################################################")
                        logger.error("#####################################################################\n")

                    makeSdistortion(pinholelist, pinhole_coordlist, overwrite)
                
                    logger.info("###################################################################")
                    logger.info("                                                                  ")
                    logger.info("     Step 4: Spatial curvature and spectral distortion (trace     ")
                    logger.info("             the spatial curvature and spectral distortion in     ") 
                    logger.info("             the pinhole flat) - COMPLETED                        ")
                    logger.info("                                                                  ")
                    logger.info("##################################################################\n")

                ############################################################################
                ##  STEP 5: Combine arc frames.                                           ##
                ##          Determine the wavelength solution and create the wavelength   ##
                ##          referenced arc image.                                         ##
                ##  OUTPUT: Combined and wavelength-calibrated arc image.                 ##
                ############################################################################

                elif valindex == 5:
                    if manualMode:
                        a = raw_input("About to enter step 5: wavelength solution.")
                    
                    if '32' in calpath and '10' in calpath:
                        waveCal_coordlist = 'gnirs$data/lowresargon.dat'
                    elif '110' in calpath:
                        waveCal_coordlist = 'gnirs$data/argon.dat'
                    else:
                        logger.error("######################################################################")
                        logger.error("######################################################################")
                        logger.error("                                                                      ")
                        logger.error("         ERROR in calibrate: unknown grating. Exiting script.         ")
                        logger.error("                                                                      ")
                        logger.error("######################################################################")
                        logger.error("#####################################################################\n")

                    makeWaveCal(masterflat, sdistrefimage, waveCal_coordlist, overwrite)
                    
                    logger.info("###################################################################")
                    logger.info("                                                                  ")
                    logger.info("   STEP 5: Wavelength solution (Combine arc frames, determine     ")
                    logger.info("           the wavelength solution, and create the wavelength     ")
                    logger.info("           referenced arc image) - COMPLETED                      ")
                    logger.info("                                                                  ")
                    logger.info("##################################################################\n")

                else:
                    logger.error("ERROR in gnirsBaselineCalibration: step ", valindex, " is not valid.\n")
                    raise SystemExit

                valindex += 1

            logger.info("###########################################################################")
            logger.info("                                                                          ")
            logger.info("  COMPLETE - Calibration reductions completed for                         ")
            logger.info("  %s", calpath)
            logger.info("                                                                          ")
            logger.info("##########################################################################\n")
        
        else:
            logger.warning("GNIRS baseline calibration turned off. Not reducing baseline calibrations in %s", calpath)

    # Return to directory script was begun from.
    os.chdir(path)
        
    return

#####################################################################################
#                                        FUNCTIONS                                  #
#####################################################################################

def prepareCalibrations(allcallist, bpmfile, overwrite):
    """
    Prepare calibration frames for further processing.

    Use NSPREPARE on all calibration frames (QHflats, IRflats, arcs, pinholes, and darks to update the raw data headers 
    and attach the mask definition file (MDF) as a binary table on all files. Note that dark frames will not have an 
    MDF attached by default. Instead, the appropriate MDF is added in NSREDUCE or NSFLAT to match the data being 
    reduced.

    Use NSREDUCE to cut the calibration frame spectra to the size specified by the MDF, placing different spectral
    orders in separate image extensions. Darks will not be used for further reductions after this point.
    """
    logger = log.getLogger('PrepareCalibrations')
    
    # Update all calibration frames with mdf offset value and generate variance and data quality extensions. This code
    # structure checks if iraf output files already exist. If output files exist and overwrite is specified, iraf
    # output is overwritten.
    for image in allcallist:
        image = image.strip()
        if os.path.exists('n'+ image):
            if overwrite:
                logger.warning('Removing old n%s', image)
                os.remove('n'+image)
            else:
                logger.warning("Output exists and -overwrite not set - skipping nsprepare for all calibration frames.")
    iraf.nsprepare(inimages='@QHflats.list,@IRflats.list,@arcs.list,@pinholes.list', rawpath='', outimages='', \
        outprefix='n', bpm=bpmfile, logfile=logger.root.handlers[0].baseFilename, fl_vardq='yes', fl_cravg='no', \
        crradius=0.0, fl_dark_mdf='no', fl_correct='no', fl_saturated='yes', fl_nonlinear='yes', fl_checkwcs='yes', \
        fl_forcewcs='yes', arraytable='gnirs$data/array.fits', configtable='gnirs$data/config.fits', specsec='[*,*]', \
        offsetsec='none', pixscale='0.15', shiftimage='', shiftx='INDEF', shifty='INDEF', obstype='FLAT', \
        fl_inter='no', verbose='yes', mode='al')
    
    # Cut all calibration frames according to the size specified by the MDFs.
    for image in allcallist:
        image = image.strip()
        if os.path.exists('rn'+ image):
            if overwrite:
                logger.warning('Removing old rn%s', image)
                os.remove('rn'+image)
            else:
                logger.info("Output exists and -overwrite not set - skipping nsreduce for all calibration frames.")
    iraf.nsreduce(inimages='n//@QHflats.list,n//@IRflats.list,n//@arcs.list,n//@pinholes.list', outimages='', \
        outprefix='r', fl_cut='yes', section='', fl_corner='yes', fl_process_cut='yes', fl_nsappwave='no', \
        nsappwavedb='gnirs$data/nsappwave.fits', crval='INDEF', cdelt='INDEF', fl_dark='no', darkimage='', \
        fl_save_dark='no', fl_sky='no', skyimages='', skysection='', combtype='median', rejtype='avsigclip', \
        masktype='goodvalue', maskvalue=0.0, scale='none', zero='median', weight='none', statsec='[*,*]', \
        lthreshold='INDEF', hthreshold='INDEF', nlow=1, nhigh=1, nkeep=0, mclip='yes', lsigma=3.0, hsigma=3.0, \
        snoise='0.0', sigscale=0.1, pclip=-0.5, grow=0.0, skyrange='INDEF', nodsize=3.0, fl_flat='no', flatimage='', \
        flatmin=0.0, fl_vardq='yes', logfile=logger.root.handlers[0].baseFilename, verbose='yes', debug='no', \
        force='no', mode='al')

#--------------------------------------------------------------------------------------------------------------------------------#

def makeFlat(QHflatlist, IRflatlist, masterflat, overwrite):
    """
    Make flat field and bad pixel masks.

    Use NSFLAT to generate a normalized flat field for the QHflats and IRflats (for each cross-dispersed order) frm the
    lamp flats. A mask (BPM) will also be generated by thresholding - this can be used to flag bad pixels in other data.

    Use FXCOPY followed by FXINSERT to produce the final flat by grouping order 3 of IRflats and orders 4-18 of QHflats.
    The output from this task is used as the flatfield image for further reduction.
    """
    logger = log.getLogger('MakeFlat')

    # Generate normalized flat field images for QHflats and IRflats.
    QHflat = 'QHflat.fits'
    QHflat_bpm = 'QHflat_bpm.pl'
    if os.path.exists(QHflat):
        if overwrite:
            logger.warning('Removing old %s', QHflat)
            os.remove(QHflat)
        else:
            logger.info('Old %s exists and -overwrite not set - nsflat might fail for QHflats.', QHflat)
    if os.path.exists(QHflat_bpm):
        if overwrite:
            logger.warning('Removing old %s', QHflat_bpm)
            os.remove(QHflat_bpm)
        else:
            logger.info('Old %s exists and -overwrite not set - nsflat might fail for QHflats.', QHflat_bpm)
    iraf.nsflat(lampson='rn//@QHflats.list', darks='', flatfile=QHflat, darkfile='', fl_corner='yes', \
        fl_save_darks='no', flattitle='default', bpmtitle='default', bpmfile=QHflat_bpm, process="fit", \
        statsec='MDF', fitsec='MDF', thr_flo=0.35, thr_fup=4.0, thr_dlo=-20, thr_dup=100, fl_inter='no', \
        fl_range='no', fl_fixbad='yes', fixvalue=1.0, function='spline3', order=5, normstat='midpt', \
        combtype='default', rejtype='ccdclip', masktype='goodvalue', maskvalue=0.0, scale='none', zero='none', \
        weight='none', lthreshold=50.0, hthreshold='INDEF', nlow=1, nhigh=1, nkeep=0, mclip='yes', lsigma=3.0, \
        hsigma=3.0, snoise='0.0', sigscale=0.1, pclip=-0.5, grow=0.0, box_width=20, box_length=1, trace='', \
        traceproc='none', threshold=100.0, aptable='gnirs$data/apertures.fits', database='', apsum=10, tr_step=10, \
        tr_nlost=3, tr_function='legendre', tr_order=5, tr_naver=1, tr_niter=0, tr_lowrej=3.0, tr_highrej=3.0, \
        tr_grow=0.0, ap_lower=-30, ap_upper=30, fl_vardq='yes', logfile=logger.root.handlers[0].baseFilename, \
        verbose='yes', mode='al')

    IRflat = 'IRflat.fits'
    IRflat_bpm = 'IRflat_bpm.pl'
    if os.path.exists(IRflat):
        if overwrite:
            logger.warning('Removing old %s', IRflat)
            os.remove(IRflat)
        else:
            logger.error('Old %s exists and -overwrite not set - nsflat might fail for IRflats.', IRflat)
            raise SystemExit
    if os.path.exists(IRflat_bpm):
        if overwrite:
            logger.warning('Removing old %s', IRflat_bpm)
            os.remove(IRflat_bpm)
        else:
            logger.info('Old %s exists and -overwrite not set - nsflat might fail for IRflats.', QHflat_bpm)        
    iraf.nsflat(lampson='rn//@IRflats.list', darks='', flatfile=IRflat, darkfile='', fl_corner='yes', \
        fl_save_darks='no', flattitle='default', bpmtitle=IRflat_bpm, bpmfile='default', process="fit", statsec='MDF',\
        fitsec='MDF', thr_flo=0.35, thr_fup=1.5, thr_dlo=-20, thr_dup=100, fl_inter='no', fl_range='no', \
        fl_fixbad='yes', fixvalue=1.0, function='spline3', order=10, normstat='midpt', combtype='default', \
        rejtype='ccdclip', masktype='goodvalue', maskvalue=0.0, scale='none', zero='none', weight='none', \
        lthreshold=50.0, hthreshold='INDEF', nlow=1, nhigh=1, nkeep=0, mclip='yes', lsigma=3.0, hsigma=3.0, \
        snoise='0.0', sigscale=0.1, pclip=-0.5, grow=0.0, box_width=20, box_length=1, trace='', traceproc='none', \
        threshold=100.0, aptable='gnirs$data/apertures.fits', database='', apsum=10, tr_step=10, tr_nlost=3, \
        tr_function='legendre', tr_order=5, tr_naver=1, tr_niter=0, tr_lowrej=3.0, tr_highrej=3.0, tr_grow=0.0, \
        ap_lower=-30, ap_upper=30, fl_vardq='yes', logfile=logger.root.handlers[0].baseFilename, verbose='yes', \
        mode='al')

    # Group order 3 of IRflats and orders 4-18 of QHflats to create the final flat field image. Output flat field file 
    # will have name "masterflat".
    if os.path.exists(masterflat):
        if overwrite:
            logger.warning('Removing old %s', masterflat)
            os.remove(masterflat)
            iraf.fxcopy(input=IRflat, output=masterflat, group="0-3", new_file='yes', verbose='no', mode='ql')
            iraf.fxinsert(input=QHflat, output=masterflat+'[3]', groups="4-18", verbose='no', mode='ql')
        else:
            logger.info("Output exists and -overwrite not set - skipping fxcopy and fxinsert to create the masterflat.")
    else:
        iraf.fxcopy(input=IRflat, output=masterflat, group="0-3", new_file='yes', verbose='no', mode='ql')
        iraf.fxinsert(input=QHflat, output=masterflat+'[3]', groups="4-18", verbose='no', mode='ql')

    # Put the name of the masterflat file into a text file called "masterflatfile.txt" to be used by the pipeline later.
    open("masterflatfile.txt", "w").write(masterflat)

#--------------------------------------------------------------------------------------------------------------------------------#

def makeSdistortion(pinholelist, pinhole_coordlist, overwrite):
    """
    Establish Spatial-distortion calibration.

    NSSDIST uses the information in the pinhole calibration images to calibrate the spatial distorion of the GNIRS 
    field. The pinhole frame is a dispersed flat field image with a slit-mask in the field so that the illumination on 
    GNIRS is in a pattern of 6 different holes that are stacked in the y-dimension on the field. Proper alignment of
    the slit across the pattern can be used for spatial rectification of the on-sky science data. The spatial solution
    determined by NSSSDIST is linked to the science data in NSFITCOORDS.
    """
    logger = log.getLogger('MakeSdistortion')

    # Determine the spatial distortion correction. Output: If -overwrite is set, overwrites "sdistrefimagefile.txt"  
    # containing the reference pinhole flat corrected for the spatial distortion and makes changes to files in the
    # /database directory.
    if os.path.exists('/database/idrn*'):
        if overwrite:
            logger.warning('Removing old idrn* files in /database.')
            os.remove('/database/idrn*')
            iraf.nssdist(inimages='rn'+pinholelist[0].strip(), outsuffix='_sdist', pixscale=1.0, dispaxis=1, \
                database='', firstcoord=0.0, coordlist=pinhole_coordlist, aptable='gnirs$data/apertures.fits', \
                fl_inter='no', fl_dbwrite='yes', section='default', nsum=30, ftype='emission', fwidth=10.0, \
                cradius=10.0, threshold=1000.0, minsep=5.0, match=-6.0, function="legendre", order=5, sample='', \
                niterate=3, low_reject=5.0, high_reject=5.0, grow=0.0, refit='yes', step=10, trace='no', nlost=0, \
                aiddebug='', logfile=logger.root.handlers[0].baseFilename, verbose='no', debug='no', force='no', \
                mode='al')
        else:
            logger.info("Output file exists and -overwrite not set - skipping pinhole spatial-distortion calculation.")
            logger.info("calibration with nssdist.\n")
    else:
        iraf.nssdist(inimages='rn'+pinholelist[0].strip(), outsuffix='_sdist', pixscale=1.0, dispaxis=1, database='', \
            firstcoord=0.0, coordlist=pinhole_coordlist, aptable='gnirs$data/apertures.fits', fl_inter='no', \
            fl_dbwrite='yes', section='default', nsum=30, ftype='emission', fwidth=10.0, cradius=10.0, \
            threshold=1000.0, minsep=5.0, match=-6.0, function="legendre", order=5, sample='', niterate=3, \
            low_reject=5.0, high_reject=5.0, grow=0.0, refit='yes', step=10, trace='no', nlost=0, aiddebug='', \
            logfile=logger.root.handlers[0].baseFilename, verbose='no', debug='no', force='no', mode='al')

    # Put the name of the spatially referenced pinhole flat "rn"+pinholeflat into a text file called "sdistrefimagefile' 
    # to be used by the pipeline later. Also associated files are in the "database/" directory.
    open("sdistrefimagefile.txt", "w").write("rn"+pinholelist[0])

#---------------------------------------------------------------------------------------------------------------------------------------#

def makeWaveCal(masterflat, sdistrefimage, waveCal_coordlist, overwrite):
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
    logger = log.getLogger('MakeWaveCal')

    # Store the name of the bad pixel mask in "sflat_bpm".
#    sflat_bpm = open("sflat_bpmfile", "r").readlines()[0].strip()

    combinedarc = 'arc_comb.fits'

    # Combine arc frames. Output combined file will have the name "arc_comb".
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

    if os.path.exists('tftf'+combinedarc):
        if overwrite:
            oldfiles = ['f'+combinedarc, 'tf'+combinedarc, 'wtf'+combinedarc, 'ftf'+combinedarc, 'tftf'+combinedarc, \
                '/database/idwtf*', '/database/*_sdist', '/database/*_lamp']
            for file in oldfiles:
                logger.warning('Removing old %s', file)
                os.remove(file)

            logger.info("Running nsfitcoords and nstransform on the combined arc to spatially rectify it before running nswavelength.")
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
            
            logger.info("Calling nsfitcoords to produce the combined arc with horizontal lines and then nstransform on ")
            logger.info("the transformed arc using the output from nswavelength for the 'lamptrans' parameter.")
            iraf.nsfitcoords(inimages='tf'+combinedarc, outspectra='', outprefix='f', lamptrans='wtf'+combinedarc, \
                sdisttransf='', dispaxis=1, database='', fl_inter='no', fl_align='no', \
                function='chebyshev', lxorder=2, lyorder=4, sxorder=4, syorder=4, pixscale=1.0, \
                logfile=logger.root.handlers[0].baseFilename, verbose='yes', debug='no', force='no', mode='al')
            iraf.nstransform(inimages='ftf'+combinedarc, outspectra='', outprefix='t', dispaxis=1, database='', \
                fl_stripe='no', interptype='poly3', xlog='no', ylog='no', pixscale=1.0, \
                logfile=logger.root.handlers[0].baseFilename, verbose='no', debug='no', mode='al')
        else:
            logger.warning("Output files exist and -overwrite not set - skipping wavelength calibration of arcs.")
    else:
        logger.info("Running nsfitcoords and nstransform on the combined arc to spatially rectify it before running nswavelength.")
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
            nsappwavedb='gnirs$data/nsappwave.fits', fl_median='no ', sdist='', sdorder=4, xorder=2, yorder=2, \
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
    '''
    # Combine arc dark frames, "n"+image+".fits". Output combined file will have the name of the first arc dark file.
    if os.path.exists("gn"+arcdark+".fits"):
        if over:
            iraf.delete("gn"+arcdark+".fits")
            if len(arcdarklist) > 1:
                iraf.gemcombine(listit(arcdarklist,"n"),output="gn"+arcdark, fl_dqpr='yes',fl_vardq='yes',masktype="none",logfile=log)
            else:
                iraf.copy('n'+arcdark+'.fits', 'gn'+arcdark+'.fits')
        else:
            logger.info("\nOutput file exists and -over not set - skipping gemcombine of arcdarks."
    else:
        if len(arcdarklist) > 1:
            iraf.gemcombine(listit(arcdarklist,"n"),output="gn"+arcdark, fl_dqpr='yes',fl_vardq='yes',masktype="none",logfile=log)
        else:
            iraf.copy('n'+arcdark+'.fits', 'gn'+arcdark+'.fits')
    
    # Set interactive mode. Default False for standard configurations (and True for non-standard wavelength configurations ).
    pauseFlag = False

    if band == "K" and central_wavelength == 2.20:
        clist=RUNTIME_DATA_PATH+"k_ar.dat"
        my_thresh = 50.0
    elif band == "J":
        clist=RUNTIME_DATA_PATH+"j_ar.dat"
        my_thresh=100.0
    elif band == "H":
        clist=RUNTIME_DATA_PATH+"h_ar.dat"
        my_thresh=100.0
    elif band == "Z":
        clist="nifs$data/ArXe_Z.dat"
        my_thresh=100.0
    else:
        # Print a warning that the pipeline is being run with non-standard grating.
        logger.info("\n#####################################################################"
        logger.info("#####################################################################"
        logger.info(""
        logger.info("   WARNING in calibrate: found a non-standard (non Z, J, H or K) "
        logger.info("                         wavelength configuration."
        logger.info("                         NSWAVELENGTH will be run interactively."
        logger.info(""
        logger.info("#####################################################################"
        logger.info("#####################################################################\n"

        clist="gnirs$data/argon.dat"
        my_thresh=100.0
        interactive = 'yes'
    '''
#--------------------------------------------------------------------------------------------------------------------------------#

if __name__ == '__main__':
    log.configure('gnirs.log', filelevel='INFO', screenlevel='DEBUG')
    start('gnirs.cfg')
