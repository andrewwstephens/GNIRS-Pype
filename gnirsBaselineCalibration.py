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

from gnirsUtils import datefmt, copyCalibration, copyCalibrationDatabase, replaceNameDatabaseFiles  
## listit, checkLists


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

    log.configure('gnirs.log', filelevel='INFO', screenlevel='INFO')
    logger = log.getLogger('BaselineCalibrations')

    # Store current working directory for later use.
    path = os.getcwd()

#    log = os.getcwd()+'/gnirs.log'

    logger.info('\n#################################################')
    logger.info('#                                                #')
    logger.info('# Start the GNIRS Baseline Calibration Reduction #')
    logger.info('#                                                #')
    logger.info('#################################################\n')

    logger.info("\nParameters read from %s", configfile + "\n")
    
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
    manualMode = bool(config.get('defaults','manualMode'))
    overwrite = bool(config.get('defaults','overwrite'))
    # Read baselineCalibrationReduction specfic config.
    start = bool(config.get('calibrationReductionConfig','Start'))
    stop = bool(config.get('calibrationReductionConfig','Stop'))
    cleanir_arcs = bool(config.get('calibrationReductionConfig','cleanir_arcs'))
    cleanir_QHflats = bool(config.get('calibrationReductionConfig','cleanir_QHflats'))
    cleanir_IRflats = bool(config.get('calibrationReductionConfig','cleanir_IRflats'))
    cleanir_pinholes = bool(config.get('calibrationReductionConfig','cleanir_pinholes'))
    calibrationDirectories = bool(config.options('CalibrationDirectories'))

    print("done for now!")

    ################################################################################
    # Define Variables, Reduction Lists AND identify/run number of reduction steps #
    ################################################################################
'''
    # Loop over the Calibrations directories and reduce the day calibrations in each one.
    for calpath in calibrationDirectoryList:
        os.chdir(calpath)
        pwdDir = os.getcwd()+"/"
        iraffunctions.chdir(pwdDir)

        # However, don't do the reduction for a Calibrations_"configuration" directory without associated telluric or
        # science data. Check that a Calibrations_"configuration" directory exists at the same level as the 
        # Calibrations_"configuration" directory. If not, skip the reduction of calibrations in that 
        # Calibrations_"configuration" directory. "configuration" should be the last directory of calpath.
        configuration = os.path.split(calpath)[-1]
        if not os.path.exists("../"+configuration):

            logger.info("\n######################################################################")
            logger.info("                                                                       ")
            logger.info("    No configuration directory (including science or telluric data)    ")
            logger.info("    found for  ", calpath)
            logger.info("    Skipping reduction of calibrations in that directory.              ")
            logger.info("                                                                       ")
            logger.info("######################################################################\n")

            continue

        # Create lists of each type of calibration from textfiles in Calibrations directory.
        QHflatlist = open('QHflatlist', "r").readlines()
        IRflatlist = open("IRflatlist", "r").readlines()
        arclist = open("arclist", "r").readlines()
        pinholelist = open("pinholelist", "r").readlines()
        
        masterflat = "masterflat.fits"  ## name of the final flat having order 3 of IRflats and orders 4-18 of 
                                        ## QHflats grouped together
        sdistrefimage = "sdistrefiamge"  ## name of the first flat in pinholelist on which S-distortion has been
                                          ## applied and will be used as the reference image for S-distortion 
                                          ## correction of all other images taken on that same day
#        combinedarc = "arc_comb"  ## name of the arc image combined using gscombine
        
        # Check start and stop values for reduction steps. Ask user for a correction if input is not valid.
        valindex = start
        while valindex > stop  or valindex < 1 or stop > 4:
            logger.info("\n#####################################################################")
            logger.info("######################################################################")
            logger.info("                                                                      ")
            logger.info("    WARNING in calibrate: invalid start/stop values of calibration    ")
            logger.info("                           reduction steps.                           ")
            logger.info("                                                                      ")
            logger.info("######################################################################")
            logger.info("#####################################################################\n")

            valindex = int(raw_input("\nPlease enter a valid start value (1 to 4, default 1): "))
            stop = int(raw_input("\nPlease enter a valid stop value (1 to 4, default 4): "))

        # Print the current directory of calibrations being processed.
        logger.info("\n#########################################################################")
        logger.info("                                                                          ")
        logger.info("  Currently working on calibrations in ", calpath)
        logger.info("                                                                          ")
        logger.info("#########################################################################\n")

        while valindex <= stop:
            
            #############################################################################
            ##  STEP 1: Clean calibration frames.                                      ##
            ##  Output: Cleaned QHflats, IRflats, arcs, and pinholes.                  ##
            #############################################################################

            if valindex == 1:
                if manualMode:
                    a = raw_input("About to enter step 1: clean calibration frames.")
                
                if cleanir_arcs:
                    cleanir(arclist)
                else:
                    logger.info("\n#####################################################################")
                    logger.info("######################################################################")
                    logger.info("                                                                      ")
                    logger.info("               WARNING in calibrate: arcs not cleaned.                ")
                    logger.info("                                                                      ")
                    logger.info("######################################################################")
                    logger.info("#####################################################################\n")
                
                if cleanir_QHflats:
                    cleanir(QHflatlist)
                else:
                    logger.info("\n#####################################################################")
                    logger.info("######################################################################")
                    logger.info("                                                                      ")
                    logger.info("              WARNING in calibrate: QHflats not cleaned.              ")
                    logger.info("                                                                      ")
                    logger.info("######################################################################")
                    logger.info("#####################################################################\n")
                
                if cleanir_IRflats:
                    cleanir(IRflatlist)
                else:
                    logger.info("\n#####################################################################")
                    logger.info("######################################################################")
                    logger.info("                                                                      ")
                    logger.info("               WARNING in calibrate: IRflats not cleaned.             ")
                    logger.info("                                                                      ")
                    logger.info("######################################################################")
                    logger.info("#####################################################################\n")

                if cleanir_pinholes:
                    cleanir(pinholelist)
                else:
                    logger.info("\n#####################################################################")
                    logger.info("######################################################################")
                    logger.info("                                                                      ")
                    logger.info("              WARNING in calibrate: pinholes not cleaned.             ")
                    logger.info("                                                                      ")
                    logger.info("######################################################################")
                    logger.info("#####################################################################\n")
        
                logger.info("\n#################################################################")
                logger.info("                                                                  ")
                logger.info("           STEP 1: Clean calibration frames - COMPLETED           ")
                logger.info("                                                                  ")
                logger.info("#################################################################\n")

            #############################################################################
            ##  STEP 2: Create flat field images and bad pixel masks for QHflats and   ##
            ##          IR flats.                                                      ##
            ##          Create the final flat by grouping order 3 of IRflats and       ## 
            ##          orders 4-18 of QHflats.                                        ##
            ##  Output: Flat Field image.                                              ##
            #############################################################################

            elif valindex == 2:
                if manualMode:
                    a = raw_input("About to enter step 2: create flat field.")
                makeFlat(QHflatlist, IRflatlist, masterflat, configuration, overwrite, log)
                logger.info("\n#################################################################")
                logger.info("                                                                  ")
                logger.info("    STEP 2: Create flat field (QHflat, IRflat, bad pixel masks    ")
                logger.info("            and the final flat) - COMPLETED                       ")
                logger.info("                                                                  ")
                logger.info("#################################################################\n")

            ########################################################################################
            ##  Step 3: Trace the spatial curvature and spectral distortion in the pinhole flat.  ##
            ##  Output: Pinhole flat for which the spatial curvature and S-distortion has been    ##
            ##          determined.                                                               ##
            ########################################################################################

            elif valindex == 3:
                if manualMode:
                    a = raw_input("About to enter step 3: spatial curvature and spectral distortion.")
                makeSdistortion(pinholelist, sdistrefimage, configuration, overwrite, log)
                logger.info("\n#################################################################")
                logger.info("                                                                  ")
                logger.info("     Step 3: Spatial curvature and spectral distortion (trace     ")
                logger.info("             the spatial curvature and spectral distortion in     ") 
                logger.info("             the pinhole flat) - COMPLETED                        ")
                logger.info("                                                                  ")
                logger.info("#################################################################\n")

            ############################################################################
            ##  STEP 4: Combine arc frames.                                           ##
            ##          Determine the wavelength solution and create the wavelength   ##
            ##          referenced arc image.                                         ##
            ##  OUTPUT: Combined and wavelength-calibrated arc image.                 ##
            ############################################################################

            elif valindex == 4:
                if manualMode:
                    a = raw_input("About to enter step 4: wavelength solution.")
                makeWaveCal(arclist, masterflat, sdistrefimage, configuration, overwrite, log, path)
                logger.info("\n#################################################################")
                logger.info("                                                                  ")
                logger.info("   STEP 4: Wavelength solution (Combine arc frames, determine     ")
                logger.info("           the wavelength solution, and create the wavelength     ")
                logger.info("           referenced arc image) - COMPLETED                      ")
                logger.info("                                                                  ")
                logger.info("#################################################################\n")

            else:
                logger.info("\nERROR in nifs_baseline_calibration: step ", valindex, " is not valid.\n")
                raise SystemExit

            valindex += 1

        logger.info("\n#########################################################################")
        logger.info("                                                                          ")
        logger.info("  COMPLETE - Calibration reductions completed for                         ")
        logger.info("  ", calpath)
        logger.info("                                                                          ")
        logger.info("#########################################################################\n")


    # Return to directory script was begun from.
    os.chdir(path)
    return

#####################################################################################
#                                        FUNCTIONS                                  #
#####################################################################################

def makeFlat(QHflatlist, IRflatlist, masterflat, configuration, overwrite, log):
    """
    Make flat field and bad pixel masks.

    Use NSPREPARE on the QHflats and IRflats to update the raw data headers and attach the mask definition file (MDF)
    as a binary table on all files. Note that dark frames will not have an MDF attached by default. Instead, the 
    appropriate MDF is added in NSREDUCE or NSFLAT to match the data being reduced.

    Use NSREDUCE to cut the calibration (flat/arc) spectra to the size specified by the MDF, placing different spectral
    orders in separate image extensions.

    Use NSFLAT to generate a normalized flat field for the QHflats and IRflats (for each cross-dispersed order) frm the
    lamp flats. A mask (BPM) will also be generated by thresholding - this can be used to flag bad pixels in other data.

    Use FXCOPY followed by FXINSERT to produce the final flat by grouping order 3 of IRflats and orders 4-18 of QHflats.
    The output from this task is used as the flatfield image for further reduction.
    """

    # Update QHlamps and IRlamps on flat frames with mdf offset value and generate variance and data quality extensions.
    # This code structure checks if iraf output files already exist. If output files exist and overwrite is specified, 
    # iraf output is overwritten.
    for image in QHflatlist:
        image = image.strip()
        if os.path.exists('n'+ image):
            if overwrite:
                os.remove('n'+image)
            else:
                logger.info("Output exists and -overwrite not set - skipping nsprepare for QHflats.")
    iraf.nsprepare(inimages='@QHflatlist', rawpath='.', outimages='',outprefix='n', bpm=bpm, logfile=log, \
        fl_vardq='yes', fl_cravg='no', crradius=0.0, fl_dark_mdf='no', fl_correct='no', fl_saturated='yes', \
            fl_nonlinear='yes', fl_checkwcs='yes', fl_forcewcs='yes', arraytable='gnirs$data/array.fits', \
                configtable='gnirs$data/config.fits', specsec='[*,*]', offsetsec='none', pixscale='0.15', \
                    shiftimage='', shiftx='INDEF', shifty='INDEF', obstype='FLAT', fl_inter='no', verbose='yes', \
                            mode='al')
    
    # Write the name of the first nsprepare'd flat in QHflats (to be used as the reference image to shift the MDFs for 
    # the rest of the images taken on that same day) into a text file called "mdfshiftrefimagefile" to be used by the 
    # pipeline later.
    open("mdfshiftrefimagefile", "w").write('n'+QHflatlist[0].strip()) 
    copyCalibration('n'+QHflatlist[0].strip(), 'mdfshiftrefimagefile.fits', configuration, overwrite)

    for image in IRflatlist:
        image = image.strip()
        if os.path.exists('n'+image):
            if overwrite:
                os.remove('n'+image)
            else:
                logger.info("Output exists and -overwrite not set - skipping nsprepare for IRflats.")
    iraf.nsprepare(inimages='@IRflatlist', rawpath='.', outimages='',outprefix='n', bpm=bpm, logfile=log, \
        fl_vardq='yes', fl_cravg='no', crradius=0.0, fl_dark_mdf='no', fl_correct='no', fl_saturated='yes', \
            fl_nonlinear='yes', fl_checkwcs='yes', fl_forcewcs='yes', arraytable='gnirs$data/array.fits', \
                configtable='gnirs$data/config.fits', specsec='[*,*]', offsetsec='none', pixscale='0.15', \
                    shiftimage='./mdfshiftrefimagefile.fits', shiftx='0.', shifty='0.', obstype='FLAT', fl_inter='no', \
                            verbose='yes', mode='al')

    # Cut the QHlamps and IRlamps on flat frames according to the size specified by the MDFs.
    for image in QHflatlist:
        image = image.strip()
        if os.path.exists('rn'+ image):
            if overwrite:
                os.remove('rn'+image)
            else:
                logger.info("Output exists and -overwrite not set - skipping nsreduce for QHflats.")
    iraf.nsreduce(inimages='n//@QHflatlist', rawpath='.', outimages='',outprefix='r', fl_cut='yes', section='', \
        fl_corner='yes', fl_process_cut='yes', fl_nsappwave='no', nsappwavedb='gnirs$data/nsappwave.fits', \
            crval='INDEF', cdelt='INDEF', fl_dark='no', darkimage='', fl_save_dark='no', fl_sky='no', skyimages='',\
                skysection='', combtype='median', rejtype='avsigclip', masktype='goodvalue', maskvalue=0.0, \
                    scale='none', zero='median', weight='none', statsec='[*,*]', lthreshold='INDEF', \
                        hthreshold='INDEF', nlow=1, nhigh=1, nkeep=0, mclip='yes', lsigma=3.0, hsigma=3.0, \
                            snoise='0.0', sigscale=0.1, pclip=-0.5, grow=0.0, skyrange='INDEF', nodsize=3.0, \
                                fl_flat='no', flatimage='', flatmin=0.0, fl_vardq='yes', logfile=log, verbose='yes', \
                                    debug='no', force='no', mode='al')
    
    for image in IRflatlist:
        image = image.strip()
        if os.path.exists('rn'+ image):
            if overwrite:
                os.remove('rn'+image)
            else:
                logger.info("Output exists and -overwrite not set - skipping nsreduce for IRflats.")
    iraf.nsreduce(inimages='n//@IRflatlist', rawpath='.', outimages='',outprefix='r', fl_cut='yes', section='', \
        fl_corner='yes', fl_process_cut='yes', fl_nsappwave='no', nsappwavedb='gnirs$data/nsappwave.fits', \
            crval='INDEF', cdelt='INDEF', fl_dark='no', darkimage='', fl_save_dark='no', fl_sky='no', skyimages='',\
                skysection='', combtype='median', rejtype='avsigclip', masktype='goodvalue', maskvalue=0.0, \
                    scale='none', zero='median', weight='none', statsec='[*,*]', lthreshold='INDEF', \
                        hthreshold='INDEF', nlow=1, nhigh=1, nkeep=0, mclip='yes', lsigma=3.0, hsigma=3.0, \
                            snoise='0.0', sigscale=0.1, pclip=-0.5, grow=0.0, skyrange='INDEF', nodsize=3.0, \
                                fl_flat='no', flatimage='', flatmin=0.0, fl_vardq='yes', logfile=log, verbose='yes', \
                                    debug='no', force='no', mode='al')

    # Generate normalized flat field images for QHflats and IRflats.
    QHflat = 'QHflat.fits'
    for image in QHflatlist:
        image = image.strip()
        if os.path.exists(QHflat):
            if overwrite:
                os.remove(QHflat)
            else:
                logger.info("Output exists and -overwrite not set - skipping nsflat for QHflats.")
    iraf.nsflat(lampson='rn//@QHflatlist', darks='', flatfile=QHflat, darkfile='', fl_corner='yes', \
        fl_save_darks='no', flattitle='default', bpmtitle='default', bpmfile='default', process="fit", statsec='MDF', \
            fitsec='MDF', thr_flo=0.35, thr_fup=4.0, thr_dlo=-20, thr_dup=100, fl_inter='no', fl_range='no', \
                fl_fixbad='yes', fixvalue=1.0, function='spline3', order=5, normstat='midpt', combtype='default', \
                    rejtype='ccdclip', masktype='goodvalue', maskvalue=0.0, scale='none', zero='none', weight='none', \
                        lthreshold=50.0, hthreshold='INDEF', nlow=1, nhigh=1, nkeep=0, mclip='yes', lsigma=3.0, \
                            hsigma=3.0, snoise='0.0', sigscale=0.1, pclip=-0.5, grow=0.0, box_width=20, box_length=1, \
                                trace='', traceproc='none', threshold=100.0, aptable='gnirs$data/apertures.fits', \
                                    database='', apsum=10, tr_step=10, tr_nlost=3, tr_function='legendre', tr_order=5,\
                                        tr_naver=1, tr_niter=0, tr_lowrej=3.0, tr_highrej=3.0, tr_grow=0.0, \
                                            ap_lower=-30, ap_upper=30, fl_vardq='yes', logfile=log, verbose='yes', \
                                                mode='al')

    IRflat = 'IRflat.fits'
    for image in IRflatlist:
        image = image.strip()
        if os.path.exists(IRflat):
            if overwrite:
                os.remove(IRflat)
            else:
                logger.info("Output exists and -overwrite not set - skipping nsflat for QHflats.")
    iraf.nsflat(lampson='rn//@IRflatlist', darks='', flatfile=IRflat, darkfile='', fl_corner='yes', \
        fl_save_darks='no', flattitle='default', bpmtitle='default', bpmfile='default', process="fit", statsec='MDF', \
            fitsec='MDF', thr_flo=0.35, thr_fup=1.5, thr_dlo=-20, thr_dup=100, fl_inter='no', fl_range='no', \
                fl_fixbad='yes', fixvalue=1.0, function='spline3', order=10, normstat='midpt', combtype='default', \
                    rejtype='ccdclip', masktype='goodvalue', maskvalue=0.0, scale='none', zero='none', weight='none', \
                        lthreshold=50.0, hthreshold='INDEF', nlow=1, nhigh=1, nkeep=0, mclip='yes', lsigma=3.0, \
                            hsigma=3.0, snoise='0.0', sigscale=0.1, pclip=-0.5, grow=0.0, box_width=20, box_length=1, \
                                trace='', traceproc='none', threshold=100.0, aptable='gnirs$data/apertures.fits', \
                                    database='', apsum=10, tr_step=10, tr_nlost=3, tr_function='legendre', tr_order=5,\
                                        tr_naver=1, tr_niter=0, tr_lowrej=3.0, tr_highrej=3.0, tr_grow=0.0, \
                                            ap_lower=-30, ap_upper=30, fl_vardq='yes', logfile=log, verbose='yes', \
                                                mode='al')

    # Group order 3 of IRflats and orders 4-18 of QHflats to create the final flat field image. Output flat field file 
    # will have name "masterflat".
    if os.path.exists(masterflat):
        if overwrite:
            os.remove(masterflat)
            iraf.fxcopy(input=IRflat, output=masterflat, group="0-3", new_file='yes', verbose='no', mode='ql')
            iraf.fxinsert(input=QHflat, output=masterflat+'[3]', groups="4-18", verbose='no', mode='ql')
        else:
            logger.info("Output exists and -overwrite not set - skipping fxcopy and fxinsert to create the masterflat.")
    else:
        iraf.fxcopy(input=IRflat, output=masterflat, group="0-3", new_file='yes', verbose='no', mode='ql')
        iraf.fxinsert(input=QHflat, output=masterflat+'[3]', groups="4-18", verbose='no', mode='ql')

    # Put the name of the masterflat file into a text file called "masterflatfile" to be used by the pipeline later.
    open("masterflatfile", "w").write(masterflat)
    copyCalibration(masterflat, 'masterflatfile.fits', configuration, overwrite)

#--------------------------------------------------------------------------------------------------------------------------------#

def makeSdistortion(pinholelist, sdistrefimage, configuration, overwrite, log):
    """
    Establish Spatial-distortion calibration.

    NSSDIST uses the information in the pinhole calibration images to calibrate the spatial distorion of the GNIRS 
    field. The pinhole frame is a dispersed flat field image with a slit-mask in the field so that the illumination on 
    GNIRS is in a pattern of 6 different holes that are stacked in the y-dimension on the field. Proper alignment of
    the slit across the pattern can be used for spatial rectification of the on-sky science data. The spatial solution
    determined by NSSSDIST is linked to the science data in NSFITCOORDS.
    """

    # Update pinhole flat frames with offset value and generate variance and data quality extensions.
    for image in pinholelist:
        image = image.strip()    
        if os.path.exists("n"+image):
            if overwrite:
 #               iraf.delete("n"+image+'.fits')
                os.remove("n"+image)
            else:
                logger.info("\nOutput file exists and -overwrite not set - skipping prepare of pinholes.")
    iraf.nsprepare(inimages='@pinholelist', rawpath='.', outimages='',outprefix='n', bpm=bpm, logfile=log, \
        fl_vardq='yes', fl_cravg='no', crradius=0.0, fl_dark_mdf='no', fl_correct='no', fl_saturated='yes', \
            fl_nonlinear='yes', fl_checkwcs='yes', fl_forcewcs='yes', arraytable='gnirs$data/array.fits', \
                configtable='gnirs$data/config.fits', specsec='[*,*]', offsetsec='none', pixscale='0.15', \
                    shiftimage='./mdfshiftrefimagefile.fits', shiftx='0.', shifty='0.', obstype='FLAT', fl_inter='no',\
                            verbose='yes', mode='al')
    
    # Cut the pinholes according to the size specified by the MDFs.
    for image in pinholelist:
        image = image.strip()
        if os.path.exists('rn'+ image):
            if overwrite:
                os.remove('rn'+image)
            else:
                logger.info("Output exists and -overwrite not set - skipping nsreduce for pinholes.")
    iraf.nsreduce(inimages='n//@pinholelist', rawpath='.', outimages='',outprefix='r', fl_cut='yes', section='', \
        fl_corner='yes', fl_process_cut='yes', fl_nsappwave='no', nsappwavedb='gnirs$data/nsappwave.fits', \
            crval='INDEF', cdelt='INDEF', fl_dark='no', darkimage='', fl_save_dark='no', fl_sky='no', skyimages='',\
                skysection='', combtype='median', rejtype='avsigclip', masktype='goodvalue', maskvalue=0.0, \
                    scale='none', zero='median', weight='none', statsec='[*,*]', lthreshold='INDEF', \
                        hthreshold='INDEF', nlow=1, nhigh=1, nkeep=0, mclip='yes', lsigma=3.0, hsigma=3.0, \
                            snoise='0.0', sigscale=0.1, pclip=-0.5, grow=0.0, skyrange='INDEF', nodsize=3.0, \
                                fl_flat='no', flatimage='', flatmin=0.0, fl_vardq='yes', logfile=log, verbose='yes', \
                                    debug='no', force='no', mode='al')
    
    # Determine the spatial distortion correction. Output: If -overwrite is set, overwrites "sdistrefimage" text file 
    # containing the reference pinhole flat corrected for the spatial distortion and makes changes to files in 
    # /database directory.
    if os.path.exists(sdistrefimage):
        if overwrite:
#            iraf.delete("ronchifile")
            os.remove(sdistrefimage)
            iraf.nssdist(inimages='rn//@pinholelist', outsuffix='_sdist', pixscale=1.0, dispaxis=1, database='', \
                firstcoord=0.0, coordlist=pinhole_file, aptable='gnirs$data/apertures.fits', fl_inter='no', \
                    fl_dbwrite='yes', section='default', nsum=30, ftype='emission', fwidth=10.0, cradius=10.0, \
                        threshold=1000.0, minsep=5.0, match=-6.0, function="legendre", order=5, sample='', niterate=3,\
                            low_reject=5.0, high_reject=5.0, grow=0.0, refit='yes', step=10, trace='no', nlost=0, \
                                aiddebug='', logfile=log, verbose='no', debug='no', force='no', mode='al')
        else:
            logger.info("\nOutput file exists and -overwrite not set - not performing pinhole spatial-distortion")
            logger.info("calibration with nssdist.\n")
    else:
        iraf.nssdist(inimages='rn//@pinholelist', outsuffix='_sdist', pixscale=1.0, dispaxis=1, database='', \
                firstcoord=0.0, coordlist=pinhole_file, aptable='gnirs$data/apertures.fits', fl_inter='no', \
                    fl_dbwrite='yes', section='default', nsum=30, ftype='emission', fwidth=10.0, cradius=10.0, \
                        threshold=1000.0, minsep=5.0, match=-6.0, function="legendre", order=5, sample='', niterate=3,\
                            low_reject=5.0, high_reject=5.0, grow=0.0, refit='yes', step=10, trace='no', nlost=0, \
                                aiddebug='', logfile=log, verbose='no', debug='no', force='no', mode='al')

    # Put the name of the spatially referenced pinhole flat "rn"+pinholeflat into a text file called "sdistrefimagefile' 
    # to be used by the pipeline later. Also associated files are in the "database/" directory.
    open("sdistrefimagefile", "w").write("rn"+pinholelist[0])
    # Copy to relevant scienceObservation/calibrations/ directories
    for item in sorted(glob.glob('database/idrn*')):
        replaceNameDatabaseFiles(item, "rn"+pinholelist[0], 'sdistrefimagefile')
    copyCalibration("rn"+pinholelist[0], 'sdistrefimagefile.fits', configuration, overwrite)
    copyCalibrationDatabase("idrn", configuration, "sdistrefimagefile", overwrite)

#---------------------------------------------------------------------------------------------------------------------------------------#

def makeWaveCal(arclist, configuration, log, overwrite, path):
    """
    Determine the wavelength solution of the combined arc.

    If the user wishes to change the coordinate file to a different one, they need only to change the "clist" variable 
    to their line list in the coordli= parameter in the nswavelength call.

    Uses NSWAVELENGTH to calibrate arc data (after cutting and optionally applying a flatfield with NSREDUCE in a 
    previous step).

    DATA REDUCTION HINT -
    For the nswavelength call, the different wavelength settings use different values for some of the parameters. For 
    optimal auto results, use:

    ????

    Note that better RMS fits can be obtained by running the wavelength calibration interactively and identifying all 
    of the lines manually. Tedious, but will give more accurate results than the automatic mode (i.e., fl_inter-). Use 
    fl_inter+ for manual mode.
    """

    # Store the name of the bad pixel mask in "sflat_bpm".
    sflat_bpm = open("sflat_bpmfile", "r").readlines()[0].strip()
    # Store the name of the final flat field frame in "flat".
    flat = open("flatfile", "r").readlines()[0].strip()

    # Update arc images with offset value and generate variance and data
    # quality extensions. Results in "n"+image+".fits"
    for image in arclist:
        image = str(image).strip()
        if os.path.exists("n"+image+".fits"):
            if over:
                iraf.delete("n"+image+".fits")
            else:
                logger.info("\nOutput file exists and -over not set - skipping nfprepare of arcs."
                continue
        iraf.nfprepare(image, rawpath=".", shiftimage=shiftima,bpm=sflat_bpm,\
                       fl_vardq="yes",fl_corr='no',fl_nonl='no',logfile=log)

    # Check that output files for all arc images exists from nfprepare; if output does not
    # exist remove corresponding arc images from arclist.
    arclist = checkLists(arclist, '.', 'n', '.fits')

    # Update arc dark frames with mdf offset value and generate variance and data
    # quality extensions. Results in "n"+image+".fits"
    for image in arcdarklist:
        image = str(image).strip()
        if os.path.exists("n"+image+".fits"):
            if over:
                iraf.delete("n"+image+".fits")
            else:
                logger.info("\nOutput file exists and -over not set - skipping nfprepare of arcdarks."
                continue
        iraf.nfprepare(image, rawpath=".", shiftimage=shiftima, bpm=sflat_bpm, \
                       fl_vardq='yes',fl_corr='no',fl_nonl='no',logfile=log)

    # Check that output files for all arc images exists from nfprepare; if output does not
    # exist remove corresponding arc images from arclist.
    arcdarklist = checkLists(arcdarklist, '.', 'n', '.fits')

    # Combine arc frames, "n"+image+".fits". Output combined file will have the name of the first arc file.
    if os.path.exists("gn"+arc+".fits"):
        if over:
            iraf.delete("gn"+arc+".fits")
            if len(arclist) > 1:
                iraf.gemcombine(listit(arclist,"n"),output="gn"+arc, fl_dqpr='yes',fl_vardq='yes',masktype="none",logfile=log)
            else:
                iraf.copy('n'+arc+'.fits', 'gn'+arc+'.fits')
        else:
            logger.info("\nOutput file exists and -over not set - skipping gemcombine of arcs."
    else:
        if len(arclist) > 1:
            iraf.gemcombine(listit(arclist,"n"),output="gn"+arc, fl_dqpr='yes',fl_vardq='yes',masktype="none",logfile=log)
        else:
            iraf.copy('n'+arc+'.fits', 'gn'+arc+'.fits')

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

    # Put the name of the combined and prepared arc dark frame "gn"+arcdark into a text
    # file called arcdarkfile to be used by the pipeline later.
    open("arcdarkfile", "w").write("gn"+arcdark)

    # NSREDUCE on arc images "gn"+arc+".fits" to extract the slices and apply an approximate
    # wavelength calibration. Results in "rgn"+image+".fits"
    if os.path.exists("rgn"+arc+".fits"):
        if over:
            iraf.delete("rgn"+arc+".fits")
            iraf.nsreduce("gn"+arc, darki="gn"+arcdark, flatimage=flat, \
                          fl_vardq="no", fl_cut="yes", fl_nsappw="yes", fl_sky="no", fl_dark="yes", fl_flat="yes", \
                          logfile=log)
        else:
            logger.info("\nOutput file exists and -over not set - skipping apply_flat_arc."
    else:
        iraf.nsreduce("gn"+arc, darki="gn"+arcdark, flatimage=flat, \
                      fl_vardq="no", fl_cut="yes", fl_nsappw="yes", fl_sky="no", fl_dark="yes", fl_flat="yes", \
                      logfile=log)
    #fl_dark = "no"
    #if arcdark != "":
    #    fl_dark = "yes"
    #hdulist = astropy.io.fits.open(arc+'.fits')
    #if 'K_Long' in hdulist[0].header['GRATING']:
    #    iraf.nsreduce("gn"+arc, darki=arcdark, fl_cut="yes", fl_nsappw="yes", crval = 23000., fl_dark="yes", fl_sky="no", fl_flat="yes", flatimage=flat, fl_vardq="no",logfile=log)

    # Put the name of the combined and prepared arc dark frame "gn"+arcdark into a text
    # file called arcdarkfile to be used by the pipeline later.
    open("arcdarkfile", "w").write("gn"+arcdark)

    # Determine the wavelength setting.
    hdulist = astropy.io.fits.open("rgn"+arc+".fits")
    band = hdulist[0].header['GRATING'][0:1]
    central_wavelength = float(hdulist[0].header['GRATWAVE'])

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

    # TODO(nat): I don't like this nesting at all
    if not pauseFlag:
        # Establish wavelength calibration for arclamp spectra. Output: A series of
        # files in a "database/" directory containing the wavelength solutions of
        # each slice and a reduced arc frame "wrgn"+ARC+".fits".
        if os.path.exists("wrgn"+arc+".fits"):
            if over:
                iraf.delete("wrgn"+arc+".fits")
                iraf.nswavelength("rgn"+arc, coordli=clist, nsum=10, thresho=my_thresh, \
                                  trace='yes', fwidth=2.0, match=-6,cradius=8.0,fl_inter='no',nfound=10,nlost=10, \
                                  logfile=log)
            else:
                logger.info("\nOutput file exists and -over not set - ",\
                "not determining wavelength solution and recreating the wavelength reference arc.\n"
        else:
            iraf.nswavelength("rgn"+arc, coordli=clist, nsum=10, thresho=my_thresh, \
                              trace='yes', fwidth=2.0, match=-6,cradius=8.0,fl_inter='no',nfound=10,nlost=10, \
                              logfile=log)
    else:
        a = raw_input("For now, interactive Z or non-standard wavelength calibrations are unsupported. " + \
        "Bugs running IRAF tasks interactively from python mean iraf.nswavelength cannot be activated automatically. " + \
        "Therefore please run iraf.nswavelength() interactively from Pyraf to do a wavelength calibration by hand.")

    # Copy to relevant science observation/calibrations/ directories
    for item in glob.glob('database/idwrgn*'):
        replaceNameDatabaseFiles(item, "wrgn"+arc, 'finalArc')
    copyCalibration("wrgn"+arc+".fits", 'finalArc.fits', grating, over)
    copyCalibrationDatabase("idwrgn", grating, "finalArc", over)
'''
#--------------------------------------------------------------------------------------------------------------------------------#

if __name__ == '__main__':
    log.configure('gnirs.log', filelevel='INFO', screenlevel='DEBUG')
    start('gnirs.cfg')
