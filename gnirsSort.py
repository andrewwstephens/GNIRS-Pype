#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ConfigParser
import datetime
import log
import os
import pkg_resources
import re
import shutil
import sys

# RECIPES_PATH = pkg_resources.resource_filename('nifty', 'recipes/')
RECIPES_PATH = 'recipes/'
# RUNTIME_DATA_PATH = pkg_resources.resource_filename('nifty', 'runtimeData/')
RUNTIME_DATA_PATH = 'runtimeData/'


def start(info, config):
    """
    Sort and copy the GNIRS raw data to subdirectories: target / date / config

    info (dictionary): Dictionary of relevant header keywords read from all FITS files in rawPath.
    config (ConfigParser object): contains the contents of the user-supplied configuration file.
    """

    logger = log.getLogger('Sort')
    logger.info('####################################')
    logger.info('#                                  #')
    logger.info('#  Start GNIRS sorting and copying #')
    logger.info('#                                  #')
    logger.info('####################################')

    path = os.getcwd()

    rawPath = config.get('dataConfig','rawPath')
    # manualMode = config.getboolean('defaults','manualMode')
    # overwrite = config.getboolean('defaults','overwrite')
    # sortTellurics = config.getboolean('sortConfig','sortTellurics')
    # telluricTimeThreshold = config.getint('sortConfig','telluricTimeThreshold')

    allfilelist = []    # List of Tellurics, science frames, and acquisitions (but not calibrations)
    QHflatlist = []     # List of QH lamps on flat frames.
    IRflatlist = []     # List of IR lamps on flat frames.
    pinholelist = []    # List of pinhole flat frames.
    arclist = []        # List of arc frames.
    darklist = []       # List of dark frames.
    objectDateConfigurationList = [] # List of unique objects, dates, and configurations
    skyFrameList = []   # List of sky frames.
    obsidDateList = []  # List of (date, obsid) pairs.
    sciDateList = []    # List of unique dates by science (including sky) frames.
    sciImageList = []   # List of science observation directories.
    
    path = os.getcwd()  # Store current working directory for later use.

    # If files were copied from Gemini Internal network raw files directory will be path+"/rawPath".
    if rawPath:
        rawPath = rawPath
    else:
        rawPath = path + '/rawPath'

    logger.info("Path to raw file directory is: %s", rawPath)
    logger.info("Creating lists of each type of file.")

    for filename in sorted(info.keys()):

        if info[filename]['OBSTYPE'] == 'OBJECT' and info[filename]['OBSCLASS'] in ['science', 'acq', 'acqCal', 'partnerCal']:

            if info[filename]['OBSCLASS'] == 'acq':
                # Use the P and Q offsets of the last acquisition image to calculate the absolute P and Q offsets.
                # The absolute P and Q offsets are then used to determine if the science target is in the slit.
                # This assumes that the acq and science frames are in order and uninterrupted!
                acq_poff = info[filename]['POFFSET']
                acq_qoff = info[filename]['QOFFSET']

            elif info[filename]['OBSCLASS'] == 'science':
                logger.debug('Adding %s to the science list', filename)
                sciImageList.append(filename)

                odc = [info[filename]['OBJECT'], info[filename]['DATE-OBS'], info[filename]['CONFIG'], info[filename]['OBSID']]
                if odc not in objectDateConfigurationList:
                    objectDateConfigurationList.append(odc)

                if info[filename]['DATE-OBS'] not in sciDateList:
                    sciDateList.append(info[filename]['DATE-OBS'])

                if 'SC' in info[filename]['DECKER'] and 'XD' in info[filename]['DECKER']:
                    slitlength = 7.
                elif 'LC' in info[filename]['DECKER'] and 'XD' in info[filename]['DECKER']:
                    slitlength = 5.
                else:
                    slitlength = None
                slitwidth = float(info[filename]['SLIT'][:info[filename]['SLIT'].find('arcsec')])
                logger.debug('Slit size: %.1f x %.1f arcsec', slitwidth, slitlength)

                if abs(info[filename]['POFFSET'] - acq_poff) > slitwidth/2.0 or \
                    abs(info[filename]['QOFFSET'] - acq_qoff) > slitlength/2.0:
                    logger.debug('Adding %s to the sky list', filename)
                    skyFrameList.append(filename)

            else:
                logger.error('Unknown class for %s', filename)

        elif info[filename]['OBSTYPE'] == 'ARC':
            logger.debug('Adding %s to the arc list', filename)
            arclist.append(filename)

        elif info[filename]['OBSTYPE'] == 'DARK':
            logger.debug('Adding %s to the dark list', filename)
            darklist.append(filename)

        elif info[filename]['OBSTYPE'] == 'FLAT':

            if info[filename]['GCALLAMP'] == 'QH' and 'Pinholes' in info[filename]['SLIT']:
                logger.debug('Adding %s to pinhole list', filename)
                pinholelist.append(filename)

            elif info[filename]['GCALLAMP'] == 'QH'and 'Pinholes' not in info[filename]['SLIT']:
                logger.debug('Adding %s to the QH flat list', filename)
                QHflatlist.append(filename)

            elif info[filename]['GCALLAMP'] == 'IRhigh':
                logger.debug('Adding %s to the IR flat list', filename)
                IRflatlist.append(filename)

            else:
                logger.error('Unknown flat configuration for %s', filename)
        else:
            logger.error('Unknown configuration for %s', filename)

    # Make list of unique [date, obsid] pairs from FLATS.
    # If flat was taken on the same day as a science frame, append that flat date.
    # If not, append an arbitrary unique date from sciDateList.
    # This is so we can sort calibrations later by date and obsid.
    # n = 0
    # for flat in QHflatlist:
    #    # Make sure no duplicate dates are being entered.
    #    if QHflatlist.index(flat) == 0 or not oldobsid == info[file]['OBSID'].replace('-','_'):
    #        if info[file]['DATE-OBS'].replace('-','') in sciDateList:
    #            list1 = [info[file]['DATE-OBS'].replace('-',''), info[file]['OBSID'].replace('-','_')]
    #        else:
    #            list1 = [sciDateList[n], info[file]['OBSID'].replace('-','_')]
    #        obsidDateList.append(list1)
    #        n += 1
    #    oldobsid = info[file]['OBSID'].replace('-','_')

    # Andy doesn't understand this last section.  How about something like this?
    for f in QHflatlist:
        date_obsid = [info[f]['DATE-OBS'].replace('-', ''), info[f]['OBSID'].replace('-', '_')]
        if date_obsid not in obsidDateList:
            obsidDateList.append(date_obsid)

    logger.info("")
    logger.info("Length allfilelist (science and telluric frames): %d", len(allfilelist))
    logger.info("Length sciImageList (science and science sky frames): %d", len(sciImageList))
    logger.info("Length arclist (arc frames): %d", len(arclist))
    logger.info("Length darklist (dark frames): %d", len(darklist))
    logger.info("Length QHflatlist (QHlamps on flat frames): %d", len(QHflatlist))
    logger.info("Length IRflatlist (IRlamps on flat frames): %d", len(IRflatlist))
    logger.info("Length pinholelist (pinhole flat frames): %d", len(pinholelist))
    logger.info("Length skyFrameList (science sky frames): %d", len(skyFrameList))
    logger.info("")
    logger.info("Science list: %s", sciImageList)
    logger.info("Science dates: %s", sciDateList)
    logger.info("QH flat list: %s", QHflatlist)
    logger.info("arclist: %s", arclist)
    logger.info("obsidDateList: %s", obsidDateList)

    logger.info('Need to copy %d science, sky, Telluric, and acquisition frames', len(allfilelist))
    logger.info('Need to copy %d flats, arcs, and darks',
                len(arclist) + len(darklist) + len(QHflatlist) + len(IRflatlist) + len(pinholelist))

    # ------------------------------------------------------------------------------------------------------------------
    # Copy the FITS files into the appropriate directories

    logger.info("Making new directories and copying files...")

    objDirList = []  # List of paths sorted by object and date. ['path/to/object1/date1', 'path/to/object1/date2'].

    # scienceDirList:
    # First part is a list of calculated times of each frame in that science directory.
    # Second part is a path to a science directory.
    # Eg: [[[5400,6500,7200], '/path/to/first/science'], [[3400,4300,5200], '/path/to/second/science']...]
    scienceDirList = []

    telDirList = []  # List of paths to telluric directories.

    path = os.getcwd()
    if rawPath:
        rawPath = rawPath
    else:
        rawPath = path+'/rawPath'

    # All data is sorted by science frame data.
    # For each science frame, create a "science_object_name/date/" directory in the current working directory.

    for filename in sorted(info.keys()):

        logger.debug('%s: OBSCLASS:%s OBSTYPE:%s OBJECT:%s', filename, info[filename]['OBSCLASS'],
                     info[filename]['OBSTYPE'], info[filename]['OBJECT'])

        if info[filename]['OBSTYPE'] == 'OBJECT' and info[filename]['OBSCLASS'] == 'science':

            newpath = path + '/' + \
                info[filename]['OBJECT'] + '/' + \
                info[filename]['DATE-OBS'] + '/' + \
                info[filename]['CONFIG'] + '/' + \
                info[filename]['OBSID']
            logger.debug('newpath: %s', newpath)

            if not os.path.exists(newpath):
                logger.debug('Creating %s', newpath)
                os.makedirs(newpath)
                objDirList.append(newpath)
            logger.debug('Copying %s to %s', filename, newpath)
            shutil.copy2(rawPath + '/' + filename, newpath)
            # TODO:  Write text files with the names of the science frames and sky frames for later use by the pipeline.

        elif info[filename]['OBSTYPE'] == 'OBJECT' and info[filename]['OBSCLASS'] == 'partnerCal':

            # Associate this Telluric with all the science targets with the same config within 1.5 hours:
            matches = []
            for f in info.keys():
                if info[f]['OBSCLASS'] == 'science' and \
                        info[f]['CONFIG'] == info[filename]['CONFIG'] and \
                        abs(info[f]['AVETIME'] - info[filename]['AVETIME']) < datetime.timedelta(hours=1.5) and \
                        info[f]['OBJECT'] not in matches:
                    matches.append(info[f]['OBJECT'])
            logger.debug('This Telluric matches: %s', matches)
            for m in matches:
                newpath = path + '/' + m + '/' + \
                    info[filename]['DATE-OBS'] + '/' + \
                    info[filename]['CONFIG'] + '/' + \
                    'Tellurics/' + \
                    info[filename]['OBSID']
                logger.debug('newpath: %s', newpath)
                if not os.path.exists(newpath):
                    logger.debug('Creating %s', newpath)
                    os.makedirs(newpath)
                logger.debug('Copying %s to %s', filename, newpath)
                shutil.copy2(rawPath + '/' + filename, newpath)

        elif info[filename]['OBSTYPE'] in ['FLAT', 'ARC', 'DARK']:

            # Associate these calibrations with all the science targets on the same night with the same config
            matches = []
            for f in info.keys():
                if info[f]['OBSCLASS'] == 'science' and \
                        info[f]['CONFIG'] == info[filename]['CONFIG'] and \
                        info[f]['DATE-OBS'] == info[filename]['DATE-OBS'] and \
                        info[f]['OBJECT'] not in matches:
                    matches.append(info[f]['OBJECT'])
            logger.debug('This calibration matches: %s', matches)
            for m in matches:
                newpath = path + '/' + m + '/' + \
                    info[filename]['DATE-OBS'] + '/' + \
                    info[filename]['CONFIG'] + '/' + \
                    'Calibrations'
                logger.debug('newpath: %s', newpath)
                if not os.path.exists(newpath):
                    logger.debug('Creating %s', newpath)
                    os.makedirs(newpath)
                logger.debug('Copying %s to %s', filename, newpath)
                shutil.copy2(rawPath + '/' + filename, newpath)

        else:
            logger.warning('Unhandlded file: %s', filename)









    raise SystemExit('WHOAH')

    calDirList = []
    filelist = ['arclist', 'arcdarklist', 'flatlist', 'ronchilist', 'flatdarklist']

    # Create Calibrations directories in each of the observation date directories, e.g: YYYYMMDD/Calibrations

    logger.info("Sorting flats:")

    for i in range(len(flatlist)):
        for objDir in objDirList:
            for item in obsidDateList:
                if obsid in item:
                    date = item[0]
                    if date in objDir:
                        for entry in objectDateGratingList:
                            if entry[1] == date and entry[2] == grating:
                                if not os.path.exists(path1+'/'+entry[0]+'/'+entry[1]+'/Calibrations_'+grating):
                                    os.mkdir(path1+'/'+entry[0]+'/'+entry[1]+'/Calibrations_'+grating)
                                    calDirList.append(path1+'/'+entry[0]+'/'+entry[1]+'/Calibrations_'+grating)
                                else:
                                    if path1+'/'+entry[0]+'/'+entry[1]+'/Calibrations_'+grating not in calDirList:
                                        calDirList.append(path1+'/'+entry[0]+'/'+entry[1]+'/Calibrations_'+grating)
                                # Copy lamps on flats to appropriate directory.
                                shutil.copy('./'+flatlist[i][0], objDir+'/Calibrations_'+grating+'/')
                                flatlist[i][1] = 0
                                logger.info(flatlist[i][0])
                                count += 1
                                path = objDir+'/Calibrations_'+grating+'/'
                                # Create a flatlist in the relevent directory.
                                # Create a text file called flatlist to store the names of the
                                # lamps on flats for later use by the pipeline.
                                writeList(flatlist[i][0], 'flatlist', path)

    # Sort lamps off flats.
    logger.info("Sorting lamps off flats:")
    for i in range(len(flatdarklist)):
        os.chdir(rawPath)
        header = astropy.io.fits.open(flatdarklist[i][0])
        obsid = header[0].header['OBSID']
        grating = header[0].header['GRATING'][0:1]
        for objDir in objDirList:
            for item in obsidDateList:
                if obsid in item:
                    date = item[0]
                    if date in objDir:
                        if not os.path.exists(objDir+'/Calibrations_'+grating):
                            os.mkdir(objDir+'/Calibrations_'+grating)
                        shutil.copy('./'+flatdarklist[i][0], objDir+'/Calibrations_'+grating+'/')
                        flatdarklist[i][1] = 0
                        logger.info(flatdarklist[i][0])
                        count += 1
                        path = objDir+'/Calibrations_'+grating+'/'
                        # Create a flatdarklist in the relevant directory.
                        writeList(flatdarklist[i][0], 'flatdarklist', path)

    # Sort pinholes

    # Sort arcs.
    logger.info("Sorting arcs:")
    for i in range(len(arclist)):
        header = astropy.io.fits.open(arclist[i][0])
        date = header[0].header['DATE'].replace('-','')
        grating = header[0].header['GRATING'][0:1]
        for objDir in objDirList:
            if date in objDir:
                if not os.path.exists(objDir+'/Calibrations_'+grating):
                    os.mkdir(objDir+'/Calibrations_'+grating)
                shutil.copy('./'+arclist[i][0], objDir+'/Calibrations_'+grating+'/')
                arclist[i][1] = 0
                logger.info(arclist[i][0])
                count += 1
                path = objDir+'/Calibrations_'+grating+'/'
                # Create an arclist in the relevant directory.
                writeList(arclist[i][0], 'arclist', path)

    # Sort arc darks.
    logger.info("\nSorting arc darks:")
    for i in range(len(arcdarklist)):
        header = astropy.io.fits.open(arcdarklist[i][0])
        obsid = header[0].header['OBSID']
        grating = header[0].header['GRATING'][0:1]
        for objDir in objDirList:
            for item in obsidDateList:
                if obsid in item:
                    date = item[0]
                    if date in objDir:
                        if not os.path.exists(objDir+'/Calibrations_'+grating):
                            os.mkdir(objDir+'/Calibrations_'+grating)
                        shutil.copy('./'+arcdarklist[i][0], objDir+'/Calibrations_'+grating+'/')
                        arcdarklist[i][1] = 0
                        logger.info(arcdarklist[i][0])
                        count += 1
                        path = objDir+'/Calibrations_'+grating+'/'
                        # Create an arcdarklist in the relevant directory.
                        writeList(arcdarklist[i][0], 'arcdarklist', path)










    # ------------------------------------------------------------------------------------------------------------------
    # Check that each science directory exists and has associated calibration data.

    # Pseudocode:
    # For each science directory, make sure that:
    # a calibrations directory is present.
    # flatlist exists and has more than one file.
    # flatdarklist exists and has more than one file.
    # arclist exists and has more than one file.
    # arcdarklist exists and has more than one file.
    # ronchilist exists and has more than one file.


    # TODO(nat): This is horrifying. Good grief, wrap these repetitive calls in a function!
    logger.info("\nChecking that each science image has required calibration data. ")
    # For each science image, read its header data and try to change to the appropriate directory.
    # Check that:
    for i in range(len(sciImageList)):
        header = astropy.io.fits.open(rawPath+'/'+sciImageList[i])

        obstype = header[0].header['OBSTYPE']
        obsid = header[0].header['OBSID'][-3:].replace('-','')
        grat = header[0].header['GRATING'][0:1]
        date = header[0].header[ 'DATE'].replace('-','')
        obsclass = header[0].header['OBSCLASS']
        obj = header[0].header['OBJECT']
        obj = re.sub('[^a-zA-Z0-9\n\.]', '', obj)

        # a science and Calibrations directory are present.
        try:
            os.chdir(path1+'/'+obj+'/'+date+'/'+grat+'/obs'+obsid+'/')
            os.chdir('../../Calibrations_'+grat+'/')
        except OSError:
            logger.info("\n#####################################################################")
            logger.info("#####################################################################")
            logger.info("")
            logger.info("     WARNING in sort: no Calibrations directory found for ")
            logger.info("                      science frame "+str(sciImageList[i]))
            logger.info("")
            logger.info("#####################################################################")
            logger.info("#####################################################################\n")
            continue

        # flatlist exists and has more than one file.
        try:
            flatListFile = open('flatlist', "r").readlines()
            if len(flatListFile) <= 1:
                logger.info("\n#####################################################################")
                logger.info("#####################################################################")
                logger.info("")
                logger.info("     WARNING in sort: only 1 lamps on flat frame found for science")
                logger.info("                      frame "+str(sciImageList[i]))
                logger.info("")
                logger.info("#####################################################################")
                logger.info("#####################################################################\n")
        except IOError:
            logger.info("\n#####################################################################")
            logger.info("#####################################################################")
            logger.info("")
            logger.info("     WARNING in sort: no flatlist found for science frame")
            logger.info("                      "+str(sciImageList[i]))
            logger.info("")
            logger.info("#####################################################################")
            logger.info("#####################################################################\n")

            if not manualMode:
                # Sometimes flats can be taken a day after the observing night. First
                # look for these, and if they are not found, ask the user to provide some flats.
                foundflatFlag = False
                # Get date after the science observation
                t=time.strptime(date,'%Y%m%d')
                newdate=datetime.date(t.tm_year,t.tm_mon,t.tm_mday)+datetime.timedelta(1)
                # Loop through flatlist and see if there is an flat taken on this date
                for i in range(len(flatlist)):
                    header = astropy.io.fits.open(rawPath+'/'+flatlist[i][0])
                    date = header[0].header[ 'DATE'].replace('-','')
                    if str(date) == newdate.strftime('%Y%m%d'):
                        # If so, copy it to the appropriate calibrations directory and write an flatlist.
                        shutil.copy(rawPath + '/' + flatlist[i][0], './')
                        writeList(flatlist[i][0], 'flatlist', path)
                        logger.info("\n#####################################################################")
                        logger.info("#####################################################################")
                        logger.info("")
                        logger.info("     WARNING in sort: found a flat taken one day after a science frame.")
                        logger.info("                      "+str(sciImageList[i]))
                        logger.info("                       using that.")
                        logger.info("")
                        logger.info("#####################################################################")
                        logger.info("#####################################################################\n")
                        foundflatFlag = True
                        flatlist[i][1] = 0
                        
                if not foundflatFlag:
                    # If that quick check fails, give user a chance to try and provide an flat file.
                    # a = raw_input("\n Please provide a textfile called flatlist in " + str(os.getcwd()))
                    raise SystemError('No textfile called flatlist found in ' + str(os.getcwd()))


        # flatdarklist exists and has more than one file.
        try:
            flatDarkListFile = open('flatdarklist', "r").readlines()
            if len(flatDarkListFile) <= 1:
                logger.info("\n#####################################################################")
                logger.info("#####################################################################")
                logger.info("")
                logger.info("     WARNING in sort: only 1 lamps off flat frame found for science")
                logger.info("                      frame "+str(sciImageList[i]))
                logger.info("")
                logger.info("#####################################################################")
                logger.info("#####################################################################\n")
        except IOError:
            logger.info("\n#####################################################################")
            logger.info("#####################################################################")
            logger.info("")
            logger.info("     WARNING in sort: no flatdarklist found for science frame")
            logger.info("                      "+str(sciImageList[i]))
            logger.info("")
            logger.info("#####################################################################")
            logger.info("#####################################################################\n")
            if not manualMode:
                # Sometimes flatdarks can be taken a day after the observing night. First
                # look for these, and if they are not found, ask the user to provide some flatdarks.
                foundflatdarkFlag = False
                # Get date after the science observation
                t=time.strptime(date,'%Y%m%d')
                newdate=datetime.date(t.tm_year,t.tm_mon,t.tm_mday)+datetime.timedelta(1)
                # Loop through flatdarklist and see if there is an flatdark taken on this date
                for i in range(len(flatdarklist)):
                    header = astropy.io.fits.open(rawPath+'/'+flatdarklist[i][0])
                    date = header[0].header[ 'DATE'].replace('-','')
                    if str(date) == newdate.strftime('%Y%m%d'):
                        # If so, copy it to the appropriate calibrations directory and write an flatdarklist.
                        shutil.copy(rawPath + '/' + flatdarklist[i][0], './')
                        writeList(flatdarklist[i][0], 'flatdarklist', path)
                        logger.info("\n#####################################################################")
                        logger.info("#####################################################################")
                        logger.info("")
                        logger.info("     WARNING in sort: found a flatdark taken one day after a science frame.")
                        logger.info("                      "+str(sciImageList[i]))
                        logger.info("                       using that.")
                        logger.info("")
                        logger.info("#####################################################################")
                        logger.info("#####################################################################\n")
                        foundflatdarkFlag = True
                        flatdarklist[i][1] = 0
                        
                if not foundflatdarkFlag:
                    # If that quick check fails, give user a chance to try and provide an flatdark file.
                    raise SystemError('No textfile called flatdarklist found in ' + str(os.getcwd()) ) 
#                     a = raw_input("\n Please provide a textfile called flatdarklist in " + str(os.getcwd()) + \
#                     " or be sure not to attempt a wavelength calibration for this directory.")
                    
        # Make sure flatlist and flatdarklist are the same length. nsflat() complains otherwise.
        checkSameLengthFlatLists()

        # arclist exists.
        try:
            arcListFile = open('arclist', "r").readlines()
        except IOError:
            logger.info("\n#####################################################################")
            logger.info("#####################################################################")
            logger.info("")
            logger.info("     WARNING in sort: no arclist found for science frame")
            logger.info("                      "+str(sciImageList[i]))
            logger.info("")
            logger.info("#####################################################################")
            logger.info("#####################################################################\n")

            if not manualMode:
                # Sometimes arcs can be taken a day after the observing night. First
                # look for these, and if they are not found, ask the user to provide some arcs.
                foundArcFlag = False
                # Get date after the science observation
                t=time.strptime(date,'%Y%m%d')
                newdate=datetime.date(t.tm_year,t.tm_mon,t.tm_mday)+datetime.timedelta(1)
                # Loop through arclist and see if there is an arc taken on this date
                for i in range(len(arclist)):
                    header = astropy.io.fits.open(rawPath+'/'+arclist[i][0])
                    date = header[0].header[ 'DATE'].replace('-','')
                    if str(date) == newdate.strftime('%Y%m%d'):
                        # If so, copy it to the appropriate calibrations directory and write an arclist.
                        shutil.copy(rawPath + '/' + arclist[i][0], './')
                        writeList(arclist[i][0], 'arclist', path)
                        logger.info("\n#####################################################################")
                        logger.info("#####################################################################")
                        logger.info("")
                        logger.info("     WARNING in sort: found an arc taken one day after a science frame.")
                        logger.info("                      "+str(sciImageList[i]))
                        logger.info("                       using that.")
                        logger.info("")
                        logger.info("#####################################################################")
                        logger.info("#####################################################################\n")
                        foundArcFlag = True
                        arclist[i][1] = 0
                        
                if not foundArcFlag:
                    # If that quick check fails, give user a chance to try and provide an arc file.
                    raise SystemError('No textfile called arclist found in ' + str(os.getcwd()) ) 
#                     a = raw_input("\n Please provide a textfile called arclist in " + str(os.getcwd()) + \
#                     " or be sure not to attempt a wavelength calibration for this directory.")

        # arcdarklist exists.
        try:
            arcDarkListFile = open('arcdarklist', "r").readlines()
        except IOError:
            logger.info("\n#####################################################################")
            logger.info("#####################################################################")
            logger.info("")
            logger.info("     WARNING in sort: no arcdarklist found for science frame")
            logger.info("                      "+str(sciImageList[i]))
            logger.info("")
            logger.info("#####################################################################")
            logger.info("#####################################################################\n")
            if not manualMode:
                # Sometimes arcdarks can be taken a day after the observing night. First
                # look for these, and if they are not found, ask the user to provide some arcdarks.
                foundarcdarkFlag = False
                # Get date after the science observation
                t=time.strptime(date,'%Y%m%d')
                newdate=datetime.date(t.tm_year,t.tm_mon,t.tm_mday)+datetime.timedelta(1)
                # Loop through arcdarklist and see if there is an arcdark taken on this date
                for i in range(len(arcdarklist)):
                    header = astropy.io.fits.open(rawPath+'/'+arcdarklist[i][0])
                    date = header[0].header[ 'DATE'].replace('-','')
                    if str(date) == newdate.strftime('%Y%m%d'):
                        # If so, copy it to the appropriate calibrations directory and write an arcdarklist.
                        shutil.copy(rawPath + '/' + arcdarklist[i][0], './')
                        writeList(arcdarklist[i][0], 'arcdarklist', path)
                        logger.info("\n#####################################################################")
                        logger.info("#####################################################################")
                        logger.info("")
                        logger.info("     WARNING in sort: found an arcdark taken one day after a science frame.")
                        logger.info("                      "+str(sciImageList[i]))
                        logger.info("                       using that.")
                        logger.info("")
                        logger.info("#####################################################################")
                        logger.info("#####################################################################\n")
                        foundarcdarkFlag = True
                        arcdarklist[i][1] = 0
                        
                if not foundarcdarkFlag:
                    # If that quick check fails, give user a chance to try and provide an arcdark file.
                    raise SystemError('No textfile called arcdarklist found in ' + str(os.getcwd()) ) 
#                     a = raw_input("\n Please provide a textfile called arcdarklist in " + str(os.getcwd()) + \
#                     " or be sure not to attempt a wavelength calibration for this directory.")
                    
        # ronchilist exists and has more than one file.
        try:
            ronchiListFile = open('ronchilist', "r").readlines()
            if len(ronchiListFile) <= 1:
                logger.info("\n#####################################################################")
                logger.info("#####################################################################")
                logger.info("")
                logger.info("     WARNING in sort: only 1 ronchi flat frame found for science frame")
                logger.info("                      "+str(sciImageList[i]))
                logger.info("")
                logger.info("#####################################################################")
                logger.info("#####################################################################\n")
        except IOError:
            logger.info("\n#####################################################################")
            logger.info("#####################################################################")
            logger.info("")
            logger.info("     WARNING in sort: no ronchilist found for science frame")
            logger.info("                      "+str(sciImageList[i]))
            logger.info("")
            logger.info("#####################################################################")
            logger.info("#####################################################################\n")
            if not manualMode:
                # Sometimes ronchis can be taken a day after the observing night. First
                # look for these, and if they are not found, ask the user to provide some ronchis.
                foundronchiFlag = False
                # Get date after the science observation
                t=time.strptime(date,'%Y%m%d')
                newdate=datetime.date(t.tm_year,t.tm_mon,t.tm_mday)+datetime.timedelta(1)
                # Loop through ronchilist and see if there is an ronchi taken on this date
                for i in range(len(ronchilist)):
                    header = astropy.io.fits.open(rawPath+'/'+ronchilist[i][0])
                    date = header[0].header[ 'DATE'].replace('-','')
                    if str(date) == newdate.strftime('%Y%m%d'):
                        # If so, copy it to the appropriate calibrations directory and write an ronchilist.
                        shutil.copy(rawPath + '/' + ronchilist[i][0], './')
                        writeList(ronchilist[i][0], 'ronchilist', path)
                        logger.info("\n#####################################################################")
                        logger.info("#####################################################################")
                        logger.info("")
                        logger.info("     WARNING in sort: found a ronchi taken one day after a science frame.")
                        logger.info("                      "+str(sciImageList[i]))
                        logger.info("                       using that.")
                        logger.info("")
                        logger.info("#####################################################################")
                        logger.info("#####################################################################\n")
                        foundronchiFlag = True
                        ronchilist[i][1] = 0
                        
                if not foundronchiFlag:
                    # If that quick check fails, give user a chance to try and provide an ronchi file.
                    raise SystemError('No textfile called ronchilist found in ' + str(os.getcwd()) ) 
#                     a = raw_input("\n Please provide a textfile called ronchilist in " + str(os.getcwd()) + \
#                     " or be sure not to attempt a wavelength calibration for this directory.")
                    
        os.chdir(path1)

    # Change back to original working directory.
    os.chdir(path1)

    return calDirList





# ----------------------------------------------------------------------------------------------------------------------

def matchTellurics(telDirList, obsDirList, telluricTimeThreshold):

    """
    Matches science images with the telluric frames that are closest in time.
    Creates a file in each telluric observation directory called scienceMatchedTellsList.
    scienceMatchedTellsList lists the obsid of the science images (ie. obs123) and then the
    science images with this obsid that match the telluric observation.

    EXAMPLE:    obs28
                N20130527S0264
                N20130527S0266
                obs30
                N201305727S0299

    """
    logger = log.getLogger('matchTellurics')
    logger.info("Matching science images with tellurics closest in time.")

    # Store current working directory for later use.
    path = os.getcwd()

    # Get a list of unique dates of telluric observations.
    dateList=[]
    for i in range(len(telDirList)):
        date = telDirList[i].split(os.sep)[-4]
        if i==0 or dateList[-1]!=date:
            dateList.append(date)

    # Make a 2D list; list of lists of telluric directory path, files in that directory, grating of those tellurics.
    # [['telluric_directory/', ['telluric1.fits, telluric2.fits, ...'], 'telluric_grating']]
    for date in dateList:
        tellist = []
        for telDir in telDirList:
            if date in telDir:

                os.chdir(telDir)
                telImageList = open(telDir + '/' + 'tellist', "r").readlines()
                telImageList = [image for image in telImageList]
                telluric_image = telImageList[0]
                telluric_header = astropy.io.fits.open(telDir +'/'+ telluric_image + '.fits')
                telluric_grating = telluric_header[0].header['GRATING'][0:1]

                timeList=[]
                if os.path.exists('./scienceMatchedTellsList'):
                    os.remove('./scienceMatchedTellsList')
                templist = []
                imageList=glob.glob('N*.fits')
                templist.append(telDir)
                templist.append(imageList)
                templist.append(telluric_grating)
                tellist.append(templist)

        # Create a list of the start and stop times for each observation called timeList.
        # timeList is of the form [[obsid1, start1, stop1], [obsid2, start2, stop2],...]
        for a in range(len(tellist)):
            templist=[]
            os.chdir(tellist[a][0])
            telheader = astropy.io.fits.open(tellist[a][1][0])
            start=timeCalc(tellist[a][1][0])
            stop=timeCalc(tellist[a][1][-1])
            templist.append(os.getcwd())
            templist.append(start)
            templist.append(stop)
            templist.append(tellist[a][2])
            timeList.append(templist)

        # Find a science image with the same date and grating.
        for obsDir in obsDirList:
            os.chdir(obsDir)
            if date in obsDir:
                try:
                    sciImageList = open('scienceFrameList', "r").readlines()
                except IOError:
                    sciImageList = open('skyFrameList', "r").readlines()
                sciImageList = [image for image in sciImageList]

                # Open image and get science image grating from header.

                science_image = sciImageList[0]
                science_header = astropy.io.fits.open('./'+ science_image + '.fits')
                science_grating = science_header[0].header['GRATING'][0:1]

                for image in sciImageList:
                    diffList=[]
                    imageTime = timeCalc(image+'.fits')
                    for b in range(len(timeList)):
                        # Check to make sure telluric grating and science grating match.
                        if timeList[b][3] == science_grating:
                            if abs(imageTime-timeList[b][1]) <= int(telluricTimeThreshold) or abs(imageTime-timeList[b][2]) <= int(telluricTimeThreshold):
                                if abs(imageTime-timeList[b][1]) < abs(imageTime-timeList[b][2]):
                                    diff = abs(imageTime-timeList[b][1])
                                else:
                                    diff = abs(imageTime-timeList[b][2])
                                diffList.append(timeList[b][0])
                                diffList.append(diff)

                    # Find and record the science observation that is closest in time to the telluric image.
                    # Store the science observation name in a textfile, scienceMatchedTellsList, for later use by the pipeline.
                    if diffList:
                        minDiff = min(diffList)
                        telobs = diffList[diffList.index(minDiff)-1]
                        sciheader = astropy.io.fits.open(image+'.fits')
                        sciObsid = 'obs'+ sciheader[0].header['OBSID'][-3:].replace('-','')
                        if not os.path.exists(telobs+'/scienceMatchedTellsList'):
                            writeList(sciObsid, 'scienceMatchedTellsList', telobs)
                        else:
                            scienceMatchedTellsList = open(telobs+'/scienceMatchedTellsList', 'r').readlines()
                            scienceMatchedTellsList = [item for item in scienceMatchedTellsList]
                            if sciObsid not in scienceMatchedTellsList:
                                writeList(sciObsid, 'scienceMatchedTellsList', telobs)
                        writeList(image, 'scienceMatchedTellsList', telobs)
    os.chdir(path)


    # ---------------------------- Tests ------------------------------------- #

    # Don't use tests if user doesn't want them
    tests = True
    if tests:
        # Check that each science observation has valid telluric data.

        # For each science observation:
        for science_directory in obsDirList:
            os.chdir(science_directory)
            # Store science observation name in science_observation_name
            science_observation_name = science_directory.split(os.sep)[-1]
            # Optional: store time of a science frame in science_time.
            try:
                sciImageList = open('scienceFrameList', "r").readlines()
            except IOError:
                logger.info("\n#####################################################################")
                logger.info("#####################################################################")
                logger.info("")
                logger.info("     WARNING in sort: science "+str(science_observation_name))
                logger.info("                      in " + str(os.getcwd()))
                logger.info("                      does not have a scienceFrameList.")
                logger.info("                      I am trying to rewrite it with zero point offsets.")
                logger.info("")
                logger.info("#####################################################################")
                logger.info("#####################################################################\n")

                rewriteSciImageList(2.0, "Science")
                try:
                    sciImageList = open('scienceFrameList', "r").readlines()
                    logger.info("\nSucceeded; a science frame list exists in " + str(os.getcwd()))
                except IOError:
                    logger.info("\nWARNING: no science frames found in " + str(os.getcwd()) + ". You may have to adjust the skyThreshold parameter.")
                    raise SystemError

            sciImageList = [image for image in sciImageList]
            for science_image in sciImageList:
                scienceDirectory = os.getcwd()
                # Open image and get science image grating from header.
                science_header = astropy.io.fits.open('./'+ science_image + '.fits')
                science_time = timeCalc(science_image+'.fits')
                science_date = science_header[0].header[ 'DATE'].replace('-','')

                # Check that directory obsname matches header obsname.
                temp_obs_name = 'obs' + science_header[0].header['OBSID'][-3:].replace('-','')
                if science_observation_name != temp_obs_name:
                    logger.info("\n#####################################################################")
                    logger.info("#####################################################################")
                    logger.info("")
                    logger.info("     WARNING in sort: science "+str(science_observation_name)+ " :")
                    logger.info("                      observation name data in headers and directory")
                    logger.info("                      do not match.")
                    logger.info("")
                    logger.info("#####################################################################")
                    logger.info("#####################################################################\n")

                # Check that a tellurics directory exists.
                if os.path.exists('../Tellurics/'):
                    os.chdir('../Tellurics/')
                else:
                    logger.info("\n#####################################################################")
                    logger.info("#####################################################################")
                    logger.info("")
                    logger.info("     WARNING in sort: telluric directory for science "+str(science_observation_name))
                    logger.info("                      does not exist.")
                    logger.info("")
                    logger.info("#####################################################################")
                    logger.info("#####################################################################\n")
                    continue
                found_telluric_flag = False

                # Iterate through tellurics observation directories.
                for directory in list(glob.glob('obs*')):
                    os.chdir('./'+directory)
                    # Check that a file, scienceMatchedTellsList exists.
                    try:
                        scienceMatchedTellsList = open('scienceMatchedTellsList', "r").readlines()
                        # Check that the science observation name is in the file.
                        # Check that immediately after is at least one telluric image name.
                        # Do this by checking for the science date in the telluric name.
                        for i in range(len(scienceMatchedTellsList)):
                            telluric_observation_name = scienceMatchedTellsList[i]
                            if telluric_observation_name == science_observation_name:
                                rest_list = [x for x in scienceMatchedTellsList]
                                if science_image in rest_list:
                                    found_telluric_flag = True
                                    break
                    except IOError:
                        pass

                    if found_telluric_flag:
                        os.chdir('../')
                        break
                    else:
                        os.chdir('../')

                if not found_telluric_flag:
                    os.chdir('../')
                    logger.info("\n#####################################################################")
                    logger.info("#####################################################################")
                    logger.info("")
                    logger.info("     WARNING in sort: no tellurics data found for science "+str(science_image))
                    logger.info("     in " + str(os.getcwd()))
                    logger.info("")
                    logger.info("#####################################################################")
                    logger.info("#####################################################################\n")
                    logger.info("\n#####################################################################")
                    logger.info("#####################################################################")
                    logger.info("")
                    logger.info("    TURNING OFF TELLURIC-RELATED STEPS IN CONFIG FILE.")
                    logger.info("")
                    logger.info("#####################################################################")
                    logger.info("#####################################################################\n")
                    # Turn off the telluric correction.
                    logger.info("\nnifsSort: no tellurics data found for a directory. Turning off telluric reduction, telluric correction and flux calibration in ./config.cfg.")
                    with open('../../../config.cfg') as config_file:
                        options = ConfigObj(config_file, unrepr=True)
                        nifsPipelineConfig = options['nifsPipelineConfig']
                    nifsPipelineConfig['telluricReduction'] = False
                    nifsPipelineConfig['telluricCorrection'] = False
                    nifsPipelineConfig['fluxCalibration'] = False
                    with open('../../../config.cfg', 'w') as config_file:
                        options.write(config_file)

                os.chdir(scienceDirectory)
                # TODO(nat):
                # Optional: open that telluric image and store time in telluric_time
                # Check that abs(telluric_time - science_time) < 1.5 hours

    os.chdir(path)
    logger.info("\nI am finished matching science images with telluric frames.")
    return


    # ------------------------------------------------------------------------------------------------------------------

    # Update config.cfg with the paths to the data:
    logger.info("Updating config file with scienceDirectoryList, calibrationDirectoryList and telluricDirectoryList")
    config.set('defaults','scienceDirectoryList', scienceDirectoryList)
    config.set('defaults','telluricDirectoryList', telluricDirectoryList)
    # config.set('defaults','calibrationDirectoryList', calibrationDirectoryList)
    with open(path + '/config.cfg', 'w') as f:
        config.write(f)

    # ------------------------------------------------------------------------------------------------------------------

    raise SystemExit('STOP')
