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

    for file in sorted(info.keys()):

        if info[file]['INSTRUME'] != 'GNIRS':
            logger.debug('Skipping %s from %s', file, info[file]['INSTRUME'])
            continue

        # PRISM will be one of: SB+SXD_G5536, LB+LXD_G5535, MIR_G5511, SXD_G5509, ...
        if 'SXD' not in info[file]['PRISM'] and 'LXD' not in info[file]['PRISM']:
            logger.debug('Skipping %s which is not cross-dispersed', file)
            continue

        # Make a list of science, telluric and acquisition frames.
        # Use the copied variable (1 not copied, 0 copied) to check later that the file was copied correctly.
        # allfilelist is a 2D list of [[filename1, copied, obsclass1], [filename2, copied, obsclass2]] pairs.

        if info[file]['OBSTYPE'] == 'OBJECT' and info[file]['OBSCLASS'] in ['science', 'acq', 'acqCal', 'partnerCal']:

            # Append a list of [filename, copied, obsclass] to the list.  obsclass is used later for checks.
            templist = [file, 1, info[file]['OBSCLASS']]
            allfilelist.append(templist)

            if info[file]['OBSCLASS'] == 'acq':
                # Use the P and Q offsets of the last acquisition image to calculate the absolute P and Q offsets.
                # The absolute P and Q offsets are then used to determine if the science target is in the slit.
                # This assumes that the acq and science frames are in order and uninterrupted!
                acq_poff = info[file]['POFFSET']
                acq_qoff = info[file]['QOFFSET']

            elif info[file]['OBSCLASS'] == 'science':
                logger.debug('Adding %s to the science list', file)
                sciImageList.append(file)

                # create a summary of this file's configuration:
                odc = [re.sub('[^a-zA-Z0-9]', '', info[file]['OBJECT']),  # remove irregular characters
                       info[file]['DATE-OBS'].replace('-',''),
                       info[file]['CAMERA'][:info[file]['CAMERA'].find('_')],
                       info[file]['PRISM'][info[file]['PRISM'].find('+')+1:info[file]['PRISM'].find('_')],
                       info[file]['GRATING'][:info[file]['GRATING'].find('/')],
                       info[file]['SLIT'][:info[file]['SLIT'].find('_')],
                       info[file]['GRATWAVE'], info[file]['OBSID'].replace('-','_')]
                logger.debug('%s: %s', file, odc)

                if odc not in objectDateConfigurationList:
                    objectDateConfigurationList.append(odc)

                d = info[file]['DATE-OBS'].replace('-','')  # DATE-OBS has a format '2019-07-20'
                if d not in sciDateList:
                    sciDateList.append(d)

                if 'SC' in info[file]['DECKER'] and 'XD' in info[file]['DECKER']:
                    slitlength = 7.
                elif 'LC' in info[file]['DECKER'] and 'XD' in info[file]['DECKER']:
                    slitlength = 5.
                else:
                    slitlength = None
                slitwidth = float(info[file]['SLIT'][:info[file]['SLIT'].find('arcsec')])
                logger.debug('Slit size: %.1f x %.1f arcsec', slitwidth, slitlength)

                if abs(info[file]['POFFSET'] - acq_poff) > slitwidth/2.0 or \
                    abs(info[file]['QOFFSET'] - acq_qoff) > slitlength/2.0:
                    logger.debug('Adding %s to the sky list', file)
                    skyFrameList.append(file)

            else:
                logger.error('Unknown class for %s', file)

        elif info[file]['OBSTYPE'] == 'ARC':
            logger.debug('Adding %s to the arc list', file)
            arclist.append(file)

        elif info[file]['OBSTYPE'] == 'DARK':
            logger.debug('Adding %s to the dark list', file)
            darklist.append(file)

        elif info[file]['OBSTYPE'] == 'FLAT':

            if info[file]['GCALLAMP'] == 'QH' and 'Pinholes' in info[file]['SLIT']:
                logger.debug('Adding %s to pinhole list', file)
                pinholelist.append(file)

            elif info[file]['GCALLAMP'] == 'QH'and 'Pinholes' not in info[file]['SLIT']:
                logger.debug('Adding %s to the QH flat list', file)
                QHflatlist.append(file)

            elif info[file]['GCALLAMP'] == 'IRhigh':
                logger.debug('Adding %s to the IR flat list', file)
                IRflatlist.append(file)

            else:
                logger.error('Unknown flat configuration for %s', file)
        else:
            logger.error('Unknown configuration for %s', file)

    # Make list of unique [date, obsid] pairs from FLATS.
    # If flat was taken on the same day as a science frame, append that flat date.
    # If not, append an arbitrary unique date from sciDateList.
    # This is so we can sort calibrations later by date and obsid.
    #n = 0
    #for flat in QHflatlist:
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
    for file in QHflatlist:
        date_obsid = [info[file]['DATE-OBS'].replace('-',''), info[file]['OBSID'].replace('-','_')]
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
    number_files_to_be_copied = len(allfilelist)
    number_files_that_were_copied = 0

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

        if info[file]['INSTRUME'] != 'GNIRS':
            continue  # skip non-GNIRS files

        if 'SXD' not in info[filename]['PRISM'] and 'LXD' not in info[filename]['PRISM']:
            continue  # skip non-XD data

        
        if info[filename]['OBSCLASS'] == 'science':

            objDir = path + '/' + info[filename]['DATE-OBS'].replace('-','') + '/' + \
                      re.sub('[^a-zA-Z0-9]', '', info[filename]['OBJECT']) + '/' + \
                      info[filename]['CONFIG'] + '_' + \
                      info[filename]['OBSID'].replace('-','_')
            logger.debug('objDir: %s', objDir)

            if not os.path.exists(objDir):
                logger.debug('Creating %s', objDir)
                os.makedirs(objDir)
                objDirList.append(objDir)
            shutil.copy2(rawPath + '/' + filename, objDir)
            # TODO:  Write text files with the names of the science frames and sky frames for later use by the pipeline.

        elif info[filename]['OBSCLASS'] == 'partnerCal':

            # Associate this Telluric with all the science observations with the ame config within 1.5 hours
            
            matches = []
            for f in info.keys():
                if info[f]['OBSCLASS'] == 'science' and \
                       info[f]['CONFIG'] == info[filename]['CONFIG'] and \
                       abs(info[f]['AVETIME'] - info[filename]['AVETIME']) < datetime.timedelta(hours=1.5) and \
                       f not in matches:
                    matches.append(f)
            logger.debug('This Telluric matches: %s', matches)
                

            # step through all the files, find the science exposures with the same config, check the ave obs time

            

            

            raise SystemExit('WHOAH')






    
            timeList = []
            for k in range(len(scienceDirList)):
                # Make sure date and gratings match.
                tempDir = scienceDirList[k][1].split(os.sep)
                if date in tempDir and grat in tempDir:
                    # Open the times of all science frames in science_directory.
                    times = scienceDirList[k][0]
                    # Find difference in each time from the telluric frame we're trying to sort.
                    diffList = []
                    for b in range(len(times)):
                        difference = abs(telluric_time-scienceDirList[k][0][b])
                        templist = []
                        templist.append(difference)
                        templist.append(scienceDirList[k][1])
                        diffList.append(templist)
                    # Find the science frame with the smallest difference.
                    minDiff = min(diffList)
                    # Pass that time and path out of the for loop.
                    timeList.append(minDiff)
            # Out of the for loop, compare min times from different directories.
            if timeList:
                closest_time = min(timeList)
                # Copy the telluric frame to the path of that science frame.
                path_to_science_dir = closest_time[1]
                path_to_tellurics = os.path.split(path_to_science_dir)[0]

                # Create a Tellurics directory in science_object_name/YYYYMMDD/grating.

                if not os.path.exists(path_to_tellurics + '/Tellurics'):
                    os.mkdir(path_to_tellurics + '/Tellurics')
                # Create an obsid (eg. obs25) directory in the Tellurics directory.
                if not os.path.exists(path_to_tellurics+'/Tellurics/obs'+obsid):
                    os.mkdir(path_to_tellurics+'/Tellurics/obs'+obsid)
                    telDirList.append(path_to_tellurics+'/Tellurics/obs'+obsid)
                elif not telDirList or not telDirList[-1]==path_to_tellurics+'/Tellurics/obs'+obsid:
                    telDirList.append(path_to_tellurics+'/Tellurics/obs'+obsid)
                shutil.copy(rawPath+'/'+allfilelist[i][0], path_to_tellurics+'/Tellurics/obs'+obsid+'/')
                number_files_that_were_copied += 1
                allfilelist[i][1] = 0
                # Create an scienceFrameList in the relevant directory.
                if allfilelist[i][0] not in telskyFrameList:
                    writeList(allfilelist[i][0], 'tellist', path_to_tellurics+'/Tellurics/obs'+obsid+'/')
                # Create a skyFrameList in the relevant directory.
                if allfilelist[i][0] in telskyFrameList:
                    writeList(allfilelist[i][0], 'skyFrameList', path_to_tellurics+'/Tellurics/obs'+obsid+'/')

    # Modify scienceDirList to a format telSort can use.
    tempList = []
    for i in range(len(scienceDirList)):
        tempList.append(scienceDirList[i][1])
    scienceDirList = tempList


    # Check to see which files were not copied.
    logger.info("\nChecking for non-copied science, tellurics and acquisitions.\n")
    for i in range(len(allfilelist)):
        # Check the copied flag. If not 0, logger.info("the entry.")
        if allfilelist[i][1] != 0:
            logger.info(str(allfilelist[i][0]) + " " + str(allfilelist[i][2]) + " was not copied.")
    logger.info("\nEnd non-copied science, tellurics and acquisitions.\n")

    # Check that all science frames and sky frames were copied.
    count_from_raw_files = len(sciImageList) + len(skyFrameList)

    count = 0
    for science_directory in scienceDirList:
        for file in os.listdir(science_directory):
            if file.endswith('.fits'):
                count += 1

    if count_from_raw_files != count:
        logger.info("\nWARNING: " + str(count_from_raw_files - count) + " science frames (or sky frames) \
        were not copied.\n")
    else:
        logger.info("\nExpected number of science and sky frames copied.\n")

    logger.info("\nDone sorting and copying science and tellurics. Moving on to Calibrations.\n")

    # Test for telluric directories with mis-identified sky frames. For now, we identify sky frames
    # based on absolute P and Q offsets, not relative to a zero point. This can cause problems.
    # TODO(nat): look into if it is worth it to use relative P and Q offsets.
    for telluric_directory in telDirList:
        if os.path.exists(telluric_directory + '/skyFrameList') and not os.path.exists(telluric_directory + '/tellist'):
            logger.info("\n#####################################################################")
            logger.info("#####################################################################")
            logger.info("")
            logger.info("     WARNING in sort: a telluric directory exists that Nifty thinks ")
            logger.info("                      contains only sky frames. Nifty uses absolute")
            logger.info("                      P and Q offsets to identify sky frames; a target")
            logger.info("                      not being at 0, 0 P and Q can cause this. If this is not")
            logger.info("                      the case you can try adjusting the skyThreshold")
            logger.info("                      parameter in Nifty's configuration.")
            logger.info("")
            logger.info("#####################################################################")
            logger.info("#####################################################################\n")
            os.chdir(telluric_directory)
            rewriteSciImageList(2.0, "Telluric")
            try:
                sciImageList = open('tellist', "r").readlines()
                logger.info("\nSucceeded; a telluric frame list exists in " + str(os.getcwd()))
            except IOError:
                logger.info("\nWARNING: no telluric frames found in " + str(os.getcwd()) + ". You may have to adjust the skyThreshold parameter.")
                
                #a = raw_input("Please make a tellist (list of telluric frames) in " + str(telluric_directory))
                raise SystemError('No list of telluric frames found in ' + + str(telluric_directory))
                
            os.chdir(path)

    os.chdir(path)

    return














# ----------------------------------------------------------------------------------------------------------------------


def sortCalibrations(arcdarklist, arclist, flatlist, flatdarklist, ronchilist, objectDateGratingList, objDirList, obsidDateList, sciImageList, rawPath, manualMode):
    """
    Sort calibrations into appropriate directories based on date.
    """
    logger = log.getLogger('SortCals')

    calDirList = []
    filelist = ['arclist', 'arcdarklist', 'flatlist', 'ronchilist', 'flatdarklist']

    # Save path for later use. The rawPath part is for Gemini North network sorting.
    path1 = os.getcwd()
    if rawPath:
        rawPath = rawPath
    else:
        rawPath = path1+'/rawPath'

    # Set up some tests and checks.
    count = 0
    expected_count = len(arcdarklist) + len(arclist) + len(flatlist)\
          + len(flatdarklist) + len(ronchilist)

    logger.info("Attempting to sort " + str(expected_count) + " files.")

    # To make sure data was copied later in the pipeline:
    # Add a small copied flag to each frame in calibration file lists.
    new_flatlist = []
    for i in range(len(flatlist)):
        # Transform 1D list into 2D list of [[filename, 'copied']]
        # "copied" is 1 for not copied and 0 for copied.
        templist = []
        templist.append(flatlist[i])
        templist.append(1)
        new_flatlist.append(templist)
    flatlist = new_flatlist

    new_flatdarklist = []
    for i in range(len(flatdarklist)):
        # Transform 1D list into 2D list.
        templist = []
        templist.append(flatdarklist[i])
        templist.append(1)
        new_flatdarklist.append(templist)
    flatdarklist = new_flatdarklist

    new_arclist = []
    for i in range(len(arclist)):
        # Transform 1D list into 2D list.
        templist = []
        templist.append(arclist[i])
        templist.append(1)
        new_arclist.append(templist)
    arclist = new_arclist

    new_arcdarklist = []
    for i in range(len(arcdarklist)):
        # Transform 1D list into 2D list.
        templist = []
        templist.append(arcdarklist[i])
        templist.append(1)
        new_arcdarklist.append(templist)
    arcdarklist = new_arcdarklist

    new_ronchilist = []
    for i in range(len(ronchilist)):
        # Transform 1D list into 2D list.
        templist = []
        templist.append(ronchilist[i])
        templist.append(1)
        new_ronchilist.append(templist)
    ronchilist = new_ronchilist

    os.chdir(rawPath)

    # Create Calibrations directories in each of the observation date directories based on existence of
    # lamps on flats. Eg: YYYYMMDD/Calibrations
    # Sort lamps on flats.
    logger.info("\nSorting flats:")
    # Create a flag so we only warn about non-standard gratings once.
    grating_warning_flag = False
    for i in range(len(flatlist)):
        header = astropy.io.fits.open(flatlist[i][0])
        obsid = header[0].header['OBSID']
        grating = header[0].header['GRATING'][0:1]
        if grating not in ["K", "J", "H", "Z"]:
            logger.info("\n#####################################################################")
            logger.info("#####################################################################")
            logger.info("")
            logger.info("     WARNING in sort: non-standard (non K, J, H, K) grating encountered. ")
            logger.info("                      NIFTY has not been tested with non-standard")
            logger.info("                      gratings!")
            logger.info("")
            logger.info("#####################################################################")
            logger.info("#####################################################################\n")
        # TODO(nat): this is horrendous. Do this in a better way.
        # "Flat is better than nested."
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
    logger.info("\nSorting lamps off flats:")
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

    # Sort ronchi flats.
    logger.info("\nSorting ronchi flats:")
    for i in range(len(ronchilist)):
        os.chdir(rawPath)
        header = astropy.io.fits.open(ronchilist[i][0])
        obsid = header[0].header['OBSID']
        grating = header[0].header['GRATING'][0:1]
        for objDir in objDirList:
            for item in obsidDateList:
                if obsid in item:
                    date = item[0]
                    if date in objDir:
                        if not os.path.exists(objDir+'/Calibrations_'+grating):
                            os.mkdir(objDir+'/Calibrations_'+grating)
                        shutil.copy('./'+ronchilist[i][0], objDir+'/Calibrations_'+grating+'/')
                        ronchilist[i][1] = 0
                        logger.info(ronchilist[i][0])
                        count += 1
                        path = objDir+'/Calibrations_'+grating+'/'
                        # create a ronchilist in the relevant directory
                        writeList(ronchilist[i][0], 'ronchilist', path)

    # Sort arcs.
    logger.info("\nSorting arcs:")
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

    # Check that each file in flatlist was copied.
    for i in range(len(flatlist)):
        if flatlist[i][1] == 1:
            logger.info(str(flatlist[i][0])+ " was not copied.")


    # ---------------------------- Tests ------------------------------------- #

    # Check to see how many calibrations were copied.
    if expected_count - count == 0:
        logger.info("\nI sorted the " + str(expected_count) + " expected calibrations.\n")
    else:
        logger.info("\nI did not copy " + str(expected_count - count) + " calibration file(s).\n")

    # Check each calibration file list to see which ones were not copied.
    # Check that each file in flatlist was copied.
    for i in range(len(flatlist)):
        if flatlist[i][1] == 1:
            logger.info(str(flatlist[i][0])+ " from flatlist was not copied.")

    # Check that each file in flatdarklist was copied.
    for i in range(len(flatdarklist)):
        if flatdarklist[i][1] == 1:
            logger.info(str(flatdarklist[i][0])+ " from flatdarklist was not copied.")

    # Check that each file in ronchilist was copied.
    for i in range(len(ronchilist)):
        if ronchilist[i][1] == 1:
            logger.info(str(ronchilist[i][0])+ " from ronchilist was not copied.")

    # Check that each file in arclist was copied.
    for i in range(len(arclist)):
        if arclist[i][1] == 1:
            logger.info(str(arclist[i][0])+ " from arclist was not copied.")

    # Check that each file in arcdarklist was copied.
    for i in range(len(arcdarklist)):
        if arcdarklist[i][1] == 1:
            logger.info(str(arcdarklist[i][0])+ " from arcdarklist was not copied.")


    # Change back to original working directory.
    os.chdir(path1)

    # Check that each science directory exists and has associated calibration data.
    # Pseudocode (repeated below with actual code):
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
