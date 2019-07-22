#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ConfigParser
import datetime
import log
import os


def start(configfile):
    """
    Check that each science and Telluric directory has associated calibration data:
    - science and telluric directories have a src.list
    - a Calibrations directory is present
    - flats.list exists and has more than one file
    - arcs.list exists and has more than one file
    - pinhole.list exists and has at least one file
    ? Should we also check that the science directories have valid Telluric standards ?
    """

    logger = log.getLogger('checkcals')
    logger.info('Checking that each observation has the required calibrations.')

    config = ConfigParser.RawConfigParser()
    config.optionxform = str  # make options case-sensitive
    config.read(configfile)
    scidirs = config.options('ScienceDirectories')
    logger.debug('Science directories: %s', scidirs)
    teldirs = config.options('TelluricDirectories')
    logger.debug('Telluric directories: %s', scidirs)


    for dir in scidirs:

        if config.getboolean('ScienceDirectories', dir):  # only check directories marked True
            logger.info('%s:', dir)

            contents = os.listdir(dir)
            logger.debug('contents: %s', contents)

            if 'src.list' not in contents:
                logger.warning('No src.list found in %s', dir)

            twoup = dir[:dir[:dir.rfind('/')].rfind('/')]
            logger.debug('twoup: %s', twoup)
            contents = os.listdir(twoup)
            logger.debug('contents: %s', contents)

            if 'Calibrations' not in contents:
                logger.warning('No Calibrations directory found for %s', dir)

            if 'Tellurics' not in contents:
                logger.warning('No Tellurics directory found for %s', dir)

            calibrations = twoup + '/Calibrations'
            contents = os.listdir(calibrations)
            logger.debug('contents: %s', contents)

            for c in ['arcs.list', 'darks.list', 'QHflats.list', 'IRflats.list', 'pinholes.list']:
                if c not in contents:
                    logger.warning('No %s found in %s', c, calibrations)


    raise SystemExit('EXIT - THIS FUNCTION IS A WORK IN PROGRESS')





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


# ----------------------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    log.configure('gnirs.log', filelevel='INFO', screenlevel='DEBUG')
    start('gnirs.cfg')
