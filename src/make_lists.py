#!/usr/bin/env python
"""
Create file lists
"""
import collections
import ConfigParser
import header
import log
import obslog
import os


# ----------------------------------------------------------------------------------------------------------------------
def start(configfile):
    """
    Create file lists
    :param configfile: configuration file [string]
    """
    logger = log.getLogger('make_lists')

    config = ConfigParser.RawConfigParser()
    config.optionxform = str  # make config file options case-sensitive
    config.read(configfile)

    rawpath = config.get('DownloadData', 'RawDataDir')
    scipath = config.items("ScienceDirectories")
    telpath = config.items("TelluricDirectories")
    calpath = config.items("CalibrationDirectories")

    for path, process in scipath + telpath:                                                        # Science & Tellurics
        logger.debug('%s = %s', path, process)

        if not process:
            logger.debug('Skipping %s', path)
            continue

        logger.info('%s', path)

        info = header.info(path + '/Intermediate')

        if not os.path.exists(path + '/Intermediate/obslog.csv'):
            pathparts = path.split('_')
            obsid = pathparts[-1]
            date = pathparts[-7]
            logger.debug('obsid: %s', obsid)
            logger.debug('date: %s', date)
            obslog.writecsv(obsid=obsid, date=date, rawpath=rawpath, output=path + '/Intermediate/obslog.csv')
        olog = obslog.readcsv(path + '/Intermediate/obslog.csv')

        all_exptimes = []
        for f in info.keys():
            all_exptimes.append(info[f]['EXPTIME'])
        logger.debug('all_exptimes: %s', all_exptimes)

        if len(list(set(all_exptimes))) == 1:
            exptime = all_exptimes[0]
        else:
            logger.warning('Multiple exposure times: %s', all_exptimes)
            freq = collections.Counter(all_exptimes)
            logger.debug('Exposure time frequency: %s', freq)
            exptime, num = freq.most_common(1)[0]
            logger.debug('Most common: %s  number: %s', exptime, num)
            logger.debug('freq.values: %s', freq.values())
            num_peaks = freq.values().count(num)
            logger.debug('npeak: %s', num_peaks)
            if num_peaks == 1:
                logger.info('Only using files with the most common exposure time: %.2f sec', exptime)
            else:
                logger.error('The most common value occurs %d times.', num_peaks)
                logger.error('You will need to fix this.')
                raise SystemExit

        with    open(path + '/Intermediate/all.list', 'w') as all, \
                open(path + '/Intermediate/src.list', 'w') as src, \
                open(path + '/Intermediate/sky.list', 'w') as sky:

            for f in sorted(info.keys()):
                if info[f]['EXPTIME'] == exptime:
                    all.write(f + '\n')
                    if inslit(slit=info[f]['SLIT'], decker=info[f]['DECKER'], p=olog[f]['P'], q=olog[f]['Q']):
                        src.write(f + '\n')
                    else:
                        sky.write(f + '\n')

        all_offsets = []
        for f in info.keys():
            all_offsets.append(float(olog[f]['Q']))
        logger.debug('all_offsets: %s', all_offsets)
        unique_offsets = sorted(list(set(all_offsets)), key=abs)  # sort by abs value so nodA will be min(abs(offset))
        logger.debug('Unique offsets: %s', unique_offsets)

        if len(unique_offsets) == 2:
            with    open(path + '/Intermediate/nodA.list', 'w') as nodA, \
                    open(path + '/Intermediate/nodB.list', 'w') as nodB:

                for f in sorted(info.keys()):
                    if info[f]['EXPTIME'] == exptime:
                        if float(olog[f]['Q']) == unique_offsets[0]:
                            nodA.write(f + '\n')
                        else:
                            nodB.write(f + '\n')
        else:
            logger.error('There are %d offset positions.  No nod lists created.', len(unique_offsets))

    # ------------------------------------------------------------------------------------------------------------------
    for path, process in calpath:                                                                         # Calibrations
        logger.debug('%s = %s', path, process)

        if not process:
            logger.debug('Skipping %s', path)
            continue

        logger.info('%s', path)

        info = header.info(path)

        with    open(path + '/all.list', 'w') as all, \
                open(path + '/arcs.list', 'w') as arcs, \
                open(path + '/IRflats.list', 'w') as IRflats, \
                open(path + '/QHflats.list', 'w') as QHflats, \
                open(path + '/pinholes.list', 'w') as pinholes:

            for f in sorted(info.keys()):

                all.write(f + '\n')

                if info[f]['OBSTYPE'] == 'ARC':
                    arcs.write(f + '\n')

                elif 'Pinholes' in info[f]['SLIT']:
                    pinholes.write(f + '\n')

                elif info[f]['OBSTYPE'] == 'FLAT' and \
                        info[f]['GCALLAMP'] == 'QH' and \
                        'Pinholes' not in info[f]['SLIT']:
                    QHflats.write(f + '\n')

                elif info[f]['OBSTYPE'] == 'FLAT' and \
                        info[f]['GCALLAMP'] == 'IRhigh':
                    IRflats.write(f + '\n')

    return


# ----------------------------------------------------------------------------------------------------------------------
def inslit(slit, decker, p, q):
    """
    Determine whether an exposure is in the slit based on the slit size and offsets.
    :param slit: Slit name [string]
    :param decker: Decker name [string]
    :param p: Absolute P offset (arcseconds) [float or string]
    :param q: Absolute Q offset (arcseconds) [float or string]
    :return: boolean
    """
    logger = log.getLogger('inslit')
    p = float(p)
    q = float(q)
    logger.debug('Offsets:  P=%.2f  Q=%.2f arcsec', p, q)
    if 'SC' in decker and 'XD' in decker:
        slitlength = 7.0
    elif 'LC' in decker and 'XD' in decker:
        slitlength = 5.0
    else:
        logger.error('cannot calculate slit length')
        raise SystemExit
    slitwidth = float(slit[:slit.find('arcsec')])
    logger.debug('Slit size: %.1f x %.1f arcsec', slitwidth, slitlength)
    if abs(p) < slitwidth / 2.0 and abs(q) < slitlength / 2.0:
        answer = True
    else:
        answer = False
    logger.debug('-> %s', answer)
    return answer


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    log.configure('gnirs-pype.log', filelevel='INFO', screenlevel='DEBUG')
    start('gnirs-pype.cfg')
