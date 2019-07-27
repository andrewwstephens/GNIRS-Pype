#!/usr/bin/env python
import logging


# ----------------------------------------------------------------------------------------------------------------------
# Convenience method

def getLogger(name):
    return logging.getLogger(name)

# ----------------------------------------------------------------------------------------------------------------------


def configure(logfile=None, filelevel='INFO', screenlevel='INFO'):
    """
    Configure file and console logging.

    DEBUG     Detailed information, typically of interest only when diagnosing problems.
    INFO      Confirmation that things are working as expected.
    WARNING   An indication that something unexpected happened, or indicative of some problem in the near future.
    ERROR     Due to a more serious problem, the software has not been able to perform some function.
    CRITICAL  A serious error, indicating that the program itself may be unable to continue running.
    """
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)  # set minimum threshold level for logger

    datefmt = '%Y-%m-%d %H:%M:%S'
    formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s', datefmt=datefmt)

    if logfile:  # Create file handler
        logfilehandler = logging.FileHandler(logfile)
        if filelevel.upper() == 'DEBUG':
            logfilehandler.setLevel(logging.DEBUG)
        elif filelevel.upper() == 'INFO':
            logfilehandler.setLevel(logging.INFO)
        elif filelevel.upper() == 'WARNING':
            logfilehandler.setLevel(logging.WARNING)
        elif filelevel.upper() == 'ERROR':
            logfilehandler.setLevel(logging.ERROR)
        elif filelevel.upper() == 'CRITICAL':
            logfilehandler.setLevel(logging.CRITICAL)
        else:
            print ('ERROR: Unknown log error level')
            logfilehandler.setLevel(logging.INFO)
        logfilehandler.setFormatter(formatter)
        logger.addHandler(logfilehandler)

    # Create console screen log handler
    consoleloghandler = logging.StreamHandler()
    formatter = logging.Formatter('%(name)-16s %(levelname)-8s %(message)s')
    if screenlevel.upper() == 'DEBUG':
        consoleloghandler.setLevel(logging.DEBUG)
    elif screenlevel.upper() == 'INFO':
        consoleloghandler.setLevel(logging.INFO)
    elif screenlevel.upper() == 'WARNING':
        consoleloghandler.setLevel(logging.WARNING)
    elif screenlevel.upper() == 'ERROR':
        consoleloghandler.setLevel(logging.ERROR)
    elif screenlevel.upper() == 'CRITICAL':
        consoleloghandler.setLevel(logging.CRITICAL)
    else:
        print ('ERROR: Unknown log error level')
        consoleloghandler.setLevel(logging.INFO)
    consoleloghandler.setFormatter(formatter)
    logger.addHandler(consoleloghandler)

    return logger

# ---------------------------------------------------------------------------------------------------------------------
