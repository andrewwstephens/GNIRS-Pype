#!/usr/bin/env python
import log
import os


# ----------------------------------------------------------------------------------------------------------------------
def boxit(text, character, center=True):
    lines = text.split('\n')
    max_length = max([len(line) for line in lines])
    border = character * (max_length + 4)
    print("%s" % border)

    for line in lines:
        dif = (max_length-len(line))
        if dif == 0:
            print(character + " %s " % line + character)
        else:
            if center:
                print(character + ' ' * (dif/2) + ' %s ' % line + ' ' * (dif/2) + character)
            else:
                print(character + " %s " % line + ' ' * dif + character)

    print("%s" % border)
    return


# ----------------------------------------------------------------------------------------------------------------------
def exists(inlist, overwrite):
    """
    Check for the existence of the files in the input list.
    If overwrite == True then delete any existing files and return False
    if overwrite == False:
    - if no files exist then return False
    - if all the files exist then return True, indicating that this reduction step has been done and can be skipped
    - if some (but not all) of the files exist then throw an error because something failed somewhere
    The goal is to be able to say something like:  if exists() then skip this step, or if not exists() do stuff.
    :param inlist: list of files
    :param overwrite: boolean
    :return: boolean
    """
    logger = log.getLogger('exists')

    exist = [os.path.exists(f) for f in inlist]
    logger.debug('Exist: %s', exist)

    if overwrite:
        for i in range(len(inlist)):
            if exist[i]:
                logger.warning('Removing the old %s', inlist[i])
                os.remove(inlist[i])
        return False

    else:
        if all(exist):
            logger.debug("All files exist")
            return True

        elif not any(exist):
            logger.debug('No files exist')
            return False

        else:
            # We're stuck.  We don't have all the files required to skip this step,
            # but we can't erase files and try again.  Give up and let the user figure it out.
            logger.error('Some (but not all) the files exist and overwrite = False')
            raise SystemExit


# ----------------------------------------------------------------------------------------------------------------------
def get_orders(path):
    logger = log.getLogger('get_orders')
    if 'LB_SXD' in path:
        orders = [3, 4, 5]
    elif 'LB_LXD' in path:
        orders = [3, 4, 5, 6, 7, 8]
    elif 'SB_SXD' in path:
        orders = [3, 4, 5, 6, 7, 8]
    else:
        raise SystemExit("Unknown GNIRS configuration.  Cannot determine the XD orders.")
    logger.debug('orders: %s', orders)
    return orders


# ----------------------------------------------------------------------------------------------------------------------
def make_list(prefix, orders=None, suffix='_MEF[1]', regions=None):
    """
    Take a filename and return a list of all the orders (filename_orderN)
    # If the regions dictionary is specified then append the appropriate order, e.g. filename_orderN[xxx:yyy]
    # The orders are only needed if the regions are not specified in which case we don't know the orders.
    :param prefix: string filename prefix
    :param orders: list of spectral orders, e.g. [3,4,5,6,7,8]
    :param suffix: file name suffix
    :param regions: dictionary of all orders and corresponding regions, e.g. {3:"123:456"}
    :return: list of files for each order with optional regions, e.g. [..., prefix_orderN_suffix[region], ...]
    """
    logger = log.getLogger('make_list')
    logger.debug('prefix: %s', prefix)
    logger.debug('orders: %s', orders)
    logger.debug('suffix: %s', suffix)
    logger.debug('regions: %s', regions)

    if not orders and not regions:
        logger.error('Must specify orders or regions')
        raise SystemExit

    output = []
    if regions:
        for o in sorted(regions.keys()):
            output.append('%s_order%d%s[%s]' % (prefix, o, suffix, regions[o]))
    else:
        for o in orders:
            output.append('%s_order%d%s' % (prefix, o, suffix))

    logger.debug('list: %s', output)
    return output


# ----------------------------------------------------------------------------------------------------------------------
def nofits(filename):
    return filename.replace('.fits', '')


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    log.configure('gnirs.log', filelevel='INFO', screenlevel='DEBUG')
    boxit('Some text that needs a star box around it', '*')
    boxit('Some text that needs a line box around it', '-')
    boxit('A long centered multi-line message\nthat needs to wrap.', '#')
