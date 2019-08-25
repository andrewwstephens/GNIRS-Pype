#!/usr/bin/env python
import log


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
# We need a function that checks if a file exists,
# and deletes it if overwrite = True


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    log.configure('gnirs.log', filelevel='INFO', screenlevel='DEBUG')
    boxit('Some text that needs a star box around it', '*')
    boxit('Some text that needs a line box around it', '-')
    boxit('A long centered multi-line message\nthat needs to wrap.', '#')
