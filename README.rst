GNIRS-Pype
============
#.. image:: https://badge.fury.io/py/nifty4gemini.svg
#   :target: https://badge.fury.io/py/nifty4gemini
.. image:: https://readthedocs.org/projects/gnirs-pype/badge/?version=latest
   :alt: GNIRS-Pype's documentation, hosted on ReadtheDocs.
   :target: http://gnirs-pype.readthedocs.io/en/latest/
.. image:: https://zenodo.org/badge/103719389.svg
   :alt: DOI
   :target: https://zenodo.org/badge/latestdoi/103719389
.. image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :alt: MIT license.
   :target: https://opensource.org/licenses/MIT
.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
   :alt: GNIRS-Pype uses Astropy! Here is a link to the project webpage:
   :target: http://www.astropy.org/

A Python Data Reduction Pipeline for the Gemini Near-InfraRed Spectrograph (GNIRS).

Full documentation: `ReadTheDocs <http://gnirs-pype.readthedocs.io/en/latest/>`_.

This is a new data reduction Python pipeline that uses Astroconda and the Gemini
IRAF Package to reduce GNIRS data. It offers a complete data reduction process from
sorting the data to producing a final flux calibrated and wavelength calibrated
combined spectrum with the full S/N for a science target.

This pipeline is open source and it is supported via the `Gemini Data Reduction User Forum <http://drforum.gemini.edu/>`_.

Any feedback and comments (astephens@gemini.edu) are welcome!

Copyright
---------

For more details, please read the LICENSE.


How to Submit Bugs and Requests
-------------------------------

Very important: **do not submit a Gemini help desk ticket!**

If you want to report a problem, use the `Gemini Data Reduction Forum thread <http://drforum.gemini.edu/topic/nifs-python-data-reduction-pipeline/>`_
or create an issue in this repo.

Installation
============

Pre-Requisites
--------------
Make sure you have the latest version of Gemini Astroconda installed, have activated an Astroconda environment and have set up PYRAF.
You can find instructions for installing Astroconda `here <https://astroconda.readthedocs.io/en/latest/>`_. PYRAF can be set up by running the mkiraf command
in your "~/iraf" directory.

Installing
----------


Installing in Editable Mode
---------------------------


Quick Start
===========


License
=======

See the LICENSE file in the current directory. Note that downloadFromGeminiPublicArchive does not use the MIT
license. Refer to the LICENCE file in the downloadFromGeminiPublicArchive directory to view the appropriate license.
