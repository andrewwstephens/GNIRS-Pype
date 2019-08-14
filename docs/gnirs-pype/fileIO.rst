Example of GNIRS-Pype File I/O
==============================

Note: This is for v1.0.0. It is moderately correct up telluric correction.

This is an example of how the GNIRS-Pype directory tree appears after each step of the
data reduction. These directory trees were created using a custom **gnirstree** bash command:

.. code-block:: text

  find . -name .git -prune -o -print | sed -e 's;[^/]*/;|____;g;s;____|; |;g'

Add the following line to your ~/.bash_profile to create the **gnirstree** alias:

.. code-block:: text

  alias gnirstree="find . -name .git -prune -o -print | sed -e 's;[^/]*/;|____;g;s;____|; |;g'"

Example Data Reduction Products:
--------------------------------

.. code-block:: text

  This is currently under developmemt.


.. placeholder
