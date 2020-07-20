.. batracker documentation master file, created by
   sphinx-quickstart on Mon Jul 13 09:03:44 2020.

Welcome to batracker's documentation!
=====================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:


`batracker` is designed to help you track sounds with a modular design. 
There are two main parts of any kind of acoustic tracking 

#. Signal detection:
    This part deals with how the signal of interest is actually detected. Each type of signal to be detected
    may need special detection algorithms (eg. simple thresholds, neural networks), and this part of the package
    provides the API for detection only. 

.. toctree::
    :maxdepth: 1
    
    other_rst/signal_localication.rst

This part deals with the actual act of localising. Localising a signal in space requires calculating the 
time-difference of arrivals (TDOA) of the same signal across channels. The TDOAs can be calculated using 
various algorithms too. This part of the package deals with the act of calculating the TDOAs and generating
the source signal positions.




Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
