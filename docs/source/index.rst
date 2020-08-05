.. batracker documentation master file, created by
   sphinx-quickstart on Mon Jul 13 09:03:44 2020.

Welcome to batracker's documentation!
=====================================

.. toctree::
   :maxdepth: 3
   :caption: Contents:


`batracker` is designed to help you track sounds with a modular design. 
There are four main parts of any kind of acoustic tracking:

`Project under construction !!`
-------------------------------
This project is currently under construction, and most of the functions are currently in prototyping stage or in `thought-experiment` stage.
This page is a place to dump my thoughts and organise the project structure as I implement things. I'm currently devoting a limited time per day 
on this project, and am looking to formally find funding to push this project forward majorly. If you're interested
in collaborating on this repo or funding this project do send me a message at thejasvib@gmail.com. 


Why this project
----------------
Acoustic localisation is a pretty widely used method in bioacoustics. It can be used to track birds migrating and calling as they fly by, bats as they echolocate, and 
even dolphins and whales in the ocean. Many of the acoustic tracking workflows available to date are either inhouse or commercial solutions, which means there aren't many open-source options to acoustically track animals. To my knowledge many of these solutions have also been built primarily with single animals in mind, where situations like call overlaps are rare. 

`batracker` has been conceived as an alternate tracking solution which is 1) built from the start to be open-source and encourage collaboration and 2) designed keeping the simplest (few channels, single animal, non-reverberant) and toughest tracking situations (multi-channe, multi-animal, reverberant) in mind.


Mature projects that may be of interest
---------------------------------------
I describe here a few open-source packages which may be of interest and fill similar needs. Please let me know if there are some more -- I'd love to build this list!

#. `PAMGUARD <https://www.pamguard.org/>`_
#. `Sound Finder <https://www.tandfonline.com/doi/abs/10.1080/09524622.2013.827588>`_


Signal detection
~~~~~~~~~~~~~~~~

This part deals with how the signal of interest is actually detected. Each type of signal to be detected
may need special detection algorithms (eg. simple thresholds, neural networks), and this part of the package
provides the API for detection only. 

.. toctree::
    :maxdepth: 1

    other_rst/signal_detection.rst

Correspondence matching
~~~~~~~~~~~~~~~~~~~~~~~ 
Having detected all signals within a channel, now each signal in a channel must be matched with another one from other channels. 
Here there are a couple of algorithms that can be used, based on how dense the sounds are across the channels. The most naive algorithm
assumes that sounds of interest are emitted with long gaps of silence. Each sound in a channel is matched with the sound in other channels
that is :math:`\pm \Delta T`, where :math:`\Delta T` is the maximum possible time difference of arrival between one channel to another. :math:`\Delta T`
corresponds to the largest inter-microphone distance in an array.

.. toctree::
    :maxdepth: 1

    other_rst/correspondence_matching.rst

Time delay estimation
~~~~~~~~~~~~~~~~~~~~~ 
All matched signals are compared and their time-difference of arrivals (TDOAs) are measured using various methods.  
The simplest method is of course to perform a cross-correlation, and the difference in the position of the peaks 
is the time-delay of arrival. 

.. toctree::
    :maxdepth: 1

    other_rst/tdoa_estimation.rst

Localisation
~~~~~~~~~~~~
With the TDOAs in hand, now the position of the source can be calculated. According the algorithm at hand the source position is calculated using 
all the available (independent) TDOAs or only a subset of them. 

.. toctree::
    :maxdepth: 1
    
    other_rst/signal_localication.rst

This part of the package deals with the act of using the TDOAs and generating
the source signal positions.


Implemented to date
~~~~~~~~~~~~~~~~~~~

#. Localisation with 4 channels when given the time-difference of arrivals 
#. 



Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
