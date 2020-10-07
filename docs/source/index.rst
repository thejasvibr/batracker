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

To check out all the stuff being prototyped check it out here - I'm always happy to hear any suggestions/advice!

.. toctree::
    :maxdepth: 4
    :caption: prototyping zone

    prototyping/index.rst

Head here to check out the current ideas and example implementations 
being tried out in `batracker`



Implemented to date
~~~~~~~~~~~~~~~~~~~

#. 2020-08-23 : The first start-to-end example in place. Sound sources, detected, matched and localised in this example.
#. Localisation with 4 channels when given the time-difference of arrivals 



Why this project
----------------
Acoustic localisation is a pretty widely used method in bioacoustics. It can be used to track birds migrating and calling as they fly by, bats as they echolocate, and 
even dolphins and whales in the ocean. Many of the acoustic tracking workflows available to date are either inhouse or commercial solutions, which means there aren't many open-source options to acoustically track animals. To my knowledge many of these solutions have also been built primarily with single animals in mind, where situations like call overlaps are rare. 

`batracker` has been conceived as an alternate tracking solution which is 1) built from the start to be open-source and encourage collaboration and 2) designed keeping the simplest (few channels, single animal, non-reverberant) and toughest tracking situations (multi-channel, multi-animal, reverberant) in mind.

Yes, this is kind of an attempt to re-invent the wheel, but to do so with as much of the techniques and code openly available. Inhouse scripts/packages work but just due to logistical/time constraints - there is no incentive to openly share the knowledge and document the methods involved. The project may be stable and work only for the specific workflow in a particular lab, and thus have a tiny community. Tiny communities aren't necessrily bad, but they also 1) increase the chance that a serious bug won't be caught 2) mean that someone else looking to do a similar thing won't really know where to start looking (and may have to start from scratch). For more advantages on being part of an open-source scientific community and why bigger collaborative software projects are the need of the moment in science, check out the experiences of `this author here <https://arxiv.org/abs/1301.7064>`_. 


Mature projects that may be of interest
---------------------------------------
I describe here a few open-source packages which may be of interest and fill similar needs. Please let me know if there are some more -- I'd love to build this list!

#. `PAMGUARD <https://www.pamguard.org/>`_
#. `Sound Finder <https://www.tandfonline.com/doi/abs/10.1080/09524622.2013.827588>`_


.. toctree::
        :maxdepth: 2
        :caption: Guide-throughs

        examples/index.rst

The example galleries above will help you understand the basic concepts and parameters used in various parts of the `batracker` package.

.. toctree::
    :maxdepth: 2
    :caption: batracker components + API

    other_rst/signal_detection.rst

    other_rst/correspondence_matching.rst

    other_rst/tdoa_estimation.rst

    other_rst/signal_localication.rst


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
