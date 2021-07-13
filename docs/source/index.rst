.. batracker documentation master file, created by
   sphinx-quickstart on Mon Jul 13 09:03:44 2020.

Welcome to batracker's documentation!
=====================================

.. toctree::
   :maxdepth: 3
   :caption: Contents:


`batracker` is designed to help you track sounds with a modular design. 

Why this project
----------------
Acoustic localisation is a pretty widely used method in bioacoustics. Many of the acoustic tracking workflows available to date are either inhouse or commercial solutions, which means there aren't many open-source options to acoustically track animals. To my knowledge many of these solutions have also been built primarily with single animals in mind, where situations like call overlaps are rare. Another common (and commonsensical) assumption in many of these packages is reverberation-free audio. Especially in recordings of groups of echolocators in the field both assumptions (no/rare overlaps, no reverberation) fail. `batracker` is an attempt at pushing the state of acoustic tracking workflows in bioacoustics.

`batracker` has been conceived as an alternate solution for acoustic tracking which is:

 1) built from the start to be open-source and encourage collaboration

 and

 2) designed keeping the simplest (few channels, single animal, non-reverberant) and toughest tracking situations (multi-channel, multi-animal, reverberant) in mind.

Yes, this is an attempt to 're-invent the wheel', but to do so with the techniques and code openly documented and tested. I stress, there are acoustic tracking solutions out there. Inhouse scripts/packages work, but due to logistical/time constraints of typical scientific software development - there is no incentive to share the knowledge and document the methods involved.

A project may be stable and work only for the specific workflow in a particular lab, and thus have a tiny community. Tiny communities aren't necessrily bad, but they also 1) increase the chance that a serious bug won't be caught 2) mean that someone else looking to do a similar thing won't really know where to start looking (and may have to start from scratch). For more advantages on being part of an open-source scientific community and why bigger collaborative software projects are the need of the moment in science, check out the experiences of `this author here <https://arxiv.org/abs/1301.7064>`_. 


Implemented to date
~~~~~~~~~~~~~~~~~~~

#. 2021-07-12 : Prototyping of DATEMM (disambiguation of TDOA estimation for multiple sources)...promising results...
#. 2020-08-23 : The first start-to-end example in place. Sound sources, detected, matched and localised in this example.
#. Localisation with 4 channels when given the time-difference of arrivals 


`Project under construction !!`
-------------------------------
This project is currently under construction.
 
For now, this page is a place to dump my thoughts and organise the project structure as I implement things. I'm currently devoting a limited time per month on this project, and am looking to formally find funding to push this project forward majorly. If you're interested in collaborating on this repo or funding this project do send me a message at thejasvib@gmail.com. 

To check out all the stuff being prototyped check it out here - I'm always happy to hear any suggestions/advice!

.. toctree::
    :maxdepth: 4
    :caption: prototyping zone

    prototyping/index.rst


Mature open-source projects of interest
---------------------------------------
I describe here a few open-source packages which may be of interest and fill similar needs. Please let me know if there are some more -- I'd love to build this list!

#. `PAMGUARD <https://www.pamguard.org/>`_
#. `Sound Finder <https://www.tandfonline.com/doi/abs/10.1080/09524622.2013.827588>`_


`batracker` components + API
----------------------------

.. toctree::
    :maxdepth: 2


    other_rst/signal_detection.rst

    other_rst/correspondence_matching.rst

    other_rst/tdoa_estimation.rst

    other_rst/signal_localication.rst


