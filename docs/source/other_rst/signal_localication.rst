Signal localisation (Detailed)
==============================

Let's begin with a simple example. Assuming signal detection has already been performed, here we'll localise a single set of detections for one signal 
across all the channels. 

.. code-block:: python

    from batracker import localiser
    
    # mic_array_positions is the array geometry data
    
    position = localiser(detections, mic_array_positions, v_sound=338)
    print(position)

TDOA based methods
------------------
Time difference of arrival (TDOA) based methods rely on calculating source position from the range differences between microphone pairs. Some of these formulations include:

#.  Spherical Intersection - `SX` method (Schau & Robinson 1987): formulated for 4 microphones only 
    Status:
        > `Implemented, tests to be written`

#. Spherical Interpolation - `SI` method (Smith & Abel 1987)

#. Friedlander 1987: formulated for >= 4 microphones
    Statue:
        > `Not yet  implemented`


Non-TDOA based methods
----------------------
There are whole bunch of methods which don't necessarily use the time difference of arrival, but instead use the direction-of-arrival, or phase differences between
sounds recorded across mics. These methods include (unsure of exact implementation details...):

#. MUSIC
#. Beamforming 





signal localisation API
~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: batracker.localisation.localiser
   :members:

.. automodule:: batracker.localisation.schau_robinson_1987
   :members:

.. automodule:: batracker.localisation.friedlander_1987
   :members:



