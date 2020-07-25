Signal localisation (Detailed)
==============================

Let's begin with a simple example. Assuming signal detection has already been performed, here we'll localise a single set of detections for one signal 
across all the channels. 

.. code-block:: python

    from batracker import localiser
    
    # mic_array_positions is the array geometry data
    
    position = localiser(detections, mic_array_positions, v_sound=338)
    print(position)





signal localisation API
~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: batracker.localisation.localiser
   :members:

.. automodule:: batracker.localisation.schau_robinson_1987
   :members:



