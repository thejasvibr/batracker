Correspondence matching (Detailed)
==================================


Identifying each candidate signal: final outputs
------------------------------------------------
Each candidate signal can be identified uniquely in a file by it's channel number, and its index position within the channel. 
(`Remember, in Python indexing starts with 0`). So, for instance, let's take this example candidate signal array:

.. code-block:: python

    [ [(0.1,0.3),   (0.5,0.8),   (0.9,1.2)],
      [(0.15,0.35), (0.23, 0.45)          ],
      [(0.2,0.9),   (.9,1.2),  (1.23,1.26)],
      [(0.1,0.7),   (1.,1.1),  (1.29,1.38)]  
    ]

Each candidate signal in the array above will have the following identities:

.. code-block:: python

    [ [ (0,0),   (0,1),   (0,2) ],
      [ (1,0),   (1,1)          ],
      [ (2,0),   (2,1),   (2,2) ],
      [ (3,0),   (3,1),    (3,2) ]  
    ]

The final outputs will be a 'set' of calls which have been matched, the set name is based on the ID of the call in the reference channel (`this may/not be a good idea...what if the ref. channel doesn't have good detections or the detection quality changes over time..`).

.. code-block:: python 

    matched_sounds = {(3,0):set([(0,0), (1,0), (2,0)]),
                      (3,1):set([(0,1), (1,1), (2,1)]),
                      (3,2):set([(0,2),       (2,2)])
                     }

Not all detections will necessarily be matched with another detection in all channels. Sometimes, this might be due to failure  of the matching algorithm and sometimes it might be due to the directionality  of the source sound. In the example above the sound :code:`(3,1)` doesn't have a match on channel 2, and this means that particular sound can't be localised even though it has been detected.


Technical details and open questions
------------------------------------
The 

Important variables
-------------------

#. Maximum distance from reference mic to other mics
#. Inter-sound interval
#. Sound duration 

Scenarios to handle
-------------------
The most important thing to check is if correspondence matching based on the 
max array distance works. It `will` work when 

What happens when the inter-sound interval drops below that of the max-array distance, 
and which rules work better then?

.. image:: ../_static/schematics_correspondence.png

Possible solutions 
------------------

#. Spectrogram cross-correlation?

#. Stereo image matching : Try and use established algorithms from image analysis/ structure-from-motion/etc. to generate matches across the different calls and channels using spectrogram representations. One of the issues is that these matching algorithms are built with rotational invariance in mind  -- which might mess things up.




