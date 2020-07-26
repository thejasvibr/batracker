Correspondence matching (Detailed)
==================================


Identifying each candidate signal 
---------------------------------
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
      [ (3,0),   (3,1,    (3,2) ]  
    ]


