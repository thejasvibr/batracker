.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download_prototyping_implementing_frielander1987.py>`     to download the full example code
    .. rst-class:: sphx-glr-example-title

    .. _sphx_glr_prototyping_implementing_frielander1987.py:


Implementing Friedlander 1987
=============================
Friedlander 1987 [1] formulates the localisation problem as the source being the  
intersection point of multiple major axis lines. Each set of three mics
lies on an ellipsoid, where the major axis passes  through the source. 
The final source position is calculated as the point where all the major axes
of the ellipsoids intersect. See Figure 2. of the paper for a better idea. 



References
----------
[1] Friedlander, B. (1987). A passive localization algorithm and its accuracy analysis.
    IEEE Journal of Oceanic engineering, 12(1), 234-245.

TODO:
    
    #. Run tests with simulated data and check accuracy


.. code-block:: default

    import numpy as np 
    import scipy.spatial.distance as distance
    from scipy.linalg import circulant

Notation
:math:`x_{i} = [x,y,z]^{T}` : the position of the ith sensor
:math:`R_{is}` : distance between source and sensor i 
:math:`R_{io}` : distance between origin and sensor i 
:math:`r_{ij}= R_{is}-R_{js}`

Equation 7b in action
`j` is the reference sensor


.. code-block:: default


    def solve_friedlander1987(array_geometry, d, **kwargs):
        '''
        Parameters
        ----------
        array_geometry: np.array
            A Nx3 array with xyz coordinates of N mics.
            The reference microphone channel number (j) used in the range difference
            estimation must be explicitly stated.
        d:  np.array
            A N-1x1 np.array with the range_differences to the source. 
            All range_differences (:math:`r_{ij}`), are calculated by
            taking :math:`R_{is}-R_{js}`, where :math:`R_{is}` is the direct
            range from mic i to the source, and :math:`j` is the reference
            microphone.
        j : int >=0
            The channel number of the reference microphone.
            The first channel is the 0th.
        use_analytical : bool, optional
            Whether to use the analytical solution to find the source position
            or not. Defaults to True. If False, a least squares routine is 
            used to find the source position. The difference is likely to be 
            negligible for most cases.
        
    
        Returns
        -------
        x_real_world : list
            If the array_geometry consists of >4 mics, then the output  is 
            a list with a 3x1 np.array with the x,y,z positions of the source.
            If the array_geometry consists of 4 mics, then the solution is 
            a list with 2 candidate xyz source positions. 
        '''
    
        j = kwargs['j']
    
        num_channels, coods = array_geometry.shape
        if num_channels == 4:
            raise NotImplementedError('Friedlander 1987 algorithm not yet implemented for 4 channels! Please raise issue or submit a pull request:)!')
    
        if coods != 3:
            raise IndexError(f'The current array geometry is {coods} dimensional.Only 3d coordinates accepted')
    
        x_real_world = friedlander1987(array_geometry, j, d, 
                                               kwargs.get('use_analytical', True))
    
        return x_real_world
    
    

    def friedlander1987(mic_posns, j, rij, use_analytical=True):
        '''
    
    
        '''
        Sj = make_Sj(j, mic_posns)    
        muj = make_muj(j,rij,  mic_posns).reshape(-1,1)
        rhoj = rij.copy().reshape(-1,1)
        Dj = np.linalg.inv(np.diag(rhoj.flatten()))
    
        array_to_shift = np.concatenate((np.zeros(rhoj.size-1),
                                         np.array([1])                                     
                                       )).flatten()
        Z = circulant(array_to_shift)
        Mj = (np.eye(rhoj.size)-Z).dot(Dj) # equation 8a
   
        MjSj = Mj.dot(Sj)
        Mjmuj = Mj.dot(muj)
        if use_analytical:
            left_portion = np.linalg.inv(Sj.T.dot(Mj.T.dot(Mj.dot(Sj))))
            right_portion = Sj.T.dot(Mj.T.dot(Mj.dot(muj)))
            xs = left_portion.dot(right_portion)
        else:
            xs,resid, _,_ = np.linalg.lstsq(MjSj, Mjmuj)
        return xs
    
    


    def make_Sj(j, mic_posns):
        '''
        Parameters
        ----------
        j : int >=0
            The row number of the reference microphone
        mic_posns: np.array
            A Nmicx3 np.array
    
        Returns 
        -------
        Sj : np.array
            An Nchannels-1 x 3 array with :math:`\Delta` x,y,z descirbing the variable 
            :math:`S_{j}` in equation 7b. 
        '''
        num_channels, coods = mic_posns.shape
        if coods !=3:
            raise ValueError(f'Expected 3 coordinates, but got {coods} coordinates')
    
        Sj = np.zeros((num_channels-1, 3))
    
        valid_i_array = make_valid_i_array(j, num_channels)
    
        for rownum, i in enumerate(valid_i_array):
            Sj[rownum,:] = mic_posns[i,:] - mic_posns[j,:]
        return Sj

    def make_muj(j, rij, mic_posns):
        '''

        Parameters
        ----------
        j : int >=0
            The row number of the reference microphone
        rij : np.array
            Difference in range distances (:math:`R_{is}-R_{js}`)
        mic_posns: np.array
            A Nmicx3 np.array    

        Returns
        -------
        mu_j : np.array
            N-1 x 1 np.array with the :math:`mu_{j}`
    
        Notes
        -----
        I think there may be a typo in eq. 7c because the first row entry should bee
        :math:`{R_{1o}}^2 - -{R_{jo}}^2 - {r_{1j}^2}` rather than just 
        :math:`{R_{1o}} - -{R_{jo}}^2 - {r_{1j}^2}`. This is because :math:`\mu_{j}`
        is obtained from eq. 6, where the term is :math:`({R_{io}^2}-{R_{jo}^2})-{R_{ij}^2}`
        Here I"m going with my own interpretation of what is correct. 
    
        The positions given in mic_posns are assumed to be centred around the origin 
        (0,0,0).
    
        '''
        num_channels = mic_posns.shape[0]                           
        valid_i = make_valid_i_array(j, num_channels)
    
        R_sq_io = np.sum(mic_posns[valid_i,:]**2.0,1)
        R_sq_jo = np.sum(mic_posns[j,:]**2.0)
        r_sq_ij = rij**2
    
        mu_j = 0.5 *( R_sq_io - R_sq_jo - r_sq_ij)
    
        return mu_j
    

    def distance_to_point(x,y):
        return distance.euclidean(x,y)


    def make_valid_i_array(j, num_channels):
        '''
        Parameters
        ----------
        j:  0>int
            Index of reference sensor
        num_channels: int
            Number of sensors
    
        Returns
        -------
        valid_i_array: np.array
            num_channels-1 array with all channel numbers apart from j.        
        '''
        valid_i_array = np.arange(num_channels)
        value_to_delete = int(np.argwhere(valid_i_array==j))
        valid_i_array = np.delete(valid_i_array, value_to_delete)
        return valid_i_array
    
    if __name__=='__main__':
        
        import tacost
        from tacost import calculate_toa as ctoa

    
        mic_posns = np.array([[0,0,1],
                              [1,0,0],
                              [-1,0,0],
                              [0,1,0],
                              [1,1,0]
                              ])
    
    
        test_pos = np.random.choice(np.arange(-100,100,0.5),3)
        toa = ctoa.calculate_mic_arrival_times(test_pos,
                                               array_geometry=mic_posns)
        j = 2
        vsound = 338.0 
        valid_i = make_valid_i_array(j, mic_posns.shape[0])
        rij = (toa[valid_i]-toa[j])*vsound
    
        xs = solve_friedlander1987(mic_posns, rij, j=2)
        print(xs, test_pos)
    
    #    Sj = make_Sj(j, mic_posns)    
    #    # make muj
    #    muj = make_muj(j,rij,  mic_posns).reshape(-1,1)
    #    # make rhoj
    #    rhoj = rij.copy().reshape(-1,1)
    #
    #    # %% 
    #    from scipy.linalg import circulant
    #
    #
    #    Dj = np.linalg.inv(np.diag(rhoj.flatten()))
    #    
    #    array_to_shift = np.concatenate((np.zeros(rhoj.size-1),
    #                                     np.array([1])                                     
    #                                   )).flatten()
    #    
    #    Z = circulant(array_to_shift)
    #    Mj = (np.eye(rhoj.size)-Z).dot(Dj) # equation 8a
    #
    #    
    #    MjSj = Mj.dot(Sj)
    #    Mjmuj = Mj.dot(muj)
    #    # %%
    #    #     
    #    xs,resid, _,_ = np.linalg.lstsq(MjSj, Mjmuj)
    #    
    #    left_portion = np.linalg.inv(Sj.T.dot(Mj.T.dot(Mj.dot(Sj))))
    #    right_portion = Sj.T.dot(Mj.T.dot(Mj.dot(muj)))
    #    xs_2 = left_portion.dot(right_portion)
    #    print( xs,test_pos, xs_2)
    #
    #
    #    xyz = friedlander1987(mic_posns, j, rij)
    #    print(xyz)


.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 0 minutes  0.000 seconds)


.. _sphx_glr_download_prototyping_implementing_frielander1987.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: implementing_frielander1987.py <implementing_frielander1987.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: implementing_frielander1987.ipynb <implementing_frielander1987.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
