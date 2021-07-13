.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download_prototyping_plot_correspondence_matching.py>`     to download the full example code
    .. rst-class:: sphx-glr-example-title

    .. _sphx_glr_prototyping_plot_correspondence_matching.py:


Correspondence matching prototyping
===================================
Implements the most basic type of correspondence matching -  which matches
signals across all channels to a single source emission. 

In this prototyping example, the simplest algorithm is used to find correspondences. 
The maximum inter-channel delay is calculated based on the longest dimension of the 
array, from the reference mic to the focal mic. 



Possible approaches



.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /prototyping/images/sphx_glr_plot_correspondence_matching_001.png
          :alt: plot correspondence matching
          :class: sphx-glr-multi-img

    *

      .. image:: /prototyping/images/sphx_glr_plot_correspondence_matching_002.png
          :alt: plot correspondence matching
          :class: sphx-glr-multi-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    4 125000
      0%|          | 0/4 [00:00<?, ?it/s]     25%|##5       | 1/4 [00:02<00:07,  2.40s/it]     50%|#####     | 2/4 [00:04<00:04,  2.40s/it]     75%|#######5  | 3/4 [00:07<00:02,  2.41s/it]    100%|##########| 4/4 [00:09<00:00,  2.43s/it]    100%|##########| 4/4 [00:09<00:00,  2.43s/it]
    /home/autumn/Documents/research-repos/batracker/batracker/tests/prototyping/plot_correspondence_matching.py:194: MatplotlibDeprecationWarning: Adding an axes using the same arguments as a previous axes currently reuses the earlier instance.  In a future version, a new instance will always be created and returned.  Meanwhile, this warning can be suppressed, and the future behavior ensured, by passing a unique label to each axes instance.
      plt.subplot(num_channels*100  + 10+each+1, sharex = ax0)






|


.. code-block:: default

    import matplotlib.pyplot as plt
    plt.rcParams['agg.path.chunksize'] = 10000
    import numpy as np 
    import pandas as pd
    import scipy.io.wavfile as wavfile
    import scipy.spatial as spatial

    from plot_threshold_detector import cross_channel_threshold_detector


    def calc_intermic_distances(array_geom):
        '''
        Parameters
        ----------
        array_geom: pd.DataFrame
            With compulsory column names, 'x', 'y', and 'z'
            Each mic's coordinates must be placed in one row, 
            in the order matching the recorded channels of
            the audio file. 
        Returns
        -------
        mic2mic_distance : np.array
            Nmic x Nmix array, a distance matrix
    
        '''
        mics_xyz = array_geom[['x','y','z']].to_numpy()
        mic2mic_distance = spatial.distance_matrix(mics_xyz, mics_xyz)
        return mic2mic_distance

    def largest_mic2mic_distance(array_geom):
        '''
        Provides the largest mic-2-mic distance in the array. 
        In the presence of multiple equal mic-2-mic distances, 
        provides the first one. 
    
        Parameters
        ----------
        array_geom: pd.DataFrame
            With compulsory column names, 'x', 'y', and 'z'
            Each mic's coordinates must be placed in one row, 
            in the order matching the recorded channels of
            the audio file. 
    
        Returns
        -------
        maxval_location : np.array
            Array with 2 ints, with the row and column indices. 
            This corresponds to the mic pair with the largest mic-2-mic distance.
        maxval : float
            The largest mic-2-mic distance found. 
    
        See Also
        --------
        calc_intermic_distances
    
        '''
        mic2mic_distances = calc_intermic_distances(array_geom)
        maxval = np.max(mic2mic_distances)
        maxval_location = np.argwhere(maxval==mic2mic_distances)[0]
        return maxval_location, maxval
    

    def calculate_crosscorr_boundary(mic_pair, max_dist, detections, **kwargs):
        '''
        Generates the common audio boundaries which will be used to perform
        the cross-correlation 
    
        This function does not inspect the audio content of the 
        segment that it outputs across channels. 
    
        Parameters
        ----------
        mic_pair : list/array-like
            With 2 entries to indicate the index numbers of the mics with 
            the largest distance.s
        max_dist : float
            Distance in meters
        detections: list
            List with tuples.
        vsound : float, optional
            Defaults to 338 m/s
    
        Returns
        -------
        crosscor_boundaries
    
        Note
        ----
        This function assumes there are an equal number of sounds detected across 
        all channels. An error will be thrown if there are inequal number of detections.
        '''
        vsound = kwargs.get('vsound', 338.0)
        num_detections = np.array([len(each) for each in detections])
        if not np.all(num_detections==num_detections[0]):
            raise ValueError(f'All channels do not have an equal number of detections {num_detections}')
    
        time_delay = max_dist/vsound

        crosscor_boundaries = []

        for detection_num, detection_times in enumerate(detections[0]):
            earliest_start = np.min([ detections[mic_pair[0]][detection_num][0],
                                      detections[mic_pair[1]][detection_num][0]
                                     ])
            latest_stop = np.max([ detections[mic_pair[0]][detection_num][1],
                                   detections[mic_pair[1]][detection_num][1]
                                     ])
            crosscor_boundaries.append( (earliest_start-time_delay,
                                         latest_stop+time_delay)
                                      )
        return crosscor_boundaries   


    def generate_crosscor_boundaries(detections, array_geom, **kwargs):
        '''Wrapper function that calculates the largest mic-2-mic distance
        and proceeds to generate the common boundaries for cross-correlation. 
    
        See Also
        --------
        calculate_crosscor_boundary
        largest_mic2mic_distance
    
        '''
        mic_pair, max_dist  = largest_mic2mic_distance(array_geom)
        crosscor_boundary = calculate_crosscorr_boundary(mic_pair, max_dist, 
                                                         detections, **kwargs)
        return crosscor_boundary



    if __name__ == '__main__':

        # %%
        # Load the simulated audio file that was previously made
    
        fs, audio = wavfile.read('simulated_audio/batracker_simple.wav')
        short_audio = audio[:int(fs*0.5),:]
        multi_detections = cross_channel_threshold_detector(short_audio, fs,
                                                            dbrms_window=0.5*10**-3,
                                                            dbrms_threshold=-56)
    
        # %% 
        # Also need some kind of visualiser to see the output of the all the detections 
        # and diagnose problems.
    
        num_channels = audio.shape[1]
        all_axes = []
        plt.figure(figsize=(10,8))
        ax0 = plt.subplot(num_channels*100  + 11)
        plt.specgram(short_audio[:, 0], Fs=fs, NFFT=256, noverlap=255)
    
        for each in range(1,num_channels):
            plt.subplot(num_channels*100  + 10+each+1, sharex = ax0)
            plt.specgram(short_audio[:, each], Fs=fs, NFFT=256, noverlap=255)
        
            for every in multi_detections[each]:
                plt.vlines(every, 0, fs*0.5, linewidth=0.5, label='Detections')
        plt.legend()
    
        # %% 
        # Matching detections to each other. 
        # The sounds are matched together according to the furthest microphone distances. 
        # eg. if the mics 1 and 4 are about 2 m apart, then we expect all sounds
        # to arrive within :math:`\pm`6 ms of each other from on
        array_geometry = pd.read_csv('simulated_audio/array_geom.csv')
    
    
        loc,     max_dist  = largest_mic2mic_distance(array_geometry)
    
        # %% 
        # Now, find all sounds within +/- max delay distance of that mic pair, and
        # put them into one boundary region for cross-correlation. This defines the 
        # cross-correlation boundaries

        crosscor_boundaries = generate_crosscor_boundaries(multi_detections,array_geometry)
    
        for each in range(1,num_channels):
            plt.subplot(num_channels*100  + 10+each+1, sharex = ax0)    
            for every in crosscor_boundaries:
                plt.vlines(every, 0, fs*0.5,'r', linewidth=0.5, )

        # show detections and cross-corr boundaries
        plt.figure(figsize=(10,8))
        ax1 = plt.subplot(num_channels*100  + 11)
        plt.plot(short_audio[:, 0])
    
        for each in range(1,num_channels):
            plt.subplot(num_channels*100  + 10+each+1, sharex = ax1, sharey=ax1)
            plt.plot(short_audio[:, each])
        
            for every_det in multi_detections[each]:
                det_samples = np.array(every_det)*fs
                plt.vlines(det_samples, np.min(audio), np.max(audio), 'k', linewidth=0.75)
            for every_cc_bound in  crosscor_boundaries:
                bound_samples = np.array(every_cc_bound)*fs
                plt.vlines(bound_samples, np.min(audio), np.max(audio),'r', linewidth=0.5)
            
    
    
















.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 0 minutes  22.277 seconds)


.. _sphx_glr_download_prototyping_plot_correspondence_matching.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: plot_correspondence_matching.py <plot_correspondence_matching.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: plot_correspondence_matching.ipynb <plot_correspondence_matching.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
