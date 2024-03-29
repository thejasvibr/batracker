.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download_prototyping_plot_threshold_detector.py>`     to download the full example code
    .. rst-class:: sphx-glr-example-title

    .. _sphx_glr_prototyping_plot_threshold_detector.py:


Threshold detector protoyping
=============================
The simplest method. Calculates the dbrms profile and gets the start+stop of the
continuous regions.




.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    2 387840
      0%|          | 0/2 [00:00<?, ?it/s]/home/autumn/Documents/research-repos/batracker/batracker/tests/prototyping/plot_threshold_detector.py:107: RuntimeWarning: divide by zero encountered in log10
      dB = lambda X : 20*np.log10(np.abs(X))
     50%|#####     | 1/2 [00:06<00:06,  6.24s/it]    100%|##########| 2/2 [00:13<00:00,  6.51s/it]    100%|##########| 2/2 [00:13<00:00,  6.69s/it]






|


.. code-block:: default


    import os
    import matplotlib.pyplot as plt
    plt.rcParams['agg.path.chunksize']=10000
    import numpy as np
    import scipy.io.wavfile as wav
    import scipy.ndimage as ndimage
    import tqdm
    # % Load the raw audio 

    audio_path = os.path.join('simulated_audio','batracker_simple.wav')
    fs, audio = wav.read(audio_path)

    ## % threshold detector
    #plt.figure()
    #plt.plot(audio[:,0])

    ## % 
    def cross_channel_threshold_detector(multichannel, fs, **kwargs):
        '''
        Parameters
        ----------
        multichannel : np.array
            Msamples x Nchannels audio data
        fs : float >0
    
    
        Returns
        -------
    
        Notes
        -----
        For further keyword arguments see the `threshold_detector` function
    
        See Also
        --------
        threshold_detector
    
        '''
        samples, channels = multichannel.shape
        print(channels, samples)
        all_detections = []
        for each in tqdm.tqdm(range(channels)):
            all_detections.append(threshold_detector(multichannel[:,each], fs, **kwargs))
        return all_detections
        


    ## % Detect audio within one channel 
    def threshold_detector(one_channel, fs, **kwargs):
        '''
        Calculates the dB rms profile of the input audio and 
        selects regions which arae above  the profile. 
    
        Parameters
        ----------
        one_channel
        fs
        dbrms_threshold: float, optional
            Defaults to  -50 dB rms
        dbrms_window: float, optional
            The window which is used to calculate the dB rms profile
            in seconds.  Defaults to 0.001 seconds.
    
        Returns
        -------
        detections : list with tuples
            Each tuple corresponds to a candidate signal region
        '''
        if one_channel.ndim > 1:
            raise IndexError(f'Input audio must be flattened, and have only 1 dimension. \
                             Current audio has {one_channel.ndim} dimensions')
        dbrms_window = kwargs.get('dbrms_window',0.001) # seconds
        dbrms_threshold = kwargs.get('dbrms_threshold', -50)
    
        window_samples = int(fs*dbrms_window)
        dBrms_profile = dB(moving_rms(one_channel, window_size=window_samples))
    
        labelled, num_regions = ndimage.label(dBrms_profile>dbrms_threshold)
        if num_regions==0:
            raise ValueError(f'No regions above threshold: {dbrms_threshold} dBrms found in this channel!')
        regions_above = ndimage.find_objects(labelled.flatten())
        regions_above_timestamps = [get_start_stop_times(each, fs) for each in regions_above]
    
        return regions_above_timestamps
    
    def get_start_stop_times(findobjects_tuple, fs):
        '''
    
        '''
        only_tuple = findobjects_tuple[0]
        start, stop = only_tuple.start/fs, only_tuple.stop/fs
        return start, stop



    dB = lambda X : 20*np.log10(np.abs(X))

    def rms(X):
        '''Root mean square of a signal '''
        return np.sqrt(np.mean(X**2.0))


    def moving_rms(X, **kwargs):
        '''Calculates moving rms of a signal with given window size. 
        Outputs np.array of *same* size as X. The rms of the 
        last few samples <= window_size away from the end are assigned
        to last full-window rms calculated
        Parameters
        ----------
        X :  np.array
            Signal of interest. 
        window_size : int, optional
                     Defaults to 125 samples. 
        Returns
        -------
        all_rms : np.array
            Moving rms of the signal. 
        '''
        window_size = kwargs.get('window_size', 125)
        starts = np.arange(0, X.size)
        stops = starts+window_size
        valid = stops<X.size
        valid_starts = np.int32(starts[valid])
        valid_stops = np.int32(stops[valid])
        all_rms = np.ones(X.size).reshape(-1,1)*999

        for i, (start, stop) in enumerate(zip(valid_starts, valid_stops)):
            rms_value = rms(X[start:stop])
            all_rms[i] = rms_value
    
        # replace all un-assigned samples with the last rms value
        all_rms[all_rms==999] = np.nan

        return all_rms


  
    if __name__ == '__main__':
        fs =  192000
    
        t = np.linspace(0,0.1,fs)
        freq = 20000
        sine = np.sin(2*np.pi*freq*t)
    
        silence  = np.zeros(int(fs*0.02))
        full_sine = np.concatenate((sine, silence, sine ))
        multi_full_sine = np.column_stack((full_sine,
                                          np.roll(full_sine, int(fs*0.02)))
                                        )
    
        multi_detections = cross_channel_threshold_detector(multi_full_sine, fs)
    
    
    
    
    
    
    
    
    

.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 0 minutes  13.951 seconds)


.. _sphx_glr_download_prototyping_plot_threshold_detector.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: plot_threshold_detector.py <plot_threshold_detector.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: plot_threshold_detector.ipynb <plot_threshold_detector.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
