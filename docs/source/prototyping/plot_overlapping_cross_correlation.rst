.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download_prototyping_plot_overlapping_cross_correlation.py>`     to download the full example code
    .. rst-class:: sphx-glr-example-title

    .. _sphx_glr_prototyping_plot_overlapping_cross_correlation.py:


When calls arrive almost/overlapping
====================================

Solution 1:
    
    #. Do the cross-correlation anyway - the multiple peaks might be picked up, and 
    the relative channel delays suggested by both the peaks could be tried out anyway
    #. Do the cross-correlation after 'intra-audio normalisation' - that is, make sure
    the RMS profile within the sound region to be cross correlated is uniform. This 
    prevents false positive peaks from forming during the cross correlation. 

Broad findings:
    
    #. Hilbert envelope based audio-normalisation could work, at least for the 
    simple case of overlapping calls (without reverb etc.)


.. code-block:: default

    import matplotlib.pyplot as plt
    plt.rcParams['agg.path.chunksize'] = 100000
    import scipy.io.wavfile as wavfile
    import scipy.signal as signal 
    import numpy as np

    def make_wav_and_spec(X, fs=192000):
        time = np.linspace(0, X.size/fs, X.size)
        plt.figure()
        a0  = plt.subplot(211)
        plt.specgram(X, Fs=fs)
        plt.subplot(212, sharex=a0)
        plt.plot(time, X)


    fs = 192000
    start_f, end_f =  90000, 30000
    durn = 0.004
    num_samples = int(fs*durn)
    t = np.linspace(0,durn,num_samples)
    call = signal.chirp(t, start_f, t[-1], end_f)
    call *= signal.hanning(call.size)
    call *= 0.5


    make_wav_and_spec(call)




.. image:: /prototyping/images/sphx_glr_plot_overlapping_cross_correlation_001.png
    :alt: plot overlapping cross correlation
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    /home/autumn/Documents/trying_out/batracker/batracker/tests/prototyping/plot_overlapping_cross_correlation.py:42: DeprecationWarning: `hanning` is deprecated, use `scipy.signal.windows.hann` instead!
      call *= signal.hanning(call.size)




Now make a channel with no overlaps, and a channel with two calls arriving very
close to  each other. 


.. code-block:: default


    ch1 = np.zeros(1920)
    ch1[100:100+num_samples] += call
    ch2 = np.zeros(1920)
    ch2[300:300+num_samples] += call
    ch2[700:700+num_samples] += call


    make_wav_and_spec(ch2)


    cc = signal.correlate(ch1, ch2, 'same')
    plt.figure()
    plt.plot(cc)
    plt.vlines(cc.size/2.0, 0, np.max(cc))






.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /prototyping/images/sphx_glr_plot_overlapping_cross_correlation_002.png
          :alt: plot overlapping cross correlation
          :class: sphx-glr-multi-img

    *

      .. image:: /prototyping/images/sphx_glr_plot_overlapping_cross_correlation_003.png
          :alt: plot overlapping cross correlation
          :class: sphx-glr-multi-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    /home/autumn/anaconda3/envs/batracker/lib/python3.7/site-packages/matplotlib/axes/_axes.py:7531: RuntimeWarning: divide by zero encountered in log10
      Z = 10. * np.log10(spec)

    <matplotlib.collections.LineCollection object at 0x7f4b51fbc7b8>



What happens if one call is louder than the other in the channel with overlaps


.. code-block:: default


    mixed_ampch2 = np.zeros(1920)
    mixed_ampch2[300:300+num_samples] += call*0.3
    mixed_ampch2[700:700+num_samples] += call

    make_wav_and_spec(mixed_ampch2)

    cc_mixedamp = signal.correlate(ch1, mixed_ampch2, 'same')
    plt.figure()
    plt.plot(cc_mixedamp)
    plt.vlines(cc_mixedamp.size/2.0, 0, np.max(cc_mixedamp))




.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /prototyping/images/sphx_glr_plot_overlapping_cross_correlation_004.png
          :alt: plot overlapping cross correlation
          :class: sphx-glr-multi-img

    *

      .. image:: /prototyping/images/sphx_glr_plot_overlapping_cross_correlation_005.png
          :alt: plot overlapping cross correlation
          :class: sphx-glr-multi-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    /home/autumn/anaconda3/envs/batracker/lib/python3.7/site-packages/matplotlib/axes/_axes.py:7531: RuntimeWarning: divide by zero encountered in log10
      Z = 10. * np.log10(spec)

    <matplotlib.collections.LineCollection object at 0x7f4b51e6a860>



The peaks are of very different heights which means that the delays will of course 
be interpreted badly. Now, let's equalise the waveform somehow. 


.. code-block:: default


    hilbert_tr = signal.hilbert(mixed_ampch2)
    hilbert_envelope = np.abs(hilbert_tr)
    plt.figure()
    plt.plot(mixed_ampch2)
    plt.plot(hilbert_envelope)

    max_env = np.max(hilbert_envelope)
    amp_factor = hilbert_envelope/max_env

    uniform_amp = mixed_ampch2*(1/amp_factor)

    plt.figure()
    plt.plot(hilbert_envelope)
    plt.plot(uniform_amp)


    processed_cc = signal.correlate(ch1, uniform_amp, 'same')
    plt.figure()
    plt.plot(processed_cc)
    plt.vlines(processed_cc.size/2.0, 0, np.max(processed_cc))





.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /prototyping/images/sphx_glr_plot_overlapping_cross_correlation_006.png
          :alt: plot overlapping cross correlation
          :class: sphx-glr-multi-img

    *

      .. image:: /prototyping/images/sphx_glr_plot_overlapping_cross_correlation_007.png
          :alt: plot overlapping cross correlation
          :class: sphx-glr-multi-img

    *

      .. image:: /prototyping/images/sphx_glr_plot_overlapping_cross_correlation_008.png
          :alt: plot overlapping cross correlation
          :class: sphx-glr-multi-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    <matplotlib.collections.LineCollection object at 0x7f4b51d1b630>



Let's compare the raw overlapping calls CC with the intra-audio normlised
The peaks from the intra-audio normalised CC is definitely much larger than 
just the raw CC - which shows that the intra-audio normalisation definitely helps!


.. code-block:: default


    plt.figure()
    plt.subplot(211)
    plt.plot(cc_mixedamp)
    plt.vlines(cc_mixedamp.size/2.0, 0, np.max(cc_mixedamp))
    plt.title('Raw mixed amplitude CC')
    plt.xticks([])
    plt.subplot(212)
    plt.plot(processed_cc)
    plt.vlines(processed_cc.size/2.0, 0, np.max(processed_cc))
    plt.title('Amplitude adjusted CC')





.. image:: /prototyping/images/sphx_glr_plot_overlapping_cross_correlation_009.png
    :alt: Raw mixed amplitude CC, Amplitude adjusted CC
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    Text(0.5, 1.0, 'Amplitude adjusted CC')



One thing to pay attention to - the SNR before actually doing the envelope 
based normalisation. If the SNR is already low - then it'd mean the noise would 
be amplified -- and mess up the CC. 


.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 0 minutes  1.305 seconds)


.. _sphx_glr_download_prototyping_plot_overlapping_cross_correlation.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: plot_overlapping_cross_correlation.py <plot_overlapping_cross_correlation.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: plot_overlapping_cross_correlation.ipynb <plot_overlapping_cross_correlation.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
