from os.path import join
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from acoular import TimeSamples, PowerSpectra
import sys
from warnings import warn
try:
    # if directory is the root directory:
    # adding '.' to the path seems to be necessary for debugging this file in VS Code
    sys.path.append('.')
    import impedancetube as imp
    warn('Run this from the file directory for the relative paths to the .h5 files to work. (or adjust soundfilepath and reference_data_path)\n')

except:
    # if directory is test directory:
    sys.path.append('..')
    import impedancetube as imp


def test_two_load_method():
    ##############################################################################
    # USER INPUT:
    ##############################################################################
    # path for .npy files with reference data
    reference_data_path = './reference_data_two_load'
    Path(reference_data_path).mkdir(parents=True, exist_ok=True)

    # ---------------- Amplitude and Phase Correction Measurements ---------------
    # relative path of the time data files (.h5 data format)
    soundfilepath = '../Resources/'

    # filename of empty measurement with direct configuration:
    filename_direct = 'empty_00_11_22_33_44_55.h5'
    # channels of switched mic and filenames of measurements with switched configurations
    filenames_switched = {1: 'empty_01_10_22_33_44_55.h5',  # <- here 2nd mic (index 1) was switched w/ ref (index 0)
                          2: 'empty_02_11_20_33_44_55.h5',
                          3: 'empty_03_11_22_30_44_55.h5',
                          4: 'empty_04_11_22_33_40_55.h5',
                          5: 'empty_05_11_22_33_44_50.h5'}

    # reference channel
    # important: The reference Channel has to be 0 for the amplitude/phase correction to work!:
    ref_channel = 0

    # Mic channels in positions 1-4 of the narrow and wide configuration
    # (if the channels are sorted in increasing ordner from next to loudspeaker
    # to far away from loudspeaker, this ordering is correct)
    mic_channels_narrow = [1, 2, 3, 4]
    mic_channels_wide = [0, 2, 3, 5]

    # Filenames of the measurements (One file in each list for each measurement):
    # (in the same directory as the other sound files):
    # First load case:
    filenames_measurement_one_load = ['measurement_one_load.h5',  # you can add files here
                                      ]
    # Second load case:
    filenames_measurement_two_load = ['measurement_two_load.h5',  # you can add files here
                                      ]

    # Parameters for frequency data handling:
    block_size = 4*2048
    window = 'Hanning'
    overlap = '50%'
    cached = False

    # Parameters for plot:
    showPlot = False

    ##############################################################################
    # CALCULATION: No user input from here on
    ##############################################################################

    # ---------------- Amplitude and Phase Correction  ---------------------------

    # get timedata of direct configuration:
    time_data = TimeSamples(name=join(soundfilepath, filename_direct))

    # get frequency data / csm of direct configuration:
    freq_data = PowerSpectra(time_data=time_data,
                             block_size=block_size,
                             window=window,
                             overlap=overlap,
                             cached=cached)

    # initialize correction transferfunction with ones so the
    # ref-ref transfer function stays as ones, which is correct
    H_c = np.ones((freq_data.csm.shape[0:2]), dtype=complex)

    # iterate over all switched configurations:
    for i in filenames_switched:
        # get timedata of switched configuration:
        time_data_switched = TimeSamples(
            name=join(soundfilepath, filenames_switched[i]))

        # get frequency data of switched configuration:
        freq_data_switched = PowerSpectra(time_data=time_data_switched,
                                          block_size=freq_data.block_size,
                                          window=freq_data.window,
                                          cached=freq_data.cached)

        # calculate amplitude/phase correction for switched channel:
        calib = imp.MicSwitchCalib_E2611(freq_data=freq_data,
                                         freq_data_switched=freq_data_switched,
                                         ref_channel=0,
                                         test_channel=i)

        # store result:
        H_c[:, i] = calib.H_c

    # ---------------- Measurement  ----------------------------------------------
    # iterate over all measurements
    for filename_measurement_one_load, filename_measurement_two_load in zip(filenames_measurement_one_load,
                                                                            filenames_measurement_two_load):
        td_one_load = TimeSamples(
            name=join(soundfilepath, filename_measurement_one_load))

        td_two_load = TimeSamples(
            name=join(soundfilepath, filename_measurement_two_load))

        # get frequency data / csm:
        freq_data_one_load = PowerSpectra(time_data=td_one_load,
                                          block_size=block_size,
                                          window=window,
                                          overlap=overlap,
                                          cached=cached)

        freq_data_two_load = PowerSpectra(time_data=td_two_load,
                                          block_size=block_size,
                                          window=window,
                                          overlap=overlap,
                                          cached=cached)
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        # use both narrow and wide microphone positions for lower and higher frequencies:
        for spacing in ['wide', 'narrow']:
            if spacing == 'narrow':
                tube = imp.Tube_Transmission(tube_shape='rect',
                                             tube_d=0.1,
                                             l1=0.3,   # distance between beginning of specimen and mic 2
                                             l2=0.8,   # distance between beginning of specimen and mic 3
                                             s1=0.085,  # Distance between mic 1 and 2
                                             s2=0.085,  # Distance between mic 3 and 4
                                             d=0.5)   # length of test specimen (test tube section is 0.7m))
                mic_channels = mic_channels_narrow  # indices of microphones #1-#4

            elif spacing == 'wide':
                tube = imp.Tube_Transmission(tube_shape='rect',
                                             tube_d=0.1,
                                             l1=0.3,   # distance between beginning of specimen and mic 2
                                             l2=0.8,   # distance between beginning of specimen and mic 3
                                             s1=0.5,  # Distance between mic 1 and 2
                                             s2=0.5,  # Distance between mic 3 and 4
                                             d=0.5)   # length of test specimen (test tube section is 0.7m))
                mic_channels = mic_channels_wide

            msm1 = imp.Measurement_E2611(freq_data=freq_data_one_load,
                                         freq_data_two_load=freq_data_two_load,
                                         method='two load',
                                         tube=tube,
                                         ref_channel=ref_channel,  # index of the reference microphone
                                         mic_channels=mic_channels,  # indices of the microphones in positions 1-4
                                         H_c=H_c)  # Amplitude/Phase Correction factors

            # switched first and second load case
            msm2 = imp.Measurement_E2611(freq_data=freq_data_two_load,
                                         freq_data_two_load=freq_data_one_load,
                                         method='two load',
                                         tube=tube,
                                         ref_channel=ref_channel,  # index of the reference microphone
                                         mic_channels=mic_channels,  # indices of the microphones in positions 1-4
                                         H_c=H_c)  # Amplitude/Phase Correction factors

            # get fft frequencies
            freqs1 = msm1.freq_data.fftfreq()
            freqs2 = msm2.freq_data.fftfreq()
            # np.save(f'{reference_data_path}/freqs_{spacing}', freqs1)
            assert(np.allclose(freqs1, freqs2, equal_nan=True))
            assert(np.allclose(freqs1, np.load(
                f'{reference_data_path}/freqs_{spacing}.npy'), equal_nan=True))

            # get transfer_matric
            T1 = msm1.transfer_matrix
            T2 = msm2.transfer_matrix
            # np.save(f'{reference_data_path}/transfer_matrix_{spacing}', T1)
            assert(np.allclose(T1, T2, equal_nan=True))
            assert(np.allclose(T1, np.load(
                f'{reference_data_path}/transfer_matrix_{spacing}.npy'), equal_nan=True))

            # get transmission factor
            t1 = msm1.transmission_coefficient
            t2 = msm2.transmission_coefficient
            # np.save(f'{reference_data_path}/transmission_coefficient_{spacing}', t1)
            assert(np.allclose(t1, t2, equal_nan=True))
            assert(np.allclose(t1, np.load(
                f'{reference_data_path}/transmission_coefficient_{spacing}.npy'), equal_nan=True))

            # calculate transmission loss
            transmission_loss1 = msm1.transmission_loss
            transmission_loss2 = msm2.transmission_loss
            # np.save(f'{reference_data_path}/transmission_loss_{spacing}', transmission_loss1)
            assert(np.allclose(transmission_loss1,
                               transmission_loss2, equal_nan=True))
            assert(np.allclose(transmission_loss1, np.load(
                f'{reference_data_path}/transmission_loss_{spacing}.npy'), equal_nan=True))

            # if needed: calculate Impedance, plotting is the same
            z1 = msm1.z
            z2 = msm2.z
            # np.save(f'{reference_data_path}/z_{spacing}', z1)
            assert(np.allclose(z1, z2, equal_nan=True))
            assert(np.allclose(z1, np.load(
                f'{reference_data_path}/z_{spacing}.npy'), equal_nan=True))

            # get frequency working range
            freqrange1 = msm1.working_frequency_range
            freqrange2 = msm2.working_frequency_range
            # np.save(f'{reference_data_path}/freqrange_{spacing}', freqrange1)
            assert(np.allclose(freqrange1, freqrange2, equal_nan=True))
            assert(np.allclose(freqrange1, np.load(
                f'{reference_data_path}/freqrange_{spacing}.npy'), equal_nan=True))

            # only use frequencies in the working range
            idx = np.logical_and(
                freqs1 >= freqrange1[0], freqs1 <= freqrange1[1])

            # plot
            ax.plot(freqs1[idx], transmission_loss1[idx])

        ax.set(title=f'{filename_measurement_one_load}\n{filename_measurement_two_load}',
               xlabel='f [Hz]',
               ylabel='Transmission loss [dB]')
        ax.legend(['wide', 'narrow'])

        # Save or show plot:
        if showPlot:
            plt.show()


if __name__ == "__main__":
    test_two_load_method()
    print("\nEverything passed!\n")
