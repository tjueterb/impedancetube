from os.path import join, isdir
from os import mkdir
import matplotlib.pyplot as plt
import numpy as np
from acoular import TimeSamples, PowerSpectra
# import impedancetube as imp
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
##############################################################################
# USER INPUT:
##############################################################################

# ---------------- Amplitude and Phase Correction Measurements ---------------
# relative path of the time data files (.h5 data format), edit this accordingly!
soundfilepath = '../Resources/'

# filename of empty measurement with direct configuration:
filename_direct = 'empty_00_11_22_33_44_55.h5' #TODO: correct the files
# channels of switched mic and filenames of measurements with switched configurations
filenames_switched = {1: 'empty_01_10_22_33_44_55.h5',  # <- here mics in positions 1 and 2 (indices 0 and 1) are switched
                      2: 'empty_02_11_20_33_44_55.h5',  # <- here mics in positions 1 and 3 (indices 0 and 2) are switched
                      }

# reference channel (microphone position 1)
# important: The reference Channel has to be 0 for the amplitude/phase correction to work!:
ref_channel = 0

# channel of microphone in position 2
# (if the channels are sorted in increasing ordner from next to loudspeaker
# to far away from loudspeaker, this ordering is correct)
mic_channel_narrow = 1
mic_channel_wide = 2

# Filenames of the measurements (One file for each measurement):
# (in the same directory as the other sound files):
filenames_measurement = ['example_ISO10534-2.h5',  # add files here
                         ]

# Parameters for frequency data handling:
block_size = 16*2048
window = 'Hanning'
overlap = '50%'
cached = False

# Parameters for plot:
savePlot = False
plotpath = './Plots'

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
H_c = np.ones((freq_data.csm.shape[0:2]), dtype=complex) #TODO: verify that this shape is correct, ideally yes

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

H_c = np.ones_like(H_c) #TODO: synthesize calibration data and disable this

# ---------------- Measurement  ----------------------------------------------
# iterate over all measurements
for filename_measurement in filenames_measurement:
    td = TimeSamples(name=join(soundfilepath, filename_measurement))

    # get frequency data / csm:
    freq_data = PowerSpectra(time_data=td,
                             block_size=block_size,
                             window=window,
                             overlap=overlap,
                             cached=cached)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    # use both narrow and wide microphone positions for lower and higher frequencies:
    for spacing in ['wide', 'narrow']:
        if spacing == 'narrow':
            '''Either define your tube explicitly:'''
            # tube = imp.Tube_Impedance(tube_shape='rect',
            #                           tube_d=0.5,
            #                           s=0.45,
            #                           x1=2.75)
            '''Or use a preset:'''
            tube = imp.Tube_Impedance_TAP_rect_50cm_narrow() 
            mic_channel = mic_channel_narrow  # indices of microphones #1-#4

        elif spacing == 'wide':
            '''Either define your tube explicitly:'''
            # tube = imp.Tube_Impedance(tube_shape='rect',
            #                           tube_d=0.5,
            #                           s=0.85,
            #                           x1=2.75)
            '''Or use a preset:'''
            tube = imp.Tube_Impedance_TAP_rect_50cm_wide()
            mic_channel = mic_channel_wide

        msm = imp.Measurement_ISO10534(freq_data=freq_data,
                                       tube=tube,
                                       ref_channel=ref_channel,  # index of nicrophone in position 1
                                       mic_channel=mic_channel,  # index of microphone in position 2
                                       H_c=H_c)  # Amplitude/Phase Correction factors

        # get fft frequencies
        freqs = msm.freq_data.fftfreq()

        # get reflection factor
        reflection_factor = msm.reflection_factor

        # absorption coefficient
        absorption_coefficient = msm.absorption_coefficient
        
        # specific acoustic impedance ratio
        specific_acoustic_impedance_ratio = msm.specific_acoustic_impedance_ratio
        
        # specific acoustic admittance ratio
        specific_acoustic_admittance_ratio = msm.specific_acoustic_admittance_ratio
        
        # get frequency working range
        freqrange = msm.working_frequency_range

        # only use frequencies in the working range
        idx = np.logical_and(freqs >= freqrange[0], freqs <= freqrange[1])

        # plot
        ax.plot(freqs[idx], absorption_coefficient[idx])

    ax.set(title=filename_measurement,
           xlabel='f [Hz]',
           ylabel='Absorption coefficient',
           ylim=(0,1))
    ax.legend(['wide', 'narrow'])

    # Save or show plot:
    if not savePlot:
        plt.show()
    else:
        # create plot directory if necessary:
        if not isdir(plotpath):
            mkdir(plotpath)

        # save figure as pdf:
        filename_plot = f'{filename_measurement[:-3]}.pdf'
        fig.savefig(join(plotpath, filename_plot))
