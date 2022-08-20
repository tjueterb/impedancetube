from os.path import join, isdir
from os import mkdir
import matplotlib.pyplot as plt
import numpy as np
from acoular import Calib, TimeSamples, PowerSpectra
import sys
sys.path.append('./src')
from measurement import Measurement_E2611, MicSwitchCalib_E2611
from tube import Tube_Transmission

##############################################################################
# USER INPUT:
##############################################################################

# ---------------- Amplitude Calibration (with regular calibrator) -----------
# (use create_calib_factor_xml_file.py to convert raw csv files to xml):
calibpath = './Resources'
calibfile = 'calib.xml'
calibration = Calib(from_file=join(calibpath, calibfile))

# ---------------- Amplitude and Phase Correction Measurements ---------------
# relative path of the time data files (.h5 data format)
soundfilepath = './Resources/'

# filename of empty measurement with direct configuration:
filename_direct = 'empty_00_11_22_33_44_55.h5'
#channels of switched mic and filenames of measurements with switched configurations
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
mic_channels_wide   = [0, 2, 3, 5]

# Filenames of the measurements:
# (in the same directory as the other sound files):
filenames_measurement = ['measurement.h5', # you can add files here
                        ]

# Parameters for frequency data handling:
block_size = 4*2048
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
time_data = TimeSamples(name=join(soundfilepath, filename_direct), calib=calibration)

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
    time_data_switched = TimeSamples(name=join(soundfilepath, filenames_switched[i]), calib=calibration)

    # get frequency data of switched configuration:
    freq_data_switched = PowerSpectra(time_data=time_data_switched,
                                      block_size=freq_data.block_size,
                                      window=freq_data.window,
                                      cached=freq_data.cached)

    # calculate amplitude/phase correction for switched channel:
    calib = MicSwitchCalib_E2611(freq_data=freq_data,
                                 freq_data_switched=freq_data_switched,
                                 ref_channel=0,
                                 test_channel=i)

    # store result:
    H_c[:, i] = calib.H_c
    

# ---------------- Measurement  ----------------------------------------------
# iterate over all measurements
for filename_measurement in filenames_measurement:
    td = TimeSamples(name=join(soundfilepath, filename_measurement), 
                     calib=calibration)

    # get frequency data / csm:
    freq_data = PowerSpectra(time_data=td,
                            block_size=block_size,
                            window=window,
                            overlap=overlap,
                            cached=cached)

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    # use both narrow and wide microphone positions for lower and higher frequencies:
    for spacing in ['wide', 'narrow']:
        if spacing == 'narrow':
            tube = Tube_Transmission(tube_shape='rect',
                                     tube_d=0.1,
                                     l1=0.3,   # distance between beginning of specimen and mic 2
                                     l2=0.8,   # distance between beginning of specimen and mic 3
                                     s1=0.085,  # Distance between mic 1 and 2
                                     s2=0.085,  # Distance between mic 3 and 4
                                     d=0.5)   # length of test specimen (test tube section is 0.7m))
            mic_channels = mic_channels_narrow  # indices of microphones #1-#4

        elif spacing == 'wide':
            tube = Tube_Transmission(tube_shape='rect',
                                     tube_d=0.1,
                                     l1=0.3,   # distance between beginning of specimen and mic 2
                                     l2=0.8,   # distance between beginning of specimen and mic 3
                                     s1=0.5,  # Distance between mic 1 and 2
                                     s2=0.5,  # Distance between mic 3 and 4
                                     d=0.5)   # length of test specimen (test tube section is 0.7m))
            mic_channels = mic_channels_wide

        msm = Measurement_E2611(freq_data=freq_data,
                                tube=tube,
                                ref_channel=ref_channel,  # index of the reference microphone
                                mic_channels=mic_channels,  # indices of the microphones in positions 1-4
                                H_c=H_c)  # Amplitude/Phase Correction factors

        # get fft frequencies
        freqs = msm.freq_data.fftfreq()

        # get transmission factor
        t = msm.transmission_coefficient

        # calculate transmission loss
        transmission_loss = msm.transmission_loss

        # if needed: calculate Impedance, plotting is the same
        z = msm.z
        
        # get frequency working range
        freqrange = msm.working_frequency_range

        # only use frequencies in the working range
        idx = np.logical_and(freqs >= freqrange[0], freqs <= freqrange[1])

        # plot
        ax.plot(freqs[idx], transmission_loss[idx])
        

    ax.set(title=filename_measurement,
        xlabel='f [Hz]', 
        ylabel='Transmission loss [dB]')
    ax.legend(['wide', 'narrow'])

    # Save or show plot:
    if not savePlot:
        plt.show()
    else:
        # create plot directory if necessary:
        if not isdir(plotpath):
            mkdir(plotpath)
            
        # save figure as pdf:
        filename_plot = '%s.pdf' % filename_measurement[:-3]
        fig.savefig(join(plotpath, filename_plot))
