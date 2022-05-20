from os.path import join, isdir
from os import mkdir
import matplotlib.pyplot as plt
import numpy as np
from acoular import TimeSamples, PowerSpectra

from pyTransmission import Measurement_Cottbus, MicSwitchCalib_Cottbus

##############################################################################
# USER INPUT:
##############################################################################

# # ---------------- Amplitude Calibration (with regular calibrator) -----------
# # (use create_calib_factor_xml_file.py to convert raw csv files to xml):
# calibpath = './Resources'
# calibfile = 'calib.xml'
# calibration = Calib(from_file=join(calibpath, calibfile))

# ---------------- Amplitude and Phase Correction Measurements ---------------
# relative path of the time data files (.h5 data format)
soundfilepath = './Resources/'

# filename of empty measurement with direct configuration:
filename_direct = 'empty_00_11_22_33_44_55.h5' 

# calibration files with switched microphone configurations (see README), DON'T CHANGE THE KEY FIELDS:
filenames_switched = {'00_12_21_34_43_55': 'empty_00_12_21_34_43_55.h5',  # <- here mic 1 was switched w/ mic 2 and mic 3 was switched w/ mic 4
                      '00_13_22_31_44_55': 'empty_00_13_22_31_44_55.h5',
                      '02_11_20_35_44_53': 'empty_02_11_20_35_44_53.h5',
                      '03_11_22_30_44_55': 'empty_03_11_22_30_44_55.h5',}

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
time_data = TimeSamples(name=join(soundfilepath, filename_direct))

# get frequency data / csm of direct configuration:
freq_data = PowerSpectra(time_data=time_data,
                         block_size=block_size,
                         window=window,
                         overlap=overlap,
                         cached=cached)

# initialize correction transferfunction with ones so the
# ref-ref transfer function stays as ones, which is correct
H_c_narrow = np.ones((freq_data.csm.shape[0],3), dtype=complex)
H_c_wide   = np.ones((freq_data.csm.shape[0],3), dtype=complex)

configs = ['00_12_21_34_43_55',
           '00_13_22_31_44_55',
           '02_11_20_35_44_53',
           '03_11_22_30_44_55']

for i in configs:
    # get timedata of switched configuration:
    time_data_switched = TimeSamples(name=join(soundfilepath, filenames_switched[i]))
    # get frequency data of switched configuration:
    freq_data_switched = PowerSpectra(time_data=time_data_switched,
                                      block_size=freq_data.block_size,
                                      window=freq_data.window,
                                      cached=freq_data.cached)
    
    if i == '00_12_21_34_43_55':
        # calculate amplitude/phase correction for switched channel:
        calib = MicSwitchCalib_Cottbus(freq_data=freq_data,
                                       freq_data_switched=freq_data_switched,
                                       ref_channel=1,
                                       test_channel=2)
        # store result:
        H_c_narrow[:, 0] = calib.H_c
        
        calib = MicSwitchCalib_Cottbus(freq_data=freq_data,
                                       freq_data_switched=freq_data_switched,
                                       ref_channel=3,
                                       test_channel=4)
        # store result:
        H_c_narrow[:, 2] = calib.H_c
    
    if i == '00_13_22_31_44_55':
        # calculate amplitude/phase correction for switched channel:
        calib = MicSwitchCalib_Cottbus(freq_data=freq_data,
                                       freq_data_switched=freq_data_switched,
                                       ref_channel=1,
                                       test_channel=3)
        # store result:
        H_c_narrow[:, 1] = calib.H_c
        
    if i == '02_11_20_35_44_53':
        # calculate amplitude/phase correction for switched channel:
        calib = MicSwitchCalib_Cottbus(freq_data=freq_data,
                                       freq_data_switched=freq_data_switched,
                                       ref_channel=0,
                                       test_channel=2)
        # store result:
        H_c_wide[:, 0] = calib.H_c
        
        calib = MicSwitchCalib_Cottbus(freq_data=freq_data,
                                       freq_data_switched=freq_data_switched,
                                       ref_channel=3,
                                       test_channel=5)
        # store result:
        H_c_wide[:, 2] = calib.H_c
    
    if i == '03_11_22_30_44_55':
        # calculate amplitude/phase correction for switched channel:
        calib = MicSwitchCalib_Cottbus(freq_data=freq_data,
                                       freq_data_switched=freq_data_switched,
                                       ref_channel=0,
                                       test_channel=3)
        # store result:
        H_c_wide[:, 1] = calib.H_c

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
    ax = fig.add_subplot(1,1,1)
    # use both narrow and wide microphone positions for lower and higher frequencies:
    for spacing in ['wide', 'narrow']:
        if spacing == 'narrow':
            s1 = s2 = 0.085  # distance between mics
            mic_channels = mic_channels_narrow  # indices of microphones #1-#4
            H_c = H_c_narrow # phase/amplitude correction transfer functions
            
        elif spacing == 'wide':
            s1 = s2 = 0.5 # distance between mics
            mic_channels = mic_channels_wide
            H_c = H_c_wide # phase/amplitude correction transfer functions
            
        msm = Measurement_Cottbus(freq_data=freq_data,
                        s1=s1,  # distance between mic #1 and #2
                        s2=s2,  # distance between mic #3 and #4
                        mic_channels=mic_channels, # indices of the microphones in positions 1-4
                        H_c=H_c) # Amplitude/Phase Correction factors  

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
