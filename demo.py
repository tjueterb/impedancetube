from os.path import join

import matplotlib.pyplot as plt
import numpy as np
from acoular import Calib, TimeSamples, PowerSpectra

from pyTransmission import Measurement

# define calibration data
# (use create_calib_factor_xml_file.py to convert the raw calib files to xml):
calibpath = './Resources'
calibfile = '2021-05-05_calibration.xml'
calibration = Calib(from_file=join(calibpath, calibfile))

# define timedata:
filepath = './Resources/'
filename = 'demo.h5'
td = TimeSamples(name=join(filepath, filename), calib=calibration)

# get frequency data / csm:
freq_data = PowerSpectra(time_data=td,
                         block_size=4*2048,
                         window='Hanning')

plt.figure()

# use both narrow and wide microphone positions for lower and higher frequencies:
for spacing in ['narrow', 'wide']:
    if spacing == 'narrow':
        s1 = 0.085  # distance between mic #1 and #2
        s2 = 0.085  # distance between mic #3 and #4
        ref_channel = 4  # index of reference microphone
        mic_channels = [4, 3, 2, 1]  # indices of microphones #1-#4

    elif spacing == 'wide':
        s1 = 0.5
        s2 = 0.5
        ref_channel = 5
        mic_channels = [5, 3, 2, 0]

    msm = Measurement(freq_data=freq_data,
                      s1=s1,  # distance between mic #1 and #2
                      s2=s2,  # distance between mic #3 and #4
                      ref_channel=ref_channel,  # index of the reference microphone
                      mic_channels=mic_channels)  # indices of the microphones in positions 1-4

    # get fft frequencies
    freqs = msm.freq_data.fftfreq()

    # get transmission factor
    t = msm.transmission_coefficient

    # calculate transmission loss
    transmission_loss = -10*np.log10(t)

    # get frequency working range
    freqrange = msm.working_frequency_range

    # only use frequencies in the working range
    idx = np.logical_and(freqs >= freqrange[0], freqs <= freqrange[1])

    # plot real part
    plt.plot(freqs[idx], transmission_loss[idx].real)

plt.legend(['narrow', 'wide'])
plt.title(filename)
plt.ylim([-10, 30])
plt.xlabel('f [Hz]')
plt.ylabel('Transmission loss [dB]')
plt.show
