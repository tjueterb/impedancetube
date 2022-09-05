import numpy as np
import pyfar as pf
import h5py


#TODO: add position 3 for a new 

fs = 44100
c = 20.047 * np.sqrt(273.15 + 20) # speed of sound at 20 °C

#input noise
noise = pf.signals.noise(n_samples=4*fs, spectrum='white',rms=0.1, sampling_rate=fs)

# incident signal at position 1 (can be omitted)
dirac = pf.signals.impulse(n_samples=1, sampling_rate=fs)
pos1_incident = pf.dsp.convolve(signal1=noise, signal2=dirac)

# delay input from position 1 to position 2
s_narrow = 0.45
dt = s_narrow / c
dN = int(dt * fs)
dirac_pos1_pos2 = pf.signals.impulse(n_samples=dN+1, delay=dN, sampling_rate=fs)

# incident signal at position 2
pos2_incident = pf.dsp.convolve(signal1=pos1_incident, signal2=dirac_pos1_pos2)

# incident signal at position 3
s = 0.85 - 0.45 # distance pos2-pos3
dt = s / c
dN = int(dt * fs)
dirac_pos2_pos3 = pf.signals.impulse(n_samples=dN+1, delay=dN, sampling_rate=fs)
pos3_incident = pf.dsp.convolve(signal1=pos2_incident, signal2=dirac_pos2_pos3)

# incident signal at back wall
x1 = 2.75 # position 1 to back wall
dt = (x1-0.85) / c # time delay position 3 to back wall
dN = int(dt * fs)
dirac_pos3_back = pf.signals.impulse(n_samples=dN+1, delay=dN, sampling_rate=fs)
back_wall = pf.dsp.convolve(signal1=pos3_incident, signal2=dirac_pos3_back)

# uniform spectrum 0.5 gain
reflection_factor = 0.5
back_wall_reflected = reflection_factor * back_wall

# reflected signal at position 3
pos3_reflected = pf.dsp.convolve(signal1=back_wall_reflected, signal2=dirac_pos3_back)


# reflected signal at position 2
pos2_reflected = pf.dsp.convolve(signal1=pos3_reflected, signal2=dirac_pos2_pos3)

# reflected signal at position 1
pos1_reflected = pf.dsp.convolve(signal1=pos2_reflected, signal2=dirac_pos1_pos2)

# crop signals and add incident and reflected signal
pos1_incident  = pf.dsp.time_window(signal=pos1_incident,  window='boxcar', interval=[fs,3*fs-1], crop='window')
pos1_reflected = pf.dsp.time_window(signal=pos1_reflected, window='boxcar', interval=[fs,3*fs-1], crop='window')
pos2_incident  = pf.dsp.time_window(signal=pos2_incident,  window='boxcar', interval=[fs,3*fs-1], crop='window')
pos2_reflected = pf.dsp.time_window(signal=pos2_reflected, window='boxcar', interval=[fs,3*fs-1], crop='window')
pos3_incident  = pf.dsp.time_window(signal=pos3_incident,  window='boxcar', interval=[fs,3*fs-1], crop='window')
pos3_reflected = pf.dsp.time_window(signal=pos3_reflected, window='boxcar', interval=[fs,3*fs-1], crop='window')

signal_pos1 = pos1_incident + pos1_reflected
signal_pos2 = pos2_incident + pos2_reflected
signal_pos3 = pos3_incident + pos3_reflected

# TODO: mix signals to one object with 3 channels (acoular excpects array of shape [N_samples x N_channels])
time_data = np.concatenate((signal_pos1.time, signal_pos2.time, signal_pos3.time), axis=0).T
time_data = time_data.astype('float32')
# time_data = pf.Signal(data=data, sampling_rate=fs, domain='time')

# Write to file
hf = h5py.File('../Resources/example_ISO10534-2.h5','w')
hf.create_dataset('time_data',data=time_data)
hf['time_data'].attrs['sample_freq'] = fs
hf.close()

