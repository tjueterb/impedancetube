# pyTransmission
A routine to calculate plane wave transmission loss in the rectangular tube at the TAP

Dependencies: [acoular](http://acoular.org/), [numpy](http://numpy.org)

Using conda, you can install acoular into your environment with this command:

```
conda install -c acoular acoular
```

### How to use
see demo.py for a simple script. The Measurement class expects a PowerSpectra object as input. See the acoular documentation for more details. Hanning window is required in the norm, a large blocksize is recommended for precision.

```python
from acoular import TimeSamples, PowerSpectra

time_data = TimeSamples(name='path/file.h5')
freq_data = PowerSpectra(time_data=time_data,
                         block_size=4*2048,
                         window='Hanning')
```  

Import pyTransmission (if you just use the source code you need to work in the same directory or add the pyTransmission folder to your path).

```python
import sys
sys.path.append('path_to/pyTransmission')
from pyTransmission import Measurement
```

Create your Measurement object and supply it with the frequency data. Set the distances between the microphone pairs `s1` and `s2` (0.085 or 0.5). Set the reference channel (usually the microphone closest to the sound source) and the channels of the microphones in the positions 1-4 (see picture below for an overview and keep in mind that python starts indexing at 0). 

```python
msm = Measurement(freq_data=freq_data,
                  s1=0.085,   
                  s2=0.085,
                  ref_channel=0,
                  mic_channels=[0, 1, 2, 3])
```

![setup](https://github.com/tjueterb/pyTransmission/blob/eac49d54ffd6a800107fd2fae0760da1ad3355f4/Resources/Measurement_setup.png?raw=true)

Get the FFT frequencies:

```python
freqs = msm.freq_data.fftfreq()
```

Get the frequency dependent plane wave incident transmission loss:

```python
TL = msm.transmission_loss
```

Get the working frequency range (lower and upper limits of where the measurement is valid, dependent on the microphone spacing)

```python
f_low, f_high = msm.working_frequency_range
```

Plot the transmission loss only for the valid frequencies:

```python
import numpy as np
import matplotlib.pyplot as plt

idx = np.logical_and(freqs >= f_low, freqs <= f_high)
plt.plot(freqs[idx], TL[idx])
```
Some other traits you can set and their default values:
`temperature=20.`, `atmospheric_pressure=101.325`, `l1=0.3`, `l2=0.8`, `d=0.5`

Get a list of all traits with `msm.trait_names()`




<!---This is how to do latex equations in markdown: <img src="https://render.githubusercontent.com/render/math?math=e^{i \pi} = -1">---> 
