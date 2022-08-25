# impedancetube
A routine to calculate plane wave transmission loss in the rectangular and circular tubes following [ASTM E2611](https://dx.doi.org/10.1520/E2611-19)

Dependencies: [acoular](http://acoular.org/), [numpy](http://numpy.org), [traits](https://docs.enthought.com/traits/traits_user_manual/intro.html). [matplotlib](https://matplotlib.org) is recommended but not required.

Using conda, you can install the needed packages into your environment with these commands:

```
conda install -c acoular acoular
conda install numpy matplotlib traits
```

I recommend doing the measurements as described here and modifying the demo.py script for evaluation. The demo.py assumes that you are measuring with 6 microphones simultaneously and that you are using the first microphone (nearest to the sound source) as the reference. At the end of this instruction you can find more detailed instructions if your measurements are done with custom configurations.

## How to measure
![Measurement Setup](https://github.com/tjueterb/pyTransmission/blob/main/Resources/Measurement_setup.png?raw=true)

### Switched Microphone Measurements
First, we will need to make measurement in the empty tube with an anechoic back end. These will be used to calculate amplitude and phase correction factors between the microphones. For the direct configuration, route the microphones as shown in the figure above. Recommended measuring time is 60s.

For the switched configurations, switch the microphone positions **without un- and replugging the cables**. For example: This would be the routing for the switched configuration of microphones 0 and 1:

* mic #0 in position 1 going into ch. 0
* mic #1 in position 0 going into ch. 1

The other switched configurations follow the same rule, just switch the positions of n-th microphone and the reference microphone. Do the following measurements:

1. Measure with the direct configuration
2. Measure with reference mic (#0) switched with mic #1
3. Measure with reference mic (#0) switched with mic #2
4. Measure with reference mic (#0) switched with mic #3
5. Measure with reference mic (#0) switched with mic #4
6. Measure with reference mic (#0) switched with mic #5

Afterwards, don't forget to bring the microphones back into their original order.

### Measurements
Insert your test specimen(s) and do your measurement(s)

## How do modify demo.py for your measurements
The easiest way to evaluate measurements is taking the demo.py file and modifying it to use your freshly done measurements. Here's how:

### Put the files where you need them:
Clone/download this repository. 

#### Audio Files
For the easiest workflow, copy all your audio (.h5) files into the Resources directory. Alternatively you can modify `soundfilepath = './Resources/'` to match the path where all your audio files are.


### General Parameters:
If your measurements were done as instructed, your reference channel is the first channel (0) and your four microphone channels for the narrow and wide microphone configurations don't need to be modified. T

```python
ref_channel = 0
mic_channels_narrow = [1, 2, 3, 4]
mic_channels_wide   = [0, 2, 3, 5]
```

### Set the filenames for the Amplitude/Phase correction
Set the filenames of your empty measurements with the direct configuration and the switched configurations (See section [Switched Microphone Measurements](#switched_microphone_measurements)). The key of the dictionary denotes the index of the microphone that was switched with the reference (keep in mind that indexing starts at 0). The values represent the corresponding filenames.

```python
filename_direct    = 'empty_direct.h5'
filenames_switched = {1: 'empty_switched_1-0.h5', 
                      2: 'empty_switched_2-0.h5',
                      3: 'empty_switched_3-0.h5',
                      4: 'empty_switched_4-0.h5',
                      5: 'empty_switched_5-0.h5'}
```
### Set the filename(s) for the Measurement
The measurement sound files have to be in the same directory as the other sound files defined in `soundfilepath`. You can add multiple filenames to the list.

```python
filenames_measurement = ['measurement.h5',
                        ]
```



### Set the parameters for the frequency data handling
In the norm a Hanning window is required. A large block size increases the frequency resolution. If you set `cached = True`, the PowerSpectra object will use caching. This makes calculations faster if run repeatedly.

```python
block_size = 4*2048
window = 'Hanning'
overlap = '50%'
cached = False
```

### Set parameters for plotting
Decide if you want the plot to be saved and set the save directory. The directory will be created if it doesn't exist yet.

```python
savePlot = True
plotpath = './Plots'
```

If you did everything so far correctly an run the script, you should get a plot of the transmission loss in dB for each measurement.


## Documentation of the Measurement class

Import pyTransmission (if you just use the source code you need to work in the same directory or add the pyTransmission folder to your path).

```python
import sys
sys.path.append('path_to/pyTransmission')
from measurement import Measurement
```

The Measurement class expects a PowerSpectra object as input for the frequency data. See the acoular documentation for more details. Hanning window is required in the norm, a large blocksize is recommended for precision.

```python
from acoular import TimeSamples, PowerSpectra

time_data = TimeSamples(name='path/file.h5')
freq_data = PowerSpectra(time_data=time_data,
                         block_size=4*2048,
                         window='Hanning')
```  

Create your Measurement object and supply it with the frequency data. Set the tube-specific settings using the `Tube_Transmission` class (more information below). Set the reference channel (usually the microphone closest to the sound source) and the channels of the microphones in the positions 1-4 (see picture below for an overview and keep in mind that python starts indexing at 0). 

```python
msm = Measurement(freq_data=freq_data,
                  tube = tube,
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

## Documentation of the Tube Class

The Tube class and it's subclasses are used to define the tube shape and dimensions. For a E2611 transmission measurement use the `Tube_Transmission` class. The properties are defined as follows:

```python
tube = Tube_Transmission(tube_shape='rect',
                                     tube_d=0.1,
                                     l1=0.3,
                                     l2=0.8,
                                     s1=0.085,
                                     s2=0.085,
                                     d=0.5))
```
See the figure up top for the length definitions. The `tube_shape` can be `'rect'` for rectangular tubes or `'circ'` for round tubes. `tube_d` defines the diameter or (if rectangular) largest section dimension of the tube.

<!---This is how to do latex equations in markdown: <img src="https://render.githubusercontent.com/render/math?math=e^{i \pi} = -1">---> 
