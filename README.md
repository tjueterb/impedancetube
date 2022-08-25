# impedancetube
A routine to calculate plane wave transmission loss in rectangular and circular tubes following [ASTM E2611](https://dx.doi.org/10.1520/E2611-19)

Dependencies: [acoular](http://acoular.org/), [numpy](http://numpy.org), [traits](https://docs.enthought.com/traits/traits_user_manual/intro.html). [matplotlib](https://matplotlib.org) is recommended but not required.

Using conda and pip, you can install the required packages into your environment with these commands:

```
conda install -c acoular acoular
conda install matplotlib
pip install impedancetube
```

I recommend doing the measurements as described here and modifying the example script for evaluation. The scripts assume that you are measuring with 6 microphones simultaneously and that you are using the first microphone (nearest to the sound source) as the reference. At the end of this instruction you can find more detailed instructions if your measurements are done with custom configurations.

## How to measure (E2611 one load method)
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
For the two-load method, use two different tube terminations as described in [ASTM E2611](https://dx.doi.org/10.1520/E2611-19).

## How do modify the example scripts for your measurements
The easiest way to evaluate measurements is taking the example scripts from the Examples folder and modifying them to use your measurements. Here, the process is explained exemplary for the E2611 one-load method

### Put the files where you need them:
Clone/download this repository. Open `Examples/example_E2611_one_load.py`.

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

When you run the script, you should get a plot of the transmission loss in dB for each measurement.


<!---This is how to do latex equations in markdown: <img src="https://render.githubusercontent.com/render/math?math=e^{i \pi} = -1">---> 
