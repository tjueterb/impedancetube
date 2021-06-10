# This defines a class for transmission factor measurement

from acoular.spectra import PowerSpectra
import numpy as np
from traits.api import HasPrivateTraits, Property, Float, Trait, Delegate, Int, List, Array


class Measurement(HasPrivateTraits):
    temperature = Float(default_value=20.)
    atmospheric_pressure = Float(default_value=101.325)
    c = Property(depends_on=['temperature'],
                 desc='Speed of sound')
    rho = Property(depends_on=['temperature', 'atmospheric_pressure'],
                   desc='air density')

    # Tube dimensions
    l1 = Float(0.3, desc='distance between beginning of speciman and mic 2')
    l2 = Float(0.8, desc='distance between beginning of specimen and mic 3')
    s1 = Trait(0.085, 0.5, desc='Distance between mic 1 and 2 in m')
    s2 = Trait(0.085, 0.5, desc='Distance between mic 3 and 4 in m')
    d = Float(0.5, desc='length of test specimen (test tube section is 0.7m)')

    # channels of the microphones in the given freq_data object
    ref_channel = Int(0, desc="Channel index of the reference mic")
    mic_channels = List([1, 2, 3, 4], minlen=4, maxlen=4,
                        desc="Channel indices of mics in positions 1-4")
    
    # Amplitude and Phase correction Transfer function
    H_c = Array()
    #TODO: either Initialize with ones (only possible after freq_data shape is known)
    #TODO: Or rewrite the calib class so the entire object can be handed to this class

    # Block size
    block_size = Delegate('freq_data')

    # Tube parameters
    tube_shape = Trait('rect', 'circ',
                       desc="Shape of the measurement tube")
    tube_d = Float(0.1,
                   desc="diameter or (if rectangular) largest section dimension of tube")

    #: The :class:`~acoular.spectra.PowerSpectra` object that provides the csm.
    freq_data = Trait(PowerSpectra,
                      desc="power spectra object")

    # wave number vector:
    k = Property(depends_on=['freq_data'])

    # The transfer function for all microphones:
    transfer_function = Property(
        desc='Transfer function between all mics and ref. mic (channel i_ref)')

    # The transfer Matrix
    transfer_matrix_one_load = Property(desc='Transfer Matrix for one load')

    # transmission coefficient:
    transmission_coefficient = Property()

    # Transmission loss:
    transmission_loss = Property()

    # Working frequency range:
    working_frequency_range = Property()
    
    # Characteristic Impedanze in material:
    z = Property()

    def _get_c(self):
        """Calculates speed of sound from temperature
        See 8.2 (eq. (4))
        Args:
            temperature (float): Temperature in °C.
        Returns:
            c (float): speed of sound in m/s
        """
        c = 20.047 * np.sqrt(273.15 + self.temperature)
        return c

    def _get_rho(self):
        """Getter function for the air density

        Returns:
            [type]: [description]
        """
        rho = 1.290 * (self.atmospheric_pressure/101.325) * \
            (273.15 / (273.15 + self.temperature))
        return rho

    def _get_transfer_function(self):
        """Calculates the transfer function for all microphones
        See 8.5.1 eq (15)

        NOTE: Phase and amplitude correction transfer functions 
              are not implemented yet
              
        Returns:
            [array]:    Transfer function for all microphones 
                        size: (f x n), f: frequencies, n: channels
        """
        csm = self.freq_data.csm

        H_n_ref = np.empty(csm.shape[0:2], dtype=complex)  # create empty array

        for n in range(self.freq_data.numchannels):
            #NOTE: correct indexing of the csm verified with Adam
            H_n_ref[:, n] = csm[:, n, self.ref_channel] / \
                csm[:, self.ref_channel, self.ref_channel]  # eq (15)
                
        # apply correction transfer function:
        H_n_ref = H_n_ref / self.H_c
        
        return H_n_ref

    def _get_k(self):
        """Calculates the wave number coefficients for all frequencies
        See 3.2: k = 2*pi*f / c
        
        Returns:
            [array]: size: (f x 1)
        """
        k = 2 * np.pi * self.freq_data.fftfreq() / self.c
        return k

    def _get_transfer_matrix_one_load(self):
        """Calculates the one load transfer matrix 
        See 8.5.4.2
        
        Returns:
            [array]: Transfer Matrix
                     size:(f x 2 x 2), f: frequencies
        """ 
        # Get transfer function:
        H = self.transfer_function

        # Disable divide by zero warning because first entry of k is always 0
        with np.errstate(divide='ignore', invalid='ignore'):
            # Decompose wave field:
            # eq (17):
            A = 1j * ((H[:, self.mic_channels[0]] * np.exp(-1j*self.k*(self.l1)) -
                       H[:, self.mic_channels[1]] * np.exp(-1j*self.k*(self.l1+self.s1))) /
                      (2 * np.sin(self.k*self.s1)))
            # eq (18):
            B = 1j * ((H[:, self.mic_channels[1]] * np.exp(+1j*self.k*(self.l1+self.s1)) -
                       H[:, self.mic_channels[0]] * np.exp(+1j*self.k*(self.l1))) /
                      (2 * np.sin(self.k*self.s1)))
            # eq (19):
            C = 1j * ((H[:, self.mic_channels[2]] * np.exp(+1j*self.k*(self.l2+self.s2)) -
                       H[:, self.mic_channels[3]] * np.exp(+1j*self.k*(self.l2))) /
                      (2 * np.sin(self.k*self.s2)))
            # eq (20):
            D = 1j * ((H[:, self.mic_channels[3]] * np.exp(-1j*self.k*(self.l2)) -
                       H[:, self.mic_channels[2]] * np.exp(-1j*self.k*(self.l2+self.s2))) /
                      (2 * np.sin(self.k*self.s2)))

            # Calculate acoustic pressure and velocity on both faces of specimen:
            p0 = A + B
            pd = C * np.exp(-1j*self.k*self.d) + D * np.exp(1j*self.k*self.d)
            u0 = (A-B) / (self.rho * self.c)
            ud = ((C * np.exp(-1j*self.k*self.d) - D * np.exp(1j*self.k*self.d)) /
                  (self.rho * self.c))

        # calculate Transfer Matrix:
        T = np.zeros(shape=(H.shape[0], 2, 2), dtype=complex)
        T[:, 0, 0] = (pd*ud + p0*u0) / (p0*ud + pd*u0)
        T[:, 1, 0] = (u0**2 - ud**2) / (p0*ud + pd*u0)
        T[:, 0, 1] = (p0**2 - pd**2) / (p0*ud + pd*u0)
        T[:, 1, 1] = T[:, 0, 0]

        return T

    def _get_transmission_coefficient(self):
        """Calculates transmission coefficient
        See 8.5.5.1, eq (25)

        Returns:
            [array]: Transmission coefficient
                     size: (f x 1), f: frequencies
        """
        # Transfer Matrix:
        T = self.transfer_matrix_one_load

        # Transmission Coefficient (anechoic backed) (eq (25)):
        t = (2 * np.exp(1j*self.k*self.d) /
                (T[:, 0, 0] + 
                 T[:, 0, 1] / (self.rho * self.c) +
                 T[:, 1, 0] * (self.rho * self.c) + 
                 T[:, 1, 1]))
        return t

    def _get_transmission_loss(self):
        """Calculates Normal Incidence Transmission Loss
        See 8.5.5.2, eq (26)

        Returns:
            [array]: Transmission loss
                     size: (f x 1), f: frequencies
        """
        TL = 20*np.log10(np.absolute(1/self.transmission_coefficient))
        return TL

    def _get_working_frequency_range(self):
        """Calculates the lower and upper frequency limit of the tube. 
        See 6.2 (eq. (1) and (2)) and 6.5.3 (eq. (3))

        Returns:
            tuple: lower and upper frequency limit
        """
        # distance between microphones:
        s = min(self.s1, self.s2)
        # Lower frequency limit:
        # NOTE: in ISO 10534-2, 5% is recommended instead of 1%
        f_lower = 0.05 * self.c / s

        # Upper frequency limit due to modes:
        if self.tube_shape == 'rect':
            K = 0.5
        elif self.tube_shape == 'circ':
            K = 0.586

        # eq. (2)
        f_upper_modes = K * self.c / self.tube_d

        # Upper frequency limit due to mic spacing (cf. 6.5.4):
        s = max(self.s1, self.s2)
        f_upper_spacing = 0.8 * self.c / (2 * s)

        f_upper = min(f_upper_modes, f_upper_spacing)

        return f_lower, f_upper
    
    def _get_z(self):
        """Calculate Characteristic Impedance in material
        See 8.5.5.6 (eq. (30))
        Returns:
            [array]: Characteristic Impedance
                     size: (f x 1), f: frequencies
        """
        T = self.transfer_matrix_one_load
        z = np.sqrt(T[:,0,1]/T[:,1,0])
        return z


class MicSwitchCalib(HasPrivateTraits):
    
    #: The :class:`~acoular.spectra.PowerSpectra` object that provides the csm.
    freq_data = Trait(PowerSpectra,
                      desc="power spectra object")
    
    #: The :class:`~acoular.spectra.PowerSpectra` object that provides the csm in the switched positions.
    freq_data_switched = Trait(PowerSpectra,
                               desc="power spectra object")
    
    ref_channel = Int(0, desc='reference channel')
    test_channel = Int(1, desc='Channel of the mic that is calibrated')

    H_c = Property()
    
    def _get_H_c(self):
        """Calculates the Amplitude/Phase correction transfer function for two 
        microphones, see 8.4.5. See README.md for illustrations.
        For example: Correction between mic 1 (pos 1, index 0) and mic 2 (pos 2, index 1)
                     Direct configuration:   mic 1 at pos 1 and cable going into 1st input
                                             mic 2 at pos 2 and cable going into 2nd input
                     Switched configuration: mic 1 at pos 2 and cable going into 1st input
                                             mic 2 at pos 1 and cable going into 2nd input
        Returns:
            [array]:    Correction Transfer function 
                        size: (f,), f: frequencies
        """
        #Get CSMs
        csm = self.freq_data.csm
        csm_switched = self.freq_data_switched.csm
        
        # Allocate space
        # H_1 = H_2 = H_c = np.zeros((csm.shape[0]), dtype=complex)  # create empty array
        
        H_1 = csm[:, self.test_channel, self.ref_channel] / \
              csm[:, self.ref_channel , self.ref_channel]
              
        H_2 = csm_switched[:, self.test_channel, self.ref_channel] / \
              csm_switched[:, self.ref_channel , self.ref_channel]
        
        H_c = np.sqrt( H_1 * H_2 )
        
        return H_c
    
