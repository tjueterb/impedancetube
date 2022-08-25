# This defines a class for transmission factor measurement

from acoular.sources import TimeSamples
from acoular.spectra import PowerSpectra
import numpy as np
from numpy import exp, sin, pi, sqrt
from traits.api import HasPrivateTraits, Property, Float, Trait, Delegate, Int, List, Array

from .tube import Tube_Transmission


class Measurement(HasPrivateTraits):
    temperature = Float(default_value=20.)
    atmospheric_pressure = Float(default_value=101.325)
    c = Property(depends_on=['temperature'],
                 desc='Speed of sound')
    rho = Property(depends_on=['temperature', 'atmospheric_pressure'],
                   desc='air density')

    # Block size
    block_size = Delegate('freq_data')

    #: The :class:`~acoular.spectra.PowerSpectra` object that provides the csm.
    freq_data = Trait(PowerSpectra,
                      desc="power spectra object")

    tube = Trait(Tube_Transmission)
    
    # wave number vector:
    k = Property(depends_on=['freq_data'])

    # Working frequency range:
    working_frequency_range = Property()


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

    def _get_k(self):
        """Calculates the wave number coefficients for all frequencies
        See 3.2: k = 2*pi*f / c
        
        Update: Damping constant in Neper/Meter estimated to get a complex wave number 
        (Makes virtually no difference)
        
        Returns:
            [array]: size: (f x 1)
        """
        k = (2 * np.pi * self.freq_data.fftfreq() +
             1j * 0.0194*(np.sqrt(self.freq_data.fftfreq()))/self.tube.tube_d) / self.c
        return k

    def _get_working_frequency_range(self):
        """Calculates the lower and upper frequency limit of the tube. 
        See 6.2 (eq. (1) and (2)) and 6.5.3 (eq. (3))

        Returns:
            tuple: lower and upper frequency limit
        """
        # distance between microphones:
        s = min(self.tube.s1, self.tube.s2)
        # Lower frequency limit:
        # NOTE: in ISO 10534-2, 5% is recommended instead of 1%
        f_lower = 0.05 * self.c / s

        # Upper frequency limit due to modes:
        if self.tube.tube_shape == 'rect':
            K = 0.5
        elif self.tube.tube_shape == 'circ':
            K = 0.586

        # eq. (2)
        f_upper_modes = K * self.c / self.tube.tube_d

        # Upper frequency limit due to mic spacing (cf. 6.5.4):
        s = max(self.tube.s1, self.tube.s2)
        f_upper_spacing = 0.8 * self.c / (2 * s)

        f_upper = min(f_upper_modes, f_upper_spacing)

        return f_lower, f_upper


class Measurement_E2611(Measurement):
    '''
    Transfer Matrix Method following the E2611 Norm, with one or two loads.
    
    For the two load case, the 'method' trait needs to be set to 'two load',
    and a frequency data of the second load case needs to be passed as freq_data_two_load.
    '''
    # Tube dimensions
    '''
    Values for the TAP ducts:
    Rect: l1 = 0.3, l2 = 0.8,  s1,s2 = 0.085 or 0.5
    circ: l1 = 0.2, l2 = 0.28, s1,s2 = 0.075 or 0.225
    '''
    
    # channels of the microphones in the given freq_data object
    ref_channel = Int(0, desc="Channel index of the reference mic")
    mic_channels = List([1, 2, 3, 4], minlen=4, maxlen=4,
                        desc="Channel indices of mics in positions 0-3")

    #: The :class:`~acoular.spectra.PowerSpectra` object that provides the csm for the second load coase.
    freq_data_two_load = Trait(PowerSpectra,
                               desc="power spectra object")

    # The transfer functions for all microphones:
    transfer_function = Property(
        desc='Transfer function between all mics and ref. mic (channel i_ref)')
    transfer_function_two_load = Property(
        desc='Transfer function between all mics and ref. mic (channel i_ref) for the second load case')

    # The transfer Matrix
    transfer_matrix = Property(desc='Transfer Matrix (dependent on load case)')
    transfer_matrix_one_load = Property(desc='Transfer Matrix for one load')
    transfer_matrix_two_load = Property(desc='Transfer Matrix for two loads')

    # transmission coefficient:
    transmission_coefficient = Property()
    
    # Transmission loss:
    transmission_loss = Property()

    # Characteristic Impedance in material:
    z = Property()
    
    # Reflection coefficient:
    reflection_coefficient = Property()

    # Reflection coefficient (hard backed):
    reflection_coefficient_hard_backed = Property()

    # Absorption coefficient (anechoic):
    absorption_coefficient = Property()

    # Absorption coefficient (hard backed):
    absorption_coefficient_hard_backed = Property()

    # Propagation wavenumber in material:
    propagation_wavenumber = Property()
    
    # Choose one load or two load method:
    method = Trait('one load', 'two load',
                   default_value='one load',
                   desc="Method for calculating the transfer matrix")

    # Amplitude and Phase correction Transfer function
    H_c = Array()
    #TODO: either Initialize with ones (only possible after freq_data shape is known)
    #TODO: Or rewrite the calib class so the entire object can be handed to this class

    def _get_transfer_function(self):
        """Calculates the transfer function for all microphones
        See 8.5.1 eq (15)
              
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

    def _get_transfer_function_two_load(self):
        """Calculates the transfer function for all microphones for the second load case
        See 8.5.1 eq (15)
              
        Returns:
            [array]:    Transfer function for all microphones 
                        size: (f x n), f: frequencies, n: channels
        """
        csm = self.freq_data_two_load.csm

        H_n_ref = np.empty(csm.shape[0:2], dtype=complex)  # create empty array

        for n in range(self.freq_data_two_load.numchannels):
            #NOTE: correct indexing of the csm verified with Adam
            H_n_ref[:, n] = csm[:, n, self.ref_channel] / \
                csm[:, self.ref_channel, self.ref_channel]  # eq (15)

        # apply correction transfer function:
        H_n_ref = H_n_ref / self.H_c

        return H_n_ref

    def _get_transfer_matrix(self):
        """gets the correct transfer matrix, 
        depending on the load case specified in self.method
        """

        if self.method == 'one load':
            return self.transfer_matrix_one_load
        elif self.method == 'two load':
            return self.transfer_matrix_two_load

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
            A = 1j * ((H[:, self.mic_channels[0]] * np.exp(-1j*self.k*(self.tube.l1)) -
                       H[:, self.mic_channels[1]] * np.exp(-1j*self.k*(self.tube.l1+self.tube.s1))) /
                      (2 * np.sin(self.k*self.tube.s1)))
            # eq (18):
            B = 1j * ((H[:, self.mic_channels[1]] * np.exp(+1j*self.k*(self.tube.l1+self.tube.s1)) -
                       H[:, self.mic_channels[0]] * np.exp(+1j*self.k*(self.tube.l1))) /
                      (2 * np.sin(self.k*self.tube.s1)))
            # eq (19):
            C = 1j * ((H[:, self.mic_channels[2]] * np.exp(+1j*self.k*(self.tube.l2+self.tube.s2)) -
                       H[:, self.mic_channels[3]] * np.exp(+1j*self.k*(self.tube.l2))) /
                      (2 * np.sin(self.k*self.tube.s2)))
            # eq (20):
            D = 1j * ((H[:, self.mic_channels[3]] * np.exp(-1j*self.k*(self.tube.l2)) -
                       H[:, self.mic_channels[2]] * np.exp(-1j*self.k*(self.tube.l2+self.tube.s2))) /
                      (2 * np.sin(self.k*self.tube.s2)))

            # Calculate acoustic pressure and velocity on both faces of specimen:
            p0 = A + B
            pd = C * np.exp(-1j*self.k*self.tube.d) + D * np.exp(1j*self.k*self.tube.d)
            u0 = (A-B) / (self.rho * self.c)
            ud = ((C * np.exp(-1j*self.k*self.tube.d) - D * np.exp(1j*self.k*self.tube.d)) /
                  (self.rho * self.c))

        # calculate Transfer Matrix:
        T = np.zeros(shape=(H.shape[0], 2, 2), dtype=complex)
        T[:, 0, 0] = (pd*ud + p0*u0) / (p0*ud + pd*u0)
        T[:, 1, 0] = (u0**2 - ud**2) / (p0*ud + pd*u0)
        T[:, 0, 1] = (p0**2 - pd**2) / (p0*ud + pd*u0)
        T[:, 1, 1] = T[:, 0, 0]

        return T

    def _get_transfer_matrix_two_load(self):
        """Calculates the one load transfer matrix 
        See 8.5.4.2
        
        Returns:
            [array]: Transfer Matrix
                     size:(f x 2 x 2), f: frequencies
        """
        # Get transfer functions:
        H_a = self.transfer_function
        H_b = self.transfer_function_two_load

        #TODO: Refactor this

        # Disable divide by zero warning because first entry of k is always 0
        with np.errstate(divide='ignore', invalid='ignore'):
            # Decompose wave field for first load case:
            # eq (17):
            A_a = 1j * ((H_a[:, self.mic_channels[0]] * np.exp(-1j*self.k*(self.tube.l1)) -
                         H_a[:, self.mic_channels[1]] * np.exp(-1j*self.k*(self.tube.l1+self.tube.s1))) /
                        (2 * np.sin(self.k*self.tube.s1)))
            # eq (18):
            B_a = 1j * ((H_a[:, self.mic_channels[1]] * np.exp(+1j*self.k*(self.tube.l1+self.tube.s1)) -
                         H_a[:, self.mic_channels[0]] * np.exp(+1j*self.k*(self.tube.l1))) /
                        (2 * np.sin(self.k*self.tube.s1)))
            # eq (19):
            C_a = 1j * ((H_a[:, self.mic_channels[2]] * np.exp(+1j*self.k*(self.tube.l2+self.tube.s2)) -
                         H_a[:, self.mic_channels[3]] * np.exp(+1j*self.k*(self.tube.l2))) /
                        (2 * np.sin(self.k*self.tube.s2)))
            # eq (20):
            D_a = 1j * ((H_a[:, self.mic_channels[3]] * np.exp(-1j*self.k*(self.tube.l2)) -
                         H_a[:, self.mic_channels[2]] * np.exp(-1j*self.k*(self.tube.l2+self.tube.s2))) /
                        (2 * np.sin(self.k*self.tube.s2)))

            # Calculate acoustic pressure and velocity on both faces of specimen:
            p0_a = A_a + B_a
            pd_a = C_a * np.exp(-1j*self.k*self.tube.d) + D_a * np.exp(1j*self.k*self.tube.d)
            u0_a = (A_a-B_a) / (self.rho * self.c)
            ud_a = ((C_a * np.exp(-1j*self.k*self.tube.d) - D_a * np.exp(1j*self.k*self.tube.d)) /
                    (self.rho * self.c))

            # Decompose wave field for second load case:
            # eq (17):
            A_b = 1j * ((H_b[:, self.mic_channels[0]] * np.exp(-1j*self.k*(self.tube.l1)) -
                         H_b[:, self.mic_channels[1]] * np.exp(-1j*self.k*(self.tube.l1+self.tube.s1))) /
                        (2 * np.sin(self.k*self.tube.s1)))
            # eq (18):
            B_b = 1j * ((H_b[:, self.mic_channels[1]] * np.exp(+1j*self.k*(self.tube.l1+self.tube.s1)) -
                         H_b[:, self.mic_channels[0]] * np.exp(+1j*self.k*(self.tube.l1))) /
                        (2 * np.sin(self.k*self.tube.s1)))
            # eq (19):
            C_b = 1j * ((H_b[:, self.mic_channels[2]] * np.exp(+1j*self.k*(self.tube.l2+self.tube.s2)) -
                         H_b[:, self.mic_channels[3]] * np.exp(+1j*self.k*(self.tube.l2))) /
                        (2 * np.sin(self.k*self.tube.s2)))
            # eq (20):
            D_b = 1j * ((H_b[:, self.mic_channels[3]] * np.exp(-1j*self.k*(self.tube.l2)) -
                         H_b[:, self.mic_channels[2]] * np.exp(-1j*self.k*(self.tube.l2+self.tube.s2))) /
                        (2 * np.sin(self.k*self.tube.s2)))

            # Calculate acoustic pressure and velocity on both faces of specimen:
            p0_b = A_b + B_b
            pd_b = C_b * np.exp(-1j*self.k*self.tube.d) + D_b * np.exp(1j*self.k*self.tube.d)
            u0_b = (A_b-B_b) / (self.rho * self.c)
            ud_b = ((C_b * np.exp(-1j*self.k*self.tube.d) - D_b * np.exp(1j*self.k*self.tube.d)) /
                    (self.rho * self.c))

            # calculate Transfer Matrix:
            T = np.zeros(shape=(H_a.shape[0], 2, 2), dtype=complex)
            T[:, 0, 0] = (p0_a*ud_b - p0_b*ud_a) / (pd_a*ud_b - pd_b*ud_a)
            T[:, 1, 0] = (u0_a*ud_b - u0_b*ud_a) / (pd_a*ud_b - pd_b*ud_a)
            T[:, 0, 1] = (p0_b*pd_a - p0_a*pd_b) / (pd_a*ud_b - pd_b*ud_a)
            T[:, 1, 1] = (pd_a*u0_b - pd_b*u0_a) / (pd_a*ud_b - pd_b*ud_a)

        return T

    def _get_transmission_coefficient(self):
        """Calculates transmission coefficient
        See 8.5.5.1, eq (25)

        Returns:
            [array]: Transmission coefficient
                     size: (f x 1), f: frequencies
        """
        # Transfer Matrix:
        T = self.transfer_matrix

        # Transmission Coefficient (anechoic backed) (eq (25)):
        # Disable divide by zero warning because first entry of k is always 0
        with np.errstate(divide='ignore', invalid='ignore'):
            t = (2 * np.exp(1j*self.k*self.tube.d) /
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

    def _get_reflection_coefficient(self):
        """Calculates reflection coefficient
        See Song & Bolton (2000), eq. (9)

        Returns:
            [array]: reflection coefficient
                     size: (f x 1), f: frequencies
        """
        # Transfer Matrix:
        T = self.transfer_matrix

        with np.errstate(divide='ignore', invalid='ignore'):
            R = (T[:, 0, 0] + T[:, 0, 1]/(self.rho*self.c) - self.rho*self.c * T[:, 1, 0] - T[:, 1, 1]) / \
                (T[:, 0, 0] + T[:, 0, 1]/(self.rho*self.c) +
                 self.rho*self.c * T[:, 1, 0] + T[:, 1, 1])
        return R

    def _get_reflection_coefficient_hard_backed(self):
        """Calculates reflection coefficient (hard backed)
        See 8.5.5.3, eq (27)

        Returns:
            [array]: reflection coefficient
                     size: (f x 1), f: frequencies
        """
        # Transfer Matrix:
        T = self.transfer_matrix

        # Reflection coefficient (hard backed) (eq (27)):
        # Disable divide by zero warning because first entry of k is always 0
        with np.errstate(divide='ignore', invalid='ignore'):
            R = ((T[:, 0, 0] - (self.rho * self.c * T[:, 1, 0])) /
                 (T[:, 0, 0] + (self.rho * self.c * T[:, 1, 0])))
        return R

    def _get_absorption_coefficient(self):
        """Calculates absorption coefficient (anechoic backed)
        See 8.5.5.4, eq (28)

        Returns:
            [array]: absorption coefficient
                     size: (f x 1), f: frequencies
        """

        # Absorption coefficient (hard backed) (eq (28)):
        alpha = 1 - np.abs(self.reflection_coefficient)**2
        return alpha

    def _get_absorption_coefficient_hard_backed(self):
        """Calculates absorption coefficient (hard backed)
        See 8.5.5.4, eq (28)

        Returns:
            [array]: absorption coefficient
                     size: (f x 1), f: frequencies
        """

        # Absorption coefficient (hard backed) (eq (28)):
        alpha = 1 - np.abs(self.reflection_coefficient_hard_backed)**2
        return alpha

    def _get_propagation_wavenumber(self):
        """Calculates propagation_wavenumber in the material
        See 8.5.5.5, eq (29)

        Returns:
            [array]: propagation wavenumber
                     size: (f x 1), f: frequencies
        """
        # Transfer Matrix:
        T = self.transfer_matrix

        # Disable divide by zero warning because first entry of k is always 0
        with np.errstate(divide='ignore', invalid='ignore'):
            return 1 / self.tube.d * np.arccos(T[:, 0, 0])

    def _get_z(self):
        """Calculate Characteristic Impedance in material
        See 8.5.5.6 (eq. (30))
        Returns:
            [array]: Characteristic Impedance
                     size: (f x 1), f: frequencies
        """
        # Transfer Matrix:
        T = self.transfer_matrix

        z = np.sqrt(T[:, 0, 1]/T[:, 1, 0])
        return z


class MicSwitchCalib_E2611(HasPrivateTraits):

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
              csm[:, self.ref_channel, self.ref_channel]

        H_2 = csm_switched[:, self.test_channel, self.ref_channel] / \
              csm_switched[:, self.ref_channel, self.ref_channel]

        H_c = np.sqrt(H_1 * H_2)

        return H_c
