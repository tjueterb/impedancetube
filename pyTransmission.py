# This defines a class for transmission factor measurement

from acoular.spectra import PowerSpectra
import numpy as np
from traits.api import HasPrivateTraits, Property, Float, Trait, Delegate, Int, List


class Measurement(HasPrivateTraits):
    temperature = Float(default_value=20.)
    atmospheric_pressure = Float(default_value=101.325)
    c = Property(depends_on=['temperature'],
                 desc='Speed of sound')
    air_density = Property(depends_on=['temperature', 'atmospheric_pressure'],
                           desc='air density')

    # Tube dimensions
    l1 = Float(0.3, desc='distance between beginning of speciman and mic 2')
    l2 = Float(0.8, desc='distance between beginning of specimen and mic 3')
    s1 = Trait(0.085, 0.5, desc='Distance between mic 1 and 2 in m')
    s2 = Trait(0.085, 0.5, desc='Distance between mic 3 and 4 in m')
    d = Float(0.5, desc='length of test specimen (test tube section is 0.7m)')

    # channels of the microphones in the given freq_data object
    ref_channel = Int(desc="Channel index of the reference mic")
    mic_channels = List(minlen=4, maxlen=4,
                        desc="Channel indices of mics in positions 1-4")

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

    def _get_air_density(self):
        """Calculates air density in kg/m^3
        See 8.3 (eq. (5))
        Args:
            P (float): Atmospheric pressure in kPa. 
            T (float): Room remperature in °C. 
        Returns:
            rho (float): Air density in kg/m^3
        """
        rho = 1.290 * (self.atmospheric_pressure/101.325) * \
            (273.15 / (273.15 + self.temperature))
        return rho

    def _get_transfer_function(self):
        """Calculates the transfer function for all microphones
        See 8.5.1 eq (15)

        Args:
            ps (PowerSpectra): acoular power spectra object of measurement
            i_ref (int): index of the reference microphone 

        Returns:
            [array]:    Transfer function for all microphones 
                        size: (f x n), f: frequencies as in ps, n: channels
        """
        csm = self.freq_data.csm

        H_n_ref = np.empty(csm.shape[0:2], dtype=complex)  # create empty array

        for n in range(self.freq_data.numchannels):
            H_n_ref[:, n] = csm[:, n, self.ref_channel] / \
                csm[:, self.ref_channel, self.ref_channel]  # eq (15)
        return H_n_ref

    def _get_k(self):
        k = 2*np.pi * self.freq_data.fftfreq() / self.c
        return k

    def _get_transfer_matrix_one_load(self):
        # Get transfer function:
        H = self.transfer_function

        # Disable divide by zero warning because first entry of k is always 0
        with np.errstate(divide='ignore', invalid='ignore'):
            # Decompose wave field:
            # eq (17):
            A = 1j * (H[:, self.mic_channels[0]] * np.exp(-1j*self.k*(self.l1)) -
                      H[:, self.mic_channels[1]] * np.exp(-1j*self.k*(self.l1+self.s1))) /   \
                (2 * np.sin(self.k*self.s1))
            # eq (18):
            B = 1j * (H[:, self.mic_channels[1]] * np.exp(+1j*self.k*(self.l1+self.s1)) -
                      H[:, self.mic_channels[0]] * np.exp(+1j*self.k*(self.l1))) /      \
                (2 * np.sin(self.k*self.s1))
            # eq (19):
            C = 1j * (H[:, self.mic_channels[2]] * np.exp(+1j*self.k*(self.l2+self.s2)) -
                      H[:, self.mic_channels[3]] * np.exp(+1j*self.k*(self.l2))) /      \
                (2 * np.sin(self.k*self.s2))
            # eq (20):
            D = 1j * (H[:, self.mic_channels[3]] * np.exp(-1j*self.k*(self.l2)) -
                      H[:, self.mic_channels[2]] * np.exp(-1j*self.k*(self.l2+self.s2))) /   \
                (2 * np.sin(self.k*self.s2))

            # Calculate acoustic pressure and velocity on both faces of specimen:
            p0 = A + B
            pd = C * np.exp(-1j*self.k*self.d) + D * np.exp(1j*self.k*self.d)
            u0 = (A-B) / (self.air_density * self.c)
            ud = (C * np.exp(-1j*self.k*self.d) - D * np.exp(1j*self.k*self.d)) / \
                (self.air_density * self.c)

        # calculate Transfer Matrix:
        T = np.zeros(shape=(H.shape[0], 2, 2), dtype=complex)
        T[:, 0, 0] = (pd*ud + p0*u0) / (p0*ud + pd*u0)
        T[:, 1, 0] = (u0**2 - ud**2) / (p0*ud + pd*u0)
        T[:, 0, 1] = (p0**2 - pd**2) / (p0*ud + pd*u0)
        T[:, 1, 1] = T[:, 0, 0]

        return T

    def _get_transmission_coefficient(self):

        # Air Density:
        rho = self.air_density
        # Transfer Matrix:
        T = self.transfer_matrix_one_load

        # Transmission Coefficient (anechoic backed) (eq (25)):
        t = 2 * np.exp(1j*self.k*self.d) / \
            (T[:, 0, 0] + T[:, 0, 1] / (rho * self.c) +
             T[:, 1, 0]*(rho*self.c) + T[:, 1, 1])
        return t

    def _get_transmission_loss(self):
        # Normal incidence Transmission loss (eq. (26)):
        TL = 20*np.log10(np.absolute(1/self.transmission_coefficient))
        return TL

    def _get_working_frequency_range(self):
        """Calculates the lower and upper frequency limit of the tube. 
        See chapter 6.2.2 (Eq 1 and 2)

        Args:
            x (float):      spacing between microphones in m
            d (float):      diameter (if circular) or 
                            largest section dimension (if rectangular)
                            of the tube (in m)
            shape (string, optional):   'circ' for circular tube
                                        'rect' for rectangular tube
            c (float, optional):        speed of sound in m/s, default = 343

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
