# This defines a class for transmission factor measurement

from acoular.sources import TimeSamples
from acoular.spectra import PowerSpectra
import numpy as np
from numpy import exp, sin, pi, sqrt
from traits.api import HasPrivateTraits, Property, Float, Trait, Delegate, Int, List, Array


class Measurement(HasPrivateTraits):
    temperature = Float(default_value=20.)
    atmospheric_pressure = Float(default_value=101.325)
    c = Property(depends_on=['temperature'],
                 desc='Speed of sound')
    rho = Property(depends_on=['temperature', 'atmospheric_pressure'],
                   desc='air density')

    # Tube dimensions
    '''
    Values for the TAP ducts:
    Rect: l1 = 0.3, l2 = 0.8,  s1,s2 = 0.085 or 0.5
    circ: l1 = 0.2, l2 = 0.28, s1,s2 = 0.075 or 0.225
    '''
    l1 = Float(0.3, desc='distance between beginning of specimen and mic 2')
    l2 = Float(0.8, desc='distance between beginning of specimen and mic 3')
    s1 = Float(0.085,desc='Distance between mic 1 and 2 in m')
    s2 = Float(0.085,desc='Distance between mic 3 and 4 in m')
    d = Float(0.5, desc='length of test specimen (test tube section is 0.7m)')

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

    # transmission coefficient:
    transmission_coefficient = Property()

    # Transmission loss:
    transmission_loss = Property()

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
    
    # Working frequency range:
    working_frequency_range = Property()
    
    # Characteristic Impedance in material:
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

    def _get_k(self):
        """Calculates the wave number coefficients for all frequencies
        See 3.2: k = 2*pi*f / c
        
        Update: Damping constant in Neper/Meter estimated to get a complex wave number 
        (Makes virtually no difference)
        
        Returns:
            [array]: size: (f x 1)
        """
        k = (2 * np.pi *         self.freq_data.fftfreq() + 
             1j* 0.0194*(np.sqrt(self.freq_data.fftfreq()))/self.tube_d ) / self.c
        return k
    
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
    
class Measurement_E2611(Measurement):
    '''
    Transfer Matrix Method following the E2611 Norm
    '''
    # channels of the microphones in the given freq_data object
    ref_channel = Int(0, desc="Channel index of the reference mic")
    mic_channels = List([1, 2, 3, 4], minlen=4, maxlen=4,
                        desc="Channel indices of mics in positions 0-3")
    
    # The transfer function for all microphones:
    transfer_function = Property(
        desc='Transfer function between all mics and ref. mic (channel i_ref)')

    # The transfer Matrix
    transfer_matrix_one_load = Property(desc='Transfer Matrix for one load')
    
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
        # Disable divide by zero warning because first entry of k is always 0
        with np.errstate(divide='ignore', invalid='ignore'):
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

    def _get_reflection_coefficient(self):
        """Calculates reflection coefficient
        See Song & Bolton (2000), eq. (9)

        Returns:
            [array]: reflection coefficient
                     size: (f x 1), f: frequencies
        """
        T = self.transfer_matrix_one_load
        with np.errstate(divide='ignore', invalid='ignore'):
            R = (T[:,0,0] + T[:,0,1]/(self.rho*self.c) - self.rho*self.c * T[:,1,0] - T[:,1,1]) / \
                (T[:,0,0] + T[:,0,1]/(self.rho*self.c) + self.rho*self.c * T[:,1,0] + T[:,1,1])
        return R
        
    def _get_reflection_coefficient_hard_backed(self):
        """Calculates reflection coefficient (hard backed)
        See 8.5.5.3, eq (27)

        Returns:
            [array]: reflection coefficient
                     size: (f x 1), f: frequencies
        """
        # Transfer Matrix:
        T = self.transfer_matrix_one_load

        # Reflection coefficient (hard backed) (eq (27)):
        # Disable divide by zero warning because first entry of k is always 0
        with np.errstate(divide='ignore', invalid='ignore'):
            R = ((T[:, 0, 0] - (self.rho * self.c * T[:,1,0])) / 
                 (T[:, 0, 0] + (self.rho * self.c * T[:,1,0])))
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
        # Disable divide by zero warning because first entry of k is always 0
        with np.errstate(divide='ignore', invalid='ignore'):
            return 1 / self.d * np.arccos( self.transfer_matrix_one_load[:,0,0] )
    
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
              csm[:, self.ref_channel , self.ref_channel]
              
        H_2 = csm_switched[:, self.test_channel, self.ref_channel] / \
              csm_switched[:, self.ref_channel , self.ref_channel]
        
        H_c = np.sqrt( H_1 * H_2 )
        
        return H_c
    
class Measurement_Cottbus(Measurement):
    ''' 
    Measurement with Transfer Matrix Method Cottbus style.
    
    Expected Microphone Correction transfer function:
    H_c: [f x 3]    f: frequencies
    1st column: Correction Microphones 0 and 1
    2nd column: Correction Microphones 0 and 2
    3rd column: Correction Microphones 2 and 3
    
    '''
    
    # Amplitude and Phase correction Transfer function
    H_c = Array(shape=(None,3),desc='Required shape: see docstring')
    
    # channels of the microphones in the given freq_data object
    mic_channels = List([1, 2, 3, 4], minlen=4, maxlen=4,
                        desc="Channel indices of mics in positions 0-3")
    
    # The transfer function for all microphones:
    transfer_function = Property(
        desc='Transfer function between Mics 0->1, 0->2, 2->3')
    
    # The transfer Matrix
    transfer_matrix_one_load = Property(desc='Transfer Matrix for one load')
    
    def _get_transfer_function(self):
        """Calculates the transfer function for all microphones
        1st column: Microphones 0 and 1
        2nd column: Microphones 0 and 2
        3rd column: Microphones 2 and 3              
        Returns:
            [array]:    Transfer function for all microphones 
                        size: (f x 3), f: frequencies 
        """
        csm = self.freq_data.csm

        # uncorrected transfer functions:
        Sxx_01 = csm[:,self.mic_channels[0],self.mic_channels[0]]         #Autospektraldichte 
        Sxy_01 = csm[:,self.mic_channels[0],self.mic_channels[1]]         #Kreuzspektraldichte
        H01_unkorr = Sxy_01/Sxx_01  #Übertragungsfunktion der Kanäle 0 und 1 aus Messung (m) --> unkorrigiertes Signal

        Sxx_02 = Sxx_01             #Autospektraldichte
        Sxy_02 = csm[:,self.mic_channels[0],self.mic_channels[2]]         #Kreuzspektraldichte
        H02_unkorr = Sxy_02/Sxx_02  #Übertragungsfunktion der Kanäle 0 und 2 aus Messung (m) --> unkorrigiertes Signal  

        Sxx_23 = csm[:,self.mic_channels[2],self.mic_channels[2]]           #Autospektraldichte
        Sxy_23 = csm[:,self.mic_channels[2],self.mic_channels[3]]           #Kreuzspektraldichte
        H23_unkorr = Sxy_23/Sxx_23    #Übertragungsfunktion der Kanäle 3 und 4 aus Messung (m) --> unkorrigiertes Signal   

        # Phase/Amplitude corrected transfer functions:
        Hkorr_01, Hkorr_02, Hkorr_23 = self.H_c[:,0], self.H_c[:,1], self.H_c[:,2]
        H01 = H01_unkorr/Hkorr_01
        H02 = H02_unkorr/Hkorr_02
        H23 = H23_unkorr/Hkorr_23
        
        H = np.empty((csm.shape[0],3), dtype=complex)  # create empty array
        H[:,0] = H01
        H[:,1] = H02
        H[:,2] = H23
        
        return H
    
    def _get_transfer_matrix_one_load(self):
        """ 
        Calculate Transfer Matrix
        The variable names are numbered from 1 to 4 in order to be identical to the old script.
        renaming to indices 0-3 will be done as soon as the code is tested. 

        Returns:
            [array]: Transfer Matrix
                     size: (f x 2 x 2), f: frequencies 
        """         
              
        H12 = self.transfer_function[:,0]
        H13 = self.transfer_function[:,1]
        H34 = self.transfer_function[:,2]
        
        x1 = self.l1 + self.s1
        x2 = self.l1
        x3 = self.l2
        x4 = self.l2 + self.s2
        
        # Disable divide by zero warning because first entry of k is always 0
        with np.errstate(divide='ignore', invalid='ignore'):
            # Ermittlung von Reflektions- und Transmissionsfaktor (Vgl. Bolton/Song und Extrablatt)
            E = (1j*(np.exp(1j*self.k*x2)-H12*np.exp(1j*self.k*x1))) / \
                (2*np.sin(self.k*(x1-x2)))
            F = (1j*(H12*np.exp(-1*1j*self.k*x1)-np.exp(-1*1j*self.k*x2))) / \
                (2*np.sin(self.k*(x1-x2)))
            G = (1j*(np.exp(1j*self.k*x4)-H34*np.exp(1j*self.k*x3)))/ \
                (2*np.sin(self.k*(x3-x4)))
            H = (1j*(H34*np.exp(-1*1j*self.k*x3)-np.exp(-1*1j*self.k*x4))) / \
                (2*np.sin(self.k*(x3-x4)))
        
        P0 = E + F
        V0 = (E-F)/(self.rho*self.c)
        Pd = G*np.exp(-1*1j*self.k*self.d)+H*np.exp(1j*self.k*self.d)
        Vd = (G*np.exp(-1*1j*self.k*self.d)-H*np.exp(1j*self.k*self.d))/(self.rho*self.c)

        T11 = ((H13*Pd*Vd)+(P0*V0/H13))/((P0*Vd)+(Pd*V0))
        T12 = ((P0*P0/H13)-(H13*Pd*Pd))/((P0*Vd)+(Pd*V0))
        T21 = ((V0*V0/H13)-(H13*Vd*Vd))/((P0*Vd)+(Pd*V0))
        T22 = T11
        
        # Transfer Matrix:
        T = np.zeros(shape=(H12.shape[0], 2, 2), dtype=complex)
        T[:,0,0] = T11
        T[:,0,1] = T12
        T[:,1,0] = T21
        T[:,1,1] = T22
        
        return T
    
    def _get_transmission_coefficient(self):
        """Calculates transmission coefficient
        Returns:
            [array]: Transmission coefficient
                     size: (f x 1), f: frequencies
        """
        # Transfer Matrix:
        T = self.transfer_matrix_one_load

        # Transmission Coefficient (anechoic backed) (eq (25)):
        # Disable divide by zero warning because first entry of k is always 0
        with np.errstate(divide='ignore', invalid='ignore'):
            t = (2 * np.exp(1j*self.k*self.d) /
                 (T[:, 0, 0] + 
                  T[:, 0, 1] / (self.rho * self.c) +
                  T[:, 1, 0] * (self.rho * self.c) + 
                  T[:, 1, 1]))
        return t

    def _get_transmission_loss(self):
        """Calculates Normal Incidence Transmission Loss

        Returns:
            [array]: Transmission loss
                     size: (f x 1), f: frequencies
        """
        TL = 20*np.log10(np.absolute(1/self.transmission_coefficient))
        return TL

    def _get_reflection_coefficient(self):
        """Calculates reflection coefficient (hard backed)
        See 8.5.5.3, eq (27)

        Returns:
            [array]: reflection coefficient
                     size: (f x 1), f: frequencies
        """
        # Transfer Matrix:
        T = self.transfer_matrix_one_load
        T11 = T[:,0,0]
        T12 = T[:,0,1]
        T21 = T[:,1,0]
        T22 = T[:,1,1]
        # Reflection coefficient (hard backed): #NOTE: THIS IS DIFFERENT!!!!!
        # Disable divide by zero warning because first entry of k is always 0
        with np.errstate(divide='ignore', invalid='ignore'):
            Ra = (T11+(T12/(self.rho*self.c))-((self.rho*self.c)*T21)-T22)/ \
                 (T11+(T12/(self.rho*self.c))+((self.rho*self.c)*T21)+T22)
        return Ra
 
    def _get_absorption_coefficient(self):
        """Calculate absorption coefficient (hard backed)
        """
        

        # Absorption coefficient (hard backed) (eq (28)):
        alpha = 1 - np.abs(self.reflection_coefficient)**2
        return alpha
    
    def _get_propagation_wavenumber(self):
        """Calculates propagation_wavenumber in the material
        NOTE: this is not included in the script and taken from E2611.
        Returns:
            [array]: propagation wavenumber
                     size: (f x 1), f: frequencies
        """
        # Disable divide by zero warning because first entry of k is always 0
        with np.errstate(divide='ignore', invalid='ignore'):
            return 1 / self.d * np.arccos( self.transfer_matrix_one_load[:,0,0] )
    
    def _get_z(self):
        """Calculate Characteristic Impedance in material
        Returns:
            [array]: Characteristic Impedance
                     size: (f x 1), f: frequencies
        """
        T = self.transfer_matrix_one_load
        z = np.sqrt(T[:,0,1]/T[:,1,0])
        return z
        
      
class MicSwitchCalib_Cottbus(HasPrivateTraits):
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
        microphones, 
        Returns:
            [array]:    Correction Transfer function 
                        size: (f,), f: frequencies
        """
        #Get CSMs
        csm = self.freq_data.csm
        csm_switched = self.freq_data_switched.csm

        H_a = csm[:, self.ref_channel, self.test_channel] / \
              csm[:, self.ref_channel , self.ref_channel]
              
        H_b = csm_switched[:, self.test_channel, self.test_channel] / \
              csm_switched[:, self.test_channel , self.ref_channel]
              
        H_c = np.sqrt( H_a * H_b )
        
        return H_c

class Measurement_Cottbus_OG(Measurement):
    ''' 
    Measurement with Transfer Matrix Method Cottbus style.
    '''
    # Calibration Data:
    freq_data_00_11_22_33 = Trait(PowerSpectra,
        desc = 'PowerSpectra Object containing empty measurement')
    freq_data_01_10_23_32 = Trait(PowerSpectra,
        desc = 'PowerSpectra Object containing empty measurement with mics 0/1, 2/3 switched')
    freq_data_02_11_20_33 = Trait(PowerSpectra,
        desc = 'PowerSpectra Object containing empty measurement with mics 0/2 switched')
    
    
    # channels of the microphones in the given freq_data object
    mic_channels = List([1, 2, 3, 4], minlen=4, maxlen=4,
                        desc="Channel indices of mics in positions 0-3")
    
    def _get_transmission_loss(self):
        '''transmission loss directly'''

        #Anpassung der Größen an Song und Bolton
        x1 = -1 * (self.l1+self.s1)                 #läuft in negative x-Richtung
        x2 = -1 * self.l1                           #läuft in negative x-Richtung
        x3 = self.l2
        x4 = self.l2 + self.s2
        d  = self.d

        c=self.c
        rho = self.rho
        k = self.k

        
        # Kalibrierungsdatei 1
        f = self.freq_data_00_11_22_33
        Sxx_12 = f.csm[:,self.mic_channels[0],self.mic_channels[0]]     #Autospektraldichte
        Sxy_12 = f.csm[:,self.mic_channels[0],self.mic_channels[1]]     #Kreuzspektraldichte
        Ha12 = Sxy_12/Sxx_12      #Übertragungsfunktion der Messung 00 und 11 --> Übertragungsfunktion des Prüflings im unkorrigierten Zustand

        Sxx_13 = Sxx_12           #Autospektraldichte
        Sxy_13 = f.csm[:,self.mic_channels[0],self.mic_channels[2]]     #Kreuzspektraldichte
        Ha13 = Sxy_13/Sxx_13      #Übertragungsfunktion der Messung 00 und 22 --> Übertragungsfunktion des Prüflings im unkorrigierten Zustand

        Sxx_34 = f.csm[:,self.mic_channels[2],self.mic_channels[2]]     #Autospektraldichte
        Sxy_34 = f.csm[:,self.mic_channels[2],self.mic_channels[3]]     #Kreuzspektraldichte
        Ha34 = Sxy_34/Sxx_34      #Übertragungsfunktion der Messung 22 und 33 --> Übertragungsfunktion des Prüflings im unkorrigierten Zustand

        # Kalibrierungsdatei 2
        f = self.freq_data_01_10_23_32
        Sxx_21 = f.csm[:,self.mic_channels[1],self.mic_channels[1]]
        Sxy_21 = f.csm[:,self.mic_channels[1],self.mic_channels[0]]
        Hb12 = Sxy_21/Sxx_21          #Übertragungsfunktion der Messung 01 und 10

        Sxx_43 = f.csm[:,self.mic_channels[3],self.mic_channels[3]]
        Sxy_43 = f.csm[:,self.mic_channels[3],self.mic_channels[2]]
        Hb34 = Sxy_43/Sxx_43          #Übertragungsfunktion der Messung 23 und 32

        # Kalibrierungsdatei 3
        f = self.freq_data_02_11_20_33
        Sxx_31 = f.csm[:,self.mic_channels[2],self.mic_channels[2]]
        Sxy_31 = f.csm[:,self.mic_channels[2],self.mic_channels[0]]
        Hb13 = Sxy_31/Sxx_31          #Übertragungsfunktion der Messung 02 und 20

        # Messdatei
        f = self.freq_data
        Sxx_12_m = f.csm[:,self.mic_channels[0],self.mic_channels[0]]         #Autospektraldichte 
        Sxy_12_m = f.csm[:,self.mic_channels[0],self.mic_channels[1]]         #Kreuzspektraldichte
        H12_unkorr = Sxy_12_m/Sxx_12_m  #Übertragungsfunktion der Kanäle 0 und 1 aus Messung (m) --> unkorrigiertes Signal

        Sxx_13_m = Sxx_12_m             #Autospektraldichte
        Sxy_13_m = f.csm[:,self.mic_channels[0],self.mic_channels[2]]         #Kreuzspektraldichte
        H13_unkorr = Sxy_13_m/Sxx_13_m  #Übertragungsfunktion der Kanäle 0 und 2 aus Messung (m) --> unkorrigiertes Signal  

        Sxx_34_m = f.csm[:,self.mic_channels[2],self.mic_channels[2]]           #Autospektraldichte
        Sxy_34_m = f.csm[:,self.mic_channels[2],self.mic_channels[3]]           #Kreuzspektraldichte
        H34_unkorr = Sxy_34_m/Sxx_34_m    #Übertragungsfunktion der Kanäle 3 und 4 aus Messung (m) --> unkorrigiertes Signal   

        # Korrekturfunktionen
        Hkorr_12 = np.sqrt(Ha12/Hb12)        #sqrt(abs(Ha/Hb))*exp(1j*(0.5*(angle(Ha)-angle(Hb))))
        Hkorr_13 = np.sqrt(Ha13/Hb13)
        Hkorr_34 = np.sqrt(Ha34/Hb34)    

        #korrigierte Übertragungsfunktionen
        H12 = H12_unkorr/Hkorr_12
        H13 = H13_unkorr/Hkorr_13
        H34 = H34_unkorr/Hkorr_34

        # Ermittlung von Reflektions- und Transmissionsfaktor (Vgl. Bolton/Song und Extrablatt)

        E = (1j*(exp(1j*k*x2)-H12*exp(1j*k*x1)))/(2*sin(k*(x1-x2)))
        F = (1j*(H12*exp(-1*1j*k*x1)-exp(-1*1j*k*x2)))/(2*sin(k*(x1-x2)))
        G = (1j*(exp(1j*k*x4)-H34*exp(1j*k*x3)))/(2*sin(k*(x3-x4)))
        H = (1j*(H34*exp(-1*1j*k*x3)-exp(-1*1j*k*x4)))/(2*sin(k*(x3-x4)))

        P0 = E + F
        V0 = (E-F)/(rho*c)
        Pd = G*exp(-1*1j*k*d)+H*exp(1j*k*d)
        Vd = (G*exp(-1*1j*k*d)-H*exp(1j*k*d))/(rho*c)

        T11 = ((H13*Pd*Vd)+(P0*V0/H13))/((P0*Vd)+(Pd*V0))
        T12 = ((P0*P0/H13)-(H13*Pd*Pd))/((P0*Vd)+(Pd*V0))
        T21 = ((V0*V0/H13)-(H13*Vd*Vd))/((P0*Vd)+(Pd*V0))
        T22 = T11

        Za = sqrt(T12/T21)

        Ta = (2*exp(1j*k*d))/(T11+(T12/(rho*c))+((rho*c)*T21)+T22)
        Ra = (T11+(T12/(rho*c))-((rho*c)*T21)-T22)/(T11+(T12/(rho*c))+((rho*c)*T21)+T22)

        # Ermittlung von Reflektions-, Transmissions-, Absorptionsgrad und Schalldämmaß sowie Test auf Leistungsgleichgewicht

        t_grad=abs(Ta)**2
        schall_dm=20*np.log10(1/t_grad)
        return schall_dm

    