from traits.api import HasPrivateTraits, Property, Float, Trait, Delegate, Int, List, Array


class Tube(HasPrivateTraits):
    # Tube parameters
    tube_shape = Trait('rect', 'circ',
                       desc="Shape of the measurement tube")
    tube_d = Float(
        desc="diameter or (if rectangular) largest section dimension of tube")
        
    


class Tube_Impedance(Tube):
    '''
    Impedance tube with 2 Microphone positions with the specimen at the end of the tube 
    following ISO 10534-2 (for reflection/impedance measurements)
    The parameters are defined in the README #TODO
    '''
    s = Float(desc='distance between the microphones')
    x1 = Float(
        desc='distance between specimen and microphone 1 (furthest from the specimen)')


class Tube_Transmission(Tube):
    '''
    Tube with 4 Microphone positions for transmission measurements following the E2611 Norm
    The parameters are defined in the README
    '''
    l1 = Float(desc='distance between beginning of specimen and mic 2')
    l2 = Float(desc='distance between beginning of specimen and mic 3')
    s1 = Float(desc='Distance between mic 1 and 2 in m')
    s2 = Float(desc='Distance between mic 3 and 4 in m')
    d = Float(desc='length of test specimen (test tube section is 0.7m)')


# Presets for tubes at TAP, TU Berlin:

class Tube_TAP_Transmission_rect_narrow(Tube_Transmission):
    '''
    Tube Preset for the rectangular transmission tube at TAP, TU Berlin
    With microphone positions 1, 2, 3, 4 for higher frequencies

    '''
    tube_shape = 'rect'
    tube_d = 0.1
    l1 = 0.3   # distance between beginning of specimen and mic 2
    l2 = 0.8   # distance between beginning of specimen and mic 3
    s1 = 0.085  # Distance between mic 1 and 2
    s2 = 0.085  # Distance between mic 3 and 4
    d = 0.5   # length of test specimen (test tube section is 0.7m)


class Tube_TAP_Transmission_rect_wide(Tube_Transmission):
    '''
    Preset for the rectangular transmission tube at TAP, TU Berlin
    With microphone positions 0, 2, 3, 5 for lower frequencies
    '''
    tube_shape = 'rect'
    tube_d = 0.1
    l1 = 0.3   # distance between beginning of specimen and mic 2
    l2 = 0.8   # distance between beginning of specimen and mic 3
    s1 = 0.5   # Distance between mic 1 and 2
    s2 = 0.5   # Distance between mic 3 and 4
    d = 0.5   # length of test specimen (test tube section is 0.7m)


class Tube_TAP_Impedance_rect_25cm_narrow(Tube_Impedance):
    '''
    Preset for the rectangular impedance tube with side length 25cm at TAP, TU Berlin
    With microphone positions 1 and 2
    '''
    tube_shape = 'rect'
    tube_d = 0.25
    s = 0.225
    x1 = 2.75


class Tube_TAP_Impedance_rect_25cm_wide(Tube_Impedance):
    '''
    Preset for the rectangular impedance tube with side length 25cm at TAP, TU Berlin
    With microphone positions 1 and 3
    '''
    tube_shape = 'rect'
    tube_d = 0.25
    s = 0.85
    x1 = 2.75


class Tube_TAP_Impedance_rect_50cm_narrow(Tube_Impedance):
    '''
    Preset for the rectangular impedance tube with side length 25cm at TAP, TU Berlin
    With microphone positions 1 and 2
    '''
    tube_shape = 'rect'
    tube_d = 0.5
    s = 0.45
    x1 = 2.75


class Tube_TAP_Impedance_rect_50cm_wide(Tube_Impedance):
    '''
    Preset for the rectangular impedance tube with side length 25cm at TAP, TU Berlin
    With microphone positions 1 and 3
    '''
    tube_shape = 'rect'
    tube_d = 0.5
    s = 0.85
    x1 = 2.75
