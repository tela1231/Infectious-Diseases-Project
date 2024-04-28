import matplotlib
matplotlib.use('Agg')
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np
import sys


sys.set_int_max_str_digits(100000)

def forward_euler(f, E_OV, E_OT, E_W, I_OV, I_OT, I_W, R_OV, R_OT, R_W, S0, EOV_0, EOT_0, EW_0, IOV_0, IOT_0, IW_0, t0, tn, h):
    ##define values for T S I and R
    t_values = [t0]
    S_values = [S0]

    EOV_values = [EOV_0]
    EOT_values = [EOT_0]
    EW_values = [EW_0]
    
    IOV_values = [IOV_0]
    IOT_values = [IOT_0]
    IW_values = [IW_0]

    ROV_values = [0]
    ROT_values = [0]
    RW_values = [0]

    while t_values[-1] < tn:
        t_n = t_values[-1] + h ##incriment by h

        S_n = S_values[-1] + h * f(t_values[-1], S_values[-1], IOV_values[-1], IOT_values[-1], IW_values[-1]) ##use eulers

        EOV_n = EOV_values[-1] + h * E_OV(t_values[-1], S_values[-1], EOV_values[-1])
        EOT_n = EOT_values[-1] + h * E_OT(t_values[-1], S_values[-1], EOT_values[-1], IOV_values[-1], IOT_values[-1])
        EW_n = EW_values[-1] + h * E_W(t_values[-1], S_values[-1], EW_values[-1], IW_values[-1])

        IOV_n = IOV_values[-1] + h * I_OV(t_values[-1], EOV_values[-1], IOV_values[-1])
        IOT_n = IOT_values[-1] + h * I_OT(t_values[-1], EOT_values[-1], IOT_values[-1], IOV_values[-1])
        IW_n = IW_values[-1] + h * I_W(t_values[-1], EW_values[-1], IW_values[-1])

        ROV_n = ROV_values[-1] + h * R_OV(t_values[-1], IOV_values[-1])
        ROT_n = ROT_values[-1] + h * R_OT(t_values[-1], IOT_values[-1])
        RW_n = RW_values[-1] + h * R_W(t_values[-1], IW_values[-1])


        t_values.append(t_n) ##store values

        S_values.append(S_n)

        IOV_values.append(IOV_n)
        IOT_values.append(IOT_n)
        IW_values.append(IW_n)

        EOV_values.append(EOV_n)
        EOT_values.append(EOT_n)
        EW_values.append(EW_n)

        ROV_values.append(ROV_n)
        ROT_values.append(ROT_n)
        RW_values.append(RW_n)

        ##print(f"t_n: {t_n}, S_n: {S_n}, I_n: {I_n}, R_n: {R_n}")

    return t_values, S_values, EOV_values, EOT_values, EW_values, IOV_values, IOT_values, IW_values, ROV_values, ROT_values, RW_values

def EW0_forward_euler(f, E_OV, E_OT, E_W, I_OV, I_OT, I_W, R_OV, R_OT, R_W, S0, EOV_0, EOT_0, EW_0, IOV_0, IOT_0, IW_0, t0, tn, h):
 ##define values for T S I and R
    t_values = [t0]
    S_values = [S0]

    EOV_values = [EOV_0]
    EOT_values = [EOT_0]
    EW_values = [EW_0]
    
    IOV_values = [IOV_0]
    IOT_values = [IOT_0]
    IW_values = [IW_0]

    ROV_values = [0]
    ROT_values = [0]
    RW_values = [0]

    while t_values[-1] < tn:

        if IW_values[-1] == 0.00434:
            value = EW_values[-1]
            print("found")
            print(value)
            return value

        t_n = t_values[-1] + h ##incriment by h

        S_n = S_values[-1] + h * f(t_values[-1], S_values[-1], IOV_values[-1], IOT_values[-1], IW_values[-1]) ##use eulers

        EOV_n = EOV_values[-1] + h * E_OV(t_values[-1], S_values[-1], EOV_values[-1])
        EOT_n = EOT_values[-1] + h * E_OT(t_values[-1], S_values[-1], EOT_values[-1], IOV_values[-1], IOT_values[-1])
        EW_n = EW_values[-1] + h * E_W(t_values[-1], S_values[-1], EW_values[-1], IW_values[-1])

        IOV_n = IOV_values[-1] + h * I_OV(t_values[-1], EOV_values[-1], IOV_values[-1])
        IOT_n = IOT_values[-1] + h * I_OT(t_values[-1], EOT_values[-1], IOT_values[-1], IOV_values[-1])
        IW_n = IW_values[-1] + h * I_W(t_values[-1], EW_values[-1], IW_values[-1])

        ROV_n = ROV_values[-1] + h * R_OV(t_values[-1], IOV_values[-1])
        ROT_n = ROT_values[-1] + h * R_OT(t_values[-1], IOT_values[-1])
        RW_n = RW_values[-1] + h * R_W(t_values[-1], IW_values[-1])


        t_values.append(t_n) ##store values

        S_values.append(S_n)

        IOV_values.append(IOV_n)
        IOT_values.append(IOT_n)
        IW_values.append(IW_n)

        EOV_values.append(EOV_n)
        EOT_values.append(EOT_n)
        EW_values.append(EW_n)

        ROV_values.append(ROV_n)
        ROT_values.append(ROT_n)
        RW_values.append(RW_n)

        print(IW_n)

        ##print(f"t_n: {t_n}, S_n: {S_n}, I_n: {I_n}, R_n: {R_n}")

    return t_values, S_values, EOV_values, EOT_values, EW_values, IOV_values, IOT_values, IW_values, ROV_values, ROT_values, RW_values
    


##S I and R functions
def S_dot(t, S, I_OV, I_OT, I_W):
    B = 0.151
    vR = ((0.8) * (0.074))
    p = 0.37
    value = -(B * S * p)* (I_OV +I_OT) - (B * S * I_W) - (vR * S)

    return value

def I_OV(t, E_OV, I_OV):
    alpha = 0.25
    gamma = 0.93
    value = (alpha * E_OV) - (gamma * I_OV) 

    return value

def I_OT(t, E_OT, I_OT, I_OV):
    alpha = 0.25
    gamma = 0.93
    value = (alpha * E_OT) - (gamma * (I_OV +I_OT)) 

    return value

def I_W(t, E_W, I_W):
    alpha = 0.25
    gamma = 0.93
    value = (alpha * E_W) - (gamma * I_W) 

    return value

def E_OV(t, S, E_OV):
    alpha = 0.25
    vR = ((0.8) * (0.074))
    value = (vR * S) - (alpha * E_OV) 

    return value

def E_OT(t, S, E_OT, I_OV, I_OT):
    alpha = 0.25
    B = 0.151
    p = 0.37
    value = (p * B * S) * (I_OV + I_OT) - (alpha * E_OT) 

    return value

def E_W(t, S, E_W, I_W):
    alpha = 0.25
    B = 0.151
    value = (B * S * I_W) - (alpha * E_W) 

    return value

def R_OV(t, I_OV):
    gamma = 0.93
    value = gamma * I_OV

    return value

def R_OT(t, I_OT):
    gamma = 0.93
    value = gamma * I_OT
    return value

def R_W(t, I_W):
    gamma = 0.93
    value = gamma * I_W
    return value

if __name__ == "__main__":

##define zero initial variables
    EOV_0 = 0
    EOT_0 = 0

    IOV_0 = 0
    IOT_0 = 0

##time and step size variables
    t0 = 0
    tn = 295
    h = .001

##define variables for finding initial EW0 value
    alpha = 0.25
    gamma = 0.93

    IW0_Star = 0.001
    EW0_Star = (gamma / alpha) * IW0_Star
    S0_Star = 1 - IW0_Star - EW0_Star

    ##EW_0 = EW0_forward_euler(S_dot, E_OV, E_OT, E_W, I_OV, I_OT, I_W, R_OV, R_OT, R_W, S0_Star, EOV_0, EOT_0, EW0_Star, IOV_0, IOT_0, IW0_Star, t0, tn, h)
    graph_values = forward_euler(S_dot, E_OV, E_OT, E_W, I_OV, I_OT, I_W, R_OV, R_OT, R_W, S0_Star, EOV_0, EOT_0, EW0_Star, IOV_0, IOT_0, IW0_Star, t0, tn, h)
    t_values, S_values, EOV_values, EOT_values, EW_values, IOV_values, IOT_values, IW_values, ROV_values, ROT_values, RW_values = graph_values
    
    target_value = 0.00434
    values = IW_values
    values_array = np.array(values)

   # Sort the array
    sorted_indices = np.argsort(np.abs(values_array - target_value))
    sorted_values = values_array[sorted_indices]

    sorted_indices = np.argsort(np.abs(values_array - target_value))
    sorted_values = values_array[sorted_indices]

    # Lower tolerance level
    tolerance = 1e-6

    # Find the index of the closest value to target_value within the specified tolerance
    index = sorted_indices[np.argmin(np.abs(sorted_values - target_value) < tolerance)]

    # Print the index
    print("Index of the closest value to 0.00434 with lower tolerance:", index)

   
"""
##define variables for main model
    IW_0 = (4.34) * (10 ** -3)

    S0 = 1 - IW_0 - EW_0

##use function
    graph_values = forward_euler(S_dot, E_OV, E_OT, E_W, I_OV, I_OT, I_W, R_OV, R_OT, R_W, S0, EOV_0, EOT_0, EW_0, IOV_0, IOT_0, IW_0, t0, tn, h)


    
    t_values, S_values, EOV_values, EOT_values, EW_values, IOV_values, IOT_values, IW_values, ROV_values, ROT_values, RW_values = graph_values



    out_file = "iw Plot"

    plt.plot(t_values, IW_values, color = 'red')

  

    plt.savefig(out_file,bbox_inches='tight')

    plt.show()
"""
