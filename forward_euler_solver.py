import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys


sys.set_int_max_str_digits(100000)

def forward_euler(f, E_OV, E_OT, E_W, I_OV, I_OT, I_W, R_OV, R_OT, R_W, S0, E0_OV, E0_OT, E0_W, I0_OV, I0_OT, I0_W, t0, tn, h):
    ##define values for T S I and R
    t_values = [t0]
    S_values = [S0]

    EOV_values = [E0_OV]
    EOT_values = [E0_OT]
    EW_values = [E0_W]
    
    IOV_values = [I0_OV]
    IOT_values = [I0_OT]
    IW_values = [I0_W]

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


##S I and R functions
def S_dot(t, S, I_OV, I_OT, I_W):
    B = 0.151
    vR = ((0.8) * (0.074))
    p = 0.37
    value = -(B * S * p)* (I_OV +I_OT) - (B * S * I_W) - (vR * S)

    return value

def I_OV(t, E_OV, I_OV):
    alp = 0.25
    gamma = 0.93
    value = (alp * E_OV) - (gamma * I_OV) 

    return value

def I_OT(t, E_OT, I_OT, I_OV):
    alp = 0.25
    gamma = 0.93
    value = (alp * E_OT) - (gamma * (I_OV +I_OT)) 

    return value

def I_W(t, E_W, I_W):
    alp = 0.25
    gamma = 0.93
    value = (alp * E_W) - (gamma * I_W) 

    return value

def E_OV(t, S, E_OV):
    alp = 0.25
    vR = ((0.8) * (0.074))
    value = (vR * S) - (alp * E_OV) 

    return value

def E_OT(t, S, E_OT, I_OV, I_OT):
    alp = 0.25
    B = 0.151
    p = 0.37
    value = (p * B * S) * (I_OV + I_OT) - (alp * E_OT) 

    return value

def E_W(t, S, E_W, I_W):
    alp = 0.25
    B = 0.151
    value = (B * S * I_W) - (alp * E_W) 

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
    ##define int variables
    S0 = 0.99

    E0_OV = 0
    E0_OT = 0
    E0_W = .001

    I0_OV = 0
    I0_OT = 0
    I0_W = (1.2161) * (10 ** -4)

    t0 = 0
    tn = 295
    h = 1

    ##use function
    graph_values = forward_euler(S_dot, E_OV, E_OT, E_W, I_OV, I_OT, I_W, R_OV, R_OT, R_W, S0, E0_OV, E0_OT, E0_W, I0_OV, I0_OT, I0_W, t0, tn, h)


    
    t_values, S_values, EOV_values, EOT_values, EW_values, IOV_values, IOT_values, IW_values, ROV_values, ROT_values, RW_values = graph_values

    print(t_values)
    print(IW_values)


    out_file = "iw Plot"

    plt.plot(t_values, IW_values, color = 'red')

  

    plt.savefig(out_file,bbox_inches='tight')

    plt.show()

