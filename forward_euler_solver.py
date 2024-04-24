import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def forward_euler(f, I, y0, I0, t0, tn, h):
    ##define values for T S I and R
    t_values = [t0]
    y_values = [y0]
    I_values = [I0]
    R_values = [0]

    while t_values[-1] < tn:
        t_n = t_values[-1] + h ##incriment by h
        y_n = y_values[-1] + h * f(t_values[-1], y_values[-1], I_values[-1]) ##use eulers
        I_n = I_values[-1] + h * I(t_values[-1], y_values[-1], I_values[-1])
        t_values.append(t_n) ##store values
        y_values.append(y_n)
        I_values.append(I_n)

        ##print(f"t_n: {t_n}, y_n: {y_n}, I_n: {I_n}, R_n: {R_n}")

    return t_values, y_values, I_values


##S I and R functions
def S_dot(t, S, I):
    B = 3
    lam = 2

    value = -(B * S * I) + lam * I

    return value

def I_dot(t, S, I):
    B = 3
    lam = 2

    value = (B * S * I) - lam * I

    return value


##Make Analytical Function
def I(t, I):
    B = 3
    lam = 2
    R0 = B/lam
    R0_inv = 1/R0

    value = (1 - R0_inv) / (1 + (((1 - R0_inv - I) / I) * np.exp(-(B - lam) * t)))

    return value


def E(h):
    S0 = 0.99
    t0 = 0
    tn = 25
    I0 = 0.01

    forward_values = forward_euler(S_dot, I_dot, S0, I0, t0, tn, h)
    t_values, S_values, I_values = forward_values

    error_values = []

    for i in range(len(t_values)):
        abs_error = abs(I_values[i] - I(t_values[i],I0)) 
        error_values.append(abs_error)

    max_error = max(error_values)


    return max_error


if __name__ == "__main__":
    ##define int variables
    S0 = 0.99
    t0 = 0
    tn = 25
    h = 4
    I0 = 0.01

    ##use function
    graph_values = forward_euler(S_dot, I_dot, S0, I0, 0, tn, h)


    I_t = np.linspace(t0, tn, 100)

    h_values = [2, 1, 1/2 , 1/4, 1/8, 1/16, 1/32]
    error_values = []

    for i in range(len(h_values)):
        error_values.append(E(h_values[i]))

    ##store values
    t_values, S_values, I_values = graph_values


    out_file = "i Plot"
    out_file2 = "E Plot"

    plt.plot(t_values, I_values, color = 'red')
    plt.plot(I_t, I(I_t, I0), color = 'black', linestyle = '--')

    plt.ylim(0, 0.5)

    plt.xlabel('Time')
    plt.ylabel('I/N')
    plt.title('Forward Euler vs Analytical soln of i')
    plt.legend(['Forward Euler', 'Analytical'])

    plt.figure()

    plt.loglog(h_values, error_values, marker = 'o')
    plt.xlabel('Step sizes')
    plt.ylabel("Maximum Absolute Error")
    plt.title('Log-log of Max Absolute Error vs Step Size')

    plt.show()

    plt.savefig(out_file,bbox_inches='tight')

