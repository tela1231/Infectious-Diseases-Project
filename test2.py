import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import math



# Define the differential OPV_equations
def OPV_equations(t, y, beta, rho, sigma, gamma):
    S, IOV, IOT, EOV, EOT, ROV, ROT= y
    
    # Define phi based on the time t
    if 144 <= t <= 175 or 206 <= t <= 216:
        phi = 0.8 * 0.074
    else:
        phi = 0
    
    S_dot = -(beta * S * rho) * (IOV + IOT)  - (phi * S)
    IOV_dot = sigma * EOV - gamma * IOV
    IOT_dot = sigma * EOT - gamma * (IOT)
    EOV_dot = (phi * S) - (sigma * EOV)
    EOT_dot = (rho * beta * S) * (IOV + IOT) - (sigma * EOT)
    ROV_dot = gamma * IOV
    ROT_dot = gamma * IOT
    
    return [S_dot, IOV_dot, IOT_dot, EOV_dot, EOT_dot, ROV_dot, ROT_dot]

def WPV_equations(t, y, beta, rho, sigma, gamma):
    S, IW, EW, RW = y
    
    S_dot = - (beta * S * IW) 
    
    IW_dot = sigma * EW - gamma * IW
   
    EW_dot = (beta * S * IW) - (sigma * EW)
    
    RW_dot = gamma * IW
    
    return [S_dot, IW_dot, EW_dot, RW_dot]



# Define initial conditions
IOV_0 = 0
IOT_0 = 0
ROV_0 = 0
ROT_0 = 0
RW_0 = 0
alpha_star = 0.25
gamma_star = 0.093
IW_0 = 0.00434
IW_00 = 0.001
EW_0 = (gamma_star / alpha_star) * IW_00
S0 = 1 - IW_00 - EW_0
EOV_0 = 0
EOT_0 = 0

# Initial conditions

y0_OPV = [S0, IOV_0, IOT_0, EOV_0, EOT_0, ROV_0, ROT_0]

y0_WPV = [S0, IW_0, EW_0, RW_0]

# Define parameter values
beta = 0.151
rho = 0.37
sigma = 0.25
gamma = 0.093
sewage_scaling_o = 2.02 * (10**7)
sewage_scaling_w = 1.07 * (10**7)

# Define the time span
t_span = [0, 295]

# Solve the differential OPV_equations
OPV_solution = solve_ivp(OPV_equations, t_span, y0_OPV, args=(beta, rho, sigma, gamma), dense_output=True)

# Solve the differential WPV_equations
WPV_solution = solve_ivp(WPV_equations, t_span, y0_WPV, args=(beta, rho, sigma, gamma), dense_output=True)

# Extract the OPV_solution
t_values = np.linspace(0, 295, 100)
y_OPV_values = OPV_solution.sol(t_values)

# Extract the WPV_solution
y_WPV_values = WPV_solution.sol(t_values)


# Extract the OPV_solution for S, I, E, and R
S_values = y_OPV_values[0]
IOV_values = y_OPV_values[1]
IOT_values = y_OPV_values[2]
EOV_values = y_OPV_values[3]
EOT_values = y_OPV_values[4]
ROV_values = y_OPV_values[5]
ROT_values = y_OPV_values[6]

# Extract the WPV_solution for S, I, E, and R
S_values = y_WPV_values[0]
IW_values = y_WPV_values[1]
EW_values = y_WPV_values[2]
RW_values = y_WPV_values[3]


# Plotting and saving each variable separately

# Plot S
plt.figure(figsize=(10, 6))
plt.plot(t_values, S_values, color='blue')
plt.xlabel('Time')
plt.ylabel('Population')
plt.title('Susceptible Population')
plt.grid(True)
plt.savefig('S_plot.png', bbox_inches='tight')
plt.close()

# Plot IW, IOT, and IOV
plt.figure(figsize=(10, 6))
plt.plot(t_values, IW_values, color='green', label='IW')
plt.plot(t_values, IOT_values, color='orange', label='IOT')
plt.plot(t_values, IOV_values, color='blue', label='IOV')
plt.xlabel('Time')
plt.ylabel('Population')
plt.title('Infected Population')
plt.legend()
plt.grid(True)
plt.savefig('Infected_plot.png', bbox_inches='tight')
plt.close()

# Plot EW, EOT, and EOV
plt.figure(figsize=(10, 6))
plt.plot(t_values, EW_values, color='green', label='EW')
plt.plot(t_values, EOT_values, color='orange', label='EOT')
plt.plot(t_values, EOV_values, color='blue', label='EOV')
plt.xlabel('Time')
plt.ylabel('Population')
plt.title('Exposed Population')
plt.legend()
plt.grid(True)
plt.savefig('Exposed_plot.png', bbox_inches='tight')
plt.close()

# Plot RW, ROT, and ROV
plt.figure(figsize=(10, 6))
plt.plot(t_values, RW_values, color='green', label='RW')
plt.plot(t_values, ROT_values, color='orange', label='ROT')
plt.plot(t_values, ROV_values, color='blue', label='ROV')
plt.xlabel('Time')
plt.ylabel('Population')
plt.title('Recovered Population')
plt.legend()
plt.grid(True)
plt.savefig('Recovered_plot.png', bbox_inches='tight')
plt.close()

