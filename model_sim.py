import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import math
import matplotlib.dates as mdates
from datetime import datetime, timedelta



# Define the differential equations
def equations(t, y, beta, rho, sigma, gamma):
    S, IOV, IOT, IW, EOV, EOT, EW, ROV, ROT, RW = y

    # Define phi based on the time t
    if 144 <= t <= 175 or 206 <= t <= 216:
        phi = (0.8) * (0.074)
    else:
        phi = 0

    S_dot = -(beta * S * rho) * (IOV + IOT) - (beta * S * IW) - (phi * S)
    IOV_dot = sigma * EOV - gamma * IOV
    IOT_dot = sigma * EOT - gamma * (IOT)
    IW_dot = sigma * EW - gamma * IW
    EOV_dot = (phi * S) - (sigma * EOV)
    EOT_dot = (rho * beta * S) * (IOV + IOT) - (sigma * EOT)
    EW_dot = (beta * S * IW) - (sigma * EW)
    ROV_dot = gamma * IOV
    ROT_dot = gamma * IOT
    RW_dot = gamma * IW

    return [S_dot, IOV_dot, IOT_dot, IW_dot, EOV_dot, EOT_dot, EW_dot, ROV_dot, ROT_dot, RW_dot]

def calculate_y(variable, I_values):
    y_values = []
    for I in I_values:
        if I != 0:  # Check if I is not zero
            y_values.append(math.log2(variable / I))
        else:
            # Handle the case where I is zero (for example, set y to 60)
            y_values.append(60)
    return y_values



def find_closest_index(numbers, target):
    closest_index = None
    min_difference = float('inf')  # Initialize with positive infinity

    for i, num in enumerate(numbers):
        difference = abs(num - target)
        if difference < min_difference:
            min_difference = difference
            closest_index = i

    return closest_index


# Define initial conditions
IOV_0 = 0
IOT_0 = 0
EOV_0 = 0
EOT_0 = 0
ROV_0 = 0
ROT_0 = 0
RW_0 = 0
IW_0 = 0.00434
EW_0 = 0 #set temp value to EW

# Define the time span
t_span = [0, 295]

# define star initial conditions for EW estimation
alpha_star = 0.25
gamma_star = 0.093
IW_0_star = 0.001
EW_0_star = (gamma_star / alpha_star) * IW_0_star
S0_star = 1 - IW_0_star - EW_0_star

# Define parameter values
beta = 0.151
rho = 0.37
sigma = 0.25
gamma = 0.093
sewage_scaling_o = 2.02 * (10**7)
sewage_scaling_w = 1.07 * (10**7)

# Star initial conditions
y0_star = [S0_star, IOV_0, IOT_0, IW_0_star, EOV_0, EOT_0, EW_0_star, ROV_0, ROT_0, RW_0]  

# Solve the differential equations for EW estimation
star_solutions = solve_ivp(equations, t_span, y0_star, args=(beta, rho, sigma, gamma), dense_output=True)

# Extract the star solutions
t_values = np.linspace(0, 295, 100)
y_star_values = star_solutions.sol(t_values)

# Extract the IW_star and EW_star soln
IW_star_values = y_star_values[3]
EW_star_values = y_star_values[6]

# Use find_closest_index to find where IW_star is equal to IW_0
index = find_closest_index(IW_star_values, IW_0)

# Assign the value of EW_0 to the index  from EW_star_values
EW_0 = EW_star_values[index]

# Determine S0 with new EW_0 value
S0 = 1 - IW_0 - EW_0

# Set new initial conditions for main modeling
y0 = [S0, IOV_0, IOT_0, IW_0, EOV_0, EOT_0, EW_0, ROV_0, ROT_0, RW_0]  

# Solve the differential equations
solution = solve_ivp(equations, t_span, y0, args=(beta, rho, sigma, gamma), dense_output=True)

# Extract the solutions
y_values = solution.sol(t_values)

# Extract the star_solutions for S, I, E, and R
S_values = y_values[0]
IW_values = y_values[3]
IOV_values = y_values[1]
IOT_values = y_values[2]
EOV_values = y_values[4]
EOT_values = y_values[5]
EW_values = y_values[6]
ROV_values = y_values[7]
ROT_values = y_values[8]
RW_values = y_values[9]
Io_values = IOT_values + IOV_values

# Add function for cycle thresholds
y_w = calculate_y(sewage_scaling_w, IW_values)
y_o = calculate_y(sewage_scaling_o, Io_values)

# Define start and end dates
start_date = datetime(2013, 3, 11)
end_date = datetime(2013, 12, 31)

# Convert numerical values of t_values to dates
num_days = 295
date_range = [start_date + timedelta(days=int(t)) for t in np.linspace(0, num_days, len(t_values))]


# Plotting and saving each variable separately

# Plot S
plt.figure(figsize=(10, 6))
plt.plot(date_range, S_values, color='blue')
plt.xlabel('Date')
plt.ylabel('Population')
plt.title('Susceptible Population')
plt.grid(True)

# Set date format on x-axis
date_format = mdates.DateFormatter('%Y-%m')
plt.gca().xaxis.set_major_formatter(date_format)
plt.gca().xaxis.set_major_locator(mdates.MonthLocator())

plt.savefig('S_plot.png', bbox_inches='tight')
plt.close()

# Plot IW, IOT, and IOV
plt.figure(figsize=(10, 6))
plt.plot(date_range, IW_values, color='green', label='IW')
plt.plot(date_range, IOT_values, color='orange', label='IOT')
plt.plot(date_range, IOV_values, color='blue', label='IOV')
plt.xlabel('Date')
plt.ylabel('Population')
plt.title('Infected Population')
plt.legend()
plt.grid(True)

# Set date format on x-axis
date_format = mdates.DateFormatter('%Y-%m')
plt.gca().xaxis.set_major_formatter(date_format)
plt.gca().xaxis.set_major_locator(mdates.MonthLocator())

plt.savefig('Infected_plot.png', bbox_inches='tight')
plt.close()

# Plot EW, EOT, and EOV
plt.figure(figsize=(10, 6))
plt.plot(date_range, EW_values, color='green', label='EW')
plt.plot(date_range, EOT_values, color='orange', label='EOT')
plt.plot(date_range, EOV_values, color='blue', label='EOV')
plt.xlabel('Date')
plt.ylabel('Population')
plt.title('Exposed Population')
plt.legend()
plt.grid(True)

# Set date format on x-axis
date_format = mdates.DateFormatter('%Y-%m')
plt.gca().xaxis.set_major_formatter(date_format)
plt.gca().xaxis.set_major_locator(mdates.MonthLocator())

plt.savefig('Exposed_plot.png', bbox_inches='tight')
plt.close()

# Plot RW, ROT, and ROV
plt.figure(figsize=(10, 6))
plt.plot(date_range, RW_values, color='green', label='RW')
plt.plot(date_range, ROT_values, color='orange', label='ROT')
plt.plot(date_range, ROV_values, color='blue', label='ROV')
plt.xlabel('Date')
plt.ylabel('Population')
plt.title('Recovered Population')
plt.legend()
plt.grid(True)

# Set date format on x-axis
date_format = mdates.DateFormatter('%Y-%m')
plt.gca().xaxis.set_major_formatter(date_format)
plt.gca().xaxis.set_major_locator(mdates.MonthLocator())

plt.savefig('Recovered_plot.png', bbox_inches='tight')
plt.close()

# Plot y_w
plt.figure(figsize=(10, 6))
plt.plot(date_range, y_w, color='blue')
plt.xlabel('Date')
plt.ylabel('PCR Cycle Threshold')
plt.title('WPV1 PCR Cycle Threshold vs Date')
plt.grid(True)
plt.gca().invert_yaxis()  # Invert the y-axis
plt.ylim(60, 20)  # Set y-axis limits from 20 to 60

# Set date format on x-axis
date_format = mdates.DateFormatter('%Y-%m')
plt.gca().xaxis.set_major_formatter(date_format)
plt.gca().xaxis.set_major_locator(mdates.MonthLocator())

plt.savefig('yw_plot.png', bbox_inches='tight')
plt.close()

# Plot y_o
plt.figure(figsize=(10, 6))
plt.plot(date_range, y_o, color='blue')
plt.xlabel('Date')
plt.ylabel('PCR Cycle Threshold')
plt.title('OPV1 PCR Cycle Threshold vs Date')
plt.grid(True)
plt.gca().invert_yaxis()  # Invert the y-axis
plt.ylim(60, 20)  # Set y-axis limits from 20 to 60

# Set date format on x-axis
date_format = mdates.DateFormatter('%Y-%m')
plt.gca().xaxis.set_major_formatter(date_format)
plt.gca().xaxis.set_major_locator(mdates.MonthLocator())

plt.savefig('yo_plot.png', bbox_inches='tight')
plt.close()