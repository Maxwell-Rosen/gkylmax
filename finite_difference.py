
def central_difference(theta_lo, theta_hi, Theta_lo, Theta_hi):
    """
    Calculate the central difference derivative.
    
    Parameters:
    theta_lo (float): Lower bound of theta
    theta_hi (float): Upper bound of theta
    Theta_lo (float): Corresponding Theta at theta_lo
    Theta_hi (float): Corresponding Theta at theta_hi
    
    Returns:
    float: The central difference derivative
    """
    h = (theta_hi - theta_lo) / 2
    dTheta_dtheta = (Theta_hi - Theta_lo) / (theta_hi - theta_lo)
    return dTheta_dtheta

# Given values
theta_curr_lo = -0.0213334
theta_curr_hi = -0.0209109
Theta_lo = -0.450867
Theta_hi = -0.450725

# Calculate the central difference derivative
dTheta_dtheta_calculated = central_difference(theta_curr_lo, theta_curr_hi, Theta_lo, Theta_hi)

# Print the result
print(f"Calculated dTheta/dtheta: {dTheta_dtheta_calculated}")
print(f"Given dTheta/dtheta: 0.334")
