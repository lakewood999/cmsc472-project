import math

def calculate_delta_r(a, delta_v, mu=3.986004418e14):
    """
    Calculate the difference in orbital radii after a Delta V maneuver.

    Makes the assumption that the initial and post-maneuver orbits are circular.
    Since the target object will typically be a menuvering spacecraft, the assumption
    is that its orbit is relatively circular. The post-maneuever spacecraft is also assumed to be 
    circular because the ESA dataset doesn't provide true anomaly, and so we can't truly know how the orbit 
    changes. We make an average approximation as the average distances between the two objects, 
    given an eccentic orbit, can be approximated by a circular orbit of same potential energy.

    Parameters:
    - a (float): Initial semi-major axis (meters)
    - delta_v (float): Delta V applied (meters per second)
    - mu (float): Gravitational parameter (m^3/s^2). Default is Earth's mu.

    Returns:
    - delta_r (float): Difference in orbital radii (meters)
    - a_new (float): New orbital radius after maneuver (meters)
    """
    v = math.sqrt(mu / a)
    
    v_new = v + delta_v
    
    # use conservation of energy
    E_new = 0.5 * v_new**2 - mu / a
    
    # New semi-major axis from total energy
    # E = -mu / (2a)
    a_new = -mu / (2 * E_new)
    
    delta_r = a_new - a
    
    return delta_r, a_new

if __name__ == "__main__":
    print("Delta V Maneuver Impact on Circular Orbit")
    print("-----------------------------------------")

    # LEO SMA
    a = 2000 * 1e3
    
    # meters per second
    delta_v = 500
    
    
    
    #  delta_r
    delta_r, a_new = calculate_delta_r(a, delta_v)
    
    a_initial_km = a / 1e3
    a_new_km = a_new / 1e3
    delta_r_km = delta_r / 1e3
    
    print("\nResults:")
    print(f"Initial orbital radius (a): {a_initial_km:,.2f} km")
    print(f"New orbital radius (a'): {a_new_km:,.2f} km")
    print(f"Difference in radii (Δr): {delta_r_km:,.2f} km")
