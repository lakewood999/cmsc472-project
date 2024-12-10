import numpy as np

def update_covariance_with_vector(cov_matrix, delta_v_vector):
    """
    Updates a 6x6 covariance matrix after a prograde or retrograde delta-V burn.
    Using a simplified model, the velocity uncertainty is increased along the
    velocity direction, with small adjustments to the position-velocity cross terms.

    Uses a process noise model based on the magnitude of the delta-V burn.

    Based on logic from paper: https://ntrs.nasa.gov/api/citations/20070018023/downloads/20070018023.pdf 

    Parameters:
    - cov_matrix (np.ndarray): Initial 6x6 covariance matrix (position and velocity).
    - delta_v_vector (np.ndarray): 3D velocity vector of the delta-V burn (e.g., [dx, dy, dz]).

    Returns:
    - np.ndarray: Updated 6x6 covariance matrix.
    """
    # conservative scaling factor for velocity uncertainty
    # since we want to yield final obs, which tend to be more certain since it's measured close to
    # TCA, use a small value of like 5% (random guess but idk)
    velocity_uncertainty_scale = 0.05  

    delta_v_magnitude = np.linalg.norm(delta_v_vector)
    
    if delta_v_magnitude > 0:
        delta_v_unit = delta_v_vector / delta_v_magnitude
    else:
        delta_v_unit = np.zeros_like(delta_v_vector)  

    # velocity covariance - noise aligned with the delta V direction
    Q_velocity = velocity_uncertainty_scale**2 * delta_v_magnitude**2 * np.outer(delta_v_unit, delta_v_unit)

    P_velocity = cov_matrix[3:, 3:]  # bottom right of covariance (3x3)
    P_velocity_updated = P_velocity + Q_velocity    # add process noise


    P_cross = cov_matrix[:3, 3:]  
    coupling_factor = 0.1  # assume some small conservative coupling factor
    P_cross_updated = P_cross + coupling_factor * Q_velocity  

    # update position covariance (minimal change)
    P_position = cov_matrix[:3, :3]  
    P_position_updated = P_position + coupling_factor * Q_velocity  
    
    # construct updated covariance matrix
    updated_cov_matrix = np.zeros((6, 6))
    updated_cov_matrix[:3, :3] = P_position_updated
    updated_cov_matrix[3:, 3:] = P_velocity_updated
    updated_cov_matrix[:3, 3:] = P_cross_updated
    updated_cov_matrix[3:, :3] = P_cross_updated.T  # Symmetry

    return updated_cov_matrix


if "__name__" == "__main__":

    # no cross-covariance between position and velocity
    initial_cov_matrix = np.eye(6)  
    delta_v_example_vector = np.array([0.1, 0.2, 0.3])  # m/s

    updated_cov_matrix_vector = update_covariance_with_vector(initial_cov_matrix, delta_v_example_vector)
    print(updated_cov_matrix_vector)
