import numpy as np
from scipy.stats import ncx2

from maneuver_covariance import update_covariance_with_vector
from maneuver_orbit import calculate_delta_r

def compute_collision_probability(relative_position, cov_rel, collision_radius):
    """
    Compute the collision probability based on relative position and covariance.

    Parameters:
    - relative_position (np.ndarray): 3D relative position vector (meters)
    - cov_rel (np.ndarray): 3x3 relative covariance matrix (meters^2)
    - collision_radius (float): Collision radius (meters)

    Returns:
    - float: Collision probability
    """
    d = relative_position
    P_rel = cov_rel
    try:
        inv_P_rel = np.linalg.inv(P_rel)
    except np.linalg.LinAlgError:
        return 1.0  # Maximum risk if covariance matrix is singular
    lambda_param = d.T @ inv_P_rel @ d
    PR = ncx2.cdf(collision_radius**2, df=3, nc=lambda_param)
    return PR

def compute_risk(original_risk, delta_v_vector, initial_cov_matrix, initial_r):
    """
    Compute the updated risk after a maneuver.

    Parameters:
    - original_risk (float): Original collision risk
    - delta_v_vector (np.ndarray): 3D velocity vector of the delta-V burn (e.g., [dx, dy, dz])
    - initial_cov_matrix (np.ndarray): Initial 6x6 covariance matrix (position and velocity)
    - initial_r (float): Initial orbital radius (meters)

    Returns:
    - float: Updated collision risk
    """

    # Update covariance matrix
    updated_cov_matrix = update_covariance_with_vector(initial_cov_matrix, delta_v_vector)

    delta_r, a_new = calculate_delta_r(initial_r, np.linalg.norm(delta_v_vector))

    relative_position = np.array([delta_r, 0, 0])
    cov_rel = updated_cov_matrix[:3, :3]

    updated_risk = compute_collision_probability(relative_position, cov_rel, 10)

    return updated_risk * original_risk