import numpy as np


def reflectance(n1, n2, inc_deg=45):
    theta_i = np.deg2rad(inc_deg)
    theta_t = np.sqrt(1 - (n1 / n2 * np.sin(theta_i)) ** 2)

    rs = (n1 * np.cos(theta_i) - n2 * theta_t) / (n1 * np.cos(theta_i) + n2 * theta_t)
    rp = (n1 * theta_t - n2 * np.cos(theta_i)) / (n1 * theta_t + n2 * np.cos(theta_i))

    refl_s = np.abs(rs) ** 2
    refl_p = np.abs(rp) ** 2

    return (refl_s + refl_p) / 2


def spectral_reflectance(mat_idx, env_idx, inc_deg=45):
    lambda_count = len(mat_idx)

    if lambda_count != len(env_idx):
        raise ValueError("Range of mat and env indices must be the same.")

    eps_mat = np.array(mat_idx)
    eps_env = np.array(env_idx)

    theta_i = np.deg2rad(inc_deg)
    theta_t = np.sqrt(1 - (np.divide(eps_env, eps_mat) * np.sin(theta_i)) ** 2)

    rs = np.divide(
        eps_env * np.cos(theta_i) - np.multiply(eps_mat, theta_t),
        eps_env * np.cos(theta_i) + np.multiply(eps_mat, theta_t)
    )
    rp = np.divide(
        np.multiply(eps_env, theta_t) - eps_mat * np.cos(theta_i),
        np.multiply(eps_env, theta_t) + eps_mat * np.cos(theta_i)
    )

    refl_s = np.abs(rs) ** 2
    refl_p = np.abs(rp) ** 2

    return (refl_s + refl_p) / 2
