"""Quaternion conventions -- stated once, tested, used everywhere.

Storage: [q1, q2, q3, q4] = [e sin(theta/2), cos(theta/2)], scalar LAST (STK).

Meaning (STK): q is the rotation FROM the reference axes (LVLH) TO the body
axes, i.e. the body triad is the reference triad rotated by theta about e.

Attitude matrix (reference -> body component transform; Shuster 1993,
Markley & Crassidis 2014 Eq. 2.125):

    A(q) = (q4^2 - |q|^2) I + 2 q q^T - 2 q4 [q x]          v_B = A(q) v_R

so a body-fixed vector (the docking-port normal) has reference components

    p_R = A(q)^T p_B = R_active(q) p_B

which is exactly ``scipy.spatial.transform.Rotation.from_quat(q).apply(p_B)``
(scipy: scalar-last, active rotation). Tested in tests/test_quaternion.py with
the canonical case q = (0, 0, sin45, cos45): body X must point along +Y_ref.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial.transform import Rotation


def attitude_matrix(q):
    """Shuster/STK attitude matrix A(q) (reference->body), q scalar-last."""
    q = np.asarray(q, dtype=float)
    v, s = q[..., :3], q[..., 3]
    I = np.eye(3)
    qq = np.einsum("...i,...j->...ij", v, v)
    cross = np.zeros(q.shape[:-1] + (3, 3))
    cross[..., 0, 1], cross[..., 0, 2] = -v[..., 2], v[..., 1]
    cross[..., 1, 0], cross[..., 1, 2] = v[..., 2], -v[..., 0]
    cross[..., 2, 0], cross[..., 2, 1] = -v[..., 1], v[..., 0]
    return (s ** 2 - np.sum(v * v, -1))[..., None, None] * I + 2 * qq - 2 * s[..., None, None] * cross


def body_to_reference(q, v_body):
    """Reference-frame components of a body-fixed vector (convention above)."""
    return Rotation.from_quat(q).apply(v_body)


def angle_between_deg(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    c = np.sum(a * b, -1) / (np.linalg.norm(a, axis=-1) * np.linalg.norm(b, axis=-1))
    return np.degrees(np.arccos(np.clip(c, -1.0, 1.0)))
