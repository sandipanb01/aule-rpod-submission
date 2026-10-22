"""Attitude: STK .a I/O, quaternion conventions, interpolation, law identification."""
from .stk_io import AttitudeProfile, read_stk_attitude  # noqa: F401
from .providers import make_attitude_provider  # noqa: F401
