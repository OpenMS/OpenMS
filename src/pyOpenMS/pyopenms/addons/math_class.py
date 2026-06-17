"""
Pure Python Math namespace wrapper for pyOpenMS compatibility.

OpenMS::Math is a C++ namespace, not a class. This provides
pyopenms.Math.getPPM() etc. as static methods on a Python class.
"""


class Math:
    """Math namespace providing mathematical utility functions for mass spectrometry."""

    @staticmethod
    def getPPM(mz_obs: float, mz_ref: float) -> float:
        """Compute parts-per-million difference between two m/z values."""
        return (mz_obs - mz_ref) / mz_ref * 1e6

    @staticmethod
    def getPPMAbs(mz_obs: float, mz_ref: float) -> float:
        """Compute absolute PPM difference."""
        return abs((mz_obs - mz_ref) / mz_ref * 1e6)

    @staticmethod
    def ppmToMass(ppm: float, mz_ref: float) -> float:
        """Convert PPM to mass difference in Da."""
        return ppm / 1e6 * mz_ref

    @staticmethod
    def ppmToMassAbs(ppm: float, mz_ref: float) -> float:
        """Convert PPM to absolute mass difference in Da."""
        return abs(ppm / 1e6 * mz_ref)
