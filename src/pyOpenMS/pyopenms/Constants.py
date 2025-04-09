
# @brief Mathematical and physical constants namespace.
# 
# This namespace contains definitions for some basic mathematical and physical constants.
# All constants are double precision.

EPSILON = 1e-6
"""EPSILON (used for comparisons)"""

PI = 3.14159265358979323846
"""PI"""

E = 2.718281828459045235
"""Euler's number - base of the natural logarithm"""

ELEMENTARY_CHARGE = 1.60217738E-19
"""Elementary charge (Coulomb)"""

e0 = ELEMENTARY_CHARGE

ELECTRON_MASS = 9.1093897E-31
"""Electron mass (kg)"""

ELECTRON_MASS_U = 1.0 / 1822.8885020477
"""Electron mass in units"""

PROTON_MASS = 1.6726230E-27
"""Proton mass (kg)"""

PROTON_MASS_U = 1.0072764667710
"""Proton mass in units"""

C13C12_MASSDIFF_U = 1.0033548378
"""Mass difference between Carbon-13 and Carbon-12 in units"""

NEUTRON_MASS = 1.6749286E-27
"""Neutron mass (kg)"""

NEUTRON_MASS_U = 1.00866491566
"""Neutron mass in units"""

AVOGADRO = 6.0221367E+23
"""Avogadro constant (1/mol)"""

# Avogadro constant (alias)
NA = AVOGADRO

# Avogadro constant (alias)
MOL = AVOGADRO

BOLTZMANN = 1.380657E-23
"""Boltzmann constant (J/K) """

# Boltzmann constant (alias)
k = BOLTZMANN

PLANCK = 6.6260754E-34
"""Planck constant J * sec"""

# Planck constant (alias)
h = PLANCK

GAS_CONSTANT = NA * k
"""Gas constant (= NA * k)"""

# Gas constant (alias)
R = GAS_CONSTANT

FARADAY = NA * e0
"""Faraday constant (= NA * e0)"""

# Faraday constant (alias)
F = FARADAY

BOHR_RADIUS = 5.29177249E-11 # m
"""Bohr radius (m)"""

# Bohr radius (alias)
a0 = BOHR_RADIUS

#  the following values from:
#  P.W.Atkins: Physical Chemistry, 5th ed., Oxford University Press, 1995

VACUUM_PERMITTIVITY = 8.85419E-12
"""Vacuum permittivity in C^2 / (J * m)"""

VACUUM_PERMEABILITY = (4 * PI * 1E-7)
"""Vacuum permeability in J s^2 / (C^2 * m)"""

SPEED_OF_LIGHT = 2.99792458E+8
"""Speed of light in m/s"""

# Speed of Light (alias)
c = SPEED_OF_LIGHT

GRAVITATIONAL_CONSTANT  = 6.67259E-11
""" Gravitational constant in N m^2 / kg^2"""

FINE_STRUCTURE_CONSTANT = 7.29735E-3
"""Fine structure constant"""

DEG_PER_RAD = 57.2957795130823209
"""Degree per rad"""

RAD_PER_DEG = 0.0174532925199432957
"""Rad per degree"""

MM_PER_INCH = 25.4
"""mm per inch"""

M_PER_FOOT = 3.048
"""m per foot"""

JOULE_PER_CAL  = 4.184
"""Joule per calorie"""

CAL_PER_JOULE = (1 / 4.184)
"""Calories per Joule"""

PRECURSOR_ERROR_PPM_USERPARAM = "precursor_mz_error_ppm"
"""User parameter name for precursor mz error in ppm"""

FRAGMENT_ANNOTATION_USERPARAM = "fragment_annotation"
"""User parameter name for precursor mz error in ppm"""

