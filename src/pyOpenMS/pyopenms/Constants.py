
# @brief Mathematical and physical constants namespace.
# 
# This namespace contains definitions for some basic mathematical and physical constants.
# All constants are double precision.

# EPSILON (used for comparisons)
EPSILON = 1e-6

# PI
PI = 3.14159265358979323846

# Euler's number - base of the natural logarithm
E = 2.718281828459045235

# Elementary charge
ELEMENTARY_CHARGE = 1.60217738E-19 # C

# Elementary charge (alias)
e0 = ELEMENTARY_CHARGE

# Electron mass
ELECTRON_MASS = 9.1093897E-31 # kg

# Electron mass in units
ELECTRON_MASS_U = 1.0 / 1822.8885020477 # u

# Proton mass
PROTON_MASS = 1.6726230E-27 # kg

# Proton mass in units
PROTON_MASS_U = 1.0072764667710 # u

# Mass difference between Carbon-13 and Carbon-12 in units
C13C12_MASSDIFF_U = 1.0033548378 # u

# Neutron mass
NEUTRON_MASS = 1.6749286E-27 # kg

# Neutron mass in units
NEUTRON_MASS_U = 1.00866491566 # u

# Avogadro constant
AVOGADRO = 6.0221367E+23 # 1 / mol

# Avogadro constant (alias)
NA = AVOGADRO

# Avogadro constant (alias)
MOL = AVOGADRO

# Boltzmann constant
BOLTZMANN = 1.380657E-23 # J / K

# Boltzmann constant (alias)
k = BOLTZMANN

# Planck constant
PLANCK = 6.6260754E-34 # J * sec

# Planck constant (alias)
h = PLANCK

# Gas constant (= NA * k)
GAS_CONSTANT = NA * k

# Gas constant (alias)
R = GAS_CONSTANT

# Faraday constant (= NA * e0)
FARADAY = NA * e0

# Faraday constant (alias)
F = FARADAY

# Bohr radius
BOHR_RADIUS = 5.29177249E-11 # m

# Bohr radius (alias)
a0 = BOHR_RADIUS

#  the following values from:
#  P.W.Atkins: Physical Chemistry, 5th ed., Oxford University Press, 1995

# Vacuum permittivity
VACUUM_PERMITTIVITY = 8.85419E-12 # C^2 / (J * m)

# Vacuum permeability
VACUUM_PERMEABILITY = (4 * PI * 1E-7) # J s^2 / (C^2 * m)

# Speed of light
SPEED_OF_LIGHT = 2.99792458E+8 # m / s

# Speed of Light (alias)
c = SPEED_OF_LIGHT

# Gravitational constant
GRAVITATIONAL_CONSTANT  = 6.67259E-11 # N m^2 / kg^2

# Fine structure constant
FINE_STRUCTURE_CONSTANT = 7.29735E-3 # 1

# Degree per rad
DEG_PER_RAD = 57.2957795130823209

# Rad per degree
RAD_PER_DEG = 0.0174532925199432957

# mm per inch
MM_PER_INCH = 25.4

# m per foot
M_PER_FOOT = 3.048

# Joule per calorie
JOULE_PER_CAL  = 4.184

# Calories per Joule
CAL_PER_JOULE = (1 / 4.184)

# User parameter name for precursor mz error in ppm
PRECURSOR_ERROR_PPM_USERPARAM = "precursor_mz_error_ppm"

# User parameter name for precursor mz error in ppm
FRAGMENT_ANNOTATION_USERPARAM = "fragment_annotation"

