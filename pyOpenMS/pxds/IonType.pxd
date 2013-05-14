from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from TheoreticalSpectrumGenerator cimport *
from SimTypes cimport *
from SVMWrapper cimport *
from Residue cimport *
from EmpiricalFormula cimport *

cdef extern from "<OpenMS/CHEMISTRY/SvmTheoreticalSpectrumGenerator.h>" namespace "OpenMS::SvmTheoreticalSpectrumGenerator":
    
    cdef cppclass IonType "OpenMS::SvmTheoreticalSpectrumGenerator::IonType":
        IonType() nogil except +
        IonType(IonType) nogil except +
        IonType(ResidueType residue, EmpiricalFormula l, Int charge) nogil except +
        bool operator<(IonType &rhs) nogil except +
        ResidueType residue
        EmpiricalFormula loss
        Int charge

