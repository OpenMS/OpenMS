from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from OpenSwathDataStructures cimport *
from ISpectrumAccess cimport *
from MSExperiment cimport *
from MzMLSqliteHandler cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessSqMass.h>" namespace "OpenMS":

  cdef cppclass SpectrumAccessSqMass(ISpectrumAccess):
        # wrap-inherits:
        #  ISpectrumAccess

        SpectrumAccessSqMass() # wrap-pass-constructor

        SpectrumAccessSqMass(SpectrumAccessSqMass &) except + nogil 
        SpectrumAccessSqMass(MzMLSqliteHandler, libcpp_vector[int] indices) except + nogil 

        # void getAllSpectra(libcpp_vector[ OpenSwath::SpectrumPtr ] & spectra, libcpp_vector< OpenSwath::SpectrumMeta > & spectra_meta) const;
        # SpectrumAccessSqMass(OpenMS::Internal::MzMLSqliteHandler handler);
        # SpectrumAccessSqMass(SpectrumAccessSqMass sp, std::vector<int> indices);
        # SpectrumAccessSqMass(OpenMS::Internal::MzMLSqliteHandler handler, std::vector<int> indices);
