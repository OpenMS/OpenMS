from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from Residue cimport *
from DefaultParamHandler cimport *
from AASequence cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/ProtonDistributionModel.h>" namespace "OpenMS":
    
    cdef cppclass ProtonDistributionModel(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        ProtonDistributionModel() nogil except +
        ProtonDistributionModel(ProtonDistributionModel) nogil except +
        void getProtonDistribution(libcpp_vector[ double ] & bb_charges, libcpp_vector[ double ] & sc_charges, AASequence & peptide, Int charge, ResidueType res_type) nogil except +
        void getChargeStateIntensities(AASequence & peptide, AASequence & n_term_ion, AASequence & c_term_ion, Int charge, ResidueType n_term_type, libcpp_vector[ double ] & n_term_intensities, libcpp_vector[ double ] & c_term_intensities, FragmentationType type_) nogil except +
        void setPeptideProtonDistribution(libcpp_vector[ double ] & bb_charge, libcpp_vector[ double ] & sc_charge) nogil except +

cdef extern from "<OpenMS/ANALYSIS/ID/ProtonDistributionModel.h>" namespace "OpenMS::ProtonDistributionModel":
    cdef enum FragmentationType "OpenMS::ProtonDistributionModel::FragmentationType":
        #wrap-attach:
        #    ProtonDistributionModel
        ChargeDirected
        ChargeRemote
        SideChain

