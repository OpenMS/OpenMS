from Types cimport *
from DefaultParamHandler cimport *
from Matrix cimport *
from String cimport *
from Peak2D cimport *
from IsobaricChannelInformation cimport *

# typedef std::vector<IsobaricChannelInformation> IsobaricChannelList;
cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>" namespace "OpenMS":

    cdef cppclass IsobaricQuantitationMethod(DefaultParamHandler) :
        # wrap-ignore
        # ABSTRACT class
        # no-pxd-import
        # wrap-doc:
        #  Abstract base class describing an isobaric quantitation method in terms of
        #  the reporter ion channels used and an isotope correction matrix. Isobaric
        #  labeling methods like TMT and iTRAQ use reporter ions for multiplexed
        #  quantitation of peptides/proteins across multiple samples

        # wrap-inherits:
        #  DefaultParamHandler
        IsobaricQuantitationMethod() except + nogil
        IsobaricQuantitationMethod(IsobaricQuantitationMethod &) except + nogil
        String getName() except + nogil  # wrap-doc:Returns the unique name or identifier of the quantitation method
        libcpp_vector[IsobaricChannelInformation]  getChannelInformation() except + nogil  # wrap-doc:Returns information on the different channels used by this quantitation method
        Size getNumberOfChannels() except + nogil  # wrap-doc:Returns the number of channels available for this quantitation method
        Matrix[ double ] getIsotopeCorrectionMatrix() except + nogil  # wrap-doc:Returns the isotope correction matrix for correcting reporter ion intensities
        Size getReferenceChannel() except + nogil  # wrap-doc:Returns the index of the reference channel used for ratio calculation
