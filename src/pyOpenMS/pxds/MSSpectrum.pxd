from libcpp.vector cimport vector as libcpp_vector
from Peak1D cimport *
from String cimport *
from RangeManager cimport *
from DataArrays cimport *
from SpectrumSettings cimport *
from IMTypes cimport *
from DataValue cimport *
from InstrumentSettings cimport *
from SourceFile cimport *
from Precursor cimport *
from Product cimport *
from AcquisitionInfo cimport *
from DataProcessing cimport *
from MetaInfoInterface cimport *

# this class has addons, see the ./addons folder (../addons/MSSpectrum.pyx)

cdef extern from "<OpenMS/KERNEL/MSSpectrum.h>" namespace "OpenMS":

    cdef cppclass MSSpectrum(RangeManagerMzInt):
        # wrap-inherits:
        #  RangeManagerMzInt
        #  
        # wrap-doc:
        #  The representation of a 1D spectrum.
        #  Raw data access is proved by `get_peaks` and `set_peaks`, which yields numpy arrays
        #  Iterations yields access to underlying peak objects but is slower
        #  Extra data arrays can be accessed through getFloatDataArrays / getIntegerDataArrays / getStringDataArrays
        #  See help(SpectrumSettings) for information about meta-information
        #  
        #  Usage:
        #  
        #  .. code-block:: python
        #  
        #    ms_level = spectrum.getMSLevel()
        #    rt = spectrum.getRT()
        #    mz, intensities = spectrum.get_peaks()
        #  
        #  
        #  Usage:
        #  
        #  .. code-block:: python
        #  
        #    from pyopenms import *
        #  
        #    spectrum = MSSpectrum()
        #    spectrum.setDriftTime(25) # 25 ms
        #    spectrum.setRT(205.2) # 205.2 s
        #    spectrum.setMSLevel(3) # MS3
        #    p = Precursor()
        #    p.setIsolationWindowLowerOffset(1.5)
        #    p.setIsolationWindowUpperOffset(1.5)
        #    p.setMZ(600) # isolation at 600 +/- 1.5 Th
        #    p.setActivationEnergy(40) # 40 eV
        #    p.setCharge(4) # 4+ ion
        #    spectrum.setPrecursors( [p] )
        #  
        #    # Add raw data to spectrum
        #    spectrum.set_peaks( ([401.5], [900]) )
        #  
        #    # Additional data arrays / peak annotations
        #    fda = FloatDataArray()
        #    fda.setName("Signal to Noise Array")
        #    fda.push_back(15)
        #    sda = StringDataArray()
        #    sda.setName("Peak annotation")
        #    sda.push_back("y15++")
        #    spectrum.setFloatDataArrays( [fda] )
        #    spectrum.setStringDataArrays( [sda] )
        #  
        #    # Add spectrum to MSExperiment
        #    exp = MSExperiment()
        #    exp.addSpectrum(spectrum)
        #  
        #    # Add second spectrum and store as mzML file
        #    spectrum2 = MSSpectrum()
        #    spectrum2.set_peaks( ([1, 2], [1, 2]) )
        #    exp.addSpectrum(spectrum2)
        #  
        #    MzMLFile().store("testfile.mzML", exp)
        #  
        #  
        MSSpectrum() except + nogil  # wrap-doc:Constructor
        MSSpectrum(MSSpectrum &) except + nogil  # wrap-doc:Copy constructor

        double getRT() except + nogil  # wrap-doc:Returns the absolute retention time (in seconds)
        void setRT(double) except + nogil   # wrap-doc:Sets the absolute retention time (in seconds)

        double getDriftTime() except + nogil  # wrap-doc:Returns the drift time (-1 if not set)
        void setDriftTime(double) except + nogil  # wrap-doc:Sets the drift time (-1 if not set)
        DriftTimeUnit getDriftTimeUnit() except + nogil  # wrap-doc:Returns the ion mobility drift time unit
        String getDriftTimeUnitAsString() except + nogil  # wrap-doc:Returns the ion mobility drift time unit as string
        void setDriftTimeUnit(DriftTimeUnit dt) except + nogil  # wrap-doc:Sets the ion mobility drift time unit

        IMFormat getIMFormat() except + nogil  # wrap-doc:Returns the ion mobility format
        void setIMFormat(IMFormat im_format) except + nogil  # wrap-doc:Sets the ion mobility format

        bool containsIMData() except + nogil  # wrap-doc:Returns whether the spectrum contains ion mobility data
        libcpp_pair[Size, DriftTimeUnit] getIMData() except + nogil  # wrap-ignore wrap-doc:Returns position of ion mobility float data array and drift time unit

        unsigned int getMSLevel() except + nogil  # wrap-doc:Returns the MS level
        void setMSLevel(unsigned int) except + nogil  # wrap-doc:Sets the MS level

        String getName() except + nogil  # wrap-doc:Returns the name of the spectrum
        void setName(String) except + nogil  # wrap-doc:Sets the name of the spectrum

        Size size() except + nogil  # wrap-doc:Returns the number of peaks in the spectrum
        void reserve(size_t n) except + nogil  # wrap-doc:Reserves space for n peaks in the underlying container
        void resize(size_t n) except + nogil  # wrap-doc:Resize the peak array

        Peak1D& operator[](size_t) except + nogil  # wrap-upper-limit:size()

        void updateRanges() except + nogil  # wrap-doc:Recalculates the m/z and intensity ranges of the spectrum
        void clear(bool clear_meta_data) except + nogil  # wrap-doc:Clears all data (and meta data if clear_meta_data is True)
        void push_back(Peak1D)  except + nogil  # wrap-doc:Append a peak

        bool isSorted() except + nogil  # wrap-doc:Returns True if the spectrum is sorted by m/z

        int findNearest(double mz) except + nogil  # wrap-doc:Returns the index of the closest peak in m/z
        int findNearest(double mz, double tolerance) except + nogil  # wrap-doc:Returns the index of the closest peak in the provided +/- m/z tolerance window (-1 if none match)
        int findNearest(double mz, double tolerance_left, double tolerance_right) except + nogil  # wrap-doc:Returns the index of the closest peak in the provided abs. m/z tolerance window to the left and right (-1 if none match)
        int findHighestInWindow(double mz, double tolerance_left, double tolerance_right) except + nogil  # wrap-doc:Returns the index of the highest peak in the provided abs. m/z tolerance window to the left and right (-1 if none match)

        MSSpectrum select(const libcpp_vector[ size_t ] & indices) except + nogil  # wrap-doc:Subset the spectrum by indices. Also applies to associated data arrays if present.

        void assign(libcpp_vector[Peak1D].iterator, libcpp_vector[Peak1D].iterator) except + nogil  # wrap-ignore
        libcpp_vector[Peak1D].iterator begin() except + nogil   # wrap-iter-begin:__iter__(Peak1D)
        libcpp_vector[Peak1D].iterator end()   except + nogil   # wrap-iter-end:__iter__(Peak1D)

        double calculateTIC() except + nogil  # wrap-doc:Returns the total ion current (=sum) of peak intensities in the spectrum

        bool operator==(MSSpectrum) except + nogil  # wrap-doc:Equality operator
        bool operator!=(MSSpectrum) except + nogil  # wrap-doc:Inequality operator

        void sortByIntensity(bool reverse) except + nogil  # wrap-doc:Sorts the peaks by intensity (ascending if reverse is False, descending if True)
        void sortByPosition() except + nogil  # wrap-doc:Sorts the peaks by m/z position

        libcpp_vector[FloatDataArray] getFloatDataArrays() except + nogil  # wrap-doc:Returns the additional float data arrays to store e.g. meta data
        libcpp_vector[IntegerDataArray] getIntegerDataArrays() except + nogil  # wrap-doc:Returns the additional int data arrays to store e.g. meta data
        libcpp_vector[StringDataArray] getStringDataArrays() except + nogil  # wrap-doc:Returns the additional string data arrays to store e.g. meta data

        void setFloatDataArrays(libcpp_vector[FloatDataArray] fda) except + nogil  # wrap-doc:Sets the additional float data arrays to store e.g. meta data
        void setIntegerDataArrays(libcpp_vector[IntegerDataArray] ida) except + nogil  # wrap-doc:Sets the additional int data arrays to store e.g. meta data
        void setStringDataArrays(libcpp_vector[StringDataArray] sda) except + nogil  # wrap-doc:Sets the additional string data arrays to store e.g. meta data

        SpectrumType getType(bool query_data) except + nogil  # wrap-doc:Returns the spectrum type (centroided, profile or unknown). If SpectrumSettings and DataProcessing information are not sufficient and query_data is True, the data will be queried (potentially expensive)

        # SpectrumSettings forwarding methods (previously inherited)
        void unify(SpectrumSettings) except + nogil
        SpectrumType getType() except + nogil  # wrap-doc:Returns the spectrum type (centroided (PEAKS) or profile data (RAW))
        void setType(SpectrumType) except + nogil  # wrap-doc:Sets the spectrum type
        String getNativeID() except + nogil  # wrap-doc:Returns the native identifier for the spectrum, used by the acquisition software
        void setNativeID(String) except + nogil  # wrap-doc:Sets the native identifier for the spectrum, used by the acquisition software
        String getComment() except + nogil  # wrap-doc:Returns the free-text comment
        void setComment(String) except + nogil  # wrap-doc:Sets the free-text comment
        InstrumentSettings getInstrumentSettings() except + nogil  # wrap-doc:Returns the instrument settings of the current spectrum
        void setInstrumentSettings(InstrumentSettings) except + nogil  # wrap-doc:Sets the instrument settings of the current spectrum
        AcquisitionInfo getAcquisitionInfo() except + nogil  # wrap-doc:Returns the acquisition info
        void setAcquisitionInfo(AcquisitionInfo) except + nogil  # wrap-doc:Sets the acquisition info
        SourceFile getSourceFile() except + nogil  # wrap-doc:Returns the source file
        void setSourceFile(SourceFile) except + nogil  # wrap-doc:Sets the source file
        libcpp_vector[Precursor] getPrecursors() except + nogil  # wrap-doc:Returns the precursors
        void setPrecursors(libcpp_vector[Precursor]) except + nogil  # wrap-doc:Sets the precursors
        libcpp_vector[Product] getProducts() except + nogil  # wrap-doc:Returns the products
        void setProducts(libcpp_vector[Product]) except + nogil  # wrap-doc:Sets the products
        libcpp_vector[ shared_ptr[DataProcessing] ] getDataProcessing() except + nogil  # wrap-doc:Returns the description of the applied processing
        void setDataProcessing(libcpp_vector[ shared_ptr[DataProcessing] ]) except + nogil  # wrap-doc:Sets the description of the applied processing

        # MetaInfoInterface forwarding methods (previously inherited via SpectrumSettings)
        bool isMetaEmpty() except + nogil  # wrap-doc:Returns if the MetaInfo is empty
        void clearMetaInfo() except + nogil  # wrap-doc:Removes all meta values
        void getKeys(libcpp_vector[String] & keys) except + nogil  # wrap-doc:Fills the given vector with a list of all keys for which a value is set
        DataValue getMetaValue(String) except + nogil  # wrap-doc:Returns the value corresponding to a string, or DataValue::EMPTY if not found
        void setMetaValue(String, DataValue) except + nogil  # wrap-doc:Sets the DataValue corresponding to a name
        bool metaValueExists(String) except + nogil  # wrap-doc:Returns whether an entry with the given name exists
        void removeMetaValue(String) except + nogil  # wrap-doc:Removes the DataValue corresponding to `name` if it exists
