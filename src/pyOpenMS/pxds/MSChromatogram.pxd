from Types cimport *
from ChromatogramSettings cimport *
from ChromatogramPeak cimport *
from String cimport *
from RangeManager cimport *
from MetaInfoDescription cimport *
from DataArrays cimport *
from DataValue cimport *
from InstrumentSettings cimport *
from SourceFile cimport *
from Precursor cimport *
from Product cimport *
from AcquisitionInfo cimport *
from DataProcessing cimport *
from MetaInfoInterface cimport *

# this class has addons, see the ./addons folder

cdef extern from "<OpenMS/KERNEL/MSChromatogram.h>" namespace "OpenMS":

    cdef cppclass MSChromatogram (RangeManagerRtInt):
        # wrap-inherits:
        #  RangeManagerRtInt
        #  
        # wrap-doc:
        #  The representation of a chromatogram.
        #  Raw data access is proved by `get_peaks` and `set_peaks`, which yields numpy arrays
        #  Iterations yields access to underlying peak objects but is slower
        #  Extra data arrays can be accessed through getFloatDataArrays / getIntegerDataArrays / getStringDataArrays
        #  See help(ChromatogramSettings) for information about meta-information
        #  
        #  Usage:
        #  
        #  .. code-block:: python
        #  
        #    precursor = chromatogram.getPrecursor()
        #    product = chromatogram.getProduct()
        #    rt, intensities = chromatogram.get_peaks()
        #  

        MSChromatogram() except + nogil  # wrap-doc:Constructor
        MSChromatogram(MSChromatogram &) except + nogil  # wrap-doc:Copy constructor
        double getMZ() except + nogil  # wrap-doc:Returns the mz of the product entry, makes sense especially for MRM scans
        # void   setMZ(double) except + nogil

        String getName() except + nogil  # wrap-doc:Returns the name
        void setName(String) except + nogil  # wrap-doc:Sets the name

        Size size() except + nogil  # wrap-doc:Returns the number of peaks in the chromatogram
        void reserve(size_t n) except + nogil  # wrap-doc:Reserves space for n peaks in the underlying container
        void resize(size_t n) except + nogil  # wrap-doc:Resize the peak array

        ChromatogramPeak & operator[](size_t) except + nogil  # wrap-upper-limit:size()

        void updateRanges() except + nogil  # wrap-doc:Recalculates the RT and intensity ranges of the chromatogram
        void clear(int) except + nogil
            # wrap-doc:
                #  Clears all data and meta data
                #  
                #  
                #  :param clear_meta_data: If true, all meta data is cleared in addition to the data

        void push_back(ChromatogramPeak)  except + nogil  # wrap-doc:Append a peak

        bool isSorted() except + nogil  # wrap-doc:Checks if all peaks are sorted with respect to ascending RT

        void sortByIntensity(bool reverse) except + nogil
            # wrap-doc:
                #  Lexicographically sorts the peaks by their intensity
                #  
                #  
                #  Sorts the peaks according to ascending intensity. Meta data arrays will be sorted accordingly

        void sortByPosition() except + nogil
            # wrap-doc:
                #  Lexicographically sorts the peaks by their position
                #  
                #  
                #  The chromatogram is sorted with respect to position. Meta data arrays will be sorted accordingly

        int findNearest(double) except + nogil
            # wrap-doc:
                #  Binary search for the peak nearest to a specific RT
                #  :note: Make sure the chromatogram is sorted with respect to RT! Otherwise the result is undefined
                #  
                #  
                #  :param rt: The searched for mass-to-charge ratio searched
                #  :return: Returns the index of the peak.
                #  :raises:
                #    Exception: Precondition is thrown if the chromatogram is empty (not only in debug mode)

        void assign(libcpp_vector[ChromatogramPeak].iterator, libcpp_vector[ChromatogramPeak].iterator) except + nogil  # wrap-ignore
        libcpp_vector[ChromatogramPeak].iterator begin() except + nogil   # wrap-iter-begin:__iter__(ChromatogramPeak)
        libcpp_vector[ChromatogramPeak].iterator end()   except + nogil   # wrap-iter-end:__iter__(ChromatogramPeak)

        libcpp_vector[FloatDataArray] getFloatDataArrays() except + nogil  # wrap-doc:Returns a reference to the float meta data arrays
        libcpp_vector[IntegerDataArray] getIntegerDataArrays() except + nogil  # wrap-doc:Returns a reference to the integer meta data arrays
        libcpp_vector[StringDataArray] getStringDataArrays() except + nogil  # wrap-doc:Returns a reference to the string meta data arrays

        void setFloatDataArrays(libcpp_vector[FloatDataArray] fda) except + nogil  # wrap-doc:Sets the float meta data arrays
        void setIntegerDataArrays(libcpp_vector[IntegerDataArray] ida) except + nogil  # wrap-doc:Sets the integer meta data arrays
        void setStringDataArrays(libcpp_vector[StringDataArray] sda) except + nogil  # wrap-doc:Sets the string meta data arrays

        # ChromatogramSettings forwarding methods (previously inherited)
        Product getProduct() except + nogil  # wrap-doc:Returns the product ion
        void setProduct(Product p) except + nogil  # wrap-doc:Sets the product ion
        String getNativeID() except + nogil  # wrap-doc:Returns the native identifier
        void setNativeID(String native_id) except + nogil  # wrap-doc:Sets the native identifier
        String getComment() except + nogil  # wrap-doc:Returns the free-text comment
        void setComment(String comment) except + nogil  # wrap-doc:Sets the free-text comment
        InstrumentSettings getInstrumentSettings() except + nogil  # wrap-doc:Returns the instrument settings
        void setInstrumentSettings(InstrumentSettings instrument_settings) except + nogil  # wrap-doc:Sets the instrument settings
        AcquisitionInfo getAcquisitionInfo() except + nogil  # wrap-doc:Returns the acquisition info
        void setAcquisitionInfo(AcquisitionInfo acquisition_info) except + nogil  # wrap-doc:Sets the acquisition info
        SourceFile getSourceFile() except + nogil  # wrap-doc:Returns the source file
        void setSourceFile(SourceFile source_file) except + nogil  # wrap-doc:Sets the source file
        Precursor getPrecursor() except + nogil  # wrap-doc:Returns the precursor
        void setPrecursor(Precursor precursor) except + nogil  # wrap-doc:Sets the precursor
        libcpp_vector[ shared_ptr[DataProcessing] ] getDataProcessing() except + nogil  # wrap-doc:Returns the description of the applied processing
        void setDataProcessing(libcpp_vector[ shared_ptr[DataProcessing] ]) except + nogil  # wrap-doc:Sets the description of the applied processing
        void setChromatogramType(ChromatogramType type) except + nogil  # wrap-doc:Sets the chromatogram type
        ChromatogramType getChromatogramType() except + nogil  # wrap-doc:Get the chromatogram type

        # MetaInfoInterface forwarding methods (previously inherited via ChromatogramSettings)
        bool isMetaEmpty() except + nogil  # wrap-doc:Returns if the MetaInfo is empty
        void clearMetaInfo() except + nogil  # wrap-doc:Removes all meta values
        void getKeys(libcpp_vector[String] & keys) except + nogil  # wrap-doc:Fills the given vector with a list of all keys for which a value is set
        DataValue getMetaValue(String) except + nogil  # wrap-doc:Returns the value corresponding to a string, or DataValue::EMPTY if not found
        void setMetaValue(String, DataValue) except + nogil  # wrap-doc:Sets the DataValue corresponding to a name
        bool metaValueExists(String) except + nogil  # wrap-doc:Returns whether an entry with the given name exists
        void removeMetaValue(String) except + nogil  # wrap-doc:Removes the DataValue corresponding to `name` if it exists
