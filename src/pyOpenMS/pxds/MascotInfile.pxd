from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from MSExperiment cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *
from String cimport *
from File cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/FORMAT/MascotInfile.h>" namespace "OpenMS":
    
    cdef cppclass MascotInfile(ProgressLogger) :
        # wrap-inherits:
        #  ProgressLogger
        MascotInfile() nogil except +
        MascotInfile(MascotInfile &) nogil except +
        void store(const String & filename, MSSpectrum & spec, double mz, double retention_time, String search_title) nogil except + # wrap-doc:Stores the peak list in a MascotInfile that can be used as input for MASCOT shell execution
        void store(const String & filename, MSExperiment & experiment, String search_title) nogil except + # wrap-doc:Stores the experiment data in a MascotInfile that can be used as input for MASCOT shell execution

        void load(const String & filename, MSExperiment & exp) nogil except +
            # wrap-doc:
                #   Loads a Mascot Generic File into a PeakMap
                #   -----
                #   :param filename: File name which the map should be read from
                #   :param exp: The map which is filled with the data from the given file
                #   :raises:
                #     Exception: FileNotFound is thrown if the given file could not be found

        String  getBoundary() nogil except + # wrap-doc:Returns the boundary used for the MIME format
        void setBoundary(const String & boundary) nogil except + # wrap-doc:Sets the boundary used for the MIME format.By default a 22 character random string is used
        String  getDB() nogil except + # wrap-doc:Returns the DB to use
        void setDB(const String & db) nogil except + # wrap-doc:Sets the DB (default MSDB). See mascot path /config/mascot.dat in "Databases" section for possible settings
        String  getSearchType() nogil except + # wrap-doc:Returns the search type
        void setSearchType(const String & search_type) nogil except + # wrap-doc:Sets the search type (default MIS). So far only MIS is supported!Valid types are "MIS" (MS/MS Ion Search), "PMF" (Peptide Mass Fingerprint) , "SQ" (Sequence Query)
        String  getHits() nogil except + # wrap-doc:Returns the number of hits to report back
        void setHits(const String & hits) nogil except + # wrap-doc:Sets the number of hits to report back (default 20)
        String  getCleavage() nogil except + # wrap-doc:Returns the enzyme used for cleavage
        void setCleavage(const String & cleavage) nogil except + # wrap-doc:Sets the enzyme used for cleavage (default Trypsin). See mascot path /config/enzymes for possible settings
        String  getMassType() nogil except + # wrap-doc:Returns the used mass type ("Monoisotopic" or "Average")
        void setMassType(const String & mass_type) nogil except + # wrap-doc:Sets the used mass type "Monoisotopic" or "Average" (default Monoisotopic)
        libcpp_vector[ String ]  getModifications() nogil except + # wrap-doc:Returns a vector containing the fixed modifications (default none)
        void setModifications(libcpp_vector[ String ] & mods) nogil except + # wrap-doc:Sets the fixed modifications (default none). See mascot path /config/mod_file for possible settings
        libcpp_vector[ String ]  getVariableModifications() nogil except + # wrap-doc:Returns a vector containing the variable modifications (default none)
        void setVariableModifications(libcpp_vector[ String ] & mods) nogil except + # wrap-doc:Sets the fixed modifications (default none). See mascot path /config/mod_file for possible settings
        String  getInstrument() nogil except + # wrap-doc:Returns the instrument type
        void setInstrument(const String & instrument) nogil except + # wrap-doc:Sets the instrument type (Default Default). Possible instruments are ESI-QUAD-TOF, MALDI-TOF-PSD, ESI-TRAP, ESI-QUAD, ESI-FTICR, MALDI-TOF-TOF, ESI-4SECTOR, FTMS-ECD, MALDI-QUAD-TOF, MALDI-QIT-TOF
        UInt getMissedCleavages() nogil except + # wrap-doc:Returns the number of allowed missed cleavages
        void setMissedCleavages(UInt missed_cleavages) nogil except + # wrap-doc:Sets the number of allowed missed cleavages (default 1)
        float getPrecursorMassTolerance() nogil except + # wrap-doc:Returns the precursor mass tolerance
        void setPrecursorMassTolerance(float precursor_mass_tolerance) nogil except + # wrap-doc:Sets the precursor mass tolerance in Da (default 2.0)
        float getPeakMassTolerance() nogil except + # wrap-doc:Returns the peak mass tolerance in Da
        void setPeakMassTolerance(float ion_mass_tolerance) nogil except + # wrap-doc:Sets the peak mass tolerance in Da (default 1.0)
        String  getTaxonomy() nogil except + # wrap-doc:Returns the taxonomy
        void setTaxonomy(const String & taxonomy) nogil except + # wrap-doc:Sets the taxonomy (default All entries). See mascot path /config/taxonomy for possible settings
        String  getFormVersion() nogil except + # wrap-doc:Returns the Mascot form version
        void setFormVersion(const String & form_version) nogil except + # wrap-doc:Sets the Mascot form version (default 1.01)
        String  getCharges() nogil except + # wrap-doc:Returns the charges
        void setCharges(libcpp_vector[ int ] & charges) nogil except + # wrap-doc:Sets the charges (default 1+, 2+ and 3+)
