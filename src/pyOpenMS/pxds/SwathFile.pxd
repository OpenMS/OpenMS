from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from MSExperiment cimport *
from MzMLFile cimport *
from MzXMLFile cimport *
from SwathMap cimport *
from StringList cimport *
from ExperimentalSettings cimport *

cdef extern from "<OpenMS/FORMAT/SwathFile.h>" namespace "OpenMS":
    
    cdef cppclass SwathFile(ProgressLogger) :
        # wrap-inherits:
        #  ProgressLogger
        SwathFile() nogil except + # compiler
        SwathFile(SwathFile &) nogil except + # compiler

        libcpp_vector[ SwathMap ] loadSplit(StringList file_list,
                                            String tmp,
                                            shared_ptr[ ExperimentalSettings ] exp_meta,
                                            String readoptions) nogil except + # wrap-doc:Loads a Swath run from a list of split mzML files

        libcpp_vector[ SwathMap ] loadMzML(String file_,
                                           String tmp,
                                           shared_ptr[ ExperimentalSettings ] exp_meta,
                                           String readoptions) nogil except +
        # wrap-doc:
                #   Loads a Swath run from a single mzML file
                #   -----
                #   Using the `plugin_consumer`, you can provide a custom consumer which will be chained
                #   into the process of loading the data and making it available (depending on `readoptions`).
                #   This is useful if you want to modify the data a priori or extract some other information using
                #   MSDataTransformingConsumer (for example). Make sure it leaves the data intact, such that the 
                #   returned SwathMaps are actually useful
                #   -----
                #   :param file: Input filename
                #   :param tmp: Temporary directory (for cached data)
                #   :param exp_meta: Experimental metadata from mzML file
                #   :param readoptions: How are spectra accessed after reading - tradeoff between memory usage and time (disk caching)
                #   :param plugin_consumer: An intermediate custom consumer
                #   :returns: Swath maps for MS2 and MS1 (unless readoptions == split, which returns no data)

        libcpp_vector[ SwathMap ] loadMzXML(String file_,
                                            String tmp,
                                            shared_ptr[ ExperimentalSettings ] exp_meta, 
                                            String readoptions) nogil except + # wrap-doc:Loads a Swath run from a single mzXML file
