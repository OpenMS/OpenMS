from libcpp.map cimport map as libcpp_map
from String cimport *
from libcpp.string cimport string as libcpp_string
# from MzTabSmallMoleculeSectionData.string cimport string as libcpp_string

cdef extern from "<OpenMS/FORMAT/MzTab.h>" namespace "OpenMS":

    cdef cppclass MzTab:
        # wrap-doc:
                #   Data model of MzTab files
                #   -----
                #   Please see the official MzTab specification at https://code.google.com/p/mztab/

        MzTab() nogil except +
        MzTab(MzTab &) nogil except + # compiler

        # autowrap cannot handle wrapped classes as keys in std::map<> (this
        # would be std::map<String, ...> and does not work)
        # TODO STL map with wrapped key

        # libcpp_map[libcpp_string, MzTabUnitIdMetaData] getMetaData() const
        # MzTabProteinSectionData getProteinSectionData() const
        # MzTabPeptideSectionData getPeptideSectionData() const
        # MzTabSmallMoleculeSectionData getSmallMoleculeSectionData() const

