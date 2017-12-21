from libcpp.map cimport map as libcpp_map
from String cimport *
from libcpp.string cimport string as libcpp_string
# from MzTabSmallMoleculeSectionData.string cimport string as libcpp_string

cdef extern from "<OpenMS/FORMAT/MzTab.h>" namespace "OpenMS":

    cdef cppclass MzTab:

        MzTab() nogil except +
        MzTab(MzTab) nogil except + # wrap-ignore

        # autowrap cannot handle wrapped classes as keys in std::map<> (this
        # would be std::map<String, ...> and does not work)
        # TODO STL map with wrapped key

        # libcpp_map[libcpp_string, MzTabUnitIdMetaData] getMetaData() const
        # MzTabProteinSectionData getProteinSectionData() const
        # MzTabPeptideSectionData getPeptideSectionData() const
        # MzTabSmallMoleculeSectionData getSmallMoleculeSectionData() const

