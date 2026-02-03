from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from String cimport *
from Types cimport *
from ProteinIdentification cimport *

cdef extern from "<OpenMS/METADATA/IdentifierMSRunMapper.h>" namespace "OpenMS":

    cdef cppclass IdentifierMSRunMapper:
        # wrap-doc:
        #  Maps protein identification identifiers to MS run file paths.
        #
        #  This class encapsulates the mapping between ProteinIdentification identifiers
        #  and their associated MS run paths. It is useful for resolving the correct
        #  source file for peptide identifications, especially in merged identification results.
        #
        #  Example usage:
        #
        #  .. code-block:: python
        #
        #     # Create mapping from protein identifications
        #     mapper = oms.IdentifierMSRunMapper(protein_ids)
        #
        #     # Get MS run paths for a peptide's identifier
        #     paths = mapper.getMSRunPaths(pep_id.getIdentifier())
        #
        #     # Build a USI using the mapping
        #     usi = pep_id.buildUSI(mapper, "PXD000561", False)

        IdentifierMSRunMapper() except + nogil
        IdentifierMSRunMapper(IdentifierMSRunMapper &) except + nogil
        IdentifierMSRunMapper(const libcpp_vector[ProteinIdentification] & prot_ids) except + nogil
            # wrap-doc:
            #  Construct mapper from a list of ProteinIdentifications.
            #
            #  :param prot_ids: List of ProteinIdentification objects

        void create(const libcpp_vector[ProteinIdentification] & prot_ids) except + nogil
            # wrap-doc:
            #  Create/update mapping from a list of ProteinIdentifications.
            #
            #  :param prot_ids: List of ProteinIdentification objects

        bool hasIdentifier(const String & identifier) except + nogil
            # wrap-doc:
            #  Check if the mapping contains an entry for the given identifier.
            #
            #  :param identifier: ProteinIdentification identifier
            #  :return: True if identifier exists in mapping

        bool empty() except + nogil
            # wrap-doc:
            #  Check if the mapping is empty.
            #
            #  :return: True if no mappings exist

        size_t size() except + nogil
            # wrap-doc:
            #  Get the number of identifier mappings.
            #
            #  :return: Number of identifiers in the mapping

        StringList getMSRunPaths(const String & identifier) except + nogil
            # wrap-doc:
            #  Get the MS run paths for a given identifier.
            #
            #  :param identifier: ProteinIdentification identifier
            #  :return: List of MS run file paths (empty if identifier not found)
