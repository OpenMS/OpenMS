from Types cimport *
from libcpp cimport bool
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector
from libcpp.set cimport set as libcpp_set
#from libcpp.map cimport map as libcpp_map
from MetaInfoInterface cimport *
from MetaData cimport *
from ScoreType cimport *
from ProcessingSoftware cimport *
from ProcessingStep cimport *
from DBSearchParam cimport *
from InputFile cimport *
from Observation cimport *
from ParentSequence cimport *
from IdentifiedSequence cimport *
from IdentifiedCompound cimport *
from Adduct cimport *
from AdductInfo cimport *
from ObservationMatch cimport *
from ObservationMatchGroup cimport *
from ParentMatch cimport *
from ParentGroup cimport *



cdef extern from "<OpenMS/METADATA/ID/IdentificationData.h>" namespace "OpenMS":
    
  cdef cppclass IdentificationData(MetaInfoInterface):
      # wrap-inherits:
      #    MetaInfoInterface

      IdentificationData() nogil except + # wrap-doc:Representation of spectrum identification results and associated data

      IdentificationData(IdentificationData&) nogil except + # wrap-doc:Copy constructor

      #IdentificationData(IdentificationData&&) nogil except + #wrap-doc:Move constructor

      InputFileRef registerInputFile(InputFile& fie) nogil except +

      ProcessingSoftwareRef registerProcessingSoftware(ProcessingSoftware& software) nogil except +

      SearchParamRef registerDBSearchParam(DBSearchParam& param) nogil except +

      ProcessingStepRef registerProcessingStep(ProcessingStep& step) nogil except +

      ProcessingStepRef registerProcessingStep(ProcessingStep& step, SearchParamRef search_ref) nogil except +

      ScoreTypeRef registerScoreType(ScoreType & score) nogil except +

      ObservationRef registerObservation(const Observation& obs) nogil except +

      ParentSequenceRef registerParentSequence(const ParentSequence& parent) nogil except +

      void registerParentGroupSet(const ParentGroupSet) nogil except +

      #IdentifiedPeptideRef registerIdentifiedPeptide(const IdentifiedPeptide& peptide) nogil except +

      IdentifiedCompoundRef registerIdentifiedCompound(const IdentifiedCompound& compound) nogil except +

      #IdentifiedOligoRef registerIdentifiedOligo(const IdentifiedOligo& oligo) nogil except +

      AdductRef registerAdduct(const AMSE_AdductInfo& adduct) nogil except +

      ObservationMatchRef registerObservationMatch(const ObservationMatch& match) nogil except +

      MatchGroupRef registerObservationMatchGroup(const ObservationMatchGroup& group)

      InputFiles& getInputFiles() nogil except +

      libcpp_set[ProcessingSoftware] getProcessingSoftwares() nogil except +

      ProcessingSteps& getProcessingSteps() nogil except +

      DBSearchParams& getDBSearchParams() nogil except +

      #DBSearchSteps& getDBSearchSteps() nogil except +

      libcpp_set[ScoreType] getScoreTypes() nogil except +

      Observations& getObservations() nogil except +

      ParentSequences& getParentSequences() nogil except +

      ParentGroupSets& getParentGroupSets() nogil except +

      #IdentifiedPeptides& getIdentifiedPeptides() nogil except +

      IdentifiedCompounds& getIdentifiedCompounds() nogil except +

      #IdentifiedOligos& getIdentifiedOligos() nogil except +

      #Adducts& getAdducts() nogil except +

      ObservationMatches& getObservationMatches() nogil except +

      ObservationMatchGroups& getObservationMatchGroups() nogil except +

  #cdef cppclass RefTranslator:
  #  libcpp_map[InputFileRef,InputFileRef] input_file_refs
  #  libcpp_map[ScoreTypeRef,ScoreTypeRef] score_type_refs
  #  libcpp_map[ProcessingSoftwareRef,ProcessingSoftwareRef] processing_software_refs
  #  libcpp_map[SearchParamRef,SearchParamRef] search_param_refs
  #  libcpp_map[ProcessingStepRef,ProcessingStepRef] processing_step_refs
  #  libcpp_map[ObservationRef,ObservationRef] observation_refs
  #  libcpp_map[IdentifiedPeptideRef,IdentifiedPeptideRef] identified_peptide_refs
  #  libcpp_map[IdentifiedOligoRef,IdentifiedOligoRef] identified_oligo_refs
  #  libcpp_map[IdentifiedCompoundRef,IdentifiedCompoundRef] identified_compound_refs
  #  libcpp_map[AdductRef,AdductRef] adduct_refs
  #  libcpp_map[ObservationMatchRef,ObservationMatchRef]observation_match_refs

  #  bool allow_missing #= false

  #  IdentifiedMolecule translate(IdentifiedMolecule old) nogil except +

  #  ObservationMatchRef translate(ObservationMatchRef old) nogil except +


