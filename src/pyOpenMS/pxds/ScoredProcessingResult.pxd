from AppliedProcessingStep cimport *
from MetaInfoInterface cimport *
from libcpp.map cimport map as libcpp_map
#from libcpp.tuple cimport tuple as libcpp_tuple #need to wrap tuple
from MetaData cimport *
from libcpp.pair cimport pair as libcpp_pair


cdef extern from "<OpenMS/METADATA/ID/ScoredProcessingResult.h>" namespace "OpenMS::IdentificationDataInternal":

  cdef cppclass ScoredProcessingResult(MetaInfoInterface):
    AppliedProcessingSteps steps_and_scores
    # AppliedProcessingSteps.nth_index<1>.type& getStepAndScoresByStep() #TODO skipping this for now as the return type is 
    void addProcessingStep(AppliedProcessingStep& step) nogil except +

    void addProcessingStep(ProcessingStepRef, libcpp_map[ScoreTypeRef, double] scores) nogil except +

    #void addScore(ScoreTypeRef score_type, double score, optional[ProcessingStepRef] & processing_step_opt) nogil except + #FIXME wrap boost optional

    ScoredProcessingResult & merge(ScoredProcessingResult & other) nogil except + 

    libcpp_pair[double,bool] getScore(ScoreTypeRef score_ref) nogil except +

    #libcpp_pair[double,bool] getScore(ScoreTypeRef score_ref, boost_optional[ProcessingStepRef] processing_step_opt ) nogil except + FIXME needs Optional

    #libcpp_tuple[double, boost_optional[ProcessingStepRef], bool] getScoreAndStep(ProcessingStepRef score_ref) nogil except + #Needs optional

    #libcpp_tuple[double, boost_optional[ProcessingStepRef], bool] getMostRecentSCore() nogil except + #needs Optional

    size_t getNumberOfScores() nogil except +

    ScoredProcessingResult(ScoredProcessingResult & other) nogil except +