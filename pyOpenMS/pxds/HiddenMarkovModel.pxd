from Types cimport *
from Map cimport *
from Types cimport *
from StringList cimport *
from String cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/HiddenMarkovModel.h>" namespace "OpenMS":
    
    cdef cppclass HiddenMarkovModel "OpenMS::HiddenMarkovModel":
        HiddenMarkovModel() nogil except +
        HiddenMarkovModel(HiddenMarkovModel) nogil except +
        void writeGraphMLFile(String & filename) nogil except +
        # NAMESPACE # void write(std::ostream & out) nogil except +
        DoubleReal getTransitionProbability(String & s1, String & s2) nogil except +
        void setTransitionProbability(String & s1, String & s2, DoubleReal prob) nogil except +
        Size getNumberOfStates() nogil except +
        # POINTER # void addNewState(HMMState * state) nogil except +
        void addNewState(String & name) nogil except +
        void addSynonymTransition(String & name1, String & name2, String & synonym1, String & synonym2) nogil except +
        void evaluate() nogil except +
        void train() nogil except +
        void setInitialTransitionProbability(String & state, DoubleReal prob) nogil except +
        void clearInitialTransitionProbabilities() nogil except +
        void setTrainingEmissionProbability(String & state, DoubleReal prob) nogil except +
        void clearTrainingEmissionProbabilities() nogil except +
        void enableTransition(String & s1, String & s2) nogil except +
        void disableTransition(String & s1, String & s2) nogil except +
        void disableTransitions() nogil except +
        # POINTER # void calculateEmissionProbabilities(Map[ HMMState *, DoubleReal ] & emission_probs) nogil except +
        void dump() nogil except +
        void forwardDump() nogil except +
        void estimateUntrainedTransitions() nogil except +
        # POINTER # HMMState * getState(String & name) nogil except +
        # POINTER # HMMState * getState(String & name) nogil except +
        void clear() nogil except +
        void setPseudoCounts(DoubleReal pseudo_counts) nogil except +
        DoubleReal getPseudoCounts() nogil except +
        void setVariableModifications(StringList & modifications) nogil except +

cdef extern from "<OpenMS/ANALYSIS/ID/HiddenMarkovModel.h>" namespace "OpenMS":
    
    cdef cppclass HMMState "OpenMS::HMMState":
        HMMState() nogil except +
        HMMState(HMMState) nogil except +
        HMMState(String & name, bool hidden) nogil except +
        void setName(String & name) nogil except +
        String  getName() nogil except +
        void setHidden(bool hidden) nogil except +
        bool isHidden() nogil except +
        # POINTER # void addPredecessorState(HMMState * state) nogil except +
        # POINTER # void deletePredecessorState(HMMState * state) nogil except +
        # POINTER # void addSuccessorState(HMMState * state) nogil except +
        # POINTER # void deleteSuccessorState(HMMState * state) nogil except +
        # POINTER # libcpp_set[ HMMState * ]  getPredecessorStates() nogil except +
        # POINTER # libcpp_set[ HMMState * ]  getSuccessorStates() nogil except +

