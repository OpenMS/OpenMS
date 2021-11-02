from Types cimport *
from Map cimport *
from Types cimport *
from StringList cimport *
from String cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/HiddenMarkovModel.h>" namespace "OpenMS":
    
    cdef cppclass HiddenMarkovModel "OpenMS::HiddenMarkovModel":

        HiddenMarkovModel() nogil except + # wrap-doc:Hidden Markov Model implementation of PILIS
        HiddenMarkovModel(HiddenMarkovModel &) nogil except +
        void writeGraphMLFile(const String & filename) nogil except + # wrap-doc:Writes the HMM into a file in GraphML format

        # NAMESPACE # void write(std::ostream & out) nogil except +
        double getTransitionProbability(const String & s1, const String & s2) nogil except + # wrap-doc:Returns the transition probability of the given state names
        void setTransitionProbability(const String & s1, const String & s2, double prob) nogil except + # wrap-doc:Sets the transition probability of the given state names to prob
        Size getNumberOfStates() nogil except + # wrap-doc:Returns the number of states
        void addNewState(HMMState * state) nogil except + # wrap-doc:Registers a new state to the HMM
        void addNewState(const String & name) nogil except + # wrap-doc:Registers a new state to the HMM
        void addSynonymTransition(const String & name1, const String & name2, const String & synonym1, const String & synonym2) nogil except + # wrap-doc:Add a new synonym transition to the given state names
        void evaluate() nogil except + # wrap-doc:Evaluate the HMM, estimates the transition probabilities from the training
        void train() nogil except + # wrap-doc:Trains the HMM. Initial probabilities and emission probabilities of the emitting states should be set
        void setInitialTransitionProbability(const String & state, double prob) nogil except + # wrap-doc:Sets the initial transition probability of the given state to prob
        void clearInitialTransitionProbabilities() nogil except + # wrap-doc:Clears the initial probabilities
        void setTrainingEmissionProbability(const String & state, double prob) nogil except + # wrap-doc:Sets the emission probability of the given state to prob
        void clearTrainingEmissionProbabilities() nogil except + # wrap-doc:Clear the emission probabilities
        void enableTransition(const String & s1, const String & s2) nogil except + # wrap-doc:Enables a transition; adds s1 to the predecessor list of s2 and s2 to the successor list of s1
        void disableTransition(const String & s1, const String & s2) nogil except + # wrap-doc:Disables the transition; deletes the nodes from the predecessor/successor list respectively
        void disableTransitions() nogil except + # wrap-doc:Disables all transitions
        # POINTER # void calculateEmissionProbabilities(Map[ HMMState *, double ] & emission_probs) nogil except +
        void dump() nogil except + # wrap-doc:Writes some stats to cerr
        void forwardDump() nogil except + # wrap-doc:Writes some info of the forward "matrix" to cerr
        void estimateUntrainedTransitions() nogil except + # wrap-doc:Estimates the transition probabilities of not trained transitions by averages similar trained ones
        HMMState * getState(const String & name) nogil except + # wrap-doc:Returns the state with the given name
        void clear() nogil except + # wrap-doc:Clears all data
        void setPseudoCounts(double pseudo_counts) nogil except + # wrap-doc:Sets the pseudo count that are added instead of zero
        double getPseudoCounts() nogil except + # wrap-doc:Returns the pseudo counts
        void setVariableModifications(StringList & modifications) nogil except +

cdef extern from "<OpenMS/ANALYSIS/ID/HiddenMarkovModel.h>" namespace "OpenMS":
    
    cdef cppclass HMMState "OpenMS::HMMState":
        HMMState() nogil except +
        HMMState(HMMState &) nogil except +
        HMMState(const String & name, bool hidden) nogil except +

        # They don't exist ...
        # bool operator==(HMMState) nogil except +
        # bool operator!=(HMMState) nogil except +

        void setName(const String & name) nogil except + # wrap-doc:Sets the name of the state
        String getName() nogil except + # wrap-doc:Returns the name of the state
        void setHidden(bool hidden) nogil except + # wrap-doc:Sets the hidden property to the state
        bool isHidden() nogil except + # wrap-doc:Returns true if the state is hidden
        void addPredecessorState(HMMState * state) nogil except + # wrap-doc:Adds the given predecessor state to the list
        void deletePredecessorState(HMMState * state) nogil except + # wrap-doc:Deletes the given predecessor state from the list
        void addSuccessorState(HMMState * state) nogil except + # wrap-doc:Add the given successor state to the list
        void deleteSuccessorState(HMMState * state) nogil except + # wrap-doc:Deletes the given successor state from the list
        libcpp_set[ HMMState * ]  getPredecessorStates() nogil except + # wrap-doc:Returns the predecessor states of the state
        libcpp_set[ HMMState * ]  getSuccessorStates() nogil except + # wrap-doc:Returns the successor states of the state
