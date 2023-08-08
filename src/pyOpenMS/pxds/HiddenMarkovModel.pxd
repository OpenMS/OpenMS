from Types cimport *
from libcpp.map cimport map as libcpp_map
from Types cimport *
from StringList cimport *
from String cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/HiddenMarkovModel.h>" namespace "OpenMS":
    
    cdef cppclass HiddenMarkovModel "OpenMS::HiddenMarkovModel":

        HiddenMarkovModel() except + nogil  # wrap-doc:Hidden Markov Model implementation of PILIS
        HiddenMarkovModel(HiddenMarkovModel &) except + nogil 
        void writeGraphMLFile(const String & filename) except + nogil  # wrap-doc:Writes the HMM into a file in GraphML format

        # NAMESPACE # void write(std::ostream & out) except + nogil 
        double getTransitionProbability(const String & s1, const String & s2) except + nogil  # wrap-doc:Returns the transition probability of the given state names
        void setTransitionProbability(const String & s1, const String & s2, double prob) except + nogil  # wrap-doc:Sets the transition probability of the given state names to prob
        Size getNumberOfStates() except + nogil  # wrap-doc:Returns the number of states
        void addNewState(HMMState * state) except + nogil  # wrap-doc:Registers a new state to the HMM
        void addNewState(const String & name) except + nogil  # wrap-doc:Registers a new state to the HMM
        void addSynonymTransition(const String & name1, const String & name2, const String & synonym1, const String & synonym2) except + nogil  # wrap-doc:Add a new synonym transition to the given state names
        void evaluate() except + nogil  # wrap-doc:Evaluate the HMM, estimates the transition probabilities from the training
        void train() except + nogil  # wrap-doc:Trains the HMM. Initial probabilities and emission probabilities of the emitting states should be set
        void setInitialTransitionProbability(const String & state, double prob) except + nogil  # wrap-doc:Sets the initial transition probability of the given state to prob
        void clearInitialTransitionProbabilities() except + nogil  # wrap-doc:Clears the initial probabilities
        void setTrainingEmissionProbability(const String & state, double prob) except + nogil  # wrap-doc:Sets the emission probability of the given state to prob
        void clearTrainingEmissionProbabilities() except + nogil  # wrap-doc:Clear the emission probabilities
        void enableTransition(const String & s1, const String & s2) except + nogil  # wrap-doc:Enables a transition; adds s1 to the predecessor list of s2 and s2 to the successor list of s1
        void disableTransition(const String & s1, const String & s2) except + nogil  # wrap-doc:Disables the transition; deletes the nodes from the predecessor/successor list respectively
        void disableTransitions() except + nogil  # wrap-doc:Disables all transitions
        # POINTER # void calculateEmissionProbabilities(Map[ HMMState *, double ] & emission_probs) except + nogil 
        void dump() except + nogil  # wrap-doc:Writes some stats to cerr
        void forwardDump() except + nogil  # wrap-doc:Writes some info of the forward "matrix" to cerr
        void estimateUntrainedTransitions() except + nogil  # wrap-doc:Estimates the transition probabilities of not trained transitions by averages similar trained ones
        HMMState * getState(const String & name) except + nogil  # wrap-doc:Returns the state with the given name
        void clear() except + nogil  # wrap-doc:Clears all data
        void setPseudoCounts(double pseudo_counts) except + nogil  # wrap-doc:Sets the pseudo count that are added instead of zero
        double getPseudoCounts() except + nogil  # wrap-doc:Returns the pseudo counts
        void setVariableModifications(StringList & modifications) except + nogil 

cdef extern from "<OpenMS/ANALYSIS/ID/HiddenMarkovModel.h>" namespace "OpenMS":
    
    cdef cppclass HMMState "OpenMS::HMMState":
        HMMState() except + nogil 
        HMMState(HMMState &) except + nogil 
        HMMState(const String & name, bool hidden) except + nogil 

        # They don't exist ...
        # bool operator==(HMMState) except + nogil 
        # bool operator!=(HMMState) except + nogil 

        void setName(const String & name) except + nogil  # wrap-doc:Sets the name of the state
        String getName() except + nogil  # wrap-doc:Returns the name of the state
        void setHidden(bool hidden) except + nogil  # wrap-doc:Sets the hidden property to the state
        bool isHidden() except + nogil  # wrap-doc:Returns true if the state is hidden
        void addPredecessorState(HMMState * state) except + nogil  # wrap-doc:Adds the given predecessor state to the list
        void deletePredecessorState(HMMState * state) except + nogil  # wrap-doc:Deletes the given predecessor state from the list
        void addSuccessorState(HMMState * state) except + nogil  # wrap-doc:Add the given successor state to the list
        void deleteSuccessorState(HMMState * state) except + nogil  # wrap-doc:Deletes the given successor state from the list
        libcpp_set[ HMMState * ]  getPredecessorStates() except + nogil  # wrap-doc:Returns the predecessor states of the state
        libcpp_set[ HMMState * ]  getSuccessorStates() except + nogil  # wrap-doc:Returns the successor states of the state
