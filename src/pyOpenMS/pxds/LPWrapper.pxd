from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from String cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/LPWrapper.h>" namespace "OpenMS":
    
    cdef cppclass LPWrapper "OpenMS::LPWrapper":
        LPWrapper() except + nogil 
        LPWrapper(LPWrapper) except + nogil  #wrap-ignore
        Int addRow(libcpp_vector[ int ] row_indices, libcpp_vector[ double ] row_values, const String & name) except + nogil  # wrap-doc:Adds a row to the LP matrix, returns index
        Int addColumn() except + nogil  # wrap-doc:Adds an empty column to the LP matrix, returns index
        Int addColumn(libcpp_vector[ int ] column_indices, libcpp_vector[ double ] column_values, const String & name) except + nogil  # wrap-doc:Adds a column to the LP matrix, returns index
        Int addRow(libcpp_vector[ int ] & row_indices, libcpp_vector[ double ] & row_values, const String & name, double lower_bound, double upper_bound, LPWrapper_Type type_) except + nogil  # wrap-doc:Adds a row with boundaries to the LP matrix, returns index
        Int addColumn(libcpp_vector[ int ] & column_indices, libcpp_vector[ double ] & column_values, const String & name, double lower_bound, double upper_bound, LPWrapper_Type type_) except + nogil  # wrap-doc:Adds a column with boundaries to the LP matrix, returns index
        void deleteRow(Int index) except + nogil  # wrap-doc:Delete index-th row
        void setColumnName(Int index, const String & name) except + nogil  # wrap-doc:Sets name of the index-th column
        String getColumnName(Int index) except + nogil  # wrap-doc:Returns name of the index-th column
        String getRowName(Int index) except + nogil  # wrap-doc:Sets name of the index-th row
        Int getRowIndex(const String & name) except + nogil  # wrap-doc:Returns index of the row with name
        Int getColumnIndex(const String & name) except + nogil  # wrap-doc:Returns index of the column with name
        double getColumnUpperBound(Int index) except + nogil  # wrap-doc:Returns column's upper bound
        double getColumnLowerBound(Int index) except + nogil  # wrap-doc:Returns column's lower bound
        double getRowUpperBound(Int index) except + nogil  # wrap-doc:Returns row's upper bound
        double getRowLowerBound(Int index) except + nogil  # wrap-doc:Returns row's lower bound
        void setRowName(Int index, const String & name) except + nogil  # wrap-doc:Sets name of the index-th row
        void setColumnBounds(Int index, double lower_bound, double upper_bound, LPWrapper_Type type_) except + nogil  # wrap-doc:Sets column bounds
        void setRowBounds(Int index, double lower_bound, double upper_bound, LPWrapper_Type type_) except + nogil  # wrap-doc:Sets row bounds
        void setColumnType(Int index, VariableType type_) except + nogil  # wrap-doc:Sets column/variable type.
        VariableType getColumnType(Int index) except + nogil  # wrap-doc:Returns column/variable type.
        void setObjective(Int index, double obj_value) except + nogil  # wrap-doc:Sets objective value for column with index
        double getObjective(Int index) except + nogil  # wrap-doc:Returns objective value for column with index
        void setObjectiveSense(Sense sense) except + nogil  # wrap-doc:Sets objective direction
        Sense getObjectiveSense() except + nogil  # wrap-doc:Returns objective sense
        Int getNumberOfColumns() except + nogil  # wrap-doc:Returns number of columns
        Int getNumberOfRows() except + nogil  # wrap-doc:Returns number of rows
        void setElement(Int row_index, Int column_index, double value) except + nogil   # wrap-doc:Sets the element
        double getElement(Int row_index, Int column_index) except + nogil  # wrap-doc:Returns the element
        void readProblem(String filename, String format_) except + nogil 
            # wrap-doc:
                #  Read LP from file
                #  
                #  
                #  :param filename: Filename where to store the LP problem
                #  :param format: LP, MPS or GLPK

        void writeProblem(const String & filename, WriteFormat format_) except + nogil 
            # wrap-doc:
                #  Write LP formulation to a file
                #  
                #  
                #  :param filename: Output filename, if the filename ends with '.gz' it will be compressed
                #  :param format: MPS-format is supported by GLPK and COIN-OR; LP and GLPK-formats only by GLPK

        Int solve(SolverParam & solver_param, Size verbose_level) except + nogil 
            # wrap-doc:
                #  Solve problems, parameters like enabled heuristics can be given via solver_param\n
                #  
                #  The verbose level (0,1,2) determines if the solver prints status messages and internals
                #  
                #  
                #  :param solver_param: Parameters of the solver introduced by SolverParam
                #  :param verbose_level: Sets verbose level
                #  :return: solver dependent 

        SolverStatus getStatus() except + nogil 
            # wrap-doc:
                #  Returns solution status
                #  
                #  
                #  :return: status: 1 - undefined, 2 - integer optimal, 3- integer feasible (no optimality proven), 4- no integer feasible solution

        double getObjectiveValue() except + nogil 
        double getColumnValue(Int index) except + nogil 
        Int getNumberOfNonZeroEntriesInRow(Int idx) except + nogil 
        void getMatrixRow(Int idx, libcpp_vector[ int ] & indexes) except + nogil 
        SOLVER getSolver() except + nogil  # wrap-doc:Returns currently active solver

cdef extern from "<OpenMS/DATASTRUCTURES/LPWrapper.h>" namespace "OpenMS::LPWrapper":
    
    cdef cppclass SolverParam "OpenMS::LPWrapper::SolverParam":
        SolverParam() except + nogil  # wrap-doc:Hold the parameters of the LP solver
        SolverParam(SolverParam) except + nogil  #wrap-ignore

        Int message_level
        Int branching_tech
        Int backtrack_tech
        Int preprocessing_tech
        bool enable_feas_pump_heuristic
        bool enable_gmi_cuts
        bool enable_mir_cuts
        bool enable_cov_cuts
        bool enable_clq_cuts
        double mip_gap
        Int time_limit
        Int output_freq
        Int output_delay
        bool enable_presolve
        bool enable_binarization


cdef extern from "<OpenMS/DATASTRUCTURES/LPWrapper.h>" namespace "OpenMS::LPWrapper":
    cdef enum LPWrapper_Type "OpenMS::LPWrapper::Type":
        #wrap-attach:
        #   LPWrapper
        UNBOUNDED
        LOWER_BOUND_ONLY
        UPPER_BOUND_ONLY
        DOUBLE_BOUNDED
        FIXED

cdef extern from "<OpenMS/DATASTRUCTURES/LPWrapper.h>" namespace "OpenMS::LPWrapper":
    cdef enum VariableType "OpenMS::LPWrapper::VariableType":
        #wrap-attach:
        #   LPWrapper
        CONTINUOUS
        INTEGER
        BINARY

cdef extern from "<OpenMS/DATASTRUCTURES/LPWrapper.h>" namespace "OpenMS::LPWrapper":
    cdef enum Sense "OpenMS::LPWrapper::Sense":
        #wrap-attach:
        #   LPWrapper
        MIN
        MAX

cdef extern from "<OpenMS/DATASTRUCTURES/LPWrapper.h>" namespace "OpenMS::LPWrapper":
    cdef enum WriteFormat "OpenMS::LPWrapper::WriteFormat":
        #wrap-attach:
        #   LPWrapper
        FORMAT_LP
        FORMAT_MPS
        FORMAT_GLPK

cdef extern from "<OpenMS/DATASTRUCTURES/LPWrapper.h>" namespace "OpenMS::LPWrapper":
    cdef enum SOLVER "OpenMS::LPWrapper::SOLVER":
        #wrap-attach:
        #   LPWrapper
        SOLVER_GLPK

cdef extern from "<OpenMS/DATASTRUCTURES/LPWrapper.h>" namespace "OpenMS::LPWrapper":
    cdef enum SolverStatus "OpenMS::LPWrapper::SolverStatus":
        #wrap-attach:
        #   LPWrapper
        UNDEFINED
        OPTIMAL
        FEASIBLE
        NO_FEASIBLE_SOL
