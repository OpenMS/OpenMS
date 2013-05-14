from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from String cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/LPWrapper.h>" namespace "OpenMS":
    
    cdef cppclass LPWrapper "OpenMS::LPWrapper":
        LPWrapper() nogil except +
        LPWrapper(LPWrapper) nogil except + #wrap-ignore
        Int addRow(libcpp_vector[ int ] row_indices, libcpp_vector[ double ] row_values, String & name) nogil except +
        Int addColumn() nogil except +
        Int addColumn(libcpp_vector[ int ] column_indices, libcpp_vector[ double ] column_values, String & name) nogil except +
        Int addRow(libcpp_vector[ int ] & row_indices, libcpp_vector[ double ] & row_values, String & name, DoubleReal lower_bound, DoubleReal upper_bound, LPWrapper_Type type_) nogil except +
        Int addColumn(libcpp_vector[ int ] & column_indices, libcpp_vector[ double ] & column_values, String & name, DoubleReal lower_bound, DoubleReal upper_bound, LPWrapper_Type type_) nogil except +
        void deleteRow(Int index) nogil except +
        void setColumnName(Int index, String & name) nogil except +
        String getColumnName(Int index) nogil except +
        String getRowName(Int index) nogil except +
        Int getRowIndex(String & name) nogil except +
        Int getColumnIndex(String & name) nogil except +
        DoubleReal getColumnUpperBound(Int index) nogil except +
        DoubleReal getColumnLowerBound(Int index) nogil except +
        DoubleReal getRowUpperBound(Int index) nogil except +
        DoubleReal getRowLowerBound(Int index) nogil except +
        void setRowName(Int index, String & name) nogil except +
        void setColumnBounds(Int index, DoubleReal lower_bound, DoubleReal upper_bound, LPWrapper_Type type_) nogil except +
        void setRowBounds(Int index, DoubleReal lower_bound, DoubleReal upper_bound, LPWrapper_Type type_) nogil except +
        void setColumnType(Int index, VariableType type_) nogil except +
        VariableType getColumnType(Int index) nogil except +
        void setObjective(Int index, DoubleReal obj_value) nogil except +
        DoubleReal getObjective(Int index) nogil except +
        void setObjectiveSense(Sense sense) nogil except +
        Sense getObjectiveSense() nogil except +
        Int getNumberOfColumns() nogil except +
        Int getNumberOfRows() nogil except +
        void setElement(Int row_index, Int column_index, DoubleReal value) nogil except +
        DoubleReal getElement(Int row_index, Int column_index) nogil except +
        void readProblem(String filename, String format_) nogil except +
        void writeProblem(String & filename, WriteFormat format_) nogil except +
        Int solve(SolverParam & solver_param, Size verbose_level) nogil except +
        SolverStatus getStatus() nogil except +
        DoubleReal getObjectiveValue() nogil except +
        DoubleReal getColumnValue(Int index) nogil except +
        Int getNumberOfNonZeroEntriesInRow(Int idx) nogil except +
        void getMatrixRow(Int idx, libcpp_vector[ int ] & indexes) nogil except +
        void setSolver(SOLVER s) nogil except +
        SOLVER getSolver() nogil except +

cdef extern from "<OpenMS/DATASTRUCTURES/LPWrapper.h>" namespace "OpenMS::LPWrapper":
    
    cdef cppclass SolverParam "OpenMS::LPWrapper::SolverParam":
        #wrap-attach:
        #    LPWrapper
        SolverParam() nogil except +
        SolverParam(SolverParam) nogil except + #wrap-ignore
        Int message_level
        Int branching_tech
        Int backtrack_tech
        Int preprocessing_tech
        bool enable_feas_pump_heuristic
        bool enable_gmi_cuts
        bool enable_mir_cuts
        bool enable_cov_cuts
        bool enable_clq_cuts
        DoubleReal mip_gap
        Int time_limit
        Int output_freq
        Int output_delay
        bool enable_presolve
        bool enable_binarization


cdef extern from "<OpenMS/DATASTRUCTURES/LPWrapper.h>" namespace "OpenMS::LPWrapper":
    cdef enum LPWrapper_Type "OpenMS::LPWrapper::Type":
        #wrap-attach:
        #    LPWrapper
        UNBOUNDED
        LOWER_BOUND_ONLY
        UPPER_BOUND_ONLY
        DOUBLE_BOUNDED
        FIXED

cdef extern from "<OpenMS/DATASTRUCTURES/LPWrapper.h>" namespace "OpenMS::LPWrapper":
    cdef enum VariableType "OpenMS::LPWrapper::VariableType":
        #wrap-attach:
        #    LPWrapper
        CONTINUOUS
        INTEGER
        BINARY

cdef extern from "<OpenMS/DATASTRUCTURES/LPWrapper.h>" namespace "OpenMS::LPWrapper":
    cdef enum Sense "OpenMS::LPWrapper::Sense":
        #wrap-attach:
        #    LPWrapper
        MIN
        MAX

cdef extern from "<OpenMS/DATASTRUCTURES/LPWrapper.h>" namespace "OpenMS::LPWrapper":
    cdef enum WriteFormat "OpenMS::LPWrapper::WriteFormat":
        #wrap-attach:
        #    LPWrapper
        FORMAT_LP
        FORMAT_MPS
        FORMAT_GLPK

cdef extern from "<OpenMS/DATASTRUCTURES/LPWrapper.h>" namespace "OpenMS::LPWrapper":
    cdef enum SOLVER "OpenMS::LPWrapper::SOLVER":
        #wrap-attach:
        #    LPWrapper
        SOLVER_GLPK

cdef extern from "<OpenMS/DATASTRUCTURES/LPWrapper.h>" namespace "OpenMS::LPWrapper":
    cdef enum SolverStatus "OpenMS::LPWrapper::SolverStatus":
        #wrap-attach:
        #    LPWrapper
        UNDEFINED
        OPTIMAL
        FEASIBLE
        NO_FEASIBLE_SOL

