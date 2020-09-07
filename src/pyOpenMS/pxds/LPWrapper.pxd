from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from String cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/LPWrapper.h>" namespace "OpenMS":
    
    cdef cppclass LPWrapper "OpenMS::LPWrapper":
        LPWrapper() nogil except +
        LPWrapper(LPWrapper) nogil except + #wrap-ignore
        Int addRow(libcpp_vector[ int ] row_indices, libcpp_vector[ double ] row_values, const String & name) nogil except +
        Int addColumn() nogil except +
        Int addColumn(libcpp_vector[ int ] column_indices, libcpp_vector[ double ] column_values, const String & name) nogil except +
        Int addRow(libcpp_vector[ int ] & row_indices, libcpp_vector[ double ] & row_values, const String & name, double lower_bound, double upper_bound, LPWrapper_Type type_) nogil except +
        Int addColumn(libcpp_vector[ int ] & column_indices, libcpp_vector[ double ] & column_values, const String & name, double lower_bound, double upper_bound, LPWrapper_Type type_) nogil except +
        void deleteRow(Int index) nogil except +
        void setColumnName(Int index, const String & name) nogil except +
        String getColumnName(Int index) nogil except +
        String getRowName(Int index) nogil except +
        Int getRowIndex(const String & name) nogil except +
        Int getColumnIndex(const String & name) nogil except +
        double getColumnUpperBound(Int index) nogil except +
        double getColumnLowerBound(Int index) nogil except +
        double getRowUpperBound(Int index) nogil except +
        double getRowLowerBound(Int index) nogil except +
        void setRowName(Int index, const String & name) nogil except +
        void setColumnBounds(Int index, double lower_bound, double upper_bound, LPWrapper_Type type_) nogil except +
        void setRowBounds(Int index, double lower_bound, double upper_bound, LPWrapper_Type type_) nogil except +
        void setColumnType(Int index, VariableType type_) nogil except +
        VariableType getColumnType(Int index) nogil except +
        void setObjective(Int index, double obj_value) nogil except +
        double getObjective(Int index) nogil except +
        void setObjectiveSense(Sense sense) nogil except +
        Sense getObjectiveSense() nogil except +
        Int getNumberOfColumns() nogil except +
        Int getNumberOfRows() nogil except +
        void setElement(Int row_index, Int column_index, double value) nogil except +
        double getElement(Int row_index, Int column_index) nogil except +
        void readProblem(String filename, String format_) nogil except +
        void writeProblem(const String & filename, WriteFormat format_) nogil except +
        Int solve(SolverParam & solver_param, Size verbose_level) nogil except +
        SolverStatus getStatus() nogil except +
        double getObjectiveValue() nogil except +
        double getColumnValue(Int index) nogil except +
        Int getNumberOfNonZeroEntriesInRow(Int idx) nogil except +
        void getMatrixRow(Int idx, libcpp_vector[ int ] & indexes) nogil except +
        SOLVER getSolver() nogil except +

cdef extern from "<OpenMS/DATASTRUCTURES/LPWrapper.h>" namespace "OpenMS::LPWrapper":
    
    cdef cppclass SolverParam "OpenMS::LPWrapper::SolverParam":
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
        double mip_gap
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

