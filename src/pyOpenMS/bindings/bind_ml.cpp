// pyOpenMS nanobind bindings
// Domain: ml

#include "all_casters.h"
#include <OpenMS/ML/CLUSTERING/ClusteringGrid.h>
#include <OpenMS/ML/CLUSTERING/GridBasedCluster.h>
#include <OpenMS/ML/INTERPOLATION/BilinearInterpolation.h>
#include <OpenMS/ML/INTERPOLATION/LinearInterpolation.h>
#include <OpenMS/ML/NNLS/NonNegativeLeastSquaresSolver.h>
#include <OpenMS/ML/RANSAC/RANSAC.h>
#include <OpenMS/ML/RANSAC/RANSACModelLinear.h>
#include <OpenMS/ML/RANSAC/RANSACModelQuadratic.h>
#include <iomanip>
#include <nanobind/make_iterator.h>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/operators.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/set.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/vector.h>
#include <sstream>

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(_pyopenms_ml, m) {
    m.doc() = "pyOpenMS ml bindings";

    // -----------------------------------------------------------------------
    // BilinearInterpolation
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Math::BilinearInterpolation<double, double>>(m, "BilinearInterpolation")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Math::BilinearInterpolation<double, double>&>())
        .def("value", [](OpenMS::Math::BilinearInterpolation<double, double>& self, double arg_pos_0, double arg_pos_1) { return self.value(arg_pos_0, arg_pos_1); }, "arg_pos_0"_a, "arg_pos_1"_a)
        .def("addValue", [](OpenMS::Math::BilinearInterpolation<double, double>& self, double arg_pos_0, double arg_pos_1, double arg_value) { self.addValue(arg_pos_0, arg_pos_1, arg_value); }, "arg_pos_0"_a, "arg_pos_1"_a, "arg_value"_a)
        .def("getData", [](OpenMS::Math::BilinearInterpolation<double, double>& self) { return self.getData(); })
        .def("setData", [](OpenMS::Math::BilinearInterpolation<double, double>& self, OpenMS::Matrix<double> data) { self.setData(data); }, "data"_a)
        .def("empty", [](OpenMS::Math::BilinearInterpolation<double, double>& self) { return self.empty(); })
        .def("key2index_0", [](OpenMS::Math::BilinearInterpolation<double, double>& self, double pos) { return self.key2index_0(pos); }, "pos"_a)
        .def("index2key_0", [](OpenMS::Math::BilinearInterpolation<double, double>& self, double pos) { return self.index2key_0(pos); }, "pos"_a)
        .def("key2index_1", [](OpenMS::Math::BilinearInterpolation<double, double>& self, double pos) { return self.key2index_1(pos); }, "pos"_a)
        .def("index2key_1", [](OpenMS::Math::BilinearInterpolation<double, double>& self, double pos) { return self.index2key_1(pos); }, "pos"_a)
        .def("getScale_0", [](OpenMS::Math::BilinearInterpolation<double, double>& self) { return self.getScale_0(); })
        .def("setScale_0", [](OpenMS::Math::BilinearInterpolation<double, double>& self, double & scale) { self.setScale_0(scale); }, "scale"_a)
        .def("getScale_1", [](OpenMS::Math::BilinearInterpolation<double, double>& self) { return self.getScale_1(); })
        .def("setScale_1", [](OpenMS::Math::BilinearInterpolation<double, double>& self, double & scale) { self.setScale_1(scale); }, "scale"_a)
        .def("getOffset_0", [](OpenMS::Math::BilinearInterpolation<double, double>& self) { return self.getOffset_0(); })
        .def("setOffset_0", [](OpenMS::Math::BilinearInterpolation<double, double>& self, double & offset) { self.setOffset_0(offset); }, "offset"_a)
        .def("getOffset_1", [](OpenMS::Math::BilinearInterpolation<double, double>& self) { return self.getOffset_1(); })
        .def("setOffset_1", [](OpenMS::Math::BilinearInterpolation<double, double>& self, double & offset) { self.setOffset_1(offset); }, "offset"_a)
        .def("setMapping_0", [](OpenMS::Math::BilinearInterpolation<double, double>& self, double & scale, double & inside, double & outside) { self.setMapping_0(scale, inside, outside); }, "scale"_a, "inside"_a, "outside"_a)
        .def("setMapping_0", [](OpenMS::Math::BilinearInterpolation<double, double>& self, double & inside_low, double & outside_low, double & inside_high, double & outside_high) { self.setMapping_0(inside_low, outside_low, inside_high, outside_high); }, "inside_low"_a, "outside_low"_a, "inside_high"_a, "outside_high"_a)
        .def("setMapping_1", [](OpenMS::Math::BilinearInterpolation<double, double>& self, double & scale, double & inside, double & outside) { self.setMapping_1(scale, inside, outside); }, "scale"_a, "inside"_a, "outside"_a)
        .def("setMapping_1", [](OpenMS::Math::BilinearInterpolation<double, double>& self, double & inside_low, double & outside_low, double & inside_high, double & outside_high) { self.setMapping_1(inside_low, outside_low, inside_high, outside_high); }, "inside_low"_a, "outside_low"_a, "inside_high"_a, "outside_high"_a)
        .def("getInsideReferencePoint_0", [](OpenMS::Math::BilinearInterpolation<double, double>& self) { return self.getInsideReferencePoint_0(); })
        .def("getInsideReferencePoint_1", [](OpenMS::Math::BilinearInterpolation<double, double>& self) { return self.getInsideReferencePoint_1(); })
        .def("getOutsideReferencePoint_0", [](OpenMS::Math::BilinearInterpolation<double, double>& self) { return self.getOutsideReferencePoint_0(); })
        .def("getOutsideReferencePoint_1", [](OpenMS::Math::BilinearInterpolation<double, double>& self) { return self.getOutsideReferencePoint_1(); })
        .def("supportMin_0", [](OpenMS::Math::BilinearInterpolation<double, double>& self) { return self.supportMin_0(); })
        .def("supportMin_1", [](OpenMS::Math::BilinearInterpolation<double, double>& self) { return self.supportMin_1(); })
        .def("supportMax_0", [](OpenMS::Math::BilinearInterpolation<double, double>& self) { return self.supportMax_0(); })
        .def("supportMax_1", [](OpenMS::Math::BilinearInterpolation<double, double>& self) { return self.supportMax_1(); })
        ;


    // -----------------------------------------------------------------------
    // ClusteringGrid
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ClusteringGrid>(m, "ClusteringGrid", 
        R"doc(
data structure to store 2D data to be clustered * e.g. (m/z,
retention time) coordinates from multiplex filtering * * @see
LocalClustering
)doc")
        .def(nb::init<std::vector<double>, std::vector<double>>())
        .def("getGridSpacingX", [](const OpenMS::ClusteringGrid& self) { return self.getGridSpacingX(); })
        .def("getGridSpacingY", [](const OpenMS::ClusteringGrid& self) { return self.getGridSpacingY(); })
        .def("addCluster", [](OpenMS::ClusteringGrid& self, const std::pair<int, int>& cell_index, const int& cluster_index) { return self.addCluster(cell_index, cluster_index); }, "cell_index"_a, "cluster_index"_a, "Adds a cluster to this grid cell")
        .def("removeCluster", [](OpenMS::ClusteringGrid& self, const std::pair<int, int>& cell_index, const int& cluster_index) { return self.removeCluster(cell_index, cluster_index); }, "cell_index"_a, "cluster_index"_a, "Removes a cluster from this grid cell and removes the cell if no other cluster left")
        .def("removeAllClusters", [](OpenMS::ClusteringGrid& self) { return self.removeAllClusters(); }, "Removes all clusters from this grid (and hence all cells)")
        .def("getIndex", [](const OpenMS::ClusteringGrid& self, const OpenMS::DPosition<2>& position) { return self.getIndex(position); }, "position"_a)
        .def("isNonEmptyCell", [](const OpenMS::ClusteringGrid& self, const std::pair<int, int>& cell_index) { return self.isNonEmptyCell(cell_index); }, "cell_index"_a, "Checks if there are clusters at this cell index")
        .def("getCellCount", [](const OpenMS::ClusteringGrid& self) { return self.getCellCount(); }, "Returns number of grid cells occupied by one or more clusters")
        ;

    // -----------------------------------------------------------------------
    // GridBasedCluster
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::GridBasedCluster>(m, "GridBasedCluster", "basic data structure for clustering")
        .def(nb::init<OpenMS::DPosition<2>, OpenMS::DBoundingBox<2>, std::vector<int>, int, std::vector<int>>())
        .def(nb::init<OpenMS::DPosition<2>, OpenMS::DBoundingBox<2>, std::vector<int>>())
        .def("getCentre", [](const OpenMS::GridBasedCluster& self) -> const OpenMS::DPosition<2> & { return self.getCentre(); }, nb::rv_policy::reference_internal, "Returns cluster centre")
        .def("getBoundingBox", [](const OpenMS::GridBasedCluster& self) -> const OpenMS::DBoundingBox<2> & { return self.getBoundingBox(); }, nb::rv_policy::reference_internal, "Returns bounding box")
        .def("getPoints", [](const OpenMS::GridBasedCluster& self) -> const std::vector<int> & { return self.getPoints(); }, nb::rv_policy::reference_internal, "Returns indices of points in cluster")
        .def("getPropertyA", [](const OpenMS::GridBasedCluster& self) { return self.getPropertyA(); }, "Returns property A")
        .def("getPropertiesB", [](const OpenMS::GridBasedCluster& self) -> const std::vector<int> & { return self.getPropertiesB(); }, nb::rv_policy::reference_internal, "Returns properties B of all points")
        ;

    // -----------------------------------------------------------------------
    // LinearInterpolation
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Math::LinearInterpolation<double, double>>(m, "LinearInterpolation")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Math::LinearInterpolation<double, double>&>())
        .def(nb::init<double, double>())
        ;


    // -----------------------------------------------------------------------
    // NonNegativeLeastSquaresSolver
    // -----------------------------------------------------------------------
    auto nonnegativeleastsquaressolver_class = nb::class_<OpenMS::NonNegativeLeastSquaresSolver>(m, "NonNegativeLeastSquaresSolver", "Wrapper for a non-negative least squares (NNLS) solver")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::NonNegativeLeastSquaresSolver &>())
        .def_static("solve", [](const OpenMS::Matrix<double>& A, const OpenMS::Matrix<double>& b, OpenMS::Matrix<double>& x) { return OpenMS::NonNegativeLeastSquaresSolver::solve(A, b, x); }, "A"_a, "b"_a, "x"_a)
        .def_static("solve", [](double * A, int A_rows, int A_cols, std::vector<double>& b, std::vector<double>& x) { return OpenMS::NonNegativeLeastSquaresSolver::solve(A, A_rows, A_cols, b, x); }, "A"_a, "A_rows"_a, "A_cols"_a, "b"_a, "x"_a)
        .def_static("solve", [](OpenMS::Matrix<double>& A, std::vector<double>& b, std::vector<double>& x) { return OpenMS::NonNegativeLeastSquaresSolver::solve(A, b, x); }, "A"_a, "b"_a, "x"_a)
        ;
    // RETURN_STATUS enum nested under NonNegativeLeastSquaresSolver
    nb::enum_<OpenMS::NonNegativeLeastSquaresSolver::RETURN_STATUS>(nonnegativeleastsquaressolver_class, "RETURN_STATUS", nb::is_arithmetic())
        .value("SOLVED", OpenMS::NonNegativeLeastSquaresSolver::RETURN_STATUS::SOLVED)
        .value("ITERATION_EXCEEDED", OpenMS::NonNegativeLeastSquaresSolver::RETURN_STATUS::ITERATION_EXCEEDED)
        .export_values();

    // -----------------------------------------------------------------------
    // RANSAC
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Math::RANSAC<OpenMS::Math::RansacModelLinear>>(m, "RANSAC")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Math::RANSAC<OpenMS::Math::RansacModelLinear>&>())
        .def(nb::init<OpenMS::UInt64>())
        ;

    // -----------------------------------------------------------------------
    // RANSACQuadratic
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Math::RANSAC<OpenMS::Math::RansacModelQuadratic>>(m, "RANSACQuadratic")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Math::RANSAC<OpenMS::Math::RansacModelQuadratic>&>())
        .def(nb::init<OpenMS::UInt64>())
        ;


    // -----------------------------------------------------------------------
    // RANSACParam
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Math::RANSACParam>(m, "RANSACParam")
        .def(nb::init<>())
        .def(nb::init<unsigned long, unsigned long, double, unsigned long, bool>())
        .def("toString", [](const OpenMS::Math::RANSACParam& self) { return self.toString(); })
        .def_rw("n", &OpenMS::Math::RANSACParam::n)
        .def_rw("k", &OpenMS::Math::RANSACParam::k)
        .def_rw("t", &OpenMS::Math::RANSACParam::t)
        .def_rw("d", &OpenMS::Math::RANSACParam::d)
        .def_rw("relative_d", &OpenMS::Math::RANSACParam::relative_d)
        ;

    // -----------------------------------------------------------------------
    // RansacModelLinear
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Math::RansacModelLinear>(m, "RansacModelLinear", "Implementation of a linear RANSAC model fit")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Math::RansacModelLinear &>())
        ;

    // -----------------------------------------------------------------------
    // RansacModelQuadratic
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Math::RansacModelQuadratic>(m, "RansacModelQuadratic", "Implementation of a quadratic RANSAC model fit")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Math::RansacModelQuadratic &>())
        ;

}