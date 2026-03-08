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
    nb::class_<OpenMS::Math::BilinearInterpolation<double, double>>(m, "BilinearInterpolation", "OpenMS class BilinearInterpolation")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Math::BilinearInterpolation<double, double>&>())
        .def("__copy__", [](const OpenMS::Math::BilinearInterpolation<double, double>& self) { return OpenMS::Math::BilinearInterpolation<double, double>(self); })
        .def("__deepcopy__", [](const OpenMS::Math::BilinearInterpolation<double, double>& self, nb::dict) { return OpenMS::Math::BilinearInterpolation<double, double>(self); }, "memo"_a)
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
    nb::class_<OpenMS::Math::LinearInterpolation<double, double>>(m, "LinearInterpolation", "OpenMS class LinearInterpolation")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Math::LinearInterpolation<double, double>&>())
        .def("__copy__", [](const OpenMS::Math::LinearInterpolation<double, double>& self) { return OpenMS::Math::LinearInterpolation<double, double>(self); })
        .def("__deepcopy__", [](const OpenMS::Math::LinearInterpolation<double, double>& self, nb::dict) { return OpenMS::Math::LinearInterpolation<double, double>(self); }, "memo"_a)
        .def(nb::init<double, double>())
        .def("value", [](const OpenMS::Math::LinearInterpolation<double, double>& self, double pos) { return self.value(pos); }, "pos"_a, "Returns interpolated value at position")
        .def("addValue", [](OpenMS::Math::LinearInterpolation<double, double>& self, double pos, double value) { self.addValue(pos, value); }, "pos"_a, "value"_a, "Adds a value at the given position")
        .def("derivative", [](const OpenMS::Math::LinearInterpolation<double, double>& self, double pos) { return self.derivative(pos); }, "pos"_a, "Returns the derivative at position")
        .def("empty", [](const OpenMS::Math::LinearInterpolation<double, double>& self) { return self.empty(); }, "Returns true if the data is empty")
        .def("getData", [](OpenMS::Math::LinearInterpolation<double, double>& self) -> std::vector<double>& { return self.getData(); }, nb::rv_policy::reference_internal, "Returns the data vector")
        .def("setData", [](OpenMS::Math::LinearInterpolation<double, double>& self, const std::vector<double>& data) { self.setData(data); }, "data"_a, "Sets the data")
        .def("getScale", [](const OpenMS::Math::LinearInterpolation<double, double>& self) { return self.getScale(); }, "Returns the scale")
        .def("setScale", [](OpenMS::Math::LinearInterpolation<double, double>& self, double scale) { self.setScale(scale); }, "scale"_a, "Sets the scale")
        .def("getOffset", [](const OpenMS::Math::LinearInterpolation<double, double>& self) { return self.getOffset(); }, "Returns the offset")
        .def("setOffset", [](OpenMS::Math::LinearInterpolation<double, double>& self, double offset) { self.setOffset(offset); }, "offset"_a, "Sets the offset")
        .def("setMapping", [](OpenMS::Math::LinearInterpolation<double, double>& self, double scale, double inside, double outside) { self.setMapping(scale, inside, outside); }, "scale"_a, "inside"_a, "outside"_a, "Sets the mapping")
        .def("setMapping", [](OpenMS::Math::LinearInterpolation<double, double>& self, double inside_low, double outside_low, double inside_high, double outside_high) { self.setMapping(inside_low, outside_low, inside_high, outside_high); }, "inside_low"_a, "outside_low"_a, "inside_high"_a, "outside_high"_a, "Sets the mapping from two pairs of reference points")
        .def("getInsideReferencePoint", [](const OpenMS::Math::LinearInterpolation<double, double>& self) { return self.getInsideReferencePoint(); }, "Returns the inside reference point")
        .def("getOutsideReferencePoint", [](const OpenMS::Math::LinearInterpolation<double, double>& self) { return self.getOutsideReferencePoint(); }, "Returns the outside reference point")
        .def("supportMin", [](const OpenMS::Math::LinearInterpolation<double, double>& self) { return self.supportMin(); }, "Returns the minimum support position")
        .def("supportMax", [](const OpenMS::Math::LinearInterpolation<double, double>& self) { return self.supportMax(); }, "Returns the maximum support position")
        .def("key2index", [](const OpenMS::Math::LinearInterpolation<double, double>& self, double pos) { return self.key2index(pos); }, "pos"_a, "Transforms key to index")
        .def("index2key", [](const OpenMS::Math::LinearInterpolation<double, double>& self, double pos) { return self.index2key(pos); }, "pos"_a, "Transforms index to key")
        ;


    // -----------------------------------------------------------------------
    // NonNegativeLeastSquaresSolver
    // -----------------------------------------------------------------------
    auto nonnegativeleastsquaressolver_class = nb::class_<OpenMS::NonNegativeLeastSquaresSolver>(m, "NonNegativeLeastSquaresSolver", "Wrapper for a non-negative least squares (NNLS) solver")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::NonNegativeLeastSquaresSolver &>())
        .def("__copy__", [](const OpenMS::NonNegativeLeastSquaresSolver& self) { return OpenMS::NonNegativeLeastSquaresSolver(self); })
        .def("__deepcopy__", [](const OpenMS::NonNegativeLeastSquaresSolver& self, nb::dict) { return OpenMS::NonNegativeLeastSquaresSolver(self); }, "memo"_a)
        .def_static("solve", [](const OpenMS::Matrix<double>& A, const OpenMS::Matrix<double>& b) {
            OpenMS::Matrix<double> x;
            auto status = OpenMS::NonNegativeLeastSquaresSolver::solve(A, b, x);
            return nb::make_tuple(status, x);
        }, "A"_a, "b"_a,
        "Solve Ax=b (x>=0). Returns tuple(status, x). Copies A and b.")
        // Note: raw pointer overload (double*, int, int, vector&, vector&) is omitted
        // because it is not safely wrappable from Python.
        .def_static("solve", [](OpenMS::Matrix<double> A, std::vector<double> b) {
            std::vector<double> x;
            auto status = OpenMS::NonNegativeLeastSquaresSolver::solve(A, b, x);
            return nb::make_tuple(status, x);
        }, "A"_a, "b"_a,
        "Solve Ax=b (x>=0) in-place. Returns tuple(status, x). Modifies A and b internally.")
        ;
    // RETURN_STATUS enum nested under NonNegativeLeastSquaresSolver
    nb::enum_<OpenMS::NonNegativeLeastSquaresSolver::RETURN_STATUS>(nonnegativeleastsquaressolver_class, "RETURN_STATUS")
        .value("SOLVED", OpenMS::NonNegativeLeastSquaresSolver::RETURN_STATUS::SOLVED)
        .value("ITERATION_EXCEEDED", OpenMS::NonNegativeLeastSquaresSolver::RETURN_STATUS::ITERATION_EXCEEDED)
        .export_values();

    // -----------------------------------------------------------------------
    // RANSAC
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Math::RANSAC<OpenMS::Math::RansacModelLinear>>(m, "RANSAC", "OpenMS class RANSAC")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Math::RANSAC<OpenMS::Math::RansacModelLinear>&>())
        .def("__copy__", [](const OpenMS::Math::RANSAC<OpenMS::Math::RansacModelLinear>& self) { return OpenMS::Math::RANSAC<OpenMS::Math::RansacModelLinear>(self); })
        .def("__deepcopy__", [](const OpenMS::Math::RANSAC<OpenMS::Math::RansacModelLinear>& self, nb::dict) { return OpenMS::Math::RANSAC<OpenMS::Math::RansacModelLinear>(self); }, "memo"_a)
        .def(nb::init<OpenMS::UInt64>())
        .def("ransac", [](OpenMS::Math::RANSAC<OpenMS::Math::RansacModelLinear>& self, const std::vector<std::pair<double, double>>& pairs, size_t n, size_t k, double t, size_t d, bool relative_d) {
            return self.ransac(pairs, n, k, t, d, relative_d);
        }, "pairs"_a, "n"_a, "k"_a, "t"_a, "d"_a, "relative_d"_a = false, "RANSAC outlier detection algorithm")
        .def("ransac", [](OpenMS::Math::RANSAC<OpenMS::Math::RansacModelLinear>& self, const std::vector<std::pair<double, double>>& pairs, const OpenMS::Math::RANSACParam& p) {
            return self.ransac(pairs, p);
        }, "pairs"_a, "p"_a, "RANSAC outlier detection with RANSACParam")
        ;

    // -----------------------------------------------------------------------
    // RANSACQuadratic
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Math::RANSAC<OpenMS::Math::RansacModelQuadratic>>(m, "RANSACQuadratic", "OpenMS class RANSACQuadratic")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Math::RANSAC<OpenMS::Math::RansacModelQuadratic>&>())
        .def("__copy__", [](const OpenMS::Math::RANSAC<OpenMS::Math::RansacModelQuadratic>& self) { return OpenMS::Math::RANSAC<OpenMS::Math::RansacModelQuadratic>(self); })
        .def("__deepcopy__", [](const OpenMS::Math::RANSAC<OpenMS::Math::RansacModelQuadratic>& self, nb::dict) { return OpenMS::Math::RANSAC<OpenMS::Math::RansacModelQuadratic>(self); }, "memo"_a)
        .def(nb::init<OpenMS::UInt64>())
        .def("ransac", [](OpenMS::Math::RANSAC<OpenMS::Math::RansacModelQuadratic>& self, const std::vector<std::pair<double, double>>& pairs, size_t n, size_t k, double t, size_t d, bool relative_d) {
            return self.ransac(pairs, n, k, t, d, relative_d);
        }, "pairs"_a, "n"_a, "k"_a, "t"_a, "d"_a, "relative_d"_a = false, "RANSAC outlier detection algorithm")
        .def("ransac", [](OpenMS::Math::RANSAC<OpenMS::Math::RansacModelQuadratic>& self, const std::vector<std::pair<double, double>>& pairs, const OpenMS::Math::RANSACParam& p) {
            return self.ransac(pairs, p);
        }, "pairs"_a, "p"_a, "RANSAC outlier detection with RANSACParam")
        ;


    // -----------------------------------------------------------------------
    // RANSACParam
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Math::RANSACParam>(m, "RANSACParam", "OpenMS class RANSACParam")
        .def(nb::init<>())
        .def(nb::init<size_t, size_t, double, size_t, bool>())
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
        .def("__copy__", [](const OpenMS::Math::RansacModelLinear& self) { return OpenMS::Math::RansacModelLinear(self); })
        .def("__deepcopy__", [](const OpenMS::Math::RansacModelLinear& self, nb::dict) { return OpenMS::Math::RansacModelLinear(self); }, "memo"_a)
        ;

    // -----------------------------------------------------------------------
    // RansacModelQuadratic
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Math::RansacModelQuadratic>(m, "RansacModelQuadratic", "Implementation of a quadratic RANSAC model fit")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Math::RansacModelQuadratic &>())
        .def("__copy__", [](const OpenMS::Math::RansacModelQuadratic& self) { return OpenMS::Math::RansacModelQuadratic(self); })
        .def("__deepcopy__", [](const OpenMS::Math::RansacModelQuadratic& self, nb::dict) { return OpenMS::Math::RansacModelQuadratic(self); }, "memo"_a)
        ;

}