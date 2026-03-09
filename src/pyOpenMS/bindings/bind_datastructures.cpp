// pyOpenMS nanobind bindings
// Domain: datastructures

#include "all_casters.h"
#include <OpenMS/COMPARISON/SpectraSTSimilarityScore.h>
#include <OpenMS/COMPARISON/SpectrumAlignmentScore.h>
#include <OpenMS/CONCEPT/LogConfigHandler.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/Adduct.h>
#include <OpenMS/DATASTRUCTURES/OSWData.h>
#include <OpenMS/DATASTRUCTURES/MassExplainer.h>
#include <OpenMS/DATASTRUCTURES/ChargePair.h>
#include <OpenMS/DATASTRUCTURES/CVMappings.h>
#include <OpenMS/DATASTRUCTURES/CVMappingTerm.h>
#include <OpenMS/DATASTRUCTURES/Compomer.h>
#include <OpenMS/DATASTRUCTURES/CVReference.h>
#include <OpenMS/DATASTRUCTURES/CVMappingRule.h>
#include <OpenMS/DATASTRUCTURES/CalibrationData.h>
#include <OpenMS/DATASTRUCTURES/ConvexHull2D.h>
#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>
#include <OpenMS/DATASTRUCTURES/Date.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/DATASTRUCTURES/DistanceMatrix.h>
#include <OpenMS/DATASTRUCTURES/IsotopeCluster.h>
#include <OpenMS/DATASTRUCTURES/LPWrapper.h>
#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/QTCluster.h>
#include <OpenMS/DATASTRUCTURES/StringView.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/MATH/MISC/BSpline2d.h>
#include <OpenMS/MATH/MISC/CubicSpline2d.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/MATH/STATISTICS/GaussFitter.h>
#include <OpenMS/MATH/STATISTICS/MultipleTesting.h>
#include <OpenMS/MATH/STATISTICS/RankData.h>
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
#include "binding_utils.h"
#include <nanobind/eigen/dense.h>
#include <Eigen/Core>

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(_pyopenms_datastructures, m) {
    m.doc() = "pyOpenMS datastructures bindings";


    // RankData::Method enum
    // -----------------------------------------------------------------------
    // Method
    // -----------------------------------------------------------------------
    nb::enum_<OpenMS::Math::RankData::Method>(m, "Method", "Method for resolving ties in rank computation", nb::is_arithmetic())
        .value("Average", OpenMS::Math::RankData::Method::Average)
        .value("Min", OpenMS::Math::RankData::Method::Min)
        .value("Max", OpenMS::Math::RankData::Method::Max)
        .value("Dense", OpenMS::Math::RankData::Method::Dense)
        .value("Ordinal", OpenMS::Math::RankData::Method::Ordinal)
        .export_values();


    // RankData::NaNPolicy enum
    // -----------------------------------------------------------------------
    // NaNPolicy
    // -----------------------------------------------------------------------
    nb::enum_<OpenMS::Math::RankData::NaNPolicy>(m, "NaNPolicy", "Policy for handling NaN values in rank computation", nb::is_arithmetic())
        .value("Propagate", OpenMS::Math::RankData::NaNPolicy::Propagate)
        .value("Omit", OpenMS::Math::RankData::NaNPolicy::Omit)
        .value("Raise", OpenMS::Math::RankData::NaNPolicy::Raise)
        .export_values();

    // -----------------------------------------------------------------------
    // Adduct
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Adduct>(m, "Adduct", "OpenMS class Adduct")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::Adduct& self) { return OpenMS::Adduct(self); })
        .def("__deepcopy__", [](const OpenMS::Adduct& self, nb::dict) { return OpenMS::Adduct(self); }, "memo"_a)
        .def(nb::init<int>())
        .def(nb::init<int, int, double, OpenMS::String, double, double, OpenMS::String>())
        .def("getCharge", [](const OpenMS::Adduct& self) { return self.getCharge(); })
        .def("setCharge", [](OpenMS::Adduct& self, const int& charge) { return self.setCharge(charge); }, "charge"_a)
        .def("getAmount", [](const OpenMS::Adduct& self) { return self.getAmount(); })
        .def("setAmount", [](OpenMS::Adduct& self, const int& amount) { return self.setAmount(amount); }, "amount"_a)
        .def("getSingleMass", [](const OpenMS::Adduct& self) { return self.getSingleMass(); })
        .def("setSingleMass", [](OpenMS::Adduct& self, const double& singleMass) { return self.setSingleMass(singleMass); }, "singleMass"_a)
        .def("getLogProb", [](const OpenMS::Adduct& self) { return self.getLogProb(); })
        .def("setLogProb", [](OpenMS::Adduct& self, const double& log_prob) { return self.setLogProb(log_prob); }, "log_prob"_a)
        .def("getFormula", [](const OpenMS::Adduct& self) { return self.getFormula(); })
        .def("setFormula", [](OpenMS::Adduct& self, const OpenMS::String& formula) { return self.setFormula(formula); }, "formula"_a)
        .def("getRTShift", [](const OpenMS::Adduct& self) { return self.getRTShift(); })
        .def("getLabel", [](const OpenMS::Adduct& self) { return self.getLabel(); })
        .def("__hash__", [](const OpenMS::Adduct& self) { return std::hash<OpenMS::Adduct>{}(self); })
        ;

    // -----------------------------------------------------------------------
    // BSpline2d
    // -----------------------------------------------------------------------
    auto bspline2d_class = nb::class_<OpenMS::BSpline2d>(m, "BSpline2d", "b spline interpolation")
        .def(nb::init<std::vector<double>, std::vector<double>, double, OpenMS::BSpline2d::BoundaryCondition, size_t>())
        .def("solve", [](OpenMS::BSpline2d& self, const std::vector<double>& y) { return self.solve(y); }, "y"_a, "Solve the spline curve for a new set of y values. Returns false if the solution fails")
        .def("eval", [](const OpenMS::BSpline2d& self, double x) { return self.eval(x); }, "x"_a, "Returns the evaluation of the smoothed curve at a particular x value. If current state is not ok(), returns zero")
        .def("derivative", [](const OpenMS::BSpline2d& self, double x) { return self.derivative(x); }, "x"_a, "Returns the first derivative of the spline curve at the given position x. Returns zero if the current state is not ok()")
        .def("ok", [](const OpenMS::BSpline2d& self) { return self.ok(); }, "Returns whether the spline fit was successful")
        .def_static("debug", [](bool enable) { return OpenMS::BSpline2d::debug(enable); }, "enable"_a, "Enable or disable debug messages from the B-spline library")
        ;
    // BoundaryCondition enum nested under BSpline2d
    nb::enum_<OpenMS::BSpline2d::BoundaryCondition>(bspline2d_class, "BoundaryCondition", nb::is_arithmetic(), "Boundary condition type for B-spline interpolation")
        .value("BC_ZERO_ENDPOINTS", OpenMS::BSpline2d::BoundaryCondition::BC_ZERO_ENDPOINTS)
        .value("BC_ZERO_FIRST", OpenMS::BSpline2d::BoundaryCondition::BC_ZERO_FIRST)
        .value("BC_ZERO_SECOND", OpenMS::BSpline2d::BoundaryCondition::BC_ZERO_SECOND)

        .export_values();

    // -----------------------------------------------------------------------
    // CVMappingTerm
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::CVMappingTerm>(m, "CVMappingTerm", "OpenMS class CVMappingTerm")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::CVMappingTerm &>())
        .def("__copy__", [](const OpenMS::CVMappingTerm& self) { return OpenMS::CVMappingTerm(self); })
        .def("__deepcopy__", [](const OpenMS::CVMappingTerm& self, nb::dict) { return OpenMS::CVMappingTerm(self); }, "memo"_a)
        .def("setAccession", [](OpenMS::CVMappingTerm& self, const OpenMS::String& accession) { return self.setAccession(accession); }, "accession"_a, "Sets the accession string of the term")
        .def("getAccession", [](const OpenMS::CVMappingTerm& self) { return self.getAccession(); }, "Returns the accession string of the term")
        .def("setUseTermName", [](OpenMS::CVMappingTerm& self, bool use_term_name) { return self.setUseTermName(use_term_name); }, "use_term_name"_a, "Sets whether the term name should be used, instead of the accession")
        .def("getUseTermName", [](const OpenMS::CVMappingTerm& self) { return self.getUseTermName(); }, "Returns whether the term name should be used, instead of the accession")
        .def("setUseTerm", [](OpenMS::CVMappingTerm& self, bool use_term) { return self.setUseTerm(use_term); }, "use_term"_a, "Sets whether the term itself can be used (or only its children)")
        .def("getUseTerm", [](const OpenMS::CVMappingTerm& self) { return self.getUseTerm(); }, "Returns true if the term can be used, false if only children are allowed")
        .def("setTermName", [](OpenMS::CVMappingTerm& self, const OpenMS::String& term_name) { return self.setTermName(term_name); }, "term_name"_a, "Sets the name of the term")
        .def("getTermName", [](const OpenMS::CVMappingTerm& self) { return self.getTermName(); }, "Returns the name of the term")
        .def("setIsRepeatable", [](OpenMS::CVMappingTerm& self, bool is_repeatable) { return self.setIsRepeatable(is_repeatable); }, "is_repeatable"_a, "Sets whether this term can be repeated")
        .def("getIsRepeatable", [](const OpenMS::CVMappingTerm& self) { return self.getIsRepeatable(); }, "Returns true if this term can be repeated, false otherwise")
        .def("setAllowChildren", [](OpenMS::CVMappingTerm& self, bool allow_children) { return self.setAllowChildren(allow_children); }, "allow_children"_a, "Sets whether children of this term are allowed")
        .def("getAllowChildren", [](const OpenMS::CVMappingTerm& self) { return self.getAllowChildren(); }, "Returns true if the children of this term are allowed to be used")
        .def("setCVIdentifierRef", [](OpenMS::CVMappingTerm& self, const OpenMS::String& cv_identifier_ref) { return self.setCVIdentifierRef(cv_identifier_ref); }, "cv_identifier_ref"_a, "Sets the CV identifier reference string, e.g. UO for unit obo")
        .def("getCVIdentifierRef", [](const OpenMS::CVMappingTerm& self) { return self.getCVIdentifierRef(); }, "Returns the CV identifier reference string")
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        ;

    // -----------------------------------------------------------------------
    // CVReference
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::CVReference>(m, "CVReference", "Controlled Vocabulary Reference")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::CVReference &>())
        .def("__copy__", [](const OpenMS::CVReference& self) { return OpenMS::CVReference(self); })
        .def("__deepcopy__", [](const OpenMS::CVReference& self, nb::dict) { return OpenMS::CVReference(self); }, "memo"_a)
        .def("setName", [](OpenMS::CVReference& self, const OpenMS::String& name) { return self.setName(name); }, "name"_a, "Sets the name of the CV reference")
        .def("getName", [](const OpenMS::CVReference& self) { return self.getName(); }, "Returns the name of the CV reference")
        .def("setIdentifier", [](OpenMS::CVReference& self, const OpenMS::String& identifier) { return self.setIdentifier(identifier); }, "identifier"_a, "Sets the CV identifier which is referenced")
        .def("getIdentifier", [](const OpenMS::CVReference& self) { return self.getIdentifier(); }, "Returns the CV identifier which is referenced")
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        ;

    // -----------------------------------------------------------------------
    // CalibrationData
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::CalibrationData>(m, "CalibrationData", "A helper class, holding all calibration points")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::CalibrationData& self) { return OpenMS::CalibrationData(self); })
        .def("__deepcopy__", [](const OpenMS::CalibrationData& self, nb::dict) { return OpenMS::CalibrationData(self); }, "memo"_a)
        .def("getMZ", [](const OpenMS::CalibrationData& self, size_t i) { return self.getMZ(i); }, "i"_a, "Retrieve the observed m/z of the i'th calibration point")
        .def("getRT", [](const OpenMS::CalibrationData& self, size_t i) { return self.getRT(i); }, "i"_a, "Retrieve the observed RT of the i'th calibration point")
        .def("getIntensity", [](const OpenMS::CalibrationData& self, size_t i) { return self.getIntensity(i); }, "i"_a, "Retrieve the intensity of the i'th calibration point")
        .def("size", [](const OpenMS::CalibrationData& self) { return self.size(); }, "Number of calibration points")
        .def("empty", [](const OpenMS::CalibrationData& self) { return self.empty(); }, "Returns `True` if there are no peaks")
        .def("clear", [](OpenMS::CalibrationData& self) { return self.clear(); }, "Remove all calibration points")
        .def("setUsePPM", [](OpenMS::CalibrationData& self, bool usePPM) { return self.setUsePPM(usePPM); }, "usePPM"_a)
        .def("usePPM", [](const OpenMS::CalibrationData& self) { return self.usePPM(); }, "Current error unit (ppm or Th)")
        .def("insertCalibrationPoint", [](OpenMS::CalibrationData& self, double rt, double mz_obs, float intensity, double mz_ref, double weight, int group) { return self.insertCalibrationPoint(rt, mz_obs, intensity, mz_ref, weight, group); }, "rt"_a, "mz_obs"_a, "intensity"_a, "mz_ref"_a, "weight"_a, "group"_a = -1)
        .def("getNrOfGroups", [](const OpenMS::CalibrationData& self) { return self.getNrOfGroups(); }, "Number of peak groups (can be 0)")
        .def("getError", [](const OpenMS::CalibrationData& self, size_t i) { return self.getError(i); }, "i"_a, "Retrieve the error for i'th calibrant in either ppm or Th (depending on usePPM())")
        .def("getRefMZ", [](const OpenMS::CalibrationData& self, size_t i) { return self.getRefMZ(i); }, "i"_a, "Retrieve the theoretical m/z of the i'th calibration point")
        .def("getWeight", [](const OpenMS::CalibrationData& self, size_t i) { return self.getWeight(i); }, "i"_a, "Retrieve the weight of the i'th calibration point")
        .def("getGroup", [](const OpenMS::CalibrationData& self, size_t i) { return self.getGroup(i); }, "i"_a, "Retrieve the group of the i'th calibration point")
        .def_static("getMetaValues", []() { return OpenMS::CalibrationData::getMetaValues(); })
        .def("median", [](const OpenMS::CalibrationData& self, double rt_left, double rt_right) { return self.median(rt_left, rt_right); }, "rt_left"_a, "rt_right"_a, "Compute the median in the given RT range for every peak group")
        .def("sortByRT", [](OpenMS::CalibrationData& self) { return self.sortByRT(); }, "Sort calibration points by RT, to allow for valid RT chunking")
        .def("__len__", [](OpenMS::CalibrationData& self) { return self.size(); })
        ;

    // -----------------------------------------------------------------------
    // ConvexHull2D
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ConvexHull2D>(m, "ConvexHull2D",
        R"doc(
A 2-dimensional hull representation in [counter]clockwise direction -
depending on axis labelling
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ConvexHull2D &>())
        .def("__copy__", [](const OpenMS::ConvexHull2D& self) { return OpenMS::ConvexHull2D(self); })
        .def("__deepcopy__", [](const OpenMS::ConvexHull2D& self, nb::dict) { return OpenMS::ConvexHull2D(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def("clear", [](OpenMS::ConvexHull2D& self) { return self.clear(); }, "Removes all points")
        .def("getHullPoints", [](const OpenMS::ConvexHull2D& self) -> const std::vector<OpenMS::DPosition<2>> & { return self.getHullPoints(); }, nb::rv_policy::reference_internal, "Accessor for the outer points")
        .def("setHullPoints", [](OpenMS::ConvexHull2D& self, const std::vector<OpenMS::DPosition<2>>& points) { return self.setHullPoints(points); }, "points"_a, "Accessor for the outer(!) points (no checking is performed if this is actually a convex hull)")
        .def("getBoundingBox", [](const OpenMS::ConvexHull2D& self) { return self.getBoundingBox(); }, "Returns the bounding box of the feature hull points")
        .def("addPoint", [](OpenMS::ConvexHull2D& self, const OpenMS::DPosition<2>& point) { return self.addPoint(point); }, "point"_a, "Adds a point to the hull if it is not already contained. Returns if the point was added. This will trigger recomputation of the outer hull points (thus points set with setHullPoints() will be lost)")
        .def("addPoints", [](OpenMS::ConvexHull2D& self, const std::vector<OpenMS::DPosition<2>>& points) { return self.addPoints(points); }, "points"_a, "Adds points to the hull if it is not already contained. This will trigger recomputation of the outer hull points (thus points set with setHullPoints() will be lost)")
        .def("compress", [](OpenMS::ConvexHull2D& self) { return self.compress(); }, "Allows to reduce the disk/memory footprint of a hull")
        .def("expandToBoundingBox", [](OpenMS::ConvexHull2D& self) { return self.expandToBoundingBox(); }, "Expand a convex hull to its bounding box.")
        .def("encloses", [](const OpenMS::ConvexHull2D& self, const OpenMS::DPosition<2>& point) { return self.encloses(point); }, "point"_a, "Returns if the `point` lies in the feature hull")
        .def("addPointXY", [](OpenMS::ConvexHull2D& self, double x, double y) {
            return self.addPoint(OpenMS::DPosition<2>(x, y));
        }, "x"_a, "y"_a, "Adds a point (x, y) to the hull")
        .def("enclosesXY", [](const OpenMS::ConvexHull2D& self, double x, double y) {
            return self.encloses(OpenMS::DPosition<2>(x, y));
        }, "x"_a, "y"_a, "Returns if the point (x, y) lies in the feature hull")
        .def("getBoundingBox2D", [](const OpenMS::ConvexHull2D& self) {
            return self.getBoundingBox();
        }, "Returns the bounding box of the feature hull points as a DBoundingBox2")
        .def("addPointsNPY", [](OpenMS::ConvexHull2D& self, nb::ndarray<double, nb::ndim<2>, nb::c_contig> arr) {
            if (arr.shape(1) != 2)
            {
                throw std::invalid_argument("Expected array with 2 columns (x, y)");
            }
            std::vector<OpenMS::DPosition<2>> points;
            points.reserve(arr.shape(0));
            const double* data = arr.data();
            for (size_t i = 0; i < arr.shape(0); ++i)
            {
                points.emplace_back(data[i * 2], data[i * 2 + 1]);
            }
            self.addPoints(points);
        }, "points"_a, "Adds points from a numpy array of shape (N, 2) to the hull")
        .def("getHullPointsNPY", [](const OpenMS::ConvexHull2D& self) {
            const auto& pts = self.getHullPoints();
            size_t n = pts.size();
            double* buf = new double[n * 2];
            for (size_t i = 0; i < n; ++i)
            {
                buf[i * 2] = pts[i][0];
                buf[i * 2 + 1] = pts[i][1];
            }
            size_t shape[2] = {n, 2};
            nb::capsule deleter(buf, [](void* p) noexcept { delete[] static_cast<double*>(p); });
            return nb::ndarray<nb::numpy, double, nb::ndim<2>>(buf, 2, shape, deleter);
        }, "Returns the hull points as a numpy array of shape (N, 2)")
        .def("setHullPointsNPY", [](OpenMS::ConvexHull2D& self, nb::ndarray<double, nb::ndim<2>, nb::c_contig> arr) {
            if (arr.shape(1) != 2)
            {
                throw std::invalid_argument("Expected array with 2 columns (x, y)");
            }
            std::vector<OpenMS::DPosition<2>> points;
            points.reserve(arr.shape(0));
            const double* data = arr.data();
            for (size_t i = 0; i < arr.shape(0); ++i)
            {
                points.emplace_back(data[i * 2], data[i * 2 + 1]);
            }
            self.setHullPoints(points);
        }, "points"_a, "Sets the hull points from a numpy array of shape (N, 2)")
        .def("__hash__", [](const OpenMS::ConvexHull2D& self) { return std::hash<OpenMS::ConvexHull2D>{}(self); })
        .def("__repr__", [](const OpenMS::ConvexHull2D& self) {
            return "ConvexHull2D(num_points=" + std::to_string(self.getHullPoints().size()) + ")";
        })
        .def("__str__", [](const OpenMS::ConvexHull2D& self) { return nb::cast(self).attr("__repr__")(); })
        ;

    // -----------------------------------------------------------------------
    // CubicSpline2d
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::CubicSpline2d>(m, "CubicSpline2d",
        R"doc(
cubic spline interpolation as described in R.L. Burden, J.D. Faires,
Numerical Analysis, 4th ed. PWS-Kent, 1989, ISBN 0-53491-585-X, pp.
126-131
)doc")
        .def("__copy__", [](const OpenMS::CubicSpline2d& self) { return OpenMS::CubicSpline2d(self); })
        .def("__deepcopy__", [](const OpenMS::CubicSpline2d& self, nb::dict) { return OpenMS::CubicSpline2d(self); }, "memo"_a)
        .def(nb::init<std::vector<double>, std::vector<double>>())
        .def(nb::init<std::map<double, double>>())
        .def("eval", [](const OpenMS::CubicSpline2d& self, double x) { return self.eval(x); }, "x"_a, "Evaluates the cubic spline")
        .def("derivatives", [](const OpenMS::CubicSpline2d& self, double x, unsigned int order) { return self.derivatives(x, order); }, "x"_a, "order"_a, "Returns first, second or third derivative of cubic spline")
        ;

    // -----------------------------------------------------------------------
    // DBoundingBox2
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::DBoundingBox<2> >(m, "DBoundingBox2", "OpenMS class DBoundingBox2")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::DBoundingBox<2>  &>())
        .def("__copy__", [](const OpenMS::DBoundingBox<2>& self) { return OpenMS::DBoundingBox<2>(self); })
        .def("__deepcopy__", [](const OpenMS::DBoundingBox<2>& self, nb::dict) { return OpenMS::DBoundingBox<2>(self); }, "memo"_a)

        .def("minPosition", [](const OpenMS::DBoundingBox<2>& self) {
            auto p = self.minPosition();
            nb::list result;
            result.append(p[0]);
            result.append(p[1]);
            return result;
        })

        .def("maxPosition", [](const OpenMS::DBoundingBox<2>& self) {
            auto p = self.maxPosition();
            nb::list result;
            result.append(p[0]);
            result.append(p[1]);
            return result;
        })
        ;

    // -----------------------------------------------------------------------
    // Date
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Date>(m, "Date", "Date class for date handling")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Date &>())
        .def("__copy__", [](const OpenMS::Date& self) { return OpenMS::Date(self); })
        .def("__deepcopy__", [](const OpenMS::Date& self, nb::dict) { return OpenMS::Date(self); }, "memo"_a)
        .def("set", static_cast<void (OpenMS::Date::*)(const OpenMS::String&)>(&OpenMS::Date::set), "date"_a, "Sets date from a string (mm/dd/yyyy, dd.mm.yyyy, or yyyy-mm-dd)")
        .def("set", static_cast<void (OpenMS::Date::*)(OpenMS::UInt, OpenMS::UInt, OpenMS::UInt)>(&OpenMS::Date::set), "month"_a, "day"_a, "year"_a, "Sets date from three integers")
        .def_static("today", []() { return OpenMS::Date::today(); }, "Returns the current date")
        .def("get", static_cast<OpenMS::String (OpenMS::Date::*)() const>(&OpenMS::Date::get), "Returns the date as string in iso/ansi format: yyyy-mm-dd")
        .def("clear", [](OpenMS::Date& self) { self.clear(); }, "Sets the undefined date: 00/00/0000")
        .def("isValid", [](const OpenMS::Date& self) { return self.isValid(); }, "Returns if the date is valid")
        .def("isNull", [](const OpenMS::Date& self) { return self.isNull(); }, "Returns if the date is null")
        .def("year", [](const OpenMS::Date& self) { return self.year(); }, "Returns the year")
        .def("month", [](const OpenMS::Date& self) { return self.month(); }, "Returns the month")
        .def("day", [](const OpenMS::Date& self) { return self.day(); }, "Returns the day")
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def(nb::self < nb::self)
        ;

    // -----------------------------------------------------------------------
    // DateTime
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::DateTime>(m, "DateTime", "DateTime Class")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::DateTime &>())
        .def("__copy__", [](const OpenMS::DateTime& self) { return OpenMS::DateTime(self); })
        .def("__deepcopy__", [](const OpenMS::DateTime& self, nb::dict) { return OpenMS::DateTime(self); }, "memo"_a)
        .def("setDate", [](OpenMS::DateTime& self, const OpenMS::String& date) { return self.setDate(date); }, "date"_a)
        .def("setTime", [](OpenMS::DateTime& self, const OpenMS::String& date) { return self.setTime(date); }, "date"_a)
        .def("setDate", [](OpenMS::DateTime& self, unsigned int month, unsigned int day, unsigned int year) { return self.setDate(month, day, year); }, "month"_a, "day"_a, "year"_a)
        .def("setTime", [](OpenMS::DateTime& self, unsigned int hour, unsigned int minute, unsigned int second) { return self.setTime(hour, minute, second); }, "hour"_a, "minute"_a, "second"_a)
        .def("set", [](OpenMS::DateTime& self, unsigned int month, unsigned int day, unsigned int year, unsigned int hour, unsigned int minute, unsigned int second) { return self.set(month, day, year, hour, minute, second); }, "month"_a, "day"_a, "year"_a, "hour"_a, "minute"_a, "second"_a,
            R"doc(
@brief Sets date and time
The following formats are supported:
- MM/dd/yyyy hh:mm:ss
- dd.MM.yyyy hh:mm:ss
- yyyy-MM-dd hh:mm:ss
- yyyy-MM-ddThh:mm:ss (ISO 8601 format)
- yyyy-MM-ddZ (ISO 8601 format)
- yyyy-MM-dd+hh:mm (ISO 8601 format)
)doc")
        .def("getDateAndTime", [](const OpenMS::DateTime& self) { unsigned int month, day, year, hour, minute, second; self.get(month, day, year, hour, minute, second); return nb::make_tuple(month, day, year, hour, minute, second); },
            R"doc(
@brief Returns date and time components as a tuple
@return Tuple of (month, day, year, hour, minute, second) as integers
)doc")
        .def("getDateComponents", [](const OpenMS::DateTime& self) { unsigned int month, day, year; self.getDate(month, day, year); return nb::make_tuple(month, day, year); }, "Returns (month, day, year) as a tuple of integers")
        .def("getDate", [](const OpenMS::DateTime& self) { return self.getDate(); })
        .def("getTimeComponents", [](const OpenMS::DateTime& self) { unsigned int hour, minute, second; self.getTime(hour, minute, second); return nb::make_tuple(hour, minute, second); }, "Returns (hour, minute, second) as a tuple of integers")
        .def("getTime", [](const OpenMS::DateTime& self) { return self.getTime(); })
        .def_static("now", []() { return OpenMS::DateTime::now(); })
        .def("clear", [](OpenMS::DateTime& self) { return self.clear(); })
        .def("get", [](const OpenMS::DateTime& self) { return self.get(); },
            R"doc(
@brief Returns a string representation of the date and time
The format of the string will be yyyy-MM-dd hh:mm:ss
)doc")
        .def("set", [](OpenMS::DateTime& self, const OpenMS::String& date) { return self.set(date); }, "date"_a,
            R"doc(
@brief Sets date and time
The following formats are supported:
- MM/dd/yyyy hh:mm:ss
- dd.MM.yyyy hh:mm:ss
- yyyy-MM-dd hh:mm:ss
- yyyy-MM-ddThh:mm:ss (ISO 8601 format)
- yyyy-MM-ddZ (ISO 8601 format)
- yyyy-MM-dd+hh:mm (ISO 8601 format)
)doc")
        .def("__hash__", [](const OpenMS::DateTime& self) { return std::hash<OpenMS::DateTime>{}(self); })
        .def("isNull", [](const OpenMS::DateTime& self) { return self.isNull(); }, "Returns true if the DateTime is null (default constructed)")
        .def("isValid", [](const OpenMS::DateTime& self) { return self.isValid(); }, "Returns true if the DateTime is valid")
        ;

    // -----------------------------------------------------------------------
    // DistanceMatrix
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::DistanceMatrix<float>>(m, "DistanceMatrix", "OpenMS class DistanceMatrix")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::DistanceMatrix<float>&>())
        .def("__copy__", [](const OpenMS::DistanceMatrix<float>& self) { return OpenMS::DistanceMatrix<float>(self); })
        .def("__deepcopy__", [](const OpenMS::DistanceMatrix<float>& self, nb::dict) { return OpenMS::DistanceMatrix<float>(self); }, "memo"_a)
        .def(nb::init<size_t, float>())
        .def("clear", [](OpenMS::DistanceMatrix<float>& self) { self.clear(); }, "Clears the matrix")
        .def("dimensionsize", [](const OpenMS::DistanceMatrix<float>& self) { return self.dimensionsize(); }, "Returns the number of rows/columns")
        .def("getMinElementCoordinates", [](const OpenMS::DistanceMatrix<float>& self) { return self.getMinElementCoordinates(); }, "Returns coordinates of minimum element")
        .def("getValue", [](const OpenMS::DistanceMatrix<float>& self, size_t i, size_t j) { return self.getValue(i, j); }, "i"_a, "j"_a, "Returns the value at (i,j)")
        .def("reduce", [](OpenMS::DistanceMatrix<float>& self, size_t j) { self.reduce(j); }, "j"_a, "Reduces the matrix by removing the j-th row and column")
        .def("resize", [](OpenMS::DistanceMatrix<float>& self, size_t dimensionsize, float value) { self.resize(dimensionsize, value); }, "dimensionsize"_a, "value"_a = 0.0f, "Resizes the matrix")
        .def("setValue", [](OpenMS::DistanceMatrix<float>& self, size_t i, size_t j, float value) { self.setValue(i, j, value); }, "i"_a, "j"_a, "value"_a, "Sets the value at (i,j)")
        .def("setValueQuick", [](OpenMS::DistanceMatrix<float>& self, size_t i, size_t j, float value) { self.setValueQuick(i, j, value); }, "i"_a, "j"_a, "value"_a, "Sets the value at (i,j) without checking bounds")
        .def("updateMinElement", [](OpenMS::DistanceMatrix<float>& self) { self.updateMinElement(); }, "Updates the minimum element")
        ;


    // -----------------------------------------------------------------------
    // GaussFitter
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Math::GaussFitter>(m, "GaussFitter", "Implements a fitter for Gaussian functions")
        .def(nb::init<>())
        .def("setInitialParameters", [](OpenMS::Math::GaussFitter& self, const OpenMS::Math::GaussFitter::GaussFitResult& result) { return self.setInitialParameters(result); }, "result"_a, "Sets the initial parameters used by the fit method as initial guess for the Gaussian")
        .def("fit", [](const OpenMS::Math::GaussFitter& self, std::vector<OpenMS::DPosition<2>>& points) { return self.fit(points); }, "points"_a, "Fits a Gaussian distribution to the given data points")
        ;

    // -----------------------------------------------------------------------
    // GaussFitResult
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Math::GaussFitter::GaussFitResult>(m, "GaussFitResult", "Result of a Gaussian fit")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::Math::GaussFitter::GaussFitResult& self) { return OpenMS::Math::GaussFitter::GaussFitResult(self); })
        .def("__deepcopy__", [](const OpenMS::Math::GaussFitter::GaussFitResult& self, nb::dict) { return OpenMS::Math::GaussFitter::GaussFitResult(self); }, "memo"_a)
        .def(nb::init<double, double, double>())
        .def_rw("A", &OpenMS::Math::GaussFitter::GaussFitResult::A)
        .def_rw("x0", &OpenMS::Math::GaussFitter::GaussFitResult::x0)
        .def_rw("sigma", &OpenMS::Math::GaussFitter::GaussFitResult::sigma)
        .def("eval", [](const OpenMS::Math::GaussFitter::GaussFitResult& self, double x) { return self.eval(x); }, "x"_a, "Evaluates the Gaussian model at the specified point")
        ;

    // -----------------------------------------------------------------------
    // ChargedIndexSet (IsotopeCluster::ChargedIndexSet)
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::IsotopeCluster::ChargedIndexSet>(m, "ChargedIndexSet",
        "Index set with associated charge estimate for isotope clusters")
        .def(nb::init<>())
        .def_rw("charge", &OpenMS::IsotopeCluster::ChargedIndexSet::charge)
        ;

    // -----------------------------------------------------------------------
    // IsotopeCluster
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::IsotopeCluster>(m, "IsotopeCluster", "OpenMS class IsotopeCluster")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::IsotopeCluster& self) { return OpenMS::IsotopeCluster(self); })
        .def("__deepcopy__", [](const OpenMS::IsotopeCluster& self, nb::dict) { return OpenMS::IsotopeCluster(self); }, "memo"_a)
        .def_rw("peaks", &OpenMS::IsotopeCluster::peaks)
        .def_rw("scans", &OpenMS::IsotopeCluster::scans)
        ;

    // -----------------------------------------------------------------------
    // SolverParam (LPWrapper::SolverParam)
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::LPWrapper::SolverParam>(m, "SolverParam",
        "Parameters for LP/MIP solver configuration")
        .def(nb::init<>())
        .def_rw("message_level", &OpenMS::LPWrapper::SolverParam::message_level)
        .def_rw("branching_tech", &OpenMS::LPWrapper::SolverParam::branching_tech)
        .def_rw("backtrack_tech", &OpenMS::LPWrapper::SolverParam::backtrack_tech)
        .def_rw("preprocessing_tech", &OpenMS::LPWrapper::SolverParam::preprocessing_tech)
        .def_rw("enable_feas_pump_heuristic", &OpenMS::LPWrapper::SolverParam::enable_feas_pump_heuristic)
        .def_rw("enable_gmi_cuts", &OpenMS::LPWrapper::SolverParam::enable_gmi_cuts)
        .def_rw("enable_mir_cuts", &OpenMS::LPWrapper::SolverParam::enable_mir_cuts)
        .def_rw("enable_cov_cuts", &OpenMS::LPWrapper::SolverParam::enable_cov_cuts)
        .def_rw("enable_clq_cuts", &OpenMS::LPWrapper::SolverParam::enable_clq_cuts)
        .def_rw("mip_gap", &OpenMS::LPWrapper::SolverParam::mip_gap)
        .def_rw("time_limit", &OpenMS::LPWrapper::SolverParam::time_limit)
        .def_rw("output_freq", &OpenMS::LPWrapper::SolverParam::output_freq)
        .def_rw("output_delay", &OpenMS::LPWrapper::SolverParam::output_delay)
        .def_rw("enable_presolve", &OpenMS::LPWrapper::SolverParam::enable_presolve)
        .def_rw("enable_binarization", &OpenMS::LPWrapper::SolverParam::enable_binarization)
        ;

    // -----------------------------------------------------------------------
    // LPWrapper
    // -----------------------------------------------------------------------
    auto lpwrapper_class = nb::class_<OpenMS::LPWrapper>(m, "LPWrapper", "A wrapper class for linear programming (LP) solvers")
        .def(nb::init<>())
        .def("addRow", [](OpenMS::LPWrapper& self, const std::vector<int>& row_indices, const std::vector<double>& row_values, const OpenMS::String& name) { return self.addRow(row_indices, row_values, name); }, "row_indices"_a, "row_values"_a, "name"_a, "Adds a row to the LP matrix, returns index")
        .def("addColumn", [](OpenMS::LPWrapper& self) { return self.addColumn(); }, "Adds an empty column to the LP matrix, returns index")
        .def("addColumn", [](OpenMS::LPWrapper& self, const std::vector<int>& column_indices, const std::vector<double>& column_values, const OpenMS::String& name) { return self.addColumn(column_indices, column_values, name); }, "column_indices"_a, "column_values"_a, "name"_a, "Adds a column to the LP matrix, returns index")
        .def("addRow", [](OpenMS::LPWrapper& self, const std::vector<int>& row_indices, const std::vector<double>& row_values, const OpenMS::String& name, double lower_bound, double upper_bound, OpenMS::LPWrapper::Type type) { return self.addRow(row_indices, row_values, name, lower_bound, upper_bound, type); }, "row_indices"_a, "row_values"_a, "name"_a, "lower_bound"_a, "upper_bound"_a, "type"_a, "Adds a row with boundaries to the LP matrix, returns index")
        .def("addColumn", [](OpenMS::LPWrapper& self, const std::vector<int>& column_indices, const std::vector<double>& column_values, const OpenMS::String& name, double lower_bound, double upper_bound, OpenMS::LPWrapper::Type type) { return self.addColumn(column_indices, column_values, name, lower_bound, upper_bound, type); }, "column_indices"_a, "column_values"_a, "name"_a, "lower_bound"_a, "upper_bound"_a, "type"_a, "Adds a column with boundaries to the LP matrix, returns index")
        .def("deleteRow", [](OpenMS::LPWrapper& self, int index) { return self.deleteRow(index); }, "index"_a, "Delete index-th row")
        .def("setColumnName", [](OpenMS::LPWrapper& self, int index, const OpenMS::String& name) { return self.setColumnName(index, name); }, "index"_a, "name"_a, "Sets name of the index-th column")
        .def("getColumnName", [](OpenMS::LPWrapper& self, int index) { return self.getColumnName(index); }, "index"_a, "Returns name of the index-th column")
        .def("getRowName", [](OpenMS::LPWrapper& self, int index) { return self.getRowName(index); }, "index"_a, "Sets name of the index-th row")
        .def("getRowIndex", [](OpenMS::LPWrapper& self, const OpenMS::String& name) { return self.getRowIndex(name); }, "name"_a, "Returns index of the row with name")
        .def("getColumnIndex", [](OpenMS::LPWrapper& self, const OpenMS::String& name) { return self.getColumnIndex(name); }, "name"_a, "Returns index of the column with name")
        .def("getColumnUpperBound", [](OpenMS::LPWrapper& self, int index) { return self.getColumnUpperBound(index); }, "index"_a, "Returns column's upper bound")
        .def("getColumnLowerBound", [](OpenMS::LPWrapper& self, int index) { return self.getColumnLowerBound(index); }, "index"_a, "Returns column's lower bound")
        .def("getRowUpperBound", [](OpenMS::LPWrapper& self, int index) { return self.getRowUpperBound(index); }, "index"_a, "Returns row's upper bound")
        .def("getRowLowerBound", [](OpenMS::LPWrapper& self, int index) { return self.getRowLowerBound(index); }, "index"_a, "Returns row's lower bound")
        .def("setRowName", [](OpenMS::LPWrapper& self, int index, const OpenMS::String& name) { return self.setRowName(index, name); }, "index"_a, "name"_a, "Sets name of the index-th row")
        .def("setColumnBounds", [](OpenMS::LPWrapper& self, int index, double lower_bound, double upper_bound, OpenMS::LPWrapper::Type type) { return self.setColumnBounds(index, lower_bound, upper_bound, type); }, "index"_a, "lower_bound"_a, "upper_bound"_a, "type"_a, "Sets column bounds")
        .def("setRowBounds", [](OpenMS::LPWrapper& self, int index, double lower_bound, double upper_bound, OpenMS::LPWrapper::Type type) { return self.setRowBounds(index, lower_bound, upper_bound, type); }, "index"_a, "lower_bound"_a, "upper_bound"_a, "type"_a, "Sets row bounds")
        .def("setColumnType", [](OpenMS::LPWrapper& self, int index, OpenMS::LPWrapper::VariableType type) { return self.setColumnType(index, type); }, "index"_a, "type"_a, "Sets column/variable type.")
        .def("getColumnType", [](OpenMS::LPWrapper& self, int index) { return self.getColumnType(index); }, "index"_a, "Returns column/variable type.")
        .def("setObjective", [](OpenMS::LPWrapper& self, int index, double obj_value) { return self.setObjective(index, obj_value); }, "index"_a, "obj_value"_a, "Sets objective value for column with index")
        .def("getObjective", [](OpenMS::LPWrapper& self, int index) { return self.getObjective(index); }, "index"_a, "Returns objective value for column with index")
        .def("setObjectiveSense", [](OpenMS::LPWrapper& self, OpenMS::LPWrapper::Sense sense) { return self.setObjectiveSense(sense); }, "sense"_a, "Sets objective direction")
        .def("getObjectiveSense", [](OpenMS::LPWrapper& self) { return self.getObjectiveSense(); }, "Returns objective sense")
        .def("getNumberOfColumns", [](OpenMS::LPWrapper& self) { return self.getNumberOfColumns(); }, "Returns number of columns")
        .def("getNumberOfRows", [](OpenMS::LPWrapper& self) { return self.getNumberOfRows(); }, "Returns number of rows")
        .def("setElement", [](OpenMS::LPWrapper& self, int row_index, int column_index, double value) { return self.setElement(row_index, column_index, value); }, "row_index"_a, "column_index"_a, "value"_a, "Sets the element")
        .def("getElement", [](OpenMS::LPWrapper& self, int row_index, int column_index) { return self.getElement(row_index, column_index); }, "row_index"_a, "column_index"_a, "Returns the element")
        .def("readProblem", [](OpenMS::LPWrapper& self, const OpenMS::String& filename, const OpenMS::String& format) { return self.readProblem(filename, format); }, "filename"_a, "format"_a)
        .def("writeProblem", [](const OpenMS::LPWrapper& self, const OpenMS::String& filename, OpenMS::LPWrapper::WriteFormat format) { return self.writeProblem(filename, format); }, "filename"_a, "format"_a,
            R"doc(
Write LP formulation to a file
:param filename: Output filename, if the filename ends with '.gz' it will be compressed
:param format: MPS-format is supported by GLPK and COIN-OR; LP and GLPK-formats only by GLPK
)doc")
        .def("solve", [](OpenMS::LPWrapper& self, OpenMS::LPWrapper::SolverParam& solver_param, size_t verbose_level) { return self.solve(solver_param, verbose_level); }, "solver_param"_a, "verbose_level"_a = 0,
            R"doc(
Solve problems, parameters like enabled heuristics can be given via solver_param
The verbose level (0,1,2) determines if the solver prints status messages and internals
:param solver_param: Parameters of the solver introduced by SolverParam
:param verbose_level: Sets verbose level
:return: solver dependent
)doc")
        .def("getStatus", [](OpenMS::LPWrapper& self) { return self.getStatus(); },
            R"doc(
Returns solution status
:return: status: 1 - undefined, 2 - integer optimal, 3- integer feasible (no optimality proven), 4- no integer feasible solution
)doc")
        .def("getObjectiveValue", [](OpenMS::LPWrapper& self) { return self.getObjectiveValue(); },
            R"doc(
Returns the objective value of the solution
:return: The optimal objective value after solving
)doc")
        .def("getColumnValue", [](OpenMS::LPWrapper& self, int index) { return self.getColumnValue(index); }, "index"_a)
        .def("getNumberOfNonZeroEntriesInRow", [](OpenMS::LPWrapper& self, int idx) { return self.getNumberOfNonZeroEntriesInRow(idx); }, "idx"_a)
        .def("getMatrixRow", [](OpenMS::LPWrapper& self, int idx) { std::vector<int> indexes; self.getMatrixRow(idx, indexes); return indexes; }, "idx"_a)
        .def("getSolver", [](const OpenMS::LPWrapper& self) { return self.getSolver(); }, "Returns currently active solver")
        ;
    // LPWrapper_Type enum nested under LPWrapper
    nb::enum_<OpenMS::LPWrapper::Type>(lpwrapper_class, "LPWrapper_Type", nb::is_arithmetic())
        .value("UNBOUNDED", OpenMS::LPWrapper::Type::UNBOUNDED)
        .value("LOWER_BOUND_ONLY", OpenMS::LPWrapper::Type::LOWER_BOUND_ONLY)
        .value("UPPER_BOUND_ONLY", OpenMS::LPWrapper::Type::UPPER_BOUND_ONLY)
        .value("DOUBLE_BOUNDED", OpenMS::LPWrapper::Type::DOUBLE_BOUNDED)
        .value("FIXED", OpenMS::LPWrapper::Type::FIXED)
        .export_values();
    // VariableType enum nested under LPWrapper
    nb::enum_<OpenMS::LPWrapper::VariableType>(lpwrapper_class, "VariableType", nb::is_arithmetic())
        .value("CONTINUOUS", OpenMS::LPWrapper::VariableType::CONTINUOUS)
        .value("INTEGER", OpenMS::LPWrapper::VariableType::INTEGER)
        .value("BINARY", OpenMS::LPWrapper::VariableType::BINARY)
        .export_values();
    // Sense enum nested under LPWrapper
    nb::enum_<OpenMS::LPWrapper::Sense>(lpwrapper_class, "Sense", nb::is_arithmetic())
        .value("MIN", OpenMS::LPWrapper::Sense::MIN)
        .value("MAX", OpenMS::LPWrapper::Sense::MAX)
        .export_values();
    // WriteFormat enum nested under LPWrapper
    nb::enum_<OpenMS::LPWrapper::WriteFormat>(lpwrapper_class, "WriteFormat", nb::is_arithmetic())
        .value("FORMAT_LP", OpenMS::LPWrapper::WriteFormat::FORMAT_LP)
        .value("FORMAT_MPS", OpenMS::LPWrapper::WriteFormat::FORMAT_MPS)
        .value("FORMAT_GLPK", OpenMS::LPWrapper::WriteFormat::FORMAT_GLPK)
        .export_values();
    // SOLVER enum nested under LPWrapper
    nb::enum_<OpenMS::LPWrapper::SOLVER>(lpwrapper_class, "SOLVER", nb::is_arithmetic())
        .value("SOLVER_GLPK", OpenMS::LPWrapper::SOLVER::SOLVER_GLPK)
        .export_values();
    // SolverStatus enum nested under LPWrapper
    nb::enum_<OpenMS::LPWrapper::SolverStatus>(lpwrapper_class, "SolverStatus", nb::is_arithmetic())
        .value("UNDEFINED", OpenMS::LPWrapper::SolverStatus::UNDEFINED)
        .value("OPTIMAL", OpenMS::LPWrapper::SolverStatus::OPTIMAL)
        .value("FEASIBLE", OpenMS::LPWrapper::SolverStatus::FEASIBLE)
        .value("NO_FEASIBLE_SOL", OpenMS::LPWrapper::SolverStatus::NO_FEASIBLE_SOL)
        .export_values();

    // -----------------------------------------------------------------------
    // LogConfigHandler
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::LogConfigHandler>(m, "LogConfigHandler",
        R"doc(
The LogConfigHandler provides the functionality to configure the
internal logging of OpenMS algorithms that use the global instances of
LogStream
)doc")
        .def("parse", [](OpenMS::LogConfigHandler& self, const std::vector<OpenMS::String>& setting) { return self.parse(setting); }, "setting"_a)
        .def("configure", [](OpenMS::LogConfigHandler& self, const OpenMS::Param& param) { return self.configure(param); }, "param"_a,
            R"doc(
Translates the given list of parameter settings into a LogStream configuration
Translates the given list of parameter settings into a LogStream configuration.
Usually this list stems from a command line call.
Each element in the stringlist should follow this naming convention
<LOG_NAME> <ACTION> <PARAMETER>
with
- LOG_NAME: DEBUG,INFO,WARNING,ERROR,FATAL_ERROR
- ACTION: add,remove,clear
- PARAMETER: for 'add'/'remove' it is the stream name (cout, cerr or a filename), 'clear' does not require any further parameter
Example:
`DEBUG add debug.log`
This function will **not** apply to settings to the log handlers. Use configure() for that.
:param setting: StringList containing the configuration options
:raises ParseError: In case of an invalid configuration.
:return: Param object containing all settings, that can be applied using the LogConfigHandler.configure() method
)doc")
        .def("setLogLevel", [](OpenMS::LogConfigHandler& self, const OpenMS::String& log_level) { return self.setLogLevel(log_level); }, "log_level"_a,
            R"doc(
Applies the given parameters (@p param) to the current configuration
<LOG_NAME> <ACTION> <PARAMETER> <STREAMTYPE>
LOG_NAME: DEBUG, INFO, WARNING, ERROR, FATAL_ERROR
ACTION: add, remove, clear
PARAMETER: for 'add'/'remove' it is the stream name ('cout', 'cerr' or a filename), 'clear' does not require any further parameter
STREAMTYPE: FILE, STRING (for a StringStream, which you can grab by this name using getStream() )
You cannot specify a file named "cout" or "cerr" even if you specify streamtype 'FILE' - the handler will mistake this for the
internal streams, but you can use "./cout" to print to a file named cout.
A classical configuration would contain a list of settings e.g.
`DEBUG add debug.log FILE`
`INFO remove cout FILE` (FILE will be ignored)
`INFO add string_stream1 STRING`
:raises ElementNotFound: If the LogStream (first argument) does not exist.
:raises FileNotWritable: If a file (or stream) should be opened as log file (or stream) that is not accessible.
:raises IllegalArgument: If a stream should be registered, that was already registered with a different type.
)doc")
        .def_static("getInstance", []() -> OpenMS::LogConfigHandler* { return OpenMS::LogConfigHandler::getInstance(); }, nb::rv_policy::reference, "Returns the singleton instance")
        ;

    // -----------------------------------------------------------------------
    // MatrixDouble
    // -----------------------------------------------------------------------
    // -----------------------------------------------------------------------
    // MatrixDouble (Optimized Zero-Copy)
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Matrix<double>>(m, "MatrixDouble", "OpenMS class MatrixDouble")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Matrix<double>&>())
        .def("__copy__", [](const OpenMS::Matrix<double>& self) { return OpenMS::Matrix<double>(self); })
        .def("__deepcopy__", [](const OpenMS::Matrix<double>& self, nb::dict) { return OpenMS::Matrix<double>(self); }, "memo"_a)
        .def(nb::init<size_t, size_t, double>())
        .def("getValue", [](OpenMS::Matrix<double>& self, size_t i, size_t j) { return self.getValue(i, j); }, "i"_a, "j"_a)
        .def("setValue", [](OpenMS::Matrix<double>& self, size_t i, size_t j, double value) { self.setValue(i, j, value); }, "i"_a, "j"_a, "value"_a)
        .def("rows", [](OpenMS::Matrix<double>& self) { return self.rows(); })
        .def("cols", [](OpenMS::Matrix<double>& self) { return self.cols(); })
        .def("size", [](OpenMS::Matrix<double>& self) { return self.size(); })
        .def("resize", [](OpenMS::Matrix<double>& self, size_t rows, size_t cols) { self.resize(rows, cols); }, "rows"_a, "cols"_a)
        .def("__len__", [](OpenMS::Matrix<double>& self) { return self.size(); })

        .def("get_matrix_view", [](nb::object self_obj) -> nb::object {
            auto& self = nb::cast<OpenMS::Matrix<double>&>(self_obj);
            size_t shape[2] = {self.rows(), self.cols()};
            int64_t strides[2] = {1, static_cast<int64_t>(self.rows())};
            double* data_ptr = (self.rows() == 0 || self.cols() == 0) ? nullptr : self.data();
            return nb::ndarray<nb::numpy, double, nb::ndim<2>>(
                data_ptr, 2, shape, self_obj, strides
            ).cast();
        }, "Returns a zero-copy numpy view of the matrix data (F-contiguous). Modifications affect the C++ object. Returns empty array if empty.")

        .def("get_matrix", [](const OpenMS::Matrix<double>& self) {
            Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>
                eigen_map(self.data(), self.rows(), self.cols());
            return nb::cast(eigen_map);
        }, "Returns a copy of the matrix as a numpy array")

        .def("set_matrix", [](OpenMS::Matrix<double>& self, nb::ndarray<double, nb::ndim<2>> arr) {
            self.resize(arr.shape(0), arr.shape(1));

            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>
                dest_map(self.data(), self.rows(), self.cols());

            if (arr.stride(1) == 1) { // C-contiguous
                Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
                    src_map(arr.data(), arr.shape(0), arr.shape(1));
                dest_map = src_map;
            } else { // F-contiguous
                Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>
                    src_map(arr.data(), arr.shape(0), arr.shape(1));
                dest_map = src_map;
            }
        }, "matrix"_a, "Fast copy from a numpy array to the OpenMS Matrix.")

        .def_static("fromNdArray", [](nb::ndarray<double, nb::ndim<2>> arr) {
            OpenMS::Matrix<double> mat(arr.shape(0), arr.shape(1), 0.0);

            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>
                dest_map(mat.data(), mat.rows(), mat.cols());

            if (arr.stride(1) == 1) { // C-contiguous
                Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
                    src_map(arr.data(), arr.shape(0), arr.shape(1));
                dest_map = src_map;
            } else { // F-contiguous
                Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>
                    src_map(arr.data(), arr.shape(0), arr.shape(1));
                dest_map = src_map;
            }
            return mat;
        }, "array"_a, "Creates a MatrixDouble from a 2D numpy array optimally.")

        .def("__str__", [](const OpenMS::Matrix<double>& self) {
            std::ostringstream oss;
            for (size_t i = 0; i < self.rows(); ++i) {
                for (size_t j = 0; j < self.cols(); ++j) {
                    if (j > 0) oss << " ";
                    oss << self.getValue(i, j);
                }
                oss << "\n";
            }
            return oss.str();
        })
        ;


    // -----------------------------------------------------------------------
    // MultipleTesting
    // -----------------------------------------------------------------------
    auto multipletesting_class = nb::class_<OpenMS::Math::MultipleTesting>(m, "MultipleTesting",
        R"doc(
Statistical functions for multiple testing correction.
Provides FDR estimation, q-value calculation, and local FDR computation
based on the Storey-Tibshirani method.
Example:
>>> import pyopenms
>>> p_values = [0.001, 0.01, 0.05, 0.1, 0.5]
>>> # Estimate pi0
>>> pi0_result = pyopenms.MultipleTesting.pi0Est(p_values, [])
>>> # Calculate q-values
>>> q_values = pyopenms.MultipleTesting.qValue(p_values, pi0_result.pi0, False)
>>> # Calculate local FDR
>>> lfdr_values = pyopenms.MultipleTesting.lfdr(
...     p_values, pi0_result.pi0, True, True,
...     pyopenms.MultipleTesting.LfdrTransform.Probit)
)doc")
        .def_static("qValue", [](const std::vector<double>& p_values, double pi0, bool pfdr) { return OpenMS::Math::MultipleTesting::qValue(p_values, pi0, pfdr); }, "p_values"_a, "pi0"_a, "pfdr"_a, "Compute q-values from p-values using the Storey-Tibshirani method")
        .def_static("pNorm", [](const std::vector<double>& stat, const std::vector<double>& stat0) { return OpenMS::Math::MultipleTesting::pNorm(stat, stat0); }, "stat"_a, "stat0"_a, "Compute p-values from observed and null statistics using the empirical distribution")

        .def_static("pi0Est", [](std::vector<double> p_values,
                                  std::vector<double> lambda,
                                  OpenMS::Math::MultipleTesting::Pi0Method method,
                                  int smooth_df,
                                  bool smooth_log_pi0) {
            return OpenMS::Math::MultipleTesting::pi0Est(p_values, lambda, method, smooth_df, smooth_log_pi0);
        }, "p_values"_a, "lambda_"_a, "method"_a, "smooth_df"_a = 3, "smooth_log_pi0"_a = false,
           "Estimate pi0 (proportion of true null hypotheses)")

        .def_static("lfdr", [](std::vector<double> p_values,
                               double pi0,
                               bool trunc,
                               bool monotone,
                               OpenMS::Math::MultipleTesting::LfdrTransform transf,
                               double adj,
                               double eps,
                               size_t gridsize,
                               double cut) {
            return OpenMS::Math::MultipleTesting::lfdr(p_values, pi0, trunc, monotone, transf, adj, eps, gridsize, cut);
        }, "p_values"_a, "pi0"_a, "trunc"_a = true, "monotone"_a = true, "transf"_a, "adj"_a = 1.5, "eps"_a = 1e-8, "gridsize"_a = 100, "cut"_a = 0.05,
           "Compute local FDR values")

        .def_static("pi0MethodToString", &OpenMS::Math::MultipleTesting::pi0MethodToString,
           "m"_a, "Convert Pi0Method enum to string")

        .def_static("toPi0Method", &OpenMS::Math::MultipleTesting::toPi0Method,
           "s"_a, "Convert string to Pi0Method enum (case insensitive)")

        .def_static("lfdrTransformToString", &OpenMS::Math::MultipleTesting::lfdrTransformToString,
           "t"_a, "Convert LfdrTransform enum to string")

        .def_static("toLfdrTransform", &OpenMS::Math::MultipleTesting::toLfdrTransform,
           "s"_a, "Convert string to LfdrTransform enum")
        ;
    // Pi0Method enum nested under MultipleTesting
    nb::enum_<OpenMS::Math::MultipleTesting::Pi0Method>(multipletesting_class, "Pi0Method", nb::is_arithmetic())
        .value("Smoother", OpenMS::Math::MultipleTesting::Pi0Method::Smoother)
        .value("Bootstrap", OpenMS::Math::MultipleTesting::Pi0Method::Bootstrap)
        ;
    // LfdrTransform enum nested under MultipleTesting
    nb::enum_<OpenMS::Math::MultipleTesting::LfdrTransform>(multipletesting_class, "LfdrTransform", nb::is_arithmetic())
        .value("Probit", OpenMS::Math::MultipleTesting::LfdrTransform::Probit)
        .value("Logit", OpenMS::Math::MultipleTesting::LfdrTransform::Logit)
        ;

    // -----------------------------------------------------------------------
    // Param
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Param>(m, "Param",
        R"doc(
Management and storage of parameters / INI files.
This class provides a means to associate string names to int/double/string/StringList values.
It allows for parameter hierarchies and to save/load the data as XML.
Hierarchy levels are separated from each other by colons (e.g., 'algorithm:common:threshold').
Each parameter and section has a description. Newline characters in the description are possible.
Each parameter can be annotated with an arbitrary number of tags (e.g., 'advanced').
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Param &>())
        .def("__copy__", [](const OpenMS::Param& self) { return OpenMS::Param(self); })
        .def("__deepcopy__", [](const OpenMS::Param& self, nb::dict) { return OpenMS::Param(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def("getValue", [](const OpenMS::Param& self, const OpenMS::String& key) { return self.getValue(key); }, "key"_a, "Returns the value of the parameter specified by key. Raises exception if not found")
        .def("getValueType", [](const OpenMS::Param& self, const OpenMS::String& key) { return self.getValueType(key); }, "key"_a, "Returns the type of the parameter specified by key. Raises exception if not found")
        .def("getEntry", [](const OpenMS::Param& self, const OpenMS::String& key) -> const OpenMS::Param::ParamEntry & { return self.getEntry(key); }, "key"_a, nb::rv_policy::reference_internal, "Returns the whole parameter entry (value, description, tags, restrictions). Raises exception if not found")
        .def("getDescription", [](const OpenMS::Param& self, const OpenMS::String& key) { return self.getDescription(key); }, "key"_a, "Returns the description of the parameter specified by key")
        .def("exists", [](const OpenMS::Param& self, const OpenMS::String& key) { return self.exists(key); }, "key"_a, "Returns True if the parameter exists, False otherwise")
        .def("addTag", [](OpenMS::Param& self, const OpenMS::String& key, const OpenMS::String& tag) { return self.addTag(key, tag); }, "key"_a, "tag"_a, "Adds a tag to the entry specified by key (e.g., 'advanced', 'required', 'input file')")
        .def("addTags", [](OpenMS::Param& self, const OpenMS::String& key, const std::vector<std::basic_string<char>>& tags) { return self.addTags(key, tags); }, "key"_a, "tags"_a, "Adds multiple tags to the entry specified by key")
        .def("hasTag", [](const OpenMS::Param& self, const OpenMS::String& key, const OpenMS::String& tag) { return self.hasTag(key, tag); }, "key"_a, "tag"_a, "Returns True if the parameter has the specified tag")
        .def("getTags", [](const OpenMS::Param& self, const OpenMS::String& key) { return self.getTags(key); }, "key"_a, "Returns the tags of the entry specified by key")
        .def("clearTags", [](OpenMS::Param& self, const OpenMS::String& key) { return self.clearTags(key); }, "key"_a, "Removes all tags from the entry specified by key")
        .def("setSectionDescription", [](OpenMS::Param& self, const OpenMS::String& key, const OpenMS::String& description) { return self.setSectionDescription(key, description); }, "key"_a, "description"_a, "Sets a description for an existing section (not for values)")
        .def("getSectionDescription", [](const OpenMS::Param& self, const OpenMS::String& key) { return self.getSectionDescription(key); }, "key"_a, "Returns the description of the section specified by key (empty string if not found)")
        .def("addSection", [](OpenMS::Param& self, const OpenMS::String& key, const OpenMS::String& description) { return self.addSection(key, description); }, "key"_a, "description"_a, "Adds a parameter section with the given key and description")
        .def("size", [](const OpenMS::Param& self) { return self.size(); }, "Returns the number of parameter entries (leaves)")
        .def("empty", [](const OpenMS::Param& self) { return self.empty(); }, "Returns True if there are no entries")
        .def("clear", [](OpenMS::Param& self) { return self.clear(); }, "Deletes all entries")
        .def("insert", [](OpenMS::Param& self, const OpenMS::String& prefix, const OpenMS::Param& param) { return self.insert(prefix, param); }, "prefix"_a, "param"_a, "Inserts all values of another Param object with the given prefix")
        .def("remove", [](OpenMS::Param& self, const OpenMS::String& key) { return self.remove(key); }, "key"_a, "Removes an entry or section (when key ends with ':') by exact name match")
        .def("removeAll", [](OpenMS::Param& self, const OpenMS::String& prefix) { return self.removeAll(prefix); }, "prefix"_a, "Removes all entries and sections that start with the given prefix")
        .def("copy", [](const OpenMS::Param& self, const OpenMS::String& prefix, bool remove_prefix) { return self.copy(prefix, remove_prefix); }, "prefix"_a, "remove_prefix"_a = false,
            R"doc(
Returns a new Param containing all entries that start with the given prefix.
If remove_prefix is True, the prefix is removed from the keys in the returned Param
)doc")
        .def("merge", [](OpenMS::Param& self, const OpenMS::Param& toMerge) { return self.merge(toMerge); }, "toMerge"_a, "Adds missing parameters from another Param object without modifying existing ones")
        .def("setDefaults", [](OpenMS::Param& self, const OpenMS::Param& defaults, const OpenMS::String& prefix, bool showMessage) { return self.setDefaults(defaults, prefix, showMessage); }, "defaults"_a, "prefix"_a = "", "showMessage"_a = false,
            R"doc(
Inserts all values from defaults that are not already set.
Optionally adds a prefix to all keys and prints a message for each default value set
)doc")
        .def("checkDefaults", [](const OpenMS::Param& self, const OpenMS::String& name, const OpenMS::Param& defaults, const OpenMS::String& prefix) { return self.checkDefaults(name, defaults, prefix); }, "name"_a, "defaults"_a, "prefix"_a = "",
            R"doc(
Checks current parameter entries against given defaults.
Validates types, string restrictions, and numeric ranges. Raises exception on invalid parameters
)doc")
        .def("setValidStrings", [](OpenMS::Param& self, const OpenMS::String& key, const std::vector<std::basic_string<char>>& strings) { return self.setValidStrings(key, strings); }, "key"_a, "strings"_a, "Sets the list of valid string values for the parameter (checked by checkDefaults)")
        .def("getValidStrings", [](const OpenMS::Param& self, const OpenMS::String& key) -> const std::vector<std::basic_string<char>> & { return self.getValidStrings(key); }, "key"_a, nb::rv_policy::reference_internal, "Returns the list of valid string values for the parameter")
        .def("setMinInt", [](OpenMS::Param& self, const OpenMS::String& key, int min) { return self.setMinInt(key, min); }, "key"_a, "min"_a, "Sets the minimum allowed value for an integer parameter")
        .def("setMaxInt", [](OpenMS::Param& self, const OpenMS::String& key, int max) { return self.setMaxInt(key, max); }, "key"_a, "max"_a, "Sets the maximum allowed value for an integer parameter")
        .def("setMinFloat", [](OpenMS::Param& self, const OpenMS::String& key, double min) { return self.setMinFloat(key, min); }, "key"_a, "min"_a, "Sets the minimum allowed value for a float parameter")
        .def("setMaxFloat", [](OpenMS::Param& self, const OpenMS::String& key, double max) { return self.setMaxFloat(key, max); }, "key"_a, "max"_a, "Sets the maximum allowed value for a float parameter")
        .def("__iter__", [](OpenMS::Param& self) { return nb::make_iterator<nb::rv_policy::reference_internal>(nb::type<OpenMS::Param>(), "Param_iter", self.begin(), self.end()); })
        .def("__len__", [](OpenMS::Param& self) { return self.size(); })
        .def("_get_all_keys", [](const OpenMS::Param& self) {
            std::vector<std::string> keys;
            for (auto it = self.begin(); it != self.end(); ++it) {
                keys.push_back(it.getName());
            }
            return keys;
        })

        .def("setValue", [](OpenMS::Param& self, const OpenMS::String& key, const OpenMS::ParamValue& value, const OpenMS::String& description, const std::vector<std::string>& tags) {
            self.setValue(key, value, description, tags);
        }, "key"_a, "value"_a, "description"_a = "", "tags"_a = std::vector<std::string>(), "Sets a value with description and tags")
        .def("setValue", [](OpenMS::Param& self, const OpenMS::String& key, const OpenMS::ParamValue& value, const OpenMS::String& description) {
            self.setValue(key, value, description);
        }, "key"_a, "value"_a, "description"_a = "", "Sets a value with description")
        .def("setValue", [](OpenMS::Param& self, const OpenMS::String& key, const OpenMS::ParamValue& value) {
            self.setValue(key, value);
        }, "key"_a, "value"_a, "Set a value for a key")

        // Dict-like API
        .def("keys", [](const OpenMS::Param& self) {
            nb::list result;
            for (auto it = self.begin(); it != self.end(); ++it) {
                result.append(nb::str(it.getName().c_str()));
            }
            return result;
        }, "Return list of parameter keys as str")
        .def("values", [](const OpenMS::Param& self) {
            nb::list result;
            for (auto it = self.begin(); it != self.end(); ++it) {
                result.append(nb::cast(self.getValue(it.getName())));
            }
            return result;
        }, "Return list of parameter values")
        .def("items", [](const OpenMS::Param& self) {
            nb::list result;
            for (auto it = self.begin(); it != self.end(); ++it) {
                std::string key = it.getName();
                result.append(nb::make_tuple(nb::str(key.c_str()), nb::cast(self.getValue(key))));
            }
            return result;
        }, "Return list of (key, value) tuples")
        .def("asDict", [](const OpenMS::Param& self) {
            nb::dict result;
            for (auto it = self.begin(); it != self.end(); ++it) {
                std::string key = it.getName();
                result[nb::str(key.c_str())] = nb::cast(self.getValue(key));
            }
            return result;
        }, "Return dict with str keys")
        .def("to_dict", [](const OpenMS::Param& self) {
            nb::dict result;
            for (auto it = self.begin(); it != self.end(); ++it) {
                std::string key = it.getName();
                result[nb::str(key.c_str())] = nb::cast(self.getValue(key));
            }
            return result;
        }, "Return dict with string keys")
        .def("get", [](const OpenMS::Param& self, const std::string& key, nb::object default_val) -> nb::object {
            if (self.exists(key)) {
                return nb::cast(self.getValue(key));
            }
            return default_val;
        }, "key"_a, "default_"_a = nb::none(), "Get a parameter value by key, returning default if not found")
        .def("__getitem__", [](const OpenMS::Param& self, const std::string& key) {
            if (!self.exists(key)) {
                throw nb::key_error(key.c_str());
            }
            return self.getValue(key);
        }, "key"_a, "Get a parameter value by key, raising KeyError if not found")
        .def("__setitem__", [](OpenMS::Param& self, const std::string& key, const OpenMS::ParamValue& value) {
            self.setValue(key, value);
        }, "key"_a, "value"_a, "Set a parameter value by key")
        .def("__contains__", [](const OpenMS::Param& self, const std::string& key) {
            return self.exists(key);
        }, "key"_a, "Check if a parameter key exists")
        .def("update", [](OpenMS::Param& self, nb::object source, nb::object flag) {
            // Check if source is a Param
            try {
                auto& param_src = nb::cast<const OpenMS::Param&>(source);
                bool filter = !flag.is_none() && nb::cast<bool>(flag);
                for (auto it = param_src.begin(); it != param_src.end(); ++it) {
                    std::string key = it.getName();
                    if (filter && !self.exists(key)) continue;
                    self.setValue(key, param_src.getValue(key));
                }
            } catch (const nb::cast_error&) {
                // Assume it's a dict
                nb::dict d = nb::cast<nb::dict>(source);
                for (auto [k, v] : d) {
                    std::string key = nb::cast<std::string>(k);
                    self.setValue(key, nb::cast<OpenMS::ParamValue>(v));
                }
            }
        }, "source"_a, "flag"_a = nb::none(), "Update parameters from a Param or dict")
        .def_static("from_dict", [](nb::dict d) {
            OpenMS::Param p;
            for (auto [k, v] : d) {
                std::string key = nb::cast<std::string>(k);
                p.setValue(key, nb::cast<OpenMS::ParamValue>(v));
            }
            return p;
        }, "d"_a, "Create a Param from a dict")
        ;


    // -----------------------------------------------------------------------
    // ParamEntry
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Param::ParamEntry>(m, "ParamEntry", "OpenMS class ParamEntry")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::Param::ParamEntry& self) { return OpenMS::Param::ParamEntry(self); })
        .def("__deepcopy__", [](const OpenMS::Param::ParamEntry& self, nb::dict) { return OpenMS::Param::ParamEntry(self); }, "memo"_a)
        .def(nb::init<const std::string&, const OpenMS::ParamValue&, const std::string&, const std::vector<std::string>&>(), "name"_a, "value"_a, "description"_a, "tags"_a = std::vector<std::string>())
        .def_rw("name", &OpenMS::Param::ParamEntry::name)
        .def_rw("description", &OpenMS::Param::ParamEntry::description)
        .def_rw("value", &OpenMS::Param::ParamEntry::value)
        .def_rw("tags", &OpenMS::Param::ParamEntry::tags)
        .def_rw("valid_strings", &OpenMS::Param::ParamEntry::valid_strings)
        .def_rw("max_float", &OpenMS::Param::ParamEntry::max_float)
        .def_rw("min_float", &OpenMS::Param::ParamEntry::min_float)
        .def_rw("max_int", &OpenMS::Param::ParamEntry::max_int)
        .def_rw("min_int", &OpenMS::Param::ParamEntry::min_int)
        .def("isValid", [](const OpenMS::Param::ParamEntry& self) {
            std::string msg;
            bool valid = self.isValid(msg);
            return nb::make_tuple(valid, msg);
        }, "Check if value fulfills restrictions. Returns (valid, message)")
        .def("__eq__", &OpenMS::Param::ParamEntry::operator==)
        ;


    // -----------------------------------------------------------------------
    // ParamNode
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Param::ParamNode>(m, "ParamNode", "OpenMS class ParamNode")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::Param::ParamNode& self) { return OpenMS::Param::ParamNode(self); })
        .def("__deepcopy__", [](const OpenMS::Param::ParamNode& self, nb::dict) { return OpenMS::Param::ParamNode(self); }, "memo"_a)
        .def(nb::init<const std::string&, const std::string&>(), "name"_a, "description"_a)
        .def_rw("name", &OpenMS::Param::ParamNode::name)
        .def_rw("description", &OpenMS::Param::ParamNode::description)
        .def_rw("entries", &OpenMS::Param::ParamNode::entries)
        .def_rw("nodes", &OpenMS::Param::ParamNode::nodes)
        .def("size", &OpenMS::Param::ParamNode::size)
        .def("suffix", &OpenMS::Param::ParamNode::suffix, "key"_a)
        .def("__eq__", &OpenMS::Param::ParamNode::operator==)
        .def("findEntryRecursive", [](OpenMS::Param::ParamNode& self, const std::string& name) -> OpenMS::Param::ParamEntry* {
            return self.findEntryRecursive(name);
        }, "name"_a, nb::rv_policy::reference_internal, "Finds an entry by name recursively")
        .def("findParentOf", [](OpenMS::Param::ParamNode& self, const std::string& name) -> OpenMS::Param::ParamNode* {
            return self.findParentOf(name);
        }, "name"_a, nb::rv_policy::reference_internal, "Finds the parent node of the entry with the given name")
        .def("insert", [](OpenMS::Param::ParamNode& self, const OpenMS::Param::ParamNode& node, const std::string& prefix) {
            self.insert(node, prefix);
        }, "node"_a, "prefix"_a = "", "Inserts a node")
        .def("insert", [](OpenMS::Param::ParamNode& self, const OpenMS::Param::ParamEntry& entry, const std::string& prefix) {
            self.insert(entry, prefix);
        }, "entry"_a, "prefix"_a = "", "Inserts an entry")
        ;

    // -----------------------------------------------------------------------
    // Pi0Result
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Math::Pi0Result>(m, "Pi0Result",
        R"doc(
Result of pi0 estimation for multiple testing correction.
Contains the estimated proportion of true null hypotheses (pi0),
which is a key parameter in FDR estimation.
Attributes:
pi0: Estimated proportion of true null hypotheses (0-1)
pi0_lambda: Vector of pi0 estimates at each lambda threshold
lambda_: Vector of lambda threshold values used
pi0_smooth: Whether smoothing was successfully applied
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Math::Pi0Result &>())
        .def("__copy__", [](const OpenMS::Math::Pi0Result& self) { return OpenMS::Math::Pi0Result(self); })
        .def("__deepcopy__", [](const OpenMS::Math::Pi0Result& self, nb::dict) { return OpenMS::Math::Pi0Result(self); }, "memo"_a)
        .def_rw("pi0", &OpenMS::Math::Pi0Result::pi0)
        .def_rw("pi0_lambda", &OpenMS::Math::Pi0Result::pi0_lambda)
        .def_rw("lambda_", &OpenMS::Math::Pi0Result::lambda_)
        .def_rw("pi0_smooth", &OpenMS::Math::Pi0Result::pi0_smooth)
        ;

    // -----------------------------------------------------------------------
    // QTCluster
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::QTCluster>(m, "QTCluster", "A representation of a QT cluster used for feature grouping")
        .def(nb::init<const OpenMS::QTCluster &>())
        .def("__copy__", [](const OpenMS::QTCluster& self) { return OpenMS::QTCluster(self); })
        .def("__deepcopy__", [](const OpenMS::QTCluster& self, nb::dict) { return OpenMS::QTCluster(self); }, "memo"_a)
        .def("getCenterRT", [](const OpenMS::QTCluster& self) { return self.getCenterRT(); }, "Returns the RT value of the cluster")
        .def("getCenterMZ", [](const OpenMS::QTCluster& self) { return self.getCenterMZ(); }, "Returns the m/z value of the cluster center")
        .def("getXCoord", [](const OpenMS::QTCluster& self) { return self.getXCoord(); }, "Returns the x coordinate in the grid")
        .def("getYCoord", [](const OpenMS::QTCluster& self) { return self.getYCoord(); }, "Returns the y coordinate in the grid")
        .def("size", [](const OpenMS::QTCluster& self) { return self.size(); }, "Returns the size of the cluster (number of elements, incl. center)")
        .def(nb::self < nb::self)
        .def("getQuality", [](OpenMS::QTCluster& self) { return self.getQuality(); }, "Returns the cluster quality and recomputes if necessary")
        .def("getAnnotations", [](OpenMS::QTCluster& self) -> const std::set<OpenMS::AASequence> & { return self.getAnnotations(); }, nb::rv_policy::reference_internal, "Returns the set of peptide sequences annotated to the cluster center")
        .def("setInvalid", [](OpenMS::QTCluster& self) { return self.setInvalid(); }, "Sets current cluster as invalid (also frees some memory)")
        .def("isInvalid", [](const OpenMS::QTCluster& self) { return self.isInvalid(); }, "Whether current cluster is invalid")
        .def("initializeCluster", [](OpenMS::QTCluster& self) { return self.initializeCluster(); }, "Has to be called before adding elements (calling QTCluster::add)")
        .def("finalizeCluster", [](OpenMS::QTCluster& self) { return self.finalizeCluster(); }, "Has to be called after adding elements (after calling QTCluster::add one or multiple times)")
        .def("__len__", [](OpenMS::QTCluster& self) { return self.size(); })
        ;

    // -----------------------------------------------------------------------
    // RankData
    // -----------------------------------------------------------------------
    auto rankdata_class = nb::class_<OpenMS::Math::RankData>(m, "RankData",
        R"doc(
Rank items (1-based) with SciPy-like tie and NaN handling.
This wrapper exposes the concrete forwarders rankdata_double/float/int
and the enum classes RankData.Method and RankData.NaNPolicy.
)doc")

        .def_static("rankdata_double", [](std::vector<double> a, OpenMS::Math::RankData::Method m, OpenMS::Math::RankData::NaNPolicy p) {
            auto result = OpenMS::Math::RankData::rankdata_double(a, m, p);
            size_t n = result.size();
            double* data = new double[n];
            std::copy(result.begin(), result.end(), data);
            nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<double*>(p); });
            return nb::ndarray<nb::numpy, double, nb::ndim<1>>(data, {n}, owner);
        }, "a"_a, "method"_a, "nan_policy"_a)

        .def_static("rankdata_float", [](std::vector<float> a, OpenMS::Math::RankData::Method m, OpenMS::Math::RankData::NaNPolicy p) {
            auto result = OpenMS::Math::RankData::rankdata_float(a, m, p);
            size_t n = result.size();
            double* data = new double[n];
            std::copy(result.begin(), result.end(), data);
            nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<double*>(p); });
            return nb::ndarray<nb::numpy, double, nb::ndim<1>>(data, {n}, owner);
        }, "a"_a, "method"_a, "nan_policy"_a)

        .def_static("rankdata_int", [](std::vector<int> a, OpenMS::Math::RankData::Method m, OpenMS::Math::RankData::NaNPolicy p) {
            auto result = OpenMS::Math::RankData::rankdata_int(a, m, p);
            size_t n = result.size();
            double* data = new double[n];
            std::copy(result.begin(), result.end(), data);
            nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<double*>(p); });
            return nb::ndarray<nb::numpy, double, nb::ndim<1>>(data, {n}, owner);
        }, "a"_a, "method"_a, "nan_policy"_a)
        ;

    // -----------------------------------------------------------------------
    // SpectraSTSimilarityScore
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SpectraSTSimilarityScore>(m, "SpectraSTSimilarityScore",
        R"doc(
Similarity score of SpectraST
PeakSpectrumCompareFunctor inheritance
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::SpectraSTSimilarityScore &>())
        .def("__copy__", [](const OpenMS::SpectraSTSimilarityScore& self) { return OpenMS::SpectraSTSimilarityScore(self); })
        .def("__deepcopy__", [](const OpenMS::SpectraSTSimilarityScore& self, nb::dict) { return OpenMS::SpectraSTSimilarityScore(self); }, "memo"_a)
        .def("preprocess", [](OpenMS::SpectraSTSimilarityScore& self, OpenMS::MSSpectrum& spec, float remove_peak_intensity_threshold, unsigned int cut_peaks_below, size_t min_peak_number, size_t max_peak_number) { return self.preprocess(spec, remove_peak_intensity_threshold, cut_peaks_below, min_peak_number, max_peak_number); }, "spec"_a, "remove_peak_intensity_threshold"_a = 2.01, "cut_peaks_below"_a = 1000, "min_peak_number"_a = 5, "max_peak_number"_a = 150)
        .def("transform", [](OpenMS::SpectraSTSimilarityScore& self, const OpenMS::MSSpectrum& spec) { return self.transform(spec); }, "spec"_a, "Spectrum is transformed into a binned spectrum with bin size 1 and spread 1 and the intensities are normalized")
        .def("dot_bias", [](const OpenMS::SpectraSTSimilarityScore& self, const OpenMS::BinnedSpectrum& bin1, const OpenMS::BinnedSpectrum& bin2, double dot_product) { return self.dot_bias(bin1, bin2, dot_product); }, "bin1"_a, "bin2"_a, "dot_product"_a = -1)
        .def("delta_D", [](OpenMS::SpectraSTSimilarityScore& self, double top_hit, double runner_up) { return self.delta_D(top_hit, runner_up); }, "top_hit"_a, "runner_up"_a,
            R"doc(
Calculates how much of the dot product is dominated by a few peaks
:param dot_product: If -1 this value will be calculated as well.
:param bin1: First spectrum in binned representation
:param bin2: Second spectrum in binned representation
)doc")
        .def("compute_F", [](OpenMS::SpectraSTSimilarityScore& self, double dot_product, double delta_D, double dot_bias) { return self.compute_F(dot_product, delta_D, dot_bias); }, "dot_product"_a, "delta_D"_a, "dot_bias"_a,
            R"doc(
Calculates the normalized distance between top_hit and runner_up
:param top_hit: Is the best score for a given match
:param runner_up: A match with a worse score than top_hit, e.g. the second best score
:returns: normalized distance
)doc")
        ;

    // -----------------------------------------------------------------------
    // SpectrumAlignmentScore
    // -----------------------------------------------------------------------
    auto spectrumalignmentscore_class = nb::class_<OpenMS::SpectrumAlignmentScore>(m, "SpectrumAlignmentScore",
        R"doc(
DefaultParamHandler

Similarity score via spectra alignment
This class implements a simple scoring based on the alignment of spectra. This alignment
is implemented in the SpectrumAlignment class and performs a dynamic programming alignment
of the peaks, minimizing the distances between the aligned peaks and maximizing the number
of peak pairs
The scoring is done via the simple formula score = sum / (sqrt(sum1 * sum2)). sum is the
product of the intensities of the aligned peaks, with the given exponent (default is 2)
sum1 and sum2 are the sum of the intensities squared for each peak of both spectra respectively
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::SpectrumAlignmentScore &>())
        .def("__copy__", [](const OpenMS::SpectrumAlignmentScore& self) { return OpenMS::SpectrumAlignmentScore(self); })
        .def("__deepcopy__", [](const OpenMS::SpectrumAlignmentScore& self, nb::dict) { return OpenMS::SpectrumAlignmentScore(self); }, "memo"_a)
        .def("__call__", [](const OpenMS::SpectrumAlignmentScore& self, const OpenMS::MSSpectrum& spec1, const OpenMS::MSSpectrum& spec2) { return self(spec1, spec2); }, "spec1"_a, "spec2"_a, "Compute the similarity score between two spectra")
        .def("__call__", [](const OpenMS::SpectrumAlignmentScore& self, const OpenMS::MSSpectrum& spec) { return self(spec); }, "spec"_a, "Compute the self-similarity score of a spectrum")
        ;
    def_DefaultParamHandler<OpenMS::SpectrumAlignmentScore>(spectrumalignmentscore_class);

    // -----------------------------------------------------------------------
    // VersionInfo
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::VersionInfo>(m, "VersionInfo", "Version information class")
        .def_static("getTime", []() { return OpenMS::VersionInfo::getTime(); })
        .def_static("getVersion", []() { return OpenMS::VersionInfo::getVersion(); })
        .def_static("getVersionStruct", []() { return OpenMS::VersionInfo::getVersionStruct(); })
        .def_static("getRevision", []() { return OpenMS::VersionInfo::getRevision(); })
        .def_static("getBranch", []() { return OpenMS::VersionInfo::getBranch(); })
        ;

    // -----------------------------------------------------------------------
    // VersionDetails
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::VersionInfo::VersionDetails>(m, "VersionDetails", "Version details struct")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::VersionInfo::VersionDetails& self) { return OpenMS::VersionInfo::VersionDetails(self); })
        .def("__deepcopy__", [](const OpenMS::VersionInfo::VersionDetails& self, nb::dict) { return OpenMS::VersionInfo::VersionDetails(self); }, "memo"_a)
        .def_rw("version_major", &OpenMS::VersionInfo::VersionDetails::version_major)
        .def_rw("version_minor", &OpenMS::VersionInfo::VersionDetails::version_minor)
        .def_rw("version_patch", &OpenMS::VersionInfo::VersionDetails::version_patch)
        .def_rw("pre_release_identifier", &OpenMS::VersionInfo::VersionDetails::pre_release_identifier)
        .def_static("create", &OpenMS::VersionInfo::VersionDetails::create, "version"_a)
        .def("__eq__", [](const OpenMS::VersionInfo::VersionDetails& a, const OpenMS::VersionInfo::VersionDetails& b) { return a == b; })
        .def("__ne__", [](const OpenMS::VersionInfo::VersionDetails& a, const OpenMS::VersionInfo::VersionDetails& b) { return a != b; })
        .def("__lt__", [](const OpenMS::VersionInfo::VersionDetails& a, const OpenMS::VersionInfo::VersionDetails& b) { return a < b; })
        .def("__gt__", [](const OpenMS::VersionInfo::VersionDetails& a, const OpenMS::VersionInfo::VersionDetails& b) { return a > b; })
        ;

    // -----------------------------------------------------------------------
    // Compomer
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Compomer>(m, "Compomer",
        "Holds information on an edge connecting two features from a (putative) charge ladder")
        .def(nb::init<>())
        .def(nb::init<OpenMS::Int, double, double>(), "net_charge"_a, "mass"_a, "log_p"_a)
        .def(nb::init<const OpenMS::Compomer &>())
        .def("__copy__", [](const OpenMS::Compomer& self) { return OpenMS::Compomer(self); })
        .def("__deepcopy__", [](const OpenMS::Compomer& self, nb::dict) { return OpenMS::Compomer(self); }, "memo"_a)
        .def("setID", [](OpenMS::Compomer& self, const OpenMS::Size& id) { self.setID(id); }, "id"_a)
        .def("getID", [](const OpenMS::Compomer& self) { return self.getID(); })
        .def("getNetCharge", [](const OpenMS::Compomer& self) { return self.getNetCharge(); })
        .def("getMass", [](const OpenMS::Compomer& self) { return self.getMass(); })
        .def("getPositiveCharges", [](const OpenMS::Compomer& self) { return self.getPositiveCharges(); })
        .def("getNegativeCharges", [](const OpenMS::Compomer& self) { return self.getNegativeCharges(); })
        .def("getLogP", [](const OpenMS::Compomer& self) { return self.getLogP(); })
        .def("getRTShift", [](const OpenMS::Compomer& self) { return self.getRTShift(); })
        .def("getAdductsAsString", [](const OpenMS::Compomer& self) { return self.getAdductsAsString(); })
        .def("getAdductsAsString", [](const OpenMS::Compomer& self, OpenMS::UInt side) { return self.getAdductsAsString(side); }, "side"_a)
        .def("add", [](OpenMS::Compomer& self, const OpenMS::Adduct& a, OpenMS::UInt side) { self.add(a, side); }, "a"_a, "side"_a)
        .def(nb::self == nb::self)
        .def("__hash__", [](const OpenMS::Compomer& self) { return std::hash<OpenMS::Compomer>{}(self); })
        .def("getLabels", [](const OpenMS::Compomer& self, OpenMS::UInt side) { return self.getLabels(side); }, "side"_a, "Returns the labels for the given side")
        .def("isConflicting", [](const OpenMS::Compomer& self, const OpenMS::Compomer& cmp, OpenMS::UInt side_this, OpenMS::UInt side_other) { return self.isConflicting(cmp, side_this, side_other); }, "cmp"_a, "side_this"_a, "side_other"_a, "Returns true if the compomers are conflicting")
        .def("isSingleAdduct", [](const OpenMS::Compomer& self, OpenMS::Adduct& a, OpenMS::UInt side) { return self.isSingleAdduct(a, side); }, "a"_a, "side"_a, "Returns true if the compomer has a single adduct on the given side")
        .def("removeAdduct", [](const OpenMS::Compomer& self, const OpenMS::Adduct& a) { return self.removeAdduct(a); }, "a"_a, "Returns a new compomer without the given adduct")
        .def("removeAdduct", [](const OpenMS::Compomer& self, const OpenMS::Adduct& a, OpenMS::UInt side) { return self.removeAdduct(a, side); }, "a"_a, "side"_a, "Returns a new compomer without the given adduct on the given side")
        ;

    // -----------------------------------------------------------------------
    // SIDE (Compomer::SIDE)
    // -----------------------------------------------------------------------
    nb::enum_<OpenMS::Compomer::SIDE>(m, "SIDE",
        "Side enum for adduct decomposition (LEFT, RIGHT, BOTH)", nb::is_arithmetic())
        .value("LEFT", OpenMS::Compomer::LEFT)
        .value("RIGHT", OpenMS::Compomer::RIGHT)
        .value("BOTH", OpenMS::Compomer::BOTH)

        ;

    // -----------------------------------------------------------------------
    // ChargePair
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ChargePair>(m, "ChargePair",
        "Representation of a (putative) link between two Features with different charge")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ChargePair &>())
        .def("__copy__", [](const OpenMS::ChargePair& self) { return OpenMS::ChargePair(self); })
        .def("__deepcopy__", [](const OpenMS::ChargePair& self, nb::dict) { return OpenMS::ChargePair(self); }, "memo"_a)
        .def("getCharge", [](const OpenMS::ChargePair& self, OpenMS::UInt pairID) { return self.getCharge(pairID); }, "pairID"_a)
        .def("setCharge", [](OpenMS::ChargePair& self, OpenMS::UInt pairID, OpenMS::Int e) { self.setCharge(pairID, e); }, "pairID"_a, "e"_a)
        .def("getElementIndex", [](const OpenMS::ChargePair& self, OpenMS::UInt pairID) { return self.getElementIndex(pairID); }, "pairID"_a)
        .def("setElementIndex", [](OpenMS::ChargePair& self, OpenMS::UInt pairID, OpenMS::Size e) { self.setElementIndex(pairID, e); }, "pairID"_a, "e"_a)
        .def("getCompomer", [](const OpenMS::ChargePair& self) -> const OpenMS::Compomer& { return self.getCompomer(); }, nb::rv_policy::reference_internal)
        .def("setCompomer", [](OpenMS::ChargePair& self, const OpenMS::Compomer& compomer) { self.setCompomer(compomer); }, "compomer"_a)
        .def("getMassDiff", [](const OpenMS::ChargePair& self) { return self.getMassDiff(); })
        .def("setMassDiff", [](OpenMS::ChargePair& self, double mass_diff) { self.setMassDiff(mass_diff); }, "mass_diff"_a)
        .def("getEdgeScore", [](const OpenMS::ChargePair& self) { return self.getEdgeScore(); })
        .def("setEdgeScore", [](OpenMS::ChargePair& self, double score) { self.setEdgeScore(score); }, "score"_a)
        .def("isActive", [](const OpenMS::ChargePair& self) { return self.isActive(); })
        .def("setActive", [](OpenMS::ChargePair& self, bool active) { self.setActive(active); }, "active"_a)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("__hash__", [](const OpenMS::ChargePair& self) { return std::hash<OpenMS::ChargePair>{}(self); })
        ;

    // -----------------------------------------------------------------------
    // CVMappings
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::CVMappings>(m, "CVMappings",
        "Representation of controlled vocabulary mapping rules (for PSI formats)")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::CVMappings &>())
        .def("__copy__", [](const OpenMS::CVMappings& self) { return OpenMS::CVMappings(self); })
        .def("__deepcopy__", [](const OpenMS::CVMappings& self, nb::dict) { return OpenMS::CVMappings(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("setCVReferences", [](OpenMS::CVMappings& self, const std::vector<OpenMS::CVReference>& cv_references) { self.setCVReferences(cv_references); }, "cv_references"_a, "Sets the CV references")
        .def("getCVReferences", [](const OpenMS::CVMappings& self) -> const std::vector<OpenMS::CVReference>& { return self.getCVReferences(); }, nb::rv_policy::reference_internal, "Returns the CV references")
        .def("addCVReference", [](OpenMS::CVMappings& self, const OpenMS::CVReference& cv_reference) { self.addCVReference(cv_reference); }, "cv_reference"_a, "Adds a CV reference")
        .def("hasCVReference", [](OpenMS::CVMappings& self, const OpenMS::String& identifier) { return self.hasCVReference(identifier); }, "identifier"_a, "Returns true if a CV reference with the given identifier exists")
        .def("setMappingRules", [](OpenMS::CVMappings& self, const std::vector<OpenMS::CVMappingRule>& cv_mapping_rules) { self.setMappingRules(cv_mapping_rules); }, "cv_mapping_rules"_a, "Sets the mapping rules")
        .def("getMappingRules", [](const OpenMS::CVMappings& self) -> const std::vector<OpenMS::CVMappingRule>& { return self.getMappingRules(); }, nb::rv_policy::reference_internal, "Returns the mapping rules")
        .def("addMappingRule", [](OpenMS::CVMappings& self, const OpenMS::CVMappingRule& cv_mapping_rule) { self.addMappingRule(cv_mapping_rule); }, "cv_mapping_rule"_a, "Adds a mapping rule")
        ;

    // -----------------------------------------------------------------------
    // QuotingMethod (String::QuotingMethod)
    // -----------------------------------------------------------------------
    nb::enum_<OpenMS::String::QuotingMethod>(m, "QuotingMethod",
        "Method for quoting strings in CSV output", nb::is_arithmetic())
        .value("NONE", OpenMS::String::NONE)
        .value("ESCAPE", OpenMS::String::ESCAPE)
        .value("DOUBLE", OpenMS::String::DOUBLE)

        ;

    // -----------------------------------------------------------------------
    // UnitType (DataValue::UnitType)
    // -----------------------------------------------------------------------
    nb::enum_<OpenMS::DataValue::UnitType>(m, "UnitType",
        "Unit ontology type for DataValue", nb::is_arithmetic())
        .value("UNIT_ONTOLOGY", OpenMS::DataValue::UNIT_ONTOLOGY)
        .value("MS_ONTOLOGY", OpenMS::DataValue::MS_ONTOLOGY)
        .value("OTHER", OpenMS::DataValue::OTHER)

        ;

    // -----------------------------------------------------------------------
    // ValueType (ParamValue::ValueType)
    // -----------------------------------------------------------------------
    nb::enum_<OpenMS::ParamValue::ValueType>(m, "ValueType",
        "Type tag for ParamValue", nb::is_arithmetic())
        .value("STRING_VALUE", OpenMS::ParamValue::STRING_VALUE)
        .value("INT_VALUE", OpenMS::ParamValue::INT_VALUE)
        .value("DOUBLE_VALUE", OpenMS::ParamValue::DOUBLE_VALUE)
        .value("STRING_LIST", OpenMS::ParamValue::STRING_LIST)
        .value("INT_LIST", OpenMS::ParamValue::INT_LIST)
        .value("DOUBLE_LIST", OpenMS::ParamValue::DOUBLE_LIST)
        .value("EMPTY_VALUE", OpenMS::ParamValue::EMPTY_VALUE)

        ;

    // -----------------------------------------------------------------------
    // CVMappingRule
    // -----------------------------------------------------------------------
    nb::enum_<OpenMS::CVMappingRule::RequirementLevel>(m, "RequirementLevel", nb::is_arithmetic(), "Requirement level for CV mapping rules (MUST, SHOULD, MAY)")
        .value("MUST", OpenMS::CVMappingRule::MUST)
        .value("SHOULD", OpenMS::CVMappingRule::SHOULD)
        .value("MAY", OpenMS::CVMappingRule::MAY)
        ;

    nb::enum_<OpenMS::CVMappingRule::CombinationsLogic>(m, "CombinationsLogic", nb::is_arithmetic(), "Logic for combining CV mapping rule terms (OR, AND, XOR)")
        .value("OR", OpenMS::CVMappingRule::OR)
        .value("AND", OpenMS::CVMappingRule::AND)
        .value("XOR", OpenMS::CVMappingRule::XOR)
        ;

    nb::class_<OpenMS::CVMappingRule>(m, "CVMappingRule",
        "Representation of a CV mapping rule")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::CVMappingRule &>())
        .def("__copy__", [](const OpenMS::CVMappingRule& self) { return OpenMS::CVMappingRule(self); })
        .def("__deepcopy__", [](const OpenMS::CVMappingRule& self, nb::dict) { return OpenMS::CVMappingRule(self); }, "memo"_a)
        .def("setIdentifier", [](OpenMS::CVMappingRule& self, const OpenMS::String& id) { self.setIdentifier(id); }, "identifier"_a)
        .def("getIdentifier", [](const OpenMS::CVMappingRule& self) { return self.getIdentifier(); })
        .def("setElementPath", [](OpenMS::CVMappingRule& self, const OpenMS::String& path) { self.setElementPath(path); }, "element_path"_a)
        .def("getElementPath", [](const OpenMS::CVMappingRule& self) { return self.getElementPath(); })
        .def("setRequirementLevel", [](OpenMS::CVMappingRule& self, OpenMS::CVMappingRule::RequirementLevel level) { self.setRequirementLevel(level); }, "level"_a)
        .def("getRequirementLevel", [](const OpenMS::CVMappingRule& self) { return self.getRequirementLevel(); })
        .def("setCombinationsLogic", [](OpenMS::CVMappingRule& self, OpenMS::CVMappingRule::CombinationsLogic logic) { self.setCombinationsLogic(logic); }, "logic"_a)
        .def("getCombinationsLogic", [](const OpenMS::CVMappingRule& self) { return self.getCombinationsLogic(); })
        .def("setScopePath", [](OpenMS::CVMappingRule& self, const OpenMS::String& path) { self.setScopePath(path); }, "path"_a)
        .def("getScopePath", [](const OpenMS::CVMappingRule& self) { return self.getScopePath(); })
        .def("setCVTerms", [](OpenMS::CVMappingRule& self, const std::vector<OpenMS::CVMappingTerm>& terms) { self.setCVTerms(terms); }, "cv_terms"_a)
        .def("getCVTerms", [](const OpenMS::CVMappingRule& self) { return self.getCVTerms(); })
        .def("addCVTerm", [](OpenMS::CVMappingRule& self, const OpenMS::CVMappingTerm& term) { self.addCVTerm(term); }, "cv_term"_a)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        ;

    // -----------------------------------------------------------------------
    // UniqueIdGenerator
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::UniqueIdGenerator>(m, "UniqueIdGenerator",
        "Generates unique 64-bit IDs")
        .def_static("getUniqueId", []() { return OpenMS::UniqueIdGenerator::getUniqueId(); },
            "Returns a new unique ID")
        .def_static("setSeed", [](OpenMS::UInt64 seed) { OpenMS::UniqueIdGenerator::setSeed(seed); },
            "seed"_a, "Set the random generator seed")
        .def_static("getSeed", []() { return OpenMS::UniqueIdGenerator::getSeed(); },
            "Get the current seed value")
        ;

    // -----------------------------------------------------------------------
    // StringView
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::StringView>(m, "StringView",
        "Lightweight non-owning view on a string")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::StringView& self) { return OpenMS::StringView(self); })
        .def("__deepcopy__", [](const OpenMS::StringView& self, nb::dict) { return OpenMS::StringView(self); }, "memo"_a)
        .def("size", [](const OpenMS::StringView& self) { return self.size(); })
        .def("getString", [](const OpenMS::StringView& self) { return self.getString(); })
        .def("__len__", [](const OpenMS::StringView& self) { return self.size(); })
        .def("__str__", [](const OpenMS::StringView& self) { return self.getString(); })
        .def(nb::self == nb::self)
        .def(nb::self < nb::self)
        .def("substr", [](const OpenMS::StringView& self, size_t start, size_t length) { return self.substr(start, length); }, "start"_a, "length"_a, "Returns a substring view")
        ;

    // -----------------------------------------------------------------------
    // Math PPM/mass functions
    // -----------------------------------------------------------------------
    m.def("getPPM", [](double mz_obs, double mz_ref) { return OpenMS::Math::getPPM(mz_obs, mz_ref); },
        "mz_obs"_a, "mz_ref"_a, "Compute parts-per-million of two m/z values (can be negative)");
    m.def("getPPMAbs", [](double mz_obs, double mz_ref) { return OpenMS::Math::getPPMAbs(mz_obs, mz_ref); },
        "mz_obs"_a, "mz_ref"_a, "Compute absolute parts-per-million of two m/z values (always >= 0)");
    m.def("ppmToMass", [](double ppm, double mz_ref) { return OpenMS::Math::ppmToMass(ppm, mz_ref); },
        "ppm"_a, "mz_ref"_a, "Compute the mass diff in Th, given a ppm value and a reference point");
    m.def("ppmToMassAbs", [](double ppm, double mz_ref) { return OpenMS::Math::ppmToMassAbs(ppm, mz_ref); },
        "ppm"_a, "mz_ref"_a, "Compute the absolute mass diff in Th, given a ppm value and a reference point");

    // -----------------------------------------------------------------------
    // __static_* module-level wrappers for Math
    // -----------------------------------------------------------------------
    m.def("__static_Math_getPPM", [](double mz_obs, double mz_ref) -> double { return OpenMS::Math::getPPM(mz_obs, mz_ref); }, "mz_obs"_a, "mz_ref"_a);
    m.def("__static_Math_getPPMAbs", [](double mz_obs, double mz_ref) -> double { return OpenMS::Math::getPPMAbs(mz_obs, mz_ref); }, "mz_obs"_a, "mz_ref"_a);
    m.def("__static_Math_ppmToMass", [](double ppm, double mz_ref) -> double { return OpenMS::Math::ppmToMass(ppm, mz_ref); }, "ppm"_a, "mz_ref"_a);
    m.def("__static_Math_ppmToMassAbs", [](double ppm, double mz_ref) -> double { return OpenMS::Math::ppmToMassAbs(ppm, mz_ref); }, "ppm"_a, "mz_ref"_a);

    // -----------------------------------------------------------------------
    // MassExplainer
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MassExplainer>(m, "MassExplainer",
        "Computes empirical formulas for given mass differences using a set of allowed elements")
        .def(nb::init<>())
        .def(nb::init<OpenMS::MassExplainer::AdductsType>(), "adduct_base"_a)
        .def(nb::init<OpenMS::Int, OpenMS::Int, OpenMS::Int, double>(),
            "q_min"_a, "q_max"_a, "max_span"_a, "thresh_logp"_a)
        .def(nb::init<OpenMS::MassExplainer::AdductsType, OpenMS::Int, OpenMS::Int, OpenMS::Int, double, OpenMS::Size>(),
            "adduct_base"_a, "q_min"_a, "q_max"_a, "max_span"_a, "thresh_logp"_a, "max_neutrals"_a)
        .def("__copy__", [](const OpenMS::MassExplainer& self) { return OpenMS::MassExplainer(self); })
        .def("__deepcopy__", [](const OpenMS::MassExplainer& self, nb::dict) { return OpenMS::MassExplainer(self); }, "memo"_a)
        .def("compute", &OpenMS::MassExplainer::compute,
            "Compute all possible mass differences and their explanations")
        .def("setAdductBase", &OpenMS::MassExplainer::setAdductBase, "adduct_base"_a,
            "Set the base set of allowed adducts")
        .def("getAdductBase", &OpenMS::MassExplainer::getAdductBase,
            "Get the current set of allowed adducts")
        .def("getCompomerById", &OpenMS::MassExplainer::getCompomerById, "id"_a,
            nb::rv_policy::reference_internal, "Get a specific compomer by its ID")
        .def("query", [](const OpenMS::MassExplainer& self, OpenMS::Int net_charge, float mass_to_explain, float mass_delta, float thresh_log_p) {
            std::vector<OpenMS::Compomer>::const_iterator first, last;
            OpenMS::SignedSize count = self.query(net_charge, mass_to_explain, mass_delta, thresh_log_p, first, last);
            std::vector<OpenMS::Compomer> results(first, last);
            return nb::make_tuple(count, results);
        }, "net_charge"_a, "mass_to_explain"_a, "mass_delta"_a, "thresh_log_p"_a,
            "Search for explanations of a given mass difference, returns (count, explanations)")
        ;

    // -----------------------------------------------------------------------
    // OSWHierarchy
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::OSWHierarchy> osw_hierarchy(m, "OSWHierarchy", "Hierarchy levels of the OSWData tree");
    nb::enum_<OpenMS::OSWHierarchy::Level>(osw_hierarchy, "Level")
        .value("PROTEIN", OpenMS::OSWHierarchy::Level::PROTEIN)
        .value("PEPTIDE", OpenMS::OSWHierarchy::Level::PEPTIDE)
        .value("FEATURE", OpenMS::OSWHierarchy::Level::FEATURE)
        .value("TRANSITION", OpenMS::OSWHierarchy::Level::TRANSITION)
        .export_values();

    // -----------------------------------------------------------------------
    // OSWIndexTrace
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::OSWIndexTrace>(m, "OSWIndexTrace", "Describes a node in the OSWData model tree")
        .def(nb::init<>())
        .def_rw("idx_prot", &OpenMS::OSWIndexTrace::idx_prot)
        .def_rw("idx_pep", &OpenMS::OSWIndexTrace::idx_pep)
        .def_rw("idx_feat", &OpenMS::OSWIndexTrace::idx_feat)
        .def_rw("idx_trans", &OpenMS::OSWIndexTrace::idx_trans)
        .def_rw("lowest", &OpenMS::OSWIndexTrace::lowest)
        .def("isSet", &OpenMS::OSWIndexTrace::isSet)
        ;

    // -----------------------------------------------------------------------
    // OSWTransition
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::OSWTransition>(m, "OSWTransition", "High-level meta data of a transition")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::OSWTransition&>())
        .def(nb::init<const OpenMS::String&, const OpenMS::UInt32, const float, const char, const bool>(),
            "annotation"_a, "id"_a, "product_mz"_a, "type"_a, "is_decoy"_a)
        .def("getAnnotation", &OpenMS::OSWTransition::getAnnotation, nb::rv_policy::reference_internal)
        .def("getID", &OpenMS::OSWTransition::getID)
        .def("getProductMZ", &OpenMS::OSWTransition::getProductMZ)
        .def("getType", [](const OpenMS::OSWTransition& self) { return std::string(1, self.getType()); })
        .def("isDecoy", &OpenMS::OSWTransition::isDecoy)
        ;

    // -----------------------------------------------------------------------
    // OSWPeakGroup
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::OSWPeakGroup>(m, "OSWPeakGroup", "A peak group (feature) defined on a small RT range")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::OSWPeakGroup&>())
        .def("getRTExperimental", &OpenMS::OSWPeakGroup::getRTExperimental)
        .def("getRTLeftWidth", &OpenMS::OSWPeakGroup::getRTLeftWidth)
        .def("getRTRightWidth", &OpenMS::OSWPeakGroup::getRTRightWidth)
        .def("getRTDelta", &OpenMS::OSWPeakGroup::getRTDelta)
        .def("getQValue", &OpenMS::OSWPeakGroup::getQValue)
        .def("getTransitionIDs", &OpenMS::OSWPeakGroup::getTransitionIDs, nb::rv_policy::reference_internal)
        ;

    // -----------------------------------------------------------------------
    // OSWPeptidePrecursor
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::OSWPeptidePrecursor>(m, "OSWPeptidePrecursor", "A peptide with a charge state")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::OSWPeptidePrecursor&>())
        .def("getSequence", &OpenMS::OSWPeptidePrecursor::getSequence, nb::rv_policy::reference_internal)
        .def("getCharge", &OpenMS::OSWPeptidePrecursor::getCharge)
        .def("isDecoy", &OpenMS::OSWPeptidePrecursor::isDecoy)
        .def("getPCMz", &OpenMS::OSWPeptidePrecursor::getPCMz)
        .def("getFeatures", &OpenMS::OSWPeptidePrecursor::getFeatures, nb::rv_policy::reference_internal)
        ;

    // -----------------------------------------------------------------------
    // OSWProtein
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::OSWProtein>(m, "OSWProtein", "A protein containing one or more peptides")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::OSWProtein&>())
        .def("getAccession", &OpenMS::OSWProtein::getAccession, nb::rv_policy::reference_internal)
        .def("getID", &OpenMS::OSWProtein::getID)
        .def("getPeptidePrecursors", &OpenMS::OSWProtein::getPeptidePrecursors, nb::rv_policy::reference_internal)
        ;

    // -----------------------------------------------------------------------
    // OSWData
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::OSWData>(m, "OSWData", "Holds all or partial information from an OSW file")
        .def(nb::init<>())
        .def("addTransition", [](OpenMS::OSWData& self, const OpenMS::OSWTransition& tr) { self.addTransition(tr); }, "tr"_a)
        .def("addProtein", [](OpenMS::OSWData& self, OpenMS::OSWProtein prot) { self.addProtein(std::move(prot)); }, "prot"_a)
        .def("getProteins", &OpenMS::OSWData::getProteins, nb::rv_policy::reference_internal)
        .def("transitionCount", &OpenMS::OSWData::transitionCount)
        .def("getTransition", &OpenMS::OSWData::getTransition, "id"_a, nb::rv_policy::reference_internal)
        .def("setSqlSourceFile", &OpenMS::OSWData::setSqlSourceFile, "filename"_a)
        .def("getSqlSourceFile", &OpenMS::OSWData::getSqlSourceFile, nb::rv_policy::reference_internal)
        .def("setRunID", &OpenMS::OSWData::setRunID, "run_id"_a)
        .def("getRunID", &OpenMS::OSWData::getRunID)
        .def("clear", &OpenMS::OSWData::clear)
        .def("clearProteins", &OpenMS::OSWData::clearProteins)
        .def("fromNativeID", &OpenMS::OSWData::fromNativeID, "transition_id"_a)
        ;

}