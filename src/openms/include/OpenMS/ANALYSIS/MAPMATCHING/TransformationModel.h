// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <tuple>

namespace OpenMS
{
  /**
    @brief Base class for transformation models

    Implements the identity (no transformation). Parameters and data are ignored.

    Note that this class and its derived classes do not allow copying/assignment, due to the need for internal memory management associated with some of the transformation models.

    @ingroup MapAlignment
  */
  class OPENMS_DLLAPI TransformationModel
  {
  public:
    /// Coordinate pair (with optional annotation)
    struct DataPoint
    {
      double first, second;
      String note;

      DataPoint(double first_ = 0.0,
                double second_ = 0.0,
                const String& note_ = "") :
        first(first_),
        second(second_),
        note(note_)
      {}

      DataPoint(const std::pair<double, double>& pair) :
        first(pair.first),
        second(pair.second),
        note("")
      {}

      bool operator<(const DataPoint& other) const
      {
        return (std::tie(first, second, note) <
                std::tie(other.first, other.second, other.note));
      }

      bool operator==(const DataPoint& other) const
      {
        return (std::tie(first, second, note) ==
                std::tie(other.first, other.second, other.note));
      }
    };

    /// Vector of coordinate pairs
    typedef std::vector<DataPoint> DataPoints;

    /// Constructor
    TransformationModel() {}

    /// Alternative constructor (derived classes should implement this one!)
    /// Both data and params must be provided, since some derived classes require both to create a model!
    TransformationModel(const TransformationModel::DataPoints&, const Param&);

    /// Destructor
    virtual ~TransformationModel();

    /// Evaluates the model at the given value
    virtual double evaluate(double value) const;
    
    /**
    @brief Weight the data by the given weight function

    Currently supported valid weighting functions include the following:
      - 1 / x.
      - 1 / x2.
      - 1 / y.
      - 1 / y2.
      - ln x.
      - ln y.
    Note that the user needs to ensure valid bounds for the data by setting
    the x_datum_min/x_datum_max and y_datum_min/y_datum_max params.
    */
    virtual void weightData(DataPoints& data);
     
    /**
    @brief Unweight the data by the given weight function
    */
    virtual void unWeightData(DataPoints& data);
    
    /**
    @brief Check for a valid weighting function string
    */
    bool checkValidWeight(const String& weight, const std::vector<String>& valid_weights) const;

    /**
    @brief Check that the datum is within the valid min and max bounds.

    The method checks if the datum is within the user specified min and max bounds.
    If the datum is below the min bounds, the min bound is returned.
    If the datum is above the max bounds, the max bound is returned.
    */
    double checkDatumRange(const double& datum, const double& datum_min, const double& datum_max);
 
    /**
    @brief Weight the data according to the weighting function
    */
    double weightDatum(const double& datum, const String& weight) const;

    /**
    @brief Apply the reverse of the weighting function to the data
    */
    double unWeightDatum(const double& datum, const String& weight) const;

    /// Gets the (actual) parameters
    const Param& getParameters() const;

    /// Returns a list of valid x weight function strings
    std::vector<String> getValidXWeights() const;

    /// Returns a list of valid y weight function strings
    std::vector<String> getValidYWeights() const;

    /// Gets the default parameters
    static void getDefaultParameters(Param& params);

  protected:
    /// Parameters
    Param params_;
    /// x weighting
    String x_weight_;
    double x_datum_min_;
    double x_datum_max_;
    /// y weighting
    String y_weight_;
    double y_datum_min_;
    double y_datum_max_;
    bool weighting_;

  private:
    /// do not allow copy
    TransformationModel(const TransformationModel&);
    /// do not allow assignment
    const TransformationModel& operator=(const TransformationModel&);

  };

} // end of namespace OpenMS

