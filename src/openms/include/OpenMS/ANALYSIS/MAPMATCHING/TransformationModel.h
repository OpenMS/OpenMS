// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_TRANSFORMATIONMODEL_H
#define OPENMS_ANALYSIS_MAPMATCHING_TRANSFORMATIONMODEL_H

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/KERNEL/StandardTypes.h>

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
    @brief Check for a valid wighting function string
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

#endif // OPENMS_ANALYSIS_MAPMATCHING_TRANSFORMATIONMODEL_H
