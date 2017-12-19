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
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche, Chris Bielow$
// --------------------------------------------------------------------------

#ifndef OPENMS_SIMULATION_RAWMSSIGNALSIMULATION_H
#define OPENMS_SIMULATION_RAWMSSIGNALSIMULATION_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <OpenMS/SIMULATION/SimTypes.h>
#include <OpenMS/SIMULATION/EGHModel.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ProductModel.h>

namespace OpenMS
{

  class IsotopeModel;

  /**
   @brief Simulates MS signals for a given set of peptides

   Simulates MS signals for a given set of peptides, with charge annotation,
   given detectabilities, predicted retention times and charge values.

   @htmlinclude OpenMS_RawMSSignalSimulation.parameters

   @ingroup Simulation
  */
  class OPENMS_DLLAPI RawMSSignalSimulation :
    public DefaultParamHandler,
    public ProgressLogger
  {

public:
    /** @name Constructors and Destructors
      */
    //@{
    /// Default constructor
    RawMSSignalSimulation();

    /// Constructor taking a random generator
    explicit RawMSSignalSimulation(SimTypes::MutableSimRandomNumberGeneratorPtr rng);

    /// Copy constructor
    RawMSSignalSimulation(const RawMSSignalSimulation& source);

    /// Destructor
    ~RawMSSignalSimulation() override;
    //@}

    RawMSSignalSimulation& operator=(const RawMSSignalSimulation& source);

    /// load the contaminants from contaminants:file param
    /// You do not have to call this function before calling generateRawSignals(), but it might
    /// be useful to check if the contaminant file is valid
    void loadContaminants();

    /// fill experiment with signals and noise
    void generateRawSignals(SimTypes::FeatureMapSim& features,
                            SimTypes::MSSimExperiment& experiment,
                            SimTypes::MSSimExperiment& experiment_ct,
                            SimTypes::FeatureMapSim& contaminants);

protected:

    enum IONIZATIONMETHOD {IM_ESI = 0, IM_MALDI = 1, IM_ALL = 2};
    enum PROFILESHAPE {RT_RECTANGULAR, RT_GAUSSIAN};
    enum RESOLUTIONMODEL {RES_CONSTANT, RES_LINEAR, RES_SQRT};

    /// Synchronize members with param class
    void updateMembers_() override;

    /// Set default parameters
    void setDefaultParams_();

    /**
     @brief Add a 1D signal for a single feature

     @param feature The feature which should be simulated
     @param experiment The experiment to which the simulated signals should be added
     @param experiment_ct Ground truth for picked peaks
     */
    void add1DSignal_(Feature& feature, SimTypes::MSSimExperiment& experiment, SimTypes::MSSimExperiment& experiment_ct);

    /**
     @brief Add a 2D signal for a single feature

     @param feature The feature which should be simulated
     @param experiment The experiment to which the simulated signals should be added
     @param experiment_ct Ground truth for picked peaks
     */
    void add2DSignal_(Feature& feature, SimTypes::MSSimExperiment& experiment, SimTypes::MSSimExperiment& experiment_ct);

    /**
     @brief Samples signals for the given 1D model

     @param iso The isotope model from which the signals will be sampled
     @param mz_start Start coordinate (in m/z dimension) of the region where the signals will be sampled
     @param mz_end End coordinate (in m/z dimension) of the region where the signals will be sampled
     @param experiment Experiment to which the sampled signals will be added
     @param experiment_ct Experiment to which the centroided Ground Truth sampled signals will be added
     @param activeFeature The current feature that is simulated
     */
    void samplePeptideModel1D_(const IsotopeModel& iso,
                               const SimTypes::SimCoordinateType mz_start,
                               const SimTypes::SimCoordinateType mz_end,
                               SimTypes::MSSimExperiment& experiment,
                               SimTypes::MSSimExperiment& experiment_ct,
                               Feature& activeFeature);

    /**
     @brief Samples signals for the given 2D model

     @param pm The product model from which the signals will be sampled
     @param mz_start Start coordinate (in m/z dimension) of the region where the signals will be sampled
     @param mz_end End coordinate (in m/z dimension) of the region where the signals will be sampled
     @param rt_start Start coordinate (in rt dimension) of the region where the signals will be sampled
     @param rt_end End coordinate (in rt dimension) of the region where the signals will be sampled
     @param experiment Experiment to which the sampled signals will be added
     @param experiment_ct Experiment to which the centroided Ground Truth sampled signals will be added
     @param activeFeature The current feature that is simulated
     */
    void samplePeptideModel2D_(const ProductModel<2>& pm,
                               const SimTypes::SimCoordinateType mz_start,
                               const SimTypes::SimCoordinateType mz_end,
                               SimTypes::SimCoordinateType rt_start,
                               SimTypes::SimCoordinateType rt_end,
                               SimTypes::MSSimExperiment& experiment,
                               SimTypes::MSSimExperiment& experiment_ct,
                               Feature& activeFeature);

    /**
     @brief Add the correct Elution profile to the passed ProductModel
     */
    void chooseElutionProfile_(EGHModel* const elutionmodel,
                               Feature& feature,
                               const double scale,
                               const double rt_sampling_rate,
                               const SimTypes::MSSimExperiment& experiment);

    /**
     @brief build contaminant feature map
    */
    void createContaminants_(SimTypes::FeatureMapSim& contaminants, SimTypes::MSSimExperiment& exp, SimTypes::MSSimExperiment& exp_ct);

    /// Add shot noise to the experiment
    void addShotNoise_(SimTypes::MSSimExperiment& experiment, SimTypes::SimCoordinateType minimal_mz_measurement_limit, SimTypes::SimCoordinateType maximal_mz_measurement_limit);

    /// Add white noise to the experiment
    void addWhiteNoise_(SimTypes::MSSimExperiment& experiment);

    /// Add detector noise to the experiment
    void addDetectorNoise_(SimTypes::MSSimExperiment& experiment);

    /// Add a base line to the experiment
    void addBaseLine_(SimTypes::MSSimExperiment& experiment, SimTypes::SimCoordinateType minimal_mz_measurement_limit);

    /// get the mz grid where all m/z values will be mapped to
    void getSamplingGrid_(std::vector<SimTypes::SimCoordinateType>& grid,
                          const SimTypes::SimCoordinateType mz_min,
                          const SimTypes::SimCoordinateType mz_max,
                          const Int step_Da);

    /// Compress signals in a single RT scan (to merge signals which were sampled overlapping)
    void compressSignals_(SimTypes::MSSimExperiment& experiment);

    /// number of points sampled per peak's FWHM
    Int sampling_points_per_FWHM_;

    /// Mean of peak m/z error
    SimTypes::SimCoordinateType mz_error_mean_;
    /// Standard deviation of peak m/z error
    SimTypes::SimCoordinateType mz_error_stddev_;

    /**
     * @brief Computes a rescaled feature intensity based on the set parameters for feature intensity scaling and the passed parameter @p natural_scaling_factor.
     *
     * @param feature_intensity Intensity of the current feature.
     * @param natural_scaling_factor Additional scaling factor used by some of the sampling models.
     *
     * @return Rescaled feature intensity.
     */
    SimTypes::SimIntensityType getFeatureScaledIntensity_(const SimTypes::SimIntensityType feature_intensity,
                                                          const SimTypes::SimIntensityType natural_scaling_factor);


    /**
      @brief Compute resolution at a given m/z given a base resolution and how it degrades with increasing m/z

      @param query_mz The m/z value where the resolution should be estimated
      @param resolution The resolution at 400 Th
      @param model The model describing how resolution behaves, i.e.
                   - RES_CONSTANT: resolution does not change with m/z (this will just return @p resolution)<br>
                   - RES_LINEAR: resolution decreases linear with m/z, i.e. at 800 Th, it will have 50% of original<br>
                   - RES_SQRT: the resolution decreases with square root of mass, i.e. at 1600 Th, it will have 50% of original (sqrt(400) = sqrt(1600)*0.5)

     */
    double getResolution_(const double query_mz, const double resolution, const RESOLUTIONMODEL model) const;

    /**
      @brief compute the peak's SD (Gaussian) at a given m/z (internally the resolution model is used)
    */
    double getPeakWidth_(const double mz, const bool is_gaussian) const;

    /// Scaling factor of peak intensities
    SimTypes::SimIntensityType intensity_scale_;
    /// Standard deviation of peak intensity scaling
    SimTypes::SimIntensityType intensity_scale_stddev_;


    /// model of how resolution behaves with increasing m/z
    RESOLUTIONMODEL res_model_;
    /// base resolution at 400 Th
    double res_base_;
    /// m/z sampling grid for all signals
    std::vector<SimTypes::SimCoordinateType> grid_;

    /// Random number generator
    SimTypes::MutableSimRandomNumberGeneratorPtr rnd_gen_;

    struct ContaminantInfo
    {
      String name;
      EmpiricalFormula sf;
      double rt_start, rt_end, intensity;
      Int q;
      PROFILESHAPE shape;
      IONIZATIONMETHOD im;
    };

    std::vector<ContaminantInfo> contaminants_;

    /**
    @p threaded_random_numbers keeps a set of random numbers for each thread simulating a feature.
      */
    std::vector<std::vector<double> > threaded_random_numbers_;

    /**
      Indicates which random numbers each thread has used already and if the random number pool
      should be rebuild.
      */
    std::vector<Size> threaded_random_numbers_index_;

    static const Size THREADED_RANDOM_NUMBER_POOL_SIZE_ = 500;

    bool contaminants_loaded_;
  };

}

#endif
