// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Author: Timo Sachsenberg, Mohammed Alhigaylan $
// $Maintainer: Timo Sachsenberg $
// -------------------------------------------------------------------------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

namespace OpenMS
{
  /**
    @brief Peak picking algorithm for ion mobility data

  @ingroup PeakPicking
      */
  class OPENMS_DLLAPI PeakPickerIM : public DefaultParamHandler
  {
  public:
    /// Default constructor initializing parameters with default values.
    PeakPickerIM();

    /// Destructor.
    ~PeakPickerIM() override = default;

    /// Picks ion mobility traces from the given spectrum.
    void pickIMTraces(MSSpectrum& spectrum);

    /// Sets the parameters for peak picking.
    using DefaultParamHandler::setParameters;
    using DefaultParamHandler::getParameters;

    /// Gets the current parameters.
    //Param getParameters() const;

        /**
     * @brief Converts an ion mobility frame to a single spectrum with averaged IM values
     *
     * This function takes an MS spectrum containing ion mobility data and reduces it to
     * a single spectrum where peaks that are close in both m/z and ion mobility space
     * are averaged together. The averaging is intensity-weighted for both m/z and ion
     * mobility values.
     *
     * The algorithm processes peaks sequentially and groups them based on two criteria:
     * 1. m/z tolerance: peaks must be within specified ppm of each other
     * 2. ion mobility tolerance: the range of IM values must not exceed the specified tolerance
     *
     * @param input Spectrum containing ion mobility data in its FloatDataArrays
     * @param ppm_tolerance_cluster Mass tolerance in parts per million (default: 50.0 ppm)
     * @param im_tolerance_cluster Ion mobility tolerance (default: 0.1 units)
     *
     * @throws Exception::MissingInformation if input spectrum lacks ion mobility data
     *
     * @note The input spectrum should contain ion mobility data in its FloatDataArrays.
     *       The output spectrum will contain averaged peaks with their corresponding
     *       intensity-weighted average ion mobility values.
     *
     * Example:
     * @code
     * MSSpectrum input;  // spectrum with IM data
     * IMFrame::toSpectrum(input_spectrum);
     * @endcode
     */
    static void pickIMCluster(MSSpectrum& spec, double ppm_tolerance_cluster = 50.0, double im_tolerance_cluster = 0.1);

    /**
     * @brief Picks ion mobility elution profiles from the given spectrum using eluting profiles.
     *
     * This function processes an MS spectrum containing ion mobility data and
     * extracts IM elution profiles based on the specified ppm tolerance.
     *
     * @param input Spectrum containing ion mobility data in its FloatDataArrays
     * @param ppm_tolerance_elution Mass tolerance in parts per million (default: 50.0 ppm)
     */
    static void pickIMElutionProfiles(MSSpectrum& input, double ppm_tolerance_elution = 50.0);

  protected:
    /// Updates internal member variables when parameters are changed.
    void updateMembers_() override;

    /// Returns the default parameters.
    Param getDefaultParameters() const;

    /// Stores the parameters for peak picking.
    Param parameters_;

  private:
    /// determine sampling rate for linear resampler
    double computeOptimalSamplingRate(const std::vector<MSSpectrum>& spectra);

    /**
     * @brief Sum up the intensity of data points with nearly identical float values.
     *
     * By default, this function assumes the tolerance provided is in parts per million.
     * But it can be adjusted to use absolute value tolerance.
     * @param input Sorted raw spectrum with duplicate peaks due to scan merging or presence of ion mobility data.
     * @param tolerance Mass tolerance between peaks
     * * @param use_ppm Whether to use parts per million tolerance. If set to False, absolute tolerance will be used.
     */
    void sumFrame_(const MSSpectrum& input_spectrum, MSSpectrum& output_spectrum, double tolerance = 0.01, bool use_ppm = true);

    /// Compute lower and upper m/z bounds based on ppm
    std::pair<double, double> ppmBounds(double mz, double ppm);

    /// Extract ion mobility traces as MSSpectra from the raw TimsTOF frame
    /// Ion mobility is temporarily written in place of m/z inside Peak1D object.
    /// raw m/z values are allocated to float data arrays with the label 'raw_mz'
    std::vector<MSSpectrum> extractIonMobilityTraces(
      const MSSpectrum& picked_spectrum,
      const MSSpectrum& raw_spectrum);

    /// compute m/z and ion mobility centers for picked traces. Returns centroided spectrum.
    MSSpectrum computeCentroids_(const std::vector<MSSpectrum>& mobilogram_traces,
                              const std::vector<MSSpectrum>& picked_traces);

    // PickIMTraces parameters
    // cached values (no Param lookups inside hot loops)
    double sum_tolerance_mz_{1.0};
    double gauss_ppm_tolerance_{5.0};
    double sum_tolerance_im_{0.0006};
    int    sgolay_frame_length_{5};
    int    sgolay_polynomial_order_{3};
    double    ppm_tolerance_cluster_{50.0};
    double    im_tolerance_cluster_{0.1};

    double    ppm_tolerance_elution_{50.0};
  };

} // namespace OpenMS