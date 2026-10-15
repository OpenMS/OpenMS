// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Patrick Boschmann, Aditya Sarna $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/IMAGING/IonImage.h>
#include <OpenMS/IMAGING/MSImagingGeometry.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <cmath>
#include <vector>

namespace OpenMS
{
namespace Internal
{

/**
  @brief Build an IonImage by summing peak intensities in an m/z window over a
         pixel subset.

  A single extraction kernel shared by @p MSImagingExperiment and
  @p OnDiscImzMLExperiment. The only substantive difference between the two
  callers is how a spectrum is obtained:

  @code
  // In-memory (no copy): lambda must declare an explicit return type to avoid
  // the silent decay-to-value footgun:
  [this](Size i) -> const MSSpectrum& { return experiment_[i]; }

  // On-disc (value, function-lifetime-extended):
  [this](Size i) -> MSSpectrum { return pimpl_->decodeSpectrum(i); }
  @endcode

  The @c const @c auto& binding inside the loop does the right thing for both
  cases: it binds the reference directly (no copy, in-memory) or
  lifetime-extends the temporary for the duration of the loop body (on-disc).
  A @c std::function<const @c MSSpectrum&(Size)> would @b not work for the
  on-disc case — it would dangle on the decoded value.

  @note The in-memory lambda @b must carry an explicit @c -> @c const
        @c MSSpectrum& return type. An unannotated @c auto-returning lambda
        decays to a by-value return and silently reintroduces a copy on the
        hot in-memory path.

  @tparam GetSpectrum  Callable with signature returning @c MSSpectrum (by
                       value, on-disc) or @c const @c MSSpectrum& (by
                       reference, in-memory). The deduced return category
                       determines the binding in @c const @c auto& @c spec.
  @param[in] geom            Imaging pixel grid (dimensions + pixel list).
  @param[in] mz              m/z center of the extraction window (>= 0).
  @param[in] tolerance_ppm   Half-window in ppm (>= 0).
  @param[in] pixel_indices   Indices into @p geom.getPixels() to process
                             (all pixels or a region subset).
  @param[in] n_spectra       Total number of spectra in the source
                             (used for bounds checking per pixel).
  @param[in] get_spectrum    Callable(spectrum_index) -> spectrum.
  @return Ion image with the same dimensions as @p geom.
  @throws Exception::InvalidValue if @p mz or @p tolerance_ppm is
          negative or non-finite, or if a pixel references an
          out-of-range spectrum_index.
*/
template <typename GetSpectrum>
IonImage extractIonImage(const MSImagingGeometry& geom,
                         double mz, double tolerance_ppm,
                         const std::vector<Size>& pixel_indices,
                         Size n_spectra,
                         GetSpectrum&& get_spectrum)
{
  if (!std::isfinite(mz) || !std::isfinite(tolerance_ppm) || mz < 0.0 || tolerance_ppm < 0.0)
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "mz and tolerance_ppm must be finite and non-negative",
                                  "mz=" + StringUtils::toStr(mz) + ", tolerance_ppm=" + StringUtils::toStr(tolerance_ppm));
  }

  const double dm    = mz * tolerance_ppm * 1e-6;
  const double mz_lo = mz - dm;
  const double mz_hi = mz + dm;

  IonImage image(geom.getWidth(), geom.getHeight());
  image.setMzRange(RangeMZ(mz_lo, mz_hi));

  const auto& pixels = geom.getPixels();
  for (Size i : pixel_indices)
  {
    const auto& p = pixels[i];
    if (p.spectrum_index >= n_spectra)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Pixel references missing spectrum",
                                    StringUtils::toStr(p.spectrum_index));
    }
    // const auto& spec: for a lambda returning const MSSpectrum& (in-memory) this
    // binds the reference with no copy; for a lambda returning MSSpectrum by value
    // (on-disc) this lifetime-extends the temporary for the duration of the loop body.
    const auto& spec = get_spectrum(p.spectrum_index);
    double sum = 0.0;
    for (auto it = spec.MZBegin(mz_lo); it != spec.MZEnd(mz_hi); ++it)
    {
      sum += static_cast<double>(it->getIntensity());
    }
    image.setIntensity(p.x, p.y, sum);
  }
  return image;
}

} // namespace Internal
} // namespace OpenMS
