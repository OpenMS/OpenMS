// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmMRM.h>


#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ProductModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgModel.h>
#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMeanIterative.h>
#include <OpenMS/FILTERING/TRANSFORMERS/LinearResampler.h>

#include <fstream>

namespace OpenMS
{
  FeatureFinderAlgorithmMRM::FeatureFinderAlgorithmMRM() :
    FeatureFinderAlgorithm()
  {
    defaults_.setValue("min_rt_distance", 10.0, "Minimal distance of MRM features in seconds.");
    defaults_.setMinFloat("min_rt_distance", 0.0);
    defaults_.setValue("min_num_peaks_per_feature", 5, "Minimal number of peaks which are needed for a single feature", ListUtils::create<String>("advanced"));
    defaults_.setMinInt("min_num_peaks_per_feature", 1);
    defaults_.setValue("min_signal_to_noise_ratio", 2.0, "Minimal S/N ratio a peak must have to be taken into account. Set to zero if the MRM-traces contains mostly signals, and no noise.");
    defaults_.setMinFloat("min_signal_to_noise_ratio", 0);
    defaults_.setValue("write_debug_files", "false", "If set to true, for each feature a plot will be created, in the subdirectory 'debug'", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("write_debug_files", ListUtils::create<String>("true,false"));

    defaults_.setValue("resample_traces", "false", "If set to true, each trace, which is in this case a part of the MRM monitoring trace with signal is resampled, using the minimal distance of two data points in RT dimension", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("resample_traces", ListUtils::create<String>("true,false"));

    defaults_.setValue("write_debuginfo", "false", "If set to true, debug messages are written, the output can be somewhat lengthy.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("write_debuginfo", ListUtils::create<String>("true,false"));

    this->defaultsToParam_();
  }

  void FeatureFinderAlgorithmMRM::run()
  {
    //-------------------------------------------------------------------------
    //General initialization
    //-------------------------------------------------------------------------

    Map<Size, Map<Size, std::vector<std::pair<double, Peak1D> > > > traces;

    SignalToNoiseEstimatorMeanIterative<PeakSpectrum> sne;
    LinearResampler resampler;

    // Split the whole map into traces (== MRM transitions)
    ff_->startProgress(0, traces.size(), "Finding features in traces.");
    Size counter(0);
    //typename Map<Size, Map<Size, std::vector<std::pair<double, Peak1D> > > >::const_iterator it1 = traces.begin();
    //typename Map<Size, std::vector<std::pair<double, Peak1D> > >::const_iterator it2;
    double min_rt_distance(param_.getValue("min_rt_distance"));
    double min_signal_to_noise_ratio(param_.getValue("min_signal_to_noise_ratio"));
    Size min_num_peaks_per_feature(param_.getValue("min_num_peaks_per_feature"));
    Size feature_id(0);
    bool write_debuginfo(param_.getValue("write_debuginfo").toBool());
    bool write_debug_files(param_.getValue("write_debug_files").toBool());
    bool resample_traces(param_.getValue("resample_traces").toBool());

    if (write_debuginfo)
    {
      std::cerr << "Starting feature finding #chromatograms=" << map_->getChromatograms().size() << ", #spectra=" << map_->size() << std::endl;
    }

    std::vector<MSChromatogram >::const_iterator first_it = map_->getChromatograms().begin();
    for (; first_it != map_->getChromatograms().end(); ++first_it)
    {
      // throw the peaks into a "spectrum" where the m/z values are RTs in reality (more a chromatogram)
      PeakSpectrum chromatogram;
      //typename std::vector<std::pair<double, Peak1D> >::const_iterator it3 = it2->second.begin();
      for (MSChromatogram::const_iterator it = first_it->begin(); it != first_it->end(); ++it)
      {
        Peak1D peak;
        peak.setMZ(it->getRT());
        peak.setIntensity(it->getIntensity());
        chromatogram.push_back(peak);
      }

      // TODO
      // find pre-section of separated RTs of peaks;
      // resampling to min distance?
      // for each of the section, try to estimate S/N
      // find core regions and fit them
      //

      if (resample_traces)
      {
        // resample the chromatogram, first find minimal distance and use this as resampling distance
        double min_distance(std::numeric_limits<double>::max()), old_rt(0);
        for (PeakSpectrum::ConstIterator it = chromatogram.begin(); it != chromatogram.end(); ++it)
        {
          if (write_debuginfo)
          {
            std::cerr << "CHROMATOGRAM: " << it->getMZ() << " " << it->getIntensity() << std::endl;
          }
          double rt_diff = it->getMZ() - old_rt;
          if (rt_diff < min_distance && rt_diff > 0)
          {
            min_distance = rt_diff;
          }
          old_rt = it->getMZ();
        }

        if (write_debuginfo)
        {
          std::cerr << "Min_distance=" << min_distance << std::endl;
        }
        if (min_distance > 50 || chromatogram.size() < min_num_peaks_per_feature)
        {
          continue;
        }
        Param resampler_param(resampler.getParameters());
        resampler_param.setValue("spacing", min_distance);
        resampler.setParameters(resampler_param);
        resampler.raster(chromatogram);
      }

      // now smooth the data
      GaussFilter filter;
      Param filter_param(filter.getParameters());
      filter.setParameters(filter_param);

      // calculate signal to noise levels
      PeakSpectrum sn_chrom;
      Param sne_param(sne.getParameters());
      // set window length to whole range, we expect only at most one signal
      if (write_debuginfo)
      {
        std::cerr << "win_len m/z: " <<  (chromatogram.end() - 1)->getMZ() << " " << chromatogram.begin()->getMZ() << std::endl;
      }
      sne_param.setValue("win_len", (chromatogram.end() - 1)->getMZ() - chromatogram.begin()->getMZ());

      if ((double)sne_param.getValue("win_len") < 10e-4)
      {
        continue;
      }
      sne.setParameters(sne_param);
      sne.init(chromatogram);

      if (write_debuginfo)
      {
        std::cerr << first_it->getPrecursor().getMZ() << " " << first_it->getProduct().getMZ() << " ";
      }

      PeakSpectrum::FloatDataArray signal_to_noise;
      for (Size i = 0; i < chromatogram.size(); ++i)
      {
        double sn(sne.getSignalToNoise(i));
        signal_to_noise.push_back(sn);
        if (write_debuginfo)
        {
          std::cerr << chromatogram[i].getMZ() << " " << chromatogram[i].getIntensity() << " " << sn << std::endl;
        }
        if (min_signal_to_noise_ratio == 0 || sn > min_signal_to_noise_ratio)
        {
          sn_chrom.push_back(chromatogram[i]);
        }
      }
      chromatogram.getFloatDataArrays().push_back(signal_to_noise);

      // now find sections in the chromatogram which have high s/n value
      double last_rt(0);
      std::vector<std::vector<DPosition<2> > > sections;
      for (PeakSpectrum::Iterator sit = sn_chrom.begin(); sit != sn_chrom.end(); ++sit)
      {
        if (write_debuginfo)
        {
          std::cerr << "SECTIONS: " << sit->getMZ() << " " << sit->getIntensity() << std::endl;
        }
        double this_rt = sit->getMZ();
        if (sections.empty() || (this_rt - last_rt) > min_rt_distance)
        {
          if (write_debuginfo)
          {
            std::cerr << "Starting new section, sections.size()=" << sections.size() << ", rt_diff=" << this_rt - last_rt << std::endl;
          }
          // new section
          std::vector<DPosition<2> > section;
          section.push_back(DPosition<2>(this_rt, sit->getIntensity()));
          sections.push_back(section);
        }
        else
        {
          sections.back().push_back(DPosition<2>(this_rt, sit->getIntensity()));
        }
        last_rt = this_rt;
      }

      // for each of the sections identify local maxima which can be used for extension
      // if needed split the sections into smaller sections
      for (Size i = 0; i != sections.size(); ++i)
      {
        if (sections[i].size() < min_num_peaks_per_feature)
        {
          continue;
        }

        std::vector<double> deltas(2, 0.0);
        double last_int(sections[i][0].getY());
        for (Size j = 1; j < sections[i].size(); ++j)
        {
          deltas.push_back((sections[i][j].getY() - last_int) / last_int);
          last_int = sections[i][j].getY();
          if (write_debuginfo)
          {
            double average_delta = std::accumulate(deltas.begin(), deltas.end(), 0.0) / (double)3.0;
            std::cerr << "AverageDelta: " << average_delta << " (" << sections[i][j].getX() << ", " << sections[i][j].getY() << ")" << std::endl;
          }
        }
      }


      // for each section estimate the rt min/max and add up the intensities
      for (Size i = 0; i != sections.size(); ++i)
      {
        if (sections[i].size() > min_num_peaks_per_feature)
        {
          // first smooth the data to prevent outliers from destroying the fit
          PeakSpectrum filter_spec;
          for (Size j = 0; j != sections[i].size(); ++j)
          {
            Peak1D p;
            p.setMZ(sections[i][j].getX());
            p.setIntensity(sections[i][j].getY());
            filter_spec.push_back(p);
          }

          // add two peaks at the beginning and at the end for better fit
          // therefore calculate average distance first
          std::vector<double> distances;
          for (Size j = 1; j < filter_spec.size(); ++j)
          {
            distances.push_back(filter_spec[j].getMZ() - filter_spec[j - 1].getMZ());
          }
          double dist_average = std::accumulate(distances.begin(), distances.end(), 0.0) / (double)distances.size();

          // append peaks
          Peak1D new_peak;
          new_peak.setIntensity(0);
          new_peak.setMZ(filter_spec.back().getMZ() + dist_average);
          filter_spec.push_back(new_peak);
          new_peak.setMZ(filter_spec.back().getMZ() + dist_average);
          filter_spec.push_back(new_peak);
          new_peak.setMZ(filter_spec.back().getMZ() + dist_average);
          filter_spec.push_back(new_peak);

          // prepend peaks
          new_peak.setMZ(filter_spec.front().getMZ() - dist_average);
          filter_spec.insert(filter_spec.begin(), new_peak);
          new_peak.setMZ(filter_spec.front().getMZ() - dist_average);
          filter_spec.insert(filter_spec.begin(), new_peak);
          new_peak.setMZ(filter_spec.front().getMZ() - dist_average);
          filter_spec.insert(filter_spec.begin(), new_peak);

          filter_param.setValue("gaussian_width", 4 * dist_average);
          filter.setParameters(filter_param);

          // smooth the data
          filter.filter(filter_spec);

          // transform the data for fitting and fit RT profile
          std::vector<Peak1D> data_to_fit;
          for (Size j = 0; j != filter_spec.size(); ++j)
          {
            Peak1D p;
            p.setPosition(filter_spec[j].getMZ());
            p.setIntensity(filter_spec[j].getIntensity());
            data_to_fit.push_back(p);
          }
          InterpolationModel* model_rt = nullptr;
          double quality = fitRT_(data_to_fit, model_rt);

          Feature f;
          f.setQuality(0, quality);
          f.setOverallQuality(quality);

          ConvexHull2D::PointArrayType hull_points(sections[i].size());
          double intensity_sum(0.0);
          for (Size j = 0; j < sections[i].size(); ++j)
          {
            hull_points[j][0] = sections[i][j].getX();
            hull_points[j][1] = first_it->getProduct().getMZ();
            intensity_sum += sections[i][j].getY();
          }

          // create the feature according to fit
          f.setRT((double)model_rt->getParameters().getValue("emg:retention"));
          f.setMZ((double)first_it->getProduct().getMZ());
          f.setIntensity(intensity_sum);
          ConvexHull2D hull;
          hull.addPoints(hull_points);
          f.getConvexHulls().push_back(hull);
          f.setMetaValue("MZ", (double)first_it->getPrecursor().getMZ());

          ++feature_id;

          // writes a feature plot using gnuplot (should be installed on computer)
          if (write_debug_files)
          {
            String base_name = "debug/" + String((double)f.getMetaValue("MZ")) + "_id" + String(feature_id) + "_RT" + String(f.getRT()) + "_Q3" + String(f.getMZ());
            std::ofstream data_out(String(base_name + "_data.dat").c_str());
            for (Size j = 0; j < sections[i].size(); ++j)
            {
              // RT intensity
              double rt = sections[i][j].getX();
              double intensity = sections[i][j].getY();
              data_out << rt << " " << intensity << std::endl;
            }
            data_out.close();

            std::ofstream smoothed_data_out(String(base_name + "_smoothed_data.dat").c_str());
            for (Size j = 0; j < filter_spec.size(); ++j)
            {
              smoothed_data_out << filter_spec[j].getMZ() + 0.5 << " " << filter_spec[j].getIntensity() << std::endl;
            }
            smoothed_data_out.close();

            std::ofstream fit_out(String(base_name + "_rt_fit.dat").c_str());
            EmgModel emg_model;
            emg_model.setParameters(model_rt->getParameters());
            emg_model.setSamples();
            double bb_min((double)emg_model.getParameters().getValue("bounding_box:min"));
            double bb_max((double)emg_model.getParameters().getValue("bounding_box:max"));
            double int_step((double)emg_model.getParameters().getValue("interpolation_step"));
            for (double pos = bb_min; pos < bb_max; pos += int_step)
            {
              // RT intensity
              fit_out << pos << " " << emg_model.getIntensity(pos) << std::endl;
            }
            fit_out.close();

            std::ofstream gnuplot_out(String(base_name + "_gnuplot.gpl").c_str());
            gnuplot_out << "set terminal png" << std::endl;
            gnuplot_out << "set output \"" << base_name << ".png\"" << std::endl;
            gnuplot_out << "plot '" << base_name << "_data.dat' w i, '" << base_name << "_smoothed_data.dat'  w i, '" << base_name << "_rt_fit.dat' w lp title 'quality=" << f.getOverallQuality() << "'" << std::endl;
            gnuplot_out.close();
            String gnuplot_call = "gnuplot " + base_name + "_gnuplot.gpl";
            int error = system(gnuplot_call.c_str());
            if (error != 0)
            {
              std::cerr << "An error occurred during the gnuplot execution" << std::endl;
            }
          }

          features_->push_back(f);
        }
      }

      ff_->setProgress(++counter);
    }
  }

  FeatureFinderAlgorithm* FeatureFinderAlgorithmMRM::create()
  {
    return new FeatureFinderAlgorithmMRM();
  }

  const String FeatureFinderAlgorithmMRM::getProductName()
  {
    return "mrm";
  }

  double FeatureFinderAlgorithmMRM::fitRT_(std::vector<Peak1D>& rt_input_data, InterpolationModel*& model) const
  {
    double quality;
    Param param;
    EmgFitter1D fitter;

    /*
param.setValue( "tolerance_stdev_bounding_box", tolerance_stdev_box_);
param.setValue( "statistics:mean", rt_stat_.mean() );
param.setValue( "statistics:variance", rt_stat_.variance() );
param.setValue( "interpolation_step", interpolation_step_rt_ );
param.setValue( "max_iteration", max_iteration_);
param.setValue( "deltaAbsError", deltaAbsError_);
param.setValue( "deltaRelError", deltaRelError_);
    */

    // Set parameter for fitter
    fitter.setParameters(param);

    // Construct model for rt
    quality = fitter.fit1d(rt_input_data, model);

    // Check quality
    if (std::isnan(quality)) quality = -1.0;

    return quality;
  }

  void FeatureFinderAlgorithmMRM::updateMembers_()
  {
  }

}
