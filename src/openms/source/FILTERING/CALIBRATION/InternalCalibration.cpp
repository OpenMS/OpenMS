// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Maintainer: Chris Bielow $
// $Authors: Alexandra Zerck, Chris Bielow $
// --------------------------------------------------------------------------


#include <OpenMS/FILTERING/CALIBRATION/InternalCalibration.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <stdio.h>

namespace OpenMS
{

  InternalCalibration::InternalCalibration() :
    DefaultParamHandler("InternalCalibration"),
    ProgressLogger(),
    trafo_(TransformationDescription::DataPoints())
  {
    defaults_.setValue("mz_tolerance", 1., "Allowed tolerance between peak and reference m/z.");
    defaults_.setMinFloat("mz_tolerance", 0.);
    defaults_.setValue("mz_tolerance_unit", "Da", "Unit for mz_tolerance.");
    defaults_.setValidStrings("mz_tolerance_unit", ListUtils::create<String>("Da,ppm"));
    defaults_.setValue("rt_tolerance", 10, "Allowed tolerance between peak and reference rt.");

    defaults_.setValue("ms_levels", ListUtils::create<Int>("1,2,3"), "MS level for which to apply the calibration. Make sure to get this correct since it depends on your instrument and the number of mass analyzers it has!");


    defaultsToParam_();
  }

  void InternalCalibration::calibrateMapGlobally(FeatureMap& feature_map,
                                                 const String& trafo_file_name)
  {
    // check if the ids
    checkReferenceIds_(feature_map);
    // first collect theoretical and observed m/z values
    std::vector<double> observed_masses;
    std::vector<double> theoretical_masses;
    for (Size f = 0; f < feature_map.size(); ++f)
    {
      // if more than one peptide id exists for this feature we can't use it as reference
      if (feature_map[f].getPeptideIdentifications().size() != 1) continue;

      Int charge = feature_map[f].getPeptideIdentifications()[0].getHits()[0].getCharge();
      double theo_mass = feature_map[f].getPeptideIdentifications()[0].getHits()[0].getSequence().getMonoWeight(Residue::Full, charge) / (double)charge;
      theoretical_masses.push_back(theo_mass);
      observed_masses.push_back(feature_map[f].getMZ());
#ifdef DEBUG_CALIBRATION
      std::cout << feature_map[f].getRT() << " " << feature_map[f].getMZ() << " " << theo_mass << std::endl;
      std::cout << feature_map[f].getPeptideIdentifications()[0].getHits().size() << std::endl;
      std::cout << feature_map[f].getPeptideIdentifications()[0].getHits()[0].getSequence() << std::endl;
      std::cout << feature_map[f].getPeptideIdentifications()[0].getHits()[0].getCharge() << std::endl;
#endif
    }
    // then make the linear regression
    makeLinearRegression_(observed_masses, theoretical_masses);
    // apply transformation
    applyTransformation_(feature_map);
    if (!trafo_file_name.empty())
    {
      TransformationXMLFile().store(trafo_file_name, trafo_);
    }
  }

  void InternalCalibration::calibrateMapGlobally(FeatureMap& feature_map, std::vector<PeptideIdentification>& ref_ids, const String& trafo_file_name)
  {
    checkReferenceIds_(ref_ids);

    FeatureMap calibrated_feature_map = feature_map;
    // clear the ids
    for (Size f = 0; f < calibrated_feature_map.size(); ++f)
    {
      calibrated_feature_map[f].getPeptideIdentifications().clear();
    }

    // map the reference ids onto the features
    IDMapper mapper;
    Param param;
    param.setValue("rt_tolerance", (double)param_.getValue("rt_tolerance"));
    param.setValue("mz_tolerance", param_.getValue("mz_tolerance"));
    param.setValue("mz_measure", param_.getValue("mz_tolerance_unit"));
    mapper.setParameters(param);
    std::vector<ProteinIdentification> vec;
    mapper.annotate(calibrated_feature_map, ref_ids, vec);

    // calibrate & store calibration function
    calibrateMapGlobally(calibrated_feature_map, trafo_file_name);

    // apply transformation again (which was just used on calibrated_feature_map)
    applyTransformation_(feature_map);
  }

  void InternalCalibration::applyTransformation_(FeatureMap& feature_map)
  {
    for (Size f = 0; f < feature_map.size(); ++f)
    {
      double f_mz = feature_map[f].getMZ();
      feature_map[f].setMZ(trafo_.apply(f_mz));

      // apply transformation to convex hulls and subordinates
      for (Size s = 0; s < feature_map[f].getSubordinates().size(); ++s)
      {
        // subordinates
        double mz = feature_map[f].getSubordinates()[s].getMZ();
        mz = trafo_.apply(mz);
        feature_map[f].getSubordinates()[s].setMZ(mz);
      }
      for (Size s = 0; s < feature_map[f].getConvexHulls().size(); ++s)
      {
        // convex hulls
        std::vector<DPosition<2> > point_vec = feature_map[f].getConvexHulls()[s].getHullPoints();
        feature_map[f].getConvexHulls()[s].clear();
        for (Size p = 0; p < point_vec.size(); ++p)
        {
          double mz = point_vec[p][1];
          mz = trafo_.apply(mz);
          point_vec[p][1] = mz;
        }
        feature_map[f].getConvexHulls()[s].setHullPoints(point_vec);
      }
    }
  }

  void InternalCalibration::checkReferenceIds_(const FeatureMap& feature_map)
  {
    Size num_ids = 0;
    for (Size f = 0; f < feature_map.size(); ++f)
    {
      if (!feature_map[f].getPeptideIdentifications().empty() && feature_map[f].getPeptideIdentifications()[0].getHits().size() > 1)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "InternalCalibration: Your feature map contains PeptideIdentifications with more than one hit, use the IDFilter to select only the best hits before you map the ids to the feature map.");
      }
      else if (!feature_map[f].getPeptideIdentifications().empty())
        ++num_ids;
    }
    if (num_ids < 2)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "InternalCalibration: Your feature map contains less than two PeptideIdentifications, can't perform a linear regression on your data.");
    }
  }

  void InternalCalibration::checkReferenceIds_(const std::vector<PeptideIdentification>& pep_ids)
  {
    for (Size p_id = 0; p_id < pep_ids.size(); ++p_id)
    {
      if (pep_ids[p_id].getHits().size() > 1)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "InternalCalibration: Your Id-file contains PeptideIdentifications with more than one hit, use the IDFilter to select only the best hits.");
      }
      if (!pep_ids[p_id].hasRT())
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, __PRETTY_FUNCTION__, "InternalCalibration: meta data value 'RT' missing for peptide identification!");
      }
      if (!pep_ids[p_id].hasMZ())
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, __PRETTY_FUNCTION__, "InternalCalibration: meta data value 'MZ' missing for peptide identification!");
      }
    }
  }

  void InternalCalibration::makeLinearRegression_(const std::vector<double>& observed_masses, const std::vector<double>& theoretical_masses)
  {
    if (observed_masses.size() != theoretical_masses.size())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Number of observed and theoretical masses must agree.");
    }
#ifdef DEBUG_CALIBRATION
    std::ofstream out("calibration_regression.txt");
    std::vector<double> rel_errors(observed_masses.size(), 0.);
    // determine rel error in ppm for the two reference masses
    for (Size ref_peak = 0; ref_peak < observed_masses.size(); ++ref_peak)
    {
      rel_errors[ref_peak] = (theoretical_masses[ref_peak] - observed_masses[ref_peak]) / theoretical_masses[ref_peak] * 1e6;

      out << observed_masses[ref_peak] << "\t" << rel_errors[ref_peak] << "\n";
      std::cout << observed_masses[ref_peak] << " " << theoretical_masses[ref_peak] << std::endl;
      // std::cout << observed_masses[ref_peak]<<"\t"<<rel_errors[ref_peak]<<std::endl;
    }
#endif

    TransformationDescription::DataPoints data;
    for (Size i = 0; i < observed_masses.size(); ++i)
    {
      data.push_back(std::make_pair(observed_masses[i], theoretical_masses[i]));
    }

    trafo_ = TransformationDescription(data);
    trafo_.fitModel("linear", Param());

#ifdef DEBUG_CALIBRATION
    //          std::cout <<"\n\n---------------------------------\n\n"<< "after calibration "<<std::endl;
    for (Size i = 0; i < observed_masses.size(); ++i)
    {
      double new_mass = trafo_.apply(observed_masses[i]);

      double rel_error = (theoretical_masses[i] - (new_mass)) / theoretical_masses[i] * 1e6;
      std::cout << observed_masses[i] << "\t" << rel_error << std::endl;
    }
#endif
  }

}
