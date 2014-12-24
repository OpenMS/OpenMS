// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/Fitter1D.h>

// include derived classes here
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BiGaussFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ExtendedIsotopeFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgFitter1D.h>

#include <OpenMS/CONCEPT/Factory.h>

namespace OpenMS
{
  Fitter1D::Fitter1D() :
    DefaultParamHandler("Fitter1D")
  {
    defaults_.setValue("interpolation_step", 0.2, "Sampling rate for the interpolation of the model function.", ListUtils::create<String>("advanced"));
    defaults_.setValue("statistics:mean", 1.0, "Centroid position of the model.", ListUtils::create<String>("advanced"));
    defaults_.setValue("statistics:variance", 1.0, "The variance of the model.", ListUtils::create<String>("advanced"));
    defaults_.setValue("tolerance_stdev_bounding_box", 3.0, "Bounding box has range [minimim of data, maximum of data] enlarged by tolerance_stdev_bounding_box times the standard deviation of the data.", ListUtils::create<String>("advanced"));

    defaultsToParam_();
  }

  Fitter1D::Fitter1D(const Fitter1D & source) :
    DefaultParamHandler(source)
  {
    setParameters(source.getParameters());
    updateMembers_();
  }

  Fitter1D & Fitter1D::operator=(const Fitter1D & source)
  {
    if (&source == this)
      return *this;

    DefaultParamHandler::operator=(source);
    setParameters(source.getParameters());
    updateMembers_();

    return *this;
  }

  void Fitter1D::registerChildren()
  {
    Factory<Fitter1D>::registerProduct(GaussFitter1D::getProductName(), &GaussFitter1D::create);
    Factory<Fitter1D>::registerProduct(BiGaussFitter1D::getProductName(), &BiGaussFitter1D::create);
    Factory<Fitter1D>::registerProduct(IsotopeFitter1D::getProductName(), &IsotopeFitter1D::create);
    Factory<Fitter1D>::registerProduct(ExtendedIsotopeFitter1D::getProductName(), &ExtendedIsotopeFitter1D::create);
    Factory<Fitter1D>::registerProduct(EmgFitter1D::getProductName(), &EmgFitter1D::create);
  }

  void Fitter1D::updateMembers_()
  {
    tolerance_stdev_box_ = param_.getValue("tolerance_stdev_bounding_box");
    interpolation_step_ = param_.getValue("interpolation_step");
    statistics_.setMean(param_.getValue("statistics:mean"));
    statistics_.setVariance(param_.getValue("statistics:variance"));
  }
  
  Fitter1D::QualityType Fitter1D::fit1d(const RawDataArrayType & /* range */, InterpolationModel * & /* model */)
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
  }

} // namespace OpenMS
