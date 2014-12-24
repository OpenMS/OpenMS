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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/BaseFeature.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

using namespace std;

namespace OpenMS
{
  const std::string BaseFeature::NamesOfAnnotationState[] = {"no ID", "single ID", "multiple IDs (identical)", "multiple IDs (divergent)"};


  BaseFeature::BaseFeature() :
    RichPeak2D(), quality_(0.0), charge_(0), width_(0), peptides_()
  {
  }

  BaseFeature::BaseFeature(const BaseFeature & rhs) :
    RichPeak2D(rhs), quality_(rhs.quality_), charge_(rhs.charge_), width_(rhs.width_),
    peptides_(rhs.peptides_)
  {
  }

  BaseFeature::BaseFeature(const RichPeak2D & point) :
    RichPeak2D(point), quality_(0.0), charge_(0), width_(0), peptides_()
  {
  }

  BaseFeature::BaseFeature(const Peak2D & point) :
    RichPeak2D(point), quality_(0.0), charge_(0), width_(0), peptides_()
  {
  }

  BaseFeature & BaseFeature::operator=(const BaseFeature & rhs)
  {
    if (&rhs == this)
      return *this;

    RichPeak2D::operator=(rhs);
    quality_ = rhs.quality_;
    charge_ = rhs.charge_;
    width_ = rhs.width_;
    peptides_ =  rhs.peptides_;

    return *this;
  }

  bool BaseFeature::operator==(const BaseFeature & rhs) const
  {
    return RichPeak2D::operator==(rhs)
           && (quality_ == rhs.quality_)
           && (charge_ == rhs.charge_)
           && (width_ == rhs.width_)
           && (peptides_ == rhs.peptides_);
  }

  bool BaseFeature::operator!=(const BaseFeature & rhs) const
  {
    return !operator==(rhs);
  }

  BaseFeature::~BaseFeature()
  {
  }

  BaseFeature::QualityType BaseFeature::getQuality() const
  {
    return quality_;
  }

  void BaseFeature::setQuality(BaseFeature::QualityType quality)
  {
    quality_ = quality;
  }

  BaseFeature::WidthType BaseFeature::getWidth() const
  {
    return width_;
  }

  void BaseFeature::setWidth(BaseFeature::WidthType fwhm)
  {
    // !!! Dirty hack: as long as featureXML doesn't support a width field,
    // we abuse the meta information for this.
    // See also FeatureXMLFile::readFeature_().
    width_ = fwhm;
    setMetaValue("FWHM", fwhm);
  }

  const BaseFeature::ChargeType & BaseFeature::getCharge() const
  {
    return charge_;
  }

  void BaseFeature::setCharge(const BaseFeature::ChargeType & charge)
  {
    charge_ = charge;
  }

  const vector<PeptideIdentification> & BaseFeature::getPeptideIdentifications()
  const
  {
    return peptides_;
  }

  vector<PeptideIdentification> & BaseFeature::getPeptideIdentifications()
  {
    return peptides_;
  }

  void BaseFeature::setPeptideIdentifications(
    const vector<PeptideIdentification> & peptides)
  {
    peptides_ = peptides;
  }


   BaseFeature::AnnotationState BaseFeature::getAnnotationState() const
   {
     if (peptides_.size() == 0) return FEATURE_ID_NONE;
     else if (peptides_.size() == 1 && peptides_[0].getHits().size() > 0) return FEATURE_ID_SINGLE;
     else
     {
       std::set<String> seqs;
       for (Size i=0; i<peptides_.size(); ++i)
       {
         if (peptides_[i].getHits().size() > 0)
         {
           PeptideIdentification id_tmp = peptides_[i];
           id_tmp.sort(); // look at best hit only - requires sorting
           seqs.insert(id_tmp.getHits()[0].getSequence().toString());
         }
       }
       if (seqs.size()==1) return FEATURE_ID_MULTIPLE_SAME; // hits have identical seqs
       else if (seqs.size()>1) return FEATURE_ID_MULTIPLE_DIVERGENT; // multiple different annotations ... probably bad mapping
       else /*if (seqs.size()==0)*/  return FEATURE_ID_NONE; // very rare case of empty hits
     }
   }


} // namespace OpenMS
