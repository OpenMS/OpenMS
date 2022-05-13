// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/VISUAL/LayerDataBase.h>
#include <OpenMS/VISUAL/INTERFACES/IPeptideIds.h>

namespace OpenMS
{

  /**
  @brief Class that stores the data for one layer of type PeptideIdentifications

  @ingroup PlotWidgets
  */
  class OPENMS_GUI_DLLAPI LayerDataIdent : public LayerDataBase, public IPeptideIds
  {
  public:
    /// Default constructor
    LayerDataIdent() :
        LayerDataBase(LayerDataBase::DT_IDENT){};
    /// no Copy-ctor (should not be needed)
    LayerDataIdent(const LayerDataIdent& ld) = delete;
    /// no assignment operator (should not be needed)
    LayerDataIdent& operator=(const LayerDataIdent& ld) = delete;
    /// move Ctor
    LayerDataIdent(LayerDataIdent&& ld) = default;
    /// move assignment
    LayerDataIdent& operator=(LayerDataIdent&& ld) = default;

    void updateRanges() override
    {
      // nothing to do...
    }

    RangeAllType getRange() const override
    {
      RangeAllType r;
      for (const PeptideIdentification& pep : peptides_)
      {
        r.extendRT(pep.getRT());
        r.extendMZ(pep.getMZ());
      }
      return r;
    }

    std::unique_ptr<LayerStatistics> getStats() const override;

    virtual const PepIds& getPeptideIds() const override
    {
      return peptides_;
    }
    virtual PepIds& getPeptideIds() override
    {
      return peptides_;
    }

    virtual void setPeptideIds(const PepIds& ids) override
    {
      peptides_ = ids;
    }
    virtual void setPeptideIds(PepIds&& ids) override
    {
      peptides_ = std::move(ids);
    }

  private:
    /// peptide identifications
    std::vector<PeptideIdentification> peptides_;
  };

}// namespace OpenMS
