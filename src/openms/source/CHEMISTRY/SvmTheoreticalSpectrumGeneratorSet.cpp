// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Authors: Sandro Andreotti $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/SvmTheoreticalSpectrumGeneratorSet.h>


namespace OpenMS
{

  // Default constructor
  SvmTheoreticalSpectrumGeneratorSet::SvmTheoreticalSpectrumGeneratorSet()
  {
  }

  // Copy constructor
  SvmTheoreticalSpectrumGeneratorSet::SvmTheoreticalSpectrumGeneratorSet(const SvmTheoreticalSpectrumGeneratorSet& source) :
    simulators_(source.simulators_)
  {
  }

  //Destructor
  SvmTheoreticalSpectrumGeneratorSet::~SvmTheoreticalSpectrumGeneratorSet()
  {
  }

  // Assignment operator
  SvmTheoreticalSpectrumGeneratorSet& SvmTheoreticalSpectrumGeneratorSet::operator=(const SvmTheoreticalSpectrumGeneratorSet& rhs)
  {
    if (this != &rhs)
    {
      simulators_ = rhs.simulators_;
    }
    return *this;
  }

  // Generate the MS/MS according to the given probabilistic model
  void SvmTheoreticalSpectrumGeneratorSet::simulate(PeakSpectrum& spectrum, const AASequence& peptide, boost::random::mt19937_64& rng, Size precursor_charge)
  {
    auto it = simulators_.find(precursor_charge);
    if (it != std::end(simulators_))
    {
      it->second.simulate(spectrum, peptide, rng, precursor_charge);
    }
    else
    {
      throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid Precursor charge, no Model available", String(precursor_charge));
    }
  }

  //Load a trained Svm and Prob. models
  void SvmTheoreticalSpectrumGeneratorSet::load(String filename)
  {
    if (!File::readable(filename)) // look in OPENMS_DATA_PATH
    {
      filename = File::find(filename);
    }

    TextFile file(filename);
    // Boost::empty()
    if (std::begin(file) == std::end(file)) { return; }// no data to load

    // skip header line
    const auto begin_from_skipped_header = std::next(std::begin(file));

    Param sim_param = SvmTheoreticalSpectrumGenerator().getDefaults();

    const auto process_content = [this, filename, &sim_param] (TextFile::Iterator::value_type i)
    {
        std::vector<String> spl;
        i.split(":", spl);
        Int precursor_charge = spl[0].toInt();

        if (spl.size() != 2 || precursor_charge < 1)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, i, " Invalid entry in SVM model File");
        }

        //load the model into the map
        sim_param.setValue("model_file_name", File::path(filename) + "/" + spl[1]);
        simulators_[precursor_charge].setParameters(sim_param);
        simulators_[precursor_charge].load();
    };

    // process content
    std::for_each(begin_from_skipped_header, std::end(file), process_content);
  }

  //Return precursor charges for which a model is contained in the set
  void SvmTheoreticalSpectrumGeneratorSet::getSupportedCharges(std::set<Size>& charges)
  {
    charges.clear();

    std::transform(std::begin(simulators_), std::end(simulators_),
                   std::inserter(charges, std::begin(charges)),
                   [] (std::pair<Size, SvmTheoreticalSpectrumGenerator> it) { return it.first; });
  }

  //return a modifiable reference to the SVM model with given charge. If charge is not supported throw exception
  SvmTheoreticalSpectrumGenerator& SvmTheoreticalSpectrumGeneratorSet::getSvmModel(Size prec_charge)
  {
    auto it = simulators_.find(prec_charge);
    if (it != std::end(simulators_))
    {
      return it->second;
    }

    throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid Precursor charge, no Model available", String(prec_charge));
  }

}
