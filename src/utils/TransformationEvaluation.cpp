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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_TransformationEvaluation TransformationEvaluation

    @brief Applies a transformation to a range of values and records the results.

    This is useful for plotting transformations for quality assessment etc.

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_TransformationEvaluation.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_TransformationEvaluation.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPTransformationEvaluation :
  public TOPPBase
{
public:

  TOPPTransformationEvaluation() :
    TOPPBase("TransformationEvaluation", "Applies a transformation to a range of values", false)
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file containing the transformation description");
    setValidFormats_("in", ListUtils::create<String>("trafoXML"));
    registerOutputFile_("out", "<file>", "", "Output file containing original and transformed values; if empty, output is written to the screen", false);
    setValidFormats_("out", ListUtils::create<String>("trafoXML"));
    registerDoubleOption_("min", "<value>", 0.0, "Minimum value to transform", false);
    registerDoubleOption_("max", "<value>", 0.0, "Maximum value to transform (if at or below 'min', select a suitable maximum based on the transformation description)", false);
    registerDoubleOption_("step", "<value>", 1.0, "Step size between 'min' and 'max'", false);
    setMinFloat_("step", 0.001);
  }

  ExitCodes main_(int, const char **) override
  {
    String in = getStringOption_("in"), out = getStringOption_("out");

    TransformationDescription trafo_in;
    TransformationXMLFile().load(in, trafo_in);
    TransformationDescription::DataPoints data;

    double min = getDoubleOption_("min"), max = getDoubleOption_("max"),
               step = getDoubleOption_("step");
    if (max <= min)
    {
      data = trafo_in.getDataPoints();
      sort(data.begin(), data.end());
      max = data.back().first;
      double magnitude = floor(log10(max));
      max = Math::ceilDecimal(max, magnitude - 1);
      if (max <= min)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "'min' must be lower than 'max'");
      }
    }

    data.clear();
    for (double value = min; value <= max; value += step)
    {
      double transformed = trafo_in.apply(value);
      if (out.empty())
      {
        cout << value << '\t' << transformed << endl;
      }
      else data.push_back(make_pair(value, transformed));
    }

    if (!out.empty())
    {
      TransformationDescription trafo_out(trafo_in);
      trafo_out.setDataPoints(data);
      TransformationXMLFile().store(out, trafo_out);
    }

    return EXECUTION_OK;
  }

};


int main(int argc, const char ** argv)
{
  TOPPTransformationEvaluation tool;
  return tool.main(argc, argv);
}

/// @endcond
