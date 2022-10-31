// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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

#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/VISUAL/LayerDataBase.h>
#include <OpenMS/VISUAL/Plot1DWidget.h>
#include <OpenMS/VISUAL/Plot2DWidget.h>
#include <OpenMS/openms_data_path.h>
#include <QApplication>

using namespace OpenMS;
using namespace std;

Int main(int argc, const char** argv)
{
  String tutorial_data_path(OPENMS_DOC_PATH + String("/code_examples/data/Tutorial_Spectrum1D.dta"));
  if (argc >= 2)
  { // the path to the data can be given on the command line
    tutorial_data_path = argv[1];
  }

  QApplication app(argc, const_cast<char**>(argv));

  PeakMap exp;
  exp.resize(1);
  DTAFile().load(tutorial_data_path, exp[0]);
  LayerDataBase::ExperimentSharedPtrType exp_sptr(new PeakMap(exp));
  LayerDataBase::ODExperimentSharedPtrType on_disc_exp_sptr(new OnDiscMSExperiment());
  auto* widget = new Plot1DWidget(Param(), DIM::Y, nullptr);
  widget->canvas()->addPeakLayer(exp_sptr, on_disc_exp_sptr);
  widget->show();

  return app.exec();
} // end of main
