// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/Spectrum2DGoToDialog.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <QtGui/QLineEdit>

using namespace std;

namespace OpenMS
{

	Spectrum2DGoToDialog::Spectrum2DGoToDialog(QWidget* parent)
		: QDialog(parent)
	{
		setupUi(this);
	}
	
	Spectrum2DGoToDialog::~Spectrum2DGoToDialog()
	{
	}
	
	void Spectrum2DGoToDialog::setRange(Real min_rt, Real max_rt, Real min_mz, Real max_mz)
	{
		min_rt_->setText(QString::number(min_rt));
		max_rt_->setText(QString::number(max_rt));
		min_mz_->setText(QString::number(min_mz));
		max_mz_->setText(QString::number(max_mz));
	}
	
	Real Spectrum2DGoToDialog::getMinRT() const
	{
	  return min_rt_->text().toFloat();
	}
	
	Real Spectrum2DGoToDialog::getMaxRT() const
	{
	  return max_rt_->text().toFloat();  
	}
	
	Real Spectrum2DGoToDialog::getMinMZ() const
	{
	  return min_mz_->text().toFloat();  
	}
	
	Real Spectrum2DGoToDialog::getMaxMZ() const
	{
	  return max_mz_->text().toFloat();  
	}

	void Spectrum2DGoToDialog::enableFeatureNumber(bool enabled)
	{
		feature_label_->setEnabled(enabled);
		nr_->setEnabled(enabled);
		feature_number_->setEnabled(enabled);
		//Reorder tab order
		if (enabled)
		{
			setTabOrder(feature_number_, ok_button_);
			setTabOrder(ok_button_, cancel_button_);
			setTabOrder(cancel_button_, min_mz_);
			setTabOrder(min_mz_, max_mz_);
			setTabOrder(max_mz_, min_rt_);
			setTabOrder(min_rt_, max_rt_);
		}
		else
		{
			setTabOrder(min_mz_, max_mz_);
			setTabOrder(max_mz_, min_rt_);
			setTabOrder(min_rt_, max_rt_);
			setTabOrder(max_rt_, ok_button_);
			setTabOrder(ok_button_, cancel_button_);
		}
	}
	
	String Spectrum2DGoToDialog::getFeatureNumber() const
	{
		return feature_number_->text(); 
	}

	bool Spectrum2DGoToDialog::showRange() const
	{
		if (feature_number_->text().trimmed()!="")
		{
			return false;
		}
		return true;
	}
	
} //namespace OpenMS
