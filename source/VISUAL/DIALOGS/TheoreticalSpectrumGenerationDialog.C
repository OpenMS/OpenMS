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
#include <OpenMS/VISUAL/DIALOGS/TheoreticalSpectrumGenerationDialog.h>

namespace OpenMS
{
	TheoreticalSpectrumGenerationDialog::TheoreticalSpectrumGenerationDialog()
	{
		setupUi(this);
		
		connect(list_widget, SIGNAL(itemChanged(QListWidgetItem*)), this, SLOT(itemChanged(QListWidgetItem*)));
		
		// select b- and y-ions as residue types by default
		list_widget->item(0)->setCheckState(Qt::Unchecked);
		list_widget->item(1)->setCheckState(Qt::Checked);
		list_widget->item(2)->setCheckState(Qt::Unchecked);
		list_widget->item(3)->setCheckState(Qt::Unchecked);
		list_widget->item(4)->setCheckState(Qt::Checked);
		list_widget->item(5)->setCheckState(Qt::Unchecked);
		list_widget->item(6)->setCheckState(Qt::Unchecked);
		list_widget->item(7)->setCheckState(Qt::Unchecked);
    list_widget->item(8)->setCheckState(Qt::Unchecked);
    list_widget->item(9)->setCheckState(Qt::Unchecked);
	}
	
	void TheoreticalSpectrumGenerationDialog::itemChanged(QListWidgetItem* item)
	{
		if (item->text() == "Isotope clusters")
		{
			if (item->checkState() == Qt::Checked)
			{
				max_iso_label->setEnabled(true);
				max_iso_spinbox->setEnabled(true);
			}
			else
			{
				max_iso_label->setEnabled(false);
				max_iso_spinbox->setEnabled(false);
			}
		}
	}

} // namespace
