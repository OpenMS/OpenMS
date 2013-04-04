// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

//OpenMS
#include <OpenMS/VISUAL/VISUALIZER/GradientVisualizer.h>
#include <OpenMS/DATASTRUCTURES/String.h>

//QT
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QValidator>
#include <QtGui/QPushButton>
#include <QtGui/QGridLayout>

//STL
#include <iostream>
#include <vector>

using namespace std;

namespace OpenMS
{

  GradientVisualizer::GradientVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<Gradient>()
  {
    nextrow_ = 0;
  }

  void GradientVisualizer::load(Gradient & g)
  {
    ptr_ = &g;
    temp_ = g;

    addLabel_("Modify Gradient information");
    addSeparator_();

    viewlayout_ = new QGridLayout();
    mainlayout_->addLayout(viewlayout_, row_, 0, 1, 3);
    row_++;
    //Get the actuall eluent, timepoint and percentage values.
    loadData_();

    addSeparator_();

    addLineEditButton_("Add new Eluent", new_eluent_, add_eluent_button_, "Add Eluent");
    addLineEditButton_("Add new Timepoint", new_timepoint_, add_timepoint_button_, "Add Timepoint");
    addLabel_("Attention: All percentage values at a certain timepoint must add up to 100.");
    addSeparator_();
    addLabel_("Remove all eluents, timepoints and percentage values.");
    addButton_(removebutton_, "Remove");

    finishAdding_();
    addSeparator_();
    connect(add_timepoint_button_, SIGNAL(clicked()), this, SLOT(addTimepoint_()));
    connect(add_eluent_button_, SIGNAL(clicked()), this, SLOT(addEluent_()));
    connect(removebutton_, SIGNAL(clicked()), this, SLOT(deleteData_()));

    //Input validator
    timepoint_vali_ = new QIntValidator(new_timepoint_);
    new_timepoint_->setValidator(timepoint_vali_);

  }

  //----------------------------------------------------------------------------
  // SLOTS
  //----------------------------------------------------------------------------
  void GradientVisualizer::deleteData_()
  {
    //Remove entries from Gradient
    temp_.clearEluents();
    temp_.clearTimepoints();
    temp_.clearPercentages();

    update_();
  }

  void GradientVisualizer::addTimepoint_()
  {
    //Check whether new timepoint is in range
    String m(new_timepoint_->text());

    if (timepoints_.empty() && m.trim().length() != 0)
    {
      temp_.addTimepoint(m.toInt());
      update_();
    }
    else
    {
      if (m.trim().length() != 0 && timepoints_.back() < m.toInt())
      {
        temp_.addTimepoint(m.toInt());
        update_();
      }
    }
  }

  void GradientVisualizer::addEluent_()
  {
    String m(new_eluent_->text());
    std::vector<String>::iterator iter;
    //check if eluent name is empty
    if (m.trim().length() != 0)
    {
      //check if eluent name already exists
      for (iter = eluents_.begin(); iter < eluents_.end(); ++iter)
      {
        if (*iter == m)
        {
          return;
        }
      }
      temp_.addEluent(m);
      update_();
    }
  }

  void GradientVisualizer::store()
  {
    //Check, whether sum is 100
    Size time_size = timepoints_.size();
    Size elu_size = eluents_.size();
    Size count = 0;
    Size elu_count = 0;
    int sum_check = 0;
    for (Size i = 0; i < time_size; ++i)
    {
      elu_count = i;
      for (Size j = 0; j < elu_size; ++j)
      {
        String value((gradientdata_[elu_count])->text());
        elu_count  = elu_count + time_size;
        sum_check = sum_check + value.toInt();
        if (j == elu_size - 1 && sum_check != 100)
        {
          cout << "The sum does not add up to 100 !" << endl;
          cout << "Please check your values." << endl;
          return;
        }
      }
      sum_check = 0;
    }

    //Store all values into the gradient object
    for (Size i = 0; i < eluents_.size(); ++i)
    {
      for (Size j = 0; j < timepoints_.size(); ++j)
      {
        String value((gradientdata_[count + j])->text());
        temp_.setPercentage(eluents_[i], timepoints_[j], value.toInt());
      }
      count = count + time_size;

    }

    //copy temporary stored data into metainfo object
    (*ptr_) = temp_;
  }

  //----------------------------------------------------------------------------
  // private
  //----------------------------------------------------------------------------

  void GradientVisualizer::loadData_()
  {
    nextrow_ = 0;

    eluents_ =    temp_.getEluents();
    timepoints_ = temp_.getTimepoints();

    //Add labels to display eluent-timepoint-percentage-triplets.
    QLabel * label = new QLabel("Eluent names\\Timepoints", this);
    viewlayout_->addWidget(label, 0, 0, 1, (int)temp_.getTimepoints().size());
    label->show();
    nextrow_++;
    gradientlabel_.push_back(label);

    //----------------------------------------------------------------------
    // Actual data
    //----------------------------------------------------------------------
    for (Size i = 0; i < timepoints_.size(); ++i)
    {
      //Add labels to display eluent-timepoint-percentage-triplets.
      QLabel * label1 = new QLabel(String(timepoints_[i]).c_str(), this);
      viewlayout_->addWidget(label1, 1, (int)(i + 1));
      label1->show();
      gradientlabel_.push_back(label1);
    }

    nextrow_++;

    //Add the percentage values for the eluents and their timepoint.
    for (Size i = 0; i < eluents_.size(); ++i)
    {

      QLabel * eluent = new QLabel(eluents_[i].c_str(), this);
      viewlayout_->addWidget(eluent, nextrow_, 0);
      eluent->show();
      gradientlabel_.push_back(eluent);

      for (Size j = 0; j < timepoints_.size(); ++j)
      {
        percentage_ = new QLineEdit(this);
        percentage_->setText(String(temp_.getPercentage(eluents_[i], timepoints_[j])).c_str());
        viewlayout_->addWidget(percentage_, nextrow_, (int)(j + 1));
        //Store pointers to the QLineEdits
        gradientdata_.push_back(percentage_);
        percentage_->show();
      }

      nextrow_++;
    }
  }

  void GradientVisualizer::removeData_()
  {
    //Remove QLineEdits
    std::vector<QLineEdit *>::iterator iter2;

    for (iter2 = gradientdata_.begin(); iter2 < gradientdata_.end(); ++iter2)
    {
      //Delete QLineEdit field from viewlayout_
      viewlayout_->removeWidget((*iter2));
      (*iter2)->hide();
      //Set pointer to 0
      (*iter2) = 0;
      //Free memory of the pointer
      delete(*iter2);
    }

    //Remove QLabels
    std::vector<QLabel *>::iterator iter_label;

    for (iter_label = gradientlabel_.begin(); iter_label < gradientlabel_.end(); ++iter_label)
    {
      viewlayout_->removeWidget((*iter_label));
      (*iter_label)->hide();
      (*iter_label) = 0;
      delete(*iter_label);
    }

    gradientdata_.clear();
    gradientlabel_.clear();
  }

  void GradientVisualizer::update_()
  {
    removeData_();
    loadData_();
  }

  void GradientVisualizer::undo_()
  {
    removeData_();
    temp_ = (*ptr_);
    loadData_();
  }

}
