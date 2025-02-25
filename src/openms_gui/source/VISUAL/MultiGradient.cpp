// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/MultiGradient.h>
#include <cstdlib>
#include <limits>
#include <sstream>
#include <utility>

using namespace std;

namespace OpenMS
{
  using namespace Exception;

  MultiGradient::MultiGradient() :
    pos_col_(),
    interpolation_mode_(IM_LINEAR),
    pre_min_(0),
    pre_size_(0),
    pre_steps_(0)
  {
    pos_col_[0] = Qt::white;
    pos_col_[100] = Qt::black;
  }

  MultiGradient::MultiGradient(const MultiGradient & multigradient) = default;

  MultiGradient & MultiGradient::operator=(const MultiGradient & rhs)
  {
    if (this == &rhs)
    {
      return *this;
    }
    pos_col_ = rhs.pos_col_,
    interpolation_mode_ = rhs.interpolation_mode_;
    pre_  = rhs.pre_;
    pre_min_ = rhs.pre_min_;
    pre_size_ = rhs.pre_size_;
    pre_steps_ = rhs.pre_steps_;
    return *this;
  }

  MultiGradient::~MultiGradient() = default;

  Size MultiGradient::size() const
  {
    return pos_col_.size();
  }

  void MultiGradient::insert(double position, QColor color)
  {
    if (position >= 0 && position <= 100)
    {
      pos_col_[position] = std::move(color);
    }
    else
    {
      throw InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
  }

  UInt MultiGradient::position(UInt index)
  {
    if (index > size() - 1)
    {
      throw IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    map<double, QColor>::iterator it = pos_col_.begin();
    for (Size i = 0; i < index; ++i)
    {
      ++it;
    }
    return it->first;
  }

  QColor MultiGradient::color(UInt index)
  {
    if (index > size() - 1)
    {
      throw IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    map<double, QColor>::iterator it = pos_col_.begin();
    for (Size i = 0; i < index; ++i)
    {
      ++it;
    }
    return it->second;
  }

  QColor MultiGradient::interpolatedColorAt(double position) const
  {
    if (position <= 0.0)
    {
      return pos_col_.begin()->second;
    }

    if (position >= 100.0)
    {
      return (--(pos_col_.end()))->second;
    }

    //linear
    if (interpolation_mode_ == IM_LINEAR)
    {
      map<double, QColor>::const_iterator it1 = pos_col_.lower_bound(position);

      if (std::abs(it1->first - position) < numeric_limits<double>::epsilon())  // compare double
      {
        return it1->second;
      }
      else
      {
        map<double, QColor>::const_iterator it0 = it1;
        --it0;
        double factor = (position - it0->first) / (it1->first - it0->first);
        return QColor(Int(factor * it1->second.red() + (1 - factor) * it0->second.red() + 0.001)
                     , Int(factor * it1->second.green() + (1 - factor) * it0->second.green() + 0.001)
                     , Int(factor * it1->second.blue() + (1 - factor) * it0->second.blue() + 0.001));
      }
    }
    //stairs
    else
    {
      map<double, QColor>::const_iterator it = pos_col_.upper_bound(position);
      --it;
      return it->second;
    }
  }

  QColor MultiGradient::interpolatedColorAt(double position, double min, double max) const
  {
    return interpolatedColorAt((position - min) / (max - min) * 100.0);
  }

  void MultiGradient::setInterpolationMode(MultiGradient::InterpolationMode mode)
  {
    interpolation_mode_ = mode;
  }

  MultiGradient::InterpolationMode MultiGradient::getInterpolationMode() const
  {
    return interpolation_mode_;
  }

  string MultiGradient::toString() const
  {
    stringstream out;

    //Interpolation Mode
    if (getInterpolationMode() == IM_LINEAR)
    {
      out << "Linear|";
    }
    else if (getInterpolationMode() == IM_STAIRS)
    {
      out << "Stairs|";
    }

    for (map<double, QColor>::const_iterator it = pos_col_.begin(); it != pos_col_.end(); ++it)
    {
      if (it != pos_col_.begin())
      {
        out << ";";
      }
      out << it->first << "," << it->second.name().toStdString();
    }
    return out.str();
  }

  void MultiGradient::fromString(const string & gradient)
  {
    pos_col_.clear();

    if (gradient.empty())
    {
      pos_col_[0] = Qt::white;
      pos_col_[100] = Qt::black;
      return;
    }

    string g(gradient);
    string::iterator tmp(g.begin());
    double tmp_pos = 0;
    for (string::iterator it = g.begin(); it != g.end(); ++it)
    {
      if (*it == '|')
      {
        //interploation mode
        if (string(tmp, it) == "Linear")
        {
          setInterpolationMode(IM_LINEAR);
        }
        else if (string(tmp, it) == "Stairs")
        {
          setInterpolationMode(IM_STAIRS);
        }
      }
      else if (*it == ';')
      {
        pos_col_[tmp_pos] = QColor(string(tmp, it).c_str());
        tmp = it + 1;
      }
      else if (*it == ',')
      {
        tmp_pos = QString(string(tmp, it).c_str()).toDouble();
        tmp = it + 1;
      }
    }
    // last entry
    pos_col_[tmp_pos] = QColor(string(tmp, g.end()).c_str());
  }

  void MultiGradient::activatePrecalculationMode(double min, double max, UInt steps)
  {
    //add security margin to range to avoid numerical problems
    pre_min_ = std::min(min, max) - 0.000005;
    pre_size_ = fabs(max - min) + 0.00001;
    pre_steps_ = steps - 1;
    pre_.clear();
    pre_.reserve(steps);
    for (Size step = 0; step < steps; ++step)
    {
      pre_.push_back(interpolatedColorAt(step, 0, pre_steps_));
      //cout << pre_.back().red() << " " << pre_.back().green() << " " << pre_.back().blue() << endl;
    }
  }

  void MultiGradient::deactivatePrecalculationMode()
  {
    pre_.clear();
  }

  bool MultiGradient::exists(double position)
  {
    return pos_col_.find(position) != pos_col_.end();
  }

  bool MultiGradient::remove(double position)
  {
    if (position < 0 + std::numeric_limits<double>::epsilon() || position > 100 - std::numeric_limits<double>::epsilon())
    {
      return false;
    }

    map<double, QColor>::iterator it = pos_col_.find(position);
    if (it != pos_col_.end())
    {
      pos_col_.erase(it);
      return true;
    }
    return false;
  }

  // static
  MultiGradient MultiGradient::getDefaultGradientLinearIntensityMode()
  {
    MultiGradient mg;
    mg.fromString("Linear|0,#eeeeee;1,#ffea00;6,#ff0000;14,#aa00ff;23,#5500ff;100,#000000");
    return mg;
  }

  // static
  MultiGradient MultiGradient::getDefaultGradientLogarithmicIntensityMode()
  {
    MultiGradient mg;
    mg.fromString("Linear|0,#ffea00;6,#ff0000;14,#aa00ff;23,#5500ff;100,#000000");
    return mg;
  }

} //namespace OpenMS
