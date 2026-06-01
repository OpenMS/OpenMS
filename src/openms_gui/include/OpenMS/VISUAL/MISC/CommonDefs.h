// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DPosition.h>


namespace OpenMS
{

  /// Macro for Qt's connect() overload resolution (in case signals/slots are overloaded and we need to tell connect what overload to pick
  /// without repeating ourselves.
  /// This can be solved in Qt 5.7 by using qOverload<>
  /// @note: provide the brackets for 'args' yourself, since there might be multiple arguments, separated by comma
  /// Example: QObject::connect(spinBox, CONNECTCAST(QSpinBox, valueChanged, (double)), slider, &QSlider::setValue);
  #define CONNECTCAST(class,func,args) static_cast<void(class::*)args>(&class::func)


  /// Enum to decide which headers(=column) names should be get/set in a table/tree widget
  enum class WidgetHeader
  {
    VISIBLE_ONLY,
    WITH_INVISIBLE,
  };

  /// Type of the Points in a 'flat' canvas (1D and 2D)
  using PointXYType = DPosition<2U>;
}
