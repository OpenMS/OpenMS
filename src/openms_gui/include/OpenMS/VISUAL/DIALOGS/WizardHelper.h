// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Jihyung Kim $
// $Authors: Jihyung Kim, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

namespace OpenMS
{
  class InputFile;
  class OutputDirectory;
  class ParamEditor;

  namespace Internal
  {
    /**
      @brief RAII class to switch to certain TabWidget, disable the GUI and go back to the orignal Tab when this class is destroyed
    */
    template<class TWidgetClass>
    class WizardGUILock
    {
    public:
        WizardGUILock(TWidgetClass* stw):
            stw_(stw),
            old_(stw->currentWidget()),
            glock_(stw)
        {
          stw->setCurrentWidget(stw->ui->tab_log);
        }

        ~WizardGUILock()
        {
          stw_->setCurrentWidget(old_);
        }

      private:
        TWidgetClass* stw_;
        QWidget* old_;
        GUIHelpers::GUILock glock_;
    };

    /// custom arguments to allow for looping calls
    struct Args
    {
      QStringList loop_arg; ///< list of arguments to insert; one for every loop
      size_t insert_pos;       ///< where to insert in the target argument list (index is 0-based)
    };

    typedef std::vector<Args> ArgLoop;

    /// Allows running an executable with arguments
    /// Multiple execution in a loop is supported by the ArgLoop argument
    /// e.g. running 'ls -la .' and 'ls -la ..'
    /// uses Command("ls", QStringList() << "-la" << "%1", ArgLoop{ Args {QStringList() << "." << "..", 1 } })
    /// All lists in loop[i].loop_arg should have the same size (i.e. same number of loops)
    struct Command
    {
      String exe;
      QStringList args;
      ArgLoop loop;

      Command(const String& e, const QStringList& a, const ArgLoop& l) :
          exe(e),
          args(a),
          loop(l) {}

      /// how many loops can we make according to the ArgLoop provided?
      /// if ArgLoop is empty, we just do a single invokation
      size_t getLoopCount() const
      {
        if (loop.empty()) return 1;
        size_t common_size = loop[0].loop_arg.size();
        for (const auto& l : loop)
        {
          if (l.loop_arg.size() != (int)common_size) throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Internal error. Not all loop arguments support the same number of loops!");
          if ((int)l.insert_pos >= args.size()) throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Internal error. Loop argument wants to insert after end of template arguments!");
        }
        return common_size;
      }
      /// for a given loop, return the substituted arguments
      /// @p loop_number of 0 is always valid, i.e. no loop args, just use the unmodified args provided
      QStringList getArgs(const int loop_number) const
      {
        if (loop_number >= (int)getLoopCount())
        {
          throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Internal error. The loop number you requested is too high!");
        }
        if (loop.empty()) return args; // no looping available

        QStringList arg_l = args;
        for (const auto& largs : loop) // replace all args for the current round
        {
          arg_l[largs.insert_pos] = args[largs.insert_pos].arg(largs.loop_arg[loop_number]);
        }
        return arg_l;
      }
    };

  }
} // ns OpenMS

// this is required to allow Ui_[tool_name]TabWidget (auto UIC'd from .ui) to have a InputFile member
using InputFile = OpenMS::InputFile;
using OutputDirectory = OpenMS::OutputDirectory;
using ParamEditor = OpenMS::ParamEditor;
using TableView = OpenMS::TableView;
