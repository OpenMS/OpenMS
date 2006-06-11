#ifndef __SPECTRUMMDIWINDOWENHANCED_H__
#define __SPECTRUMMDIWINDOWENHANCED_H__

#include "OpenMS/VISUAL/SpectrumMDIWindow.h"
#include "SpecAnnotate.h"

namespace OpenMS
  {

  class SpectrumMDIWindowEnhanced : public OpenMS::SpectrumMDIWindow
    {
      Q_OBJECT

    public:

      static SpectrumMDIWindowEnhanced* getInstance();

    private:
       static SpectrumMDIWindowEnhanced* singleton_instance_;

      SpectrumMDIWindowEnhanced(QWidget* parent=0, const char* name="SpectrumMDIWindow", WFlags f=0);

      ~SpectrumMDIWindowEnhanced();

    public slots:

      void runAnnotate();
    };
}
#endif

