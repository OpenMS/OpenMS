
#include "SpectrumMDIWindowEnhanced.h"


using namespace OpenMS;


//singleton instance
SpectrumMDIWindowEnhanced* SpectrumMDIWindowEnhanced::singleton_instance_ = NULL;



SpectrumMDIWindowEnhanced* SpectrumMDIWindowEnhanced::getInstance()
{
  if (singleton_instance_ == NULL)
    {
      singleton_instance_ = new SpectrumMDIWindowEnhanced();
    }
  return singleton_instance_;
	
//	if (!instance_)
//	{
//		instance_ = new SpectrumMDIWindowEnhanced();
//	}
//	return static_cast<SpectrumMDIWindowEnhanced*>(instance_);
}



SpectrumMDIWindowEnhanced::SpectrumMDIWindowEnhanced(QWidget* parent, const char* name, WFlags f)
    : SpectrumMDIWindow(parent, name, f)
{
   tools_menu_->insertItem(tr("Spec&Annotate (Annotate Peaks)"), this, SLOT(runAnnotate()));

}


SpectrumMDIWindowEnhanced::~SpectrumMDIWindowEnhanced()
{
   delete (singleton_instance_);
}



void SpectrumMDIWindowEnhanced::runAnnotate()
{
  SpecAnnotate * spec_annotate = new SpecAnnotate();
  spec_annotate->setCaption( "SpecAnnotate" );
  spec_annotate->show();
}
