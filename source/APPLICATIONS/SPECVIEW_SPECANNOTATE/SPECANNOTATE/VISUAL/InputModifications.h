/****************************************************************************
** Form interface generated from reading ui file 'SPECANNOTATE/VISUAL/InputModifications.ui'
**
** Created: Do Nov 4 11:38:39 2004
**      by: The User Interface Compiler ($Id: InputModifications.h,v 1.1 2005/11/11 06:29:42 marc_sturm Exp $)
**
** WARNING! All changes made in this file will be lost!
****************************************************************************/

#ifndef INPUTMODIFICATIONS_H
#define INPUTMODIFICATIONS_H

#include <qvariant.h>
#include <qdialog.h>
#include "SampleDialog.h"
#include "SpecAnnotate.h"

class QVBoxLayout;
class QHBoxLayout;
class QGridLayout;
class QLabel;
class QTextBrowser;
class QListBox;
class QListBoxItem;
class QPushButton;

namespace OpenMS
  {

  class InputModifications : public QDialog
    {
      Q_OBJECT

    public:
      InputModifications( QWidget* parent = 0, const char* name = 0, bool modal = FALSE, WFlags fl = 0 );
      ~InputModifications();

      QLabel* textLabel3;
      QTextBrowser* textBrowser1;
      QLabel* textLabel1;
      QLabel* textLabel2;
      QListBox* listBox1;
      QListBox* listBox1_2;
      QPushButton* pushButton1;
      QPushButton* pushButton1_2;
      QPushButton* pushButton4;
      QPushButton* pushButton5;
      QPushButton* pushButton3;

    public slots:
      virtual void addGroup();
      virtual void resetSelection();
      virtual void resetString();

    protected:
      SampleDialog* sd;
      SpecAnnotate* msa;
      QMap<QString, QString>* settings_;

      QVBoxLayout* InputModificationsLayout;
      QVBoxLayout* layout7;
      QHBoxLayout* layout4;
      QHBoxLayout* layout1;
      QGridLayout* layout5;

    protected slots:
      virtual void languageChange();

    private:
      virtual void init();

    private slots:
      virtual void done();

    };

}
#endif // INPUTMODIFICATIONS_H
