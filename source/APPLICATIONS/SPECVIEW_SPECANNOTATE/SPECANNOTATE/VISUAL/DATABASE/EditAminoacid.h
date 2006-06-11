/****************************************************************************
** Form interface generated from reading ui file 'SPECANNOTATE/VISUAL/DATABASE/EditAminoacid.ui'
**
** Created: Do Nov 4 11:38:39 2004
**      by: The User Interface Compiler ($Id: EditAminoacid.h,v 1.1 2005/11/11 06:29:42 marc_sturm Exp $)
**
** WARNING! All changes made in this file will be lost!
****************************************************************************/

#ifndef EDITAMINOACID_H
#define EDITAMINOACID_H

#include <qvariant.h>
#include <qdialog.h>

class QVBoxLayout;
class QHBoxLayout;
class QGridLayout;
class QSqlDatabase;
class QSqlCursor;
class QSqlForm;
class QDataBrowser;
class QLabel;
class QLineEdit;
class QPushButton;


namespace OpenMS
  {

  class EditAminoacid : public QDialog
    {
      Q_OBJECT

    public:
      EditAminoacid( QWidget* parent = 0, const char* name = 0, bool modal = FALSE, WFlags fl = 0 );
      ~EditAminoacid();

      QDataBrowser* dataBrowser3;
      QLabel* labelOne_letter_code;
      QLabel* labelAminoacid_name;
      QLineEdit* QLineEditC_term_average_mass;
      QLabel* labelMiddle_formula;
      QLineEdit* QLineEditThree_letter_code;
      QLineEdit* QLineEditAminoacid_name;
      QLabel* labelSingle_average_mass;
      QLabel* labelC_term_average_mass;
      QLabel* labelC_term_formula;
      QLineEdit* QLineEditMiddle_average_mass;
      QLineEdit* QLineEditMiddle_mono_mass;
      QLabel* labelMiddle_mono_mass;
      QLineEdit* QLineEditMiddle_formula;
      QLabel* labelN_term_formula;
      QLineEdit* QLineEditC_term_formula;
      QLineEdit* QLineEditOne_letter_code;
      QLineEdit* QLineEditSingle_average_mass;
      QLineEdit* QLineEditSingle_mono_mass;
      QLabel* labelSingle_mono_mass;
      QLineEdit* QLineEditN_term_mono_mass;
      QLabel* labelN_term_mono_mass;
      QLineEdit* QLineEditN_term_formula;
      QLabel* labelN_term_average_mass;
      QLineEdit* QLineEditSingle_formula;
      QLabel* labelMiddle_average_mass;
      QLineEdit* QLineEditC_term_mono_mass;
      QLabel* labelC_term_mono_mass;
      QLabel* labelSingle_formula;
      QLineEdit* QLineEditN_term_average_mass;
      QLabel* labelThree_letter_code;
      QPushButton* PushButtonFirst_2;
      QPushButton* PushButtonPrev_2;
      QPushButton* PushButtonNext_2;
      QPushButton* PushButtonLast_2;
      QPushButton* PushButtonInsert_2;
      QPushButton* PushButtonUpdate_2;
      QPushButton* PushButtonDelete_2;
      QPushButton* Done;
      QPushButton* Help;

    public slots:
      virtual void polish();

      virtual void help();

      virtual void insertIntoXML();


    protected:
      QGridLayout* dataBrowser3Layout;
      QGridLayout* layout7;
      QHBoxLayout* layout8;
      QHBoxLayout* layout9;

    protected slots:
      virtual void languageChange();


    };
}
#endif // EDITAMINOACID_H
