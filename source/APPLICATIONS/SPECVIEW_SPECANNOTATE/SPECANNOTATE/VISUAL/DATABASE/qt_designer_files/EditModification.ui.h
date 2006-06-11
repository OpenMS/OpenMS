/****************************************************************************
** ui.h extension file, included from the uic-generated form implementation.
**
** If you wish to add, delete or rename functions or slots use
** Qt Designer which will update this file, preserving your code. Create an
** init() function in place of a constructor, and a destroy() function in
** place of a destructor.
*****************************************************************************/

#include <qmessagebox.h>

void EditModification::help()
{
    QMessageBox::information(this, tr("Database Help: modification"), 
			     tr("To get some information about the different entries of the table, just place the cursor above their names, or use the \"What's this?\" function! \n" "To add an entry into the database, first click \"Insert\", enter your data, then press \"Update\"!"), 1);
}
