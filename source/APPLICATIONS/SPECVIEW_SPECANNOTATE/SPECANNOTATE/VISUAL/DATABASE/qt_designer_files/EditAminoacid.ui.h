/****************************************************************************
** ui.h extension file, included from the uic-generated form implementation.
**
** If you wish to add, delete or rename functions or slots use
** Qt Designer which will update this file, preserving your code. Create an
** init() function in place of a constructor, and a destroy() function in
** place of a destructor.
*****************************************************************************/


#include <qmessagebox.h>


void EditAminoacid::help()
{
QMessageBox::information(this, tr("Database Help: aminoacid"), 
			 tr("Information of different positions of amino acids are stored in database: \n \"Middle\": amino acid is in a chain between two other amino acids \n \"N-term.\": amino acid is at N-terminal position \n \"C-term.\": amino acid is on C-terminal position \n \"Single\": amino acid is in no chain at all"), 1);

}
