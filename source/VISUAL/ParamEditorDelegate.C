#include <QtGui/QItemDelegate>
#include <QtCore/QModelIndex>
#include <QtGui/QComboBox>
#include <QtCore/QAbstractItemModel>
#include <QtCore/QString>
#include <QtCore/QStringList>
#include <OpenMS/VISUAL/ParamEditorDelegate.h>

namespace OpenMS
{
	

	ParamEditorDelegate::ParamEditorDelegate(QObject *parent)
	     : QItemDelegate(parent)
	 {
	 }
	 
	 QWidget *ParamEditorDelegate::createEditor(QWidget *parent,
	     const QStyleOptionViewItem & option,
	     const QModelIndex & index ) const
	 {
		 QString str = index.model()->data(index, Qt::DisplayRole).toString();
		 if(index.column()==1 &&str.isEmpty()) return 0;
		 else if (index.column() == 2 && !str.isEmpty())
		{
			//  QVariant originalValue = index.model()->data(index, Qt::UserRole);
			//if (!isSupportedType(originalValue.type()))
			//  return 0;
			QComboBox *editor = new QComboBox(parent);
			QStringList list;
			list<<"int"<<"float"<<"string";
			editor->addItems(list);
			return editor;
		}
		else if(index.column()==2 && str.isEmpty()) return 0;
		 else return QItemDelegate::createEditor(parent,option,index);

	 }
	 
	 void ParamEditorDelegate::setEditorData(QWidget *editor,
					     const QModelIndex &index) const
	 {
		 if(index.column()==2)
		 {
	     QString str = index.model()->data(index, Qt::DisplayRole).toString();
		int pos = combo_list_.indexOf(str);
		 QComboBox *combo = static_cast<QComboBox*>(editor);
				if (pos!=-1)
				{
					combo->setCurrentIndex(pos);
				}
			}
			else QItemDelegate::setEditorData(editor,index);
	 }
	 
	 void ParamEditorDelegate::setModelData(QWidget *editor, QAbstractItemModel *model,
					    const QModelIndex &index) const
	 {
		 if(index.column()==2)
		 {
	     QComboBox *combo= static_cast<QComboBox*>(editor);
	     QString str = combo->currentText();

	     model->setData(index, str);
		 }
		 else QItemDelegate::setModelData(editor,model,index);
	 }
	 
	 void ParamEditorDelegate::updateEditorGeometry(QWidget *editor,
	     const QStyleOptionViewItem &option, const QModelIndex & index) const
	 {
		 if(index.column()==2)
		 {
	     editor->setGeometry(option.rect);
		 }
		 else QItemDelegate::updateEditorGeometry(editor,option,index);
	 }
	 
 }