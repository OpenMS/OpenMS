/****************************************************************************
** OpenMS::IdentificationVisualizer meta object code from reading C++ file 'IdentificationVisualizer.h'
**
** Created: Thu Jan 25 15:32:43 2007
**      by: The Qt MOC ($Id: qt/moc_yacc.cpp   3.3.6   edited Mar 8 17:43 $)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#undef QT_NO_COMPAT
#include "../../../include/OpenMS/VISUAL/VISUALIZER/IdentificationVisualizer.h"
#include <qmetaobject.h>
#include <qapplication.h>

#include <private/qucomextra_p.h>
#if !defined(Q_MOC_OUTPUT_REVISION) || (Q_MOC_OUTPUT_REVISION != 26)
#error "This file was generated using the moc from 3.3.6. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

const char *OpenMS::IdentificationVisualizer::className() const
{
    return "OpenMS::IdentificationVisualizer";
}

QMetaObject *OpenMS::IdentificationVisualizer::metaObj = 0;
static QMetaObjectCleanUp cleanUp_OpenMS__IdentificationVisualizer( "OpenMS::IdentificationVisualizer", &OpenMS::IdentificationVisualizer::staticMetaObject );

#ifndef QT_NO_TRANSLATION
QString OpenMS::IdentificationVisualizer::tr( const char *s, const char *c )
{
    if ( qApp )
	return qApp->translate( "OpenMS::IdentificationVisualizer", s, c, QApplication::DefaultCodec );
    else
	return QString::fromLatin1( s );
}
#ifndef QT_NO_TRANSLATION_UTF8
QString OpenMS::IdentificationVisualizer::trUtf8( const char *s, const char *c )
{
    if ( qApp )
	return qApp->translate( "OpenMS::IdentificationVisualizer", s, c, QApplication::UnicodeUTF8 );
    else
	return QString::fromUtf8( s );
}
#endif // QT_NO_TRANSLATION_UTF8

#endif // QT_NO_TRANSLATION

QMetaObject* OpenMS::IdentificationVisualizer::staticMetaObject()
{
    if ( metaObj )
	return metaObj;
    QMetaObject* parentObject = BaseVisualizer::staticMetaObject();
    static const QUMethod slot_0 = {"store", 0, 0 };
    static const QUMethod slot_1 = {"reject", 0, 0 };
    static const QUMethod slot_2 = {"updateTree", 0, 0 };
    static const QUMethod slot_3 = {"searchRefPeptides", 0, 0 };
    static const QUMethod slot_4 = {"searchNonRefPeptides", 0, 0 };
    static const QMetaData slot_tbl[] = {
	{ "store()", &slot_0, QMetaData::Private },
	{ "reject()", &slot_1, QMetaData::Private },
	{ "updateTree()", &slot_2, QMetaData::Private },
	{ "searchRefPeptides()", &slot_3, QMetaData::Private },
	{ "searchNonRefPeptides()", &slot_4, QMetaData::Private }
    };
    metaObj = QMetaObject::new_metaobject(
	"OpenMS::IdentificationVisualizer", parentObject,
	slot_tbl, 5,
	0, 0,
#ifndef QT_NO_PROPERTIES
	0, 0,
	0, 0,
#endif // QT_NO_PROPERTIES
	0, 0 );
    cleanUp_OpenMS__IdentificationVisualizer.setMetaObject( metaObj );
    return metaObj;
}

void* OpenMS::IdentificationVisualizer::qt_cast( const char* clname )
{
    if ( !qstrcmp( clname, "OpenMS::IdentificationVisualizer" ) )
	return this;
    return BaseVisualizer::qt_cast( clname );
}

bool OpenMS::IdentificationVisualizer::qt_invoke( int _id, QUObject* _o )
{
    switch ( _id - staticMetaObject()->slotOffset() ) {
    case 0: store(); break;
    case 1: reject(); break;
    case 2: updateTree(); break;
    case 3: searchRefPeptides(); break;
    case 4: searchNonRefPeptides(); break;
    default:
	return BaseVisualizer::qt_invoke( _id, _o );
    }
    return TRUE;
}

bool OpenMS::IdentificationVisualizer::qt_emit( int _id, QUObject* _o )
{
    return BaseVisualizer::qt_emit(_id,_o);
}
#ifndef QT_NO_PROPERTIES

bool OpenMS::IdentificationVisualizer::qt_property( int id, int f, QVariant* v)
{
    return BaseVisualizer::qt_property( id, f, v);
}

bool OpenMS::IdentificationVisualizer::qt_static_property( QObject* , int , int , QVariant* ){ return FALSE; }
#endif // QT_NO_PROPERTIES
