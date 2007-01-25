/****************************************************************************
** OpenMS::MetaInfoVisualizer meta object code from reading C++ file 'MetaInfoVisualizer.h'
**
** Created: Thu Jan 25 15:16:40 2007
**      by: The Qt MOC ($Id: qt/moc_yacc.cpp   3.3.6   edited Mar 8 17:43 $)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#undef QT_NO_COMPAT
#include "../../../include/OpenMS/VISUAL/VISUALIZER/MetaInfoVisualizer.h"
#include <qmetaobject.h>
#include <qapplication.h>

#include <private/qucomextra_p.h>
#if !defined(Q_MOC_OUTPUT_REVISION) || (Q_MOC_OUTPUT_REVISION != 26)
#error "This file was generated using the moc from 3.3.6. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

const char *OpenMS::MetaInfoVisualizer::className() const
{
    return "OpenMS::MetaInfoVisualizer";
}

QMetaObject *OpenMS::MetaInfoVisualizer::metaObj = 0;
static QMetaObjectCleanUp cleanUp_OpenMS__MetaInfoVisualizer( "OpenMS::MetaInfoVisualizer", &OpenMS::MetaInfoVisualizer::staticMetaObject );

#ifndef QT_NO_TRANSLATION
QString OpenMS::MetaInfoVisualizer::tr( const char *s, const char *c )
{
    if ( qApp )
	return qApp->translate( "OpenMS::MetaInfoVisualizer", s, c, QApplication::DefaultCodec );
    else
	return QString::fromLatin1( s );
}
#ifndef QT_NO_TRANSLATION_UTF8
QString OpenMS::MetaInfoVisualizer::trUtf8( const char *s, const char *c )
{
    if ( qApp )
	return qApp->translate( "OpenMS::MetaInfoVisualizer", s, c, QApplication::UnicodeUTF8 );
    else
	return QString::fromUtf8( s );
}
#endif // QT_NO_TRANSLATION_UTF8

#endif // QT_NO_TRANSLATION

QMetaObject* OpenMS::MetaInfoVisualizer::staticMetaObject()
{
    if ( metaObj )
	return metaObj;
    QMetaObject* parentObject = BaseVisualizer::staticMetaObject();
    static const QUMethod slot_0 = {"add", 0, 0 };
    static const QUMethod slot_1 = {"clear", 0, 0 };
    static const QUParameter param_slot_2[] = {
	{ 0, &static_QUType_int, 0, QUParameter::In }
    };
    static const QUMethod slot_2 = {"remove", 1, param_slot_2 };
    static const QUMethod slot_3 = {"store", 0, 0 };
    static const QUMethod slot_4 = {"reject", 0, 0 };
    static const QMetaData slot_tbl[] = {
	{ "add()", &slot_0, QMetaData::Private },
	{ "clear()", &slot_1, QMetaData::Private },
	{ "remove(int)", &slot_2, QMetaData::Private },
	{ "store()", &slot_3, QMetaData::Private },
	{ "reject()", &slot_4, QMetaData::Private }
    };
    metaObj = QMetaObject::new_metaobject(
	"OpenMS::MetaInfoVisualizer", parentObject,
	slot_tbl, 5,
	0, 0,
#ifndef QT_NO_PROPERTIES
	0, 0,
	0, 0,
#endif // QT_NO_PROPERTIES
	0, 0 );
    cleanUp_OpenMS__MetaInfoVisualizer.setMetaObject( metaObj );
    return metaObj;
}

void* OpenMS::MetaInfoVisualizer::qt_cast( const char* clname )
{
    if ( !qstrcmp( clname, "OpenMS::MetaInfoVisualizer" ) )
	return this;
    return BaseVisualizer::qt_cast( clname );
}

bool OpenMS::MetaInfoVisualizer::qt_invoke( int _id, QUObject* _o )
{
    switch ( _id - staticMetaObject()->slotOffset() ) {
    case 0: add(); break;
    case 1: clear(); break;
    case 2: remove((int)static_QUType_int.get(_o+1)); break;
    case 3: store(); break;
    case 4: reject(); break;
    default:
	return BaseVisualizer::qt_invoke( _id, _o );
    }
    return TRUE;
}

bool OpenMS::MetaInfoVisualizer::qt_emit( int _id, QUObject* _o )
{
    return BaseVisualizer::qt_emit(_id,_o);
}
#ifndef QT_NO_PROPERTIES

bool OpenMS::MetaInfoVisualizer::qt_property( int id, int f, QVariant* v)
{
    return BaseVisualizer::qt_property( id, f, v);
}

bool OpenMS::MetaInfoVisualizer::qt_static_property( QObject* , int , int , QVariant* ){ return FALSE; }
#endif // QT_NO_PROPERTIES
