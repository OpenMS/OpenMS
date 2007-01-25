/****************************************************************************
** OpenMS::InstrumentSettingsVisualizer meta object code from reading C++ file 'InstrumentSettingsVisualizer.h'
**
** Created: Thu Jan 25 15:16:41 2007
**      by: The Qt MOC ($Id: qt/moc_yacc.cpp   3.3.6   edited Mar 8 17:43 $)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#undef QT_NO_COMPAT
#include "../../../include/OpenMS/VISUAL/VISUALIZER/InstrumentSettingsVisualizer.h"
#include <qmetaobject.h>
#include <qapplication.h>

#include <private/qucomextra_p.h>
#if !defined(Q_MOC_OUTPUT_REVISION) || (Q_MOC_OUTPUT_REVISION != 26)
#error "This file was generated using the moc from 3.3.6. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

const char *OpenMS::InstrumentSettingsVisualizer::className() const
{
    return "OpenMS::InstrumentSettingsVisualizer";
}

QMetaObject *OpenMS::InstrumentSettingsVisualizer::metaObj = 0;
static QMetaObjectCleanUp cleanUp_OpenMS__InstrumentSettingsVisualizer( "OpenMS::InstrumentSettingsVisualizer", &OpenMS::InstrumentSettingsVisualizer::staticMetaObject );

#ifndef QT_NO_TRANSLATION
QString OpenMS::InstrumentSettingsVisualizer::tr( const char *s, const char *c )
{
    if ( qApp )
	return qApp->translate( "OpenMS::InstrumentSettingsVisualizer", s, c, QApplication::DefaultCodec );
    else
	return QString::fromLatin1( s );
}
#ifndef QT_NO_TRANSLATION_UTF8
QString OpenMS::InstrumentSettingsVisualizer::trUtf8( const char *s, const char *c )
{
    if ( qApp )
	return qApp->translate( "OpenMS::InstrumentSettingsVisualizer", s, c, QApplication::UnicodeUTF8 );
    else
	return QString::fromUtf8( s );
}
#endif // QT_NO_TRANSLATION_UTF8

#endif // QT_NO_TRANSLATION

QMetaObject* OpenMS::InstrumentSettingsVisualizer::staticMetaObject()
{
    if ( metaObj )
	return metaObj;
    QMetaObject* parentObject = BaseVisualizer::staticMetaObject();
    static const QUMethod slot_0 = {"store", 0, 0 };
    static const QUMethod slot_1 = {"reject", 0, 0 };
    static const QMetaData slot_tbl[] = {
	{ "store()", &slot_0, QMetaData::Private },
	{ "reject()", &slot_1, QMetaData::Private }
    };
    metaObj = QMetaObject::new_metaobject(
	"OpenMS::InstrumentSettingsVisualizer", parentObject,
	slot_tbl, 2,
	0, 0,
#ifndef QT_NO_PROPERTIES
	0, 0,
	0, 0,
#endif // QT_NO_PROPERTIES
	0, 0 );
    cleanUp_OpenMS__InstrumentSettingsVisualizer.setMetaObject( metaObj );
    return metaObj;
}

void* OpenMS::InstrumentSettingsVisualizer::qt_cast( const char* clname )
{
    if ( !qstrcmp( clname, "OpenMS::InstrumentSettingsVisualizer" ) )
	return this;
    return BaseVisualizer::qt_cast( clname );
}

bool OpenMS::InstrumentSettingsVisualizer::qt_invoke( int _id, QUObject* _o )
{
    switch ( _id - staticMetaObject()->slotOffset() ) {
    case 0: store(); break;
    case 1: reject(); break;
    default:
	return BaseVisualizer::qt_invoke( _id, _o );
    }
    return TRUE;
}

bool OpenMS::InstrumentSettingsVisualizer::qt_emit( int _id, QUObject* _o )
{
    return BaseVisualizer::qt_emit(_id,_o);
}
#ifndef QT_NO_PROPERTIES

bool OpenMS::InstrumentSettingsVisualizer::qt_property( int id, int f, QVariant* v)
{
    return BaseVisualizer::qt_property( id, f, v);
}

bool OpenMS::InstrumentSettingsVisualizer::qt_static_property( QObject* , int , int , QVariant* ){ return FALSE; }
#endif // QT_NO_PROPERTIES
