/****************************************************************************
** Meta object code from reading C++ file 'MascotRemoteQuery.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.11.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../../../../../openms/OpenMS/src/openms/include/OpenMS/FORMAT/MascotRemoteQuery.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MascotRemoteQuery.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.11.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_OpenMS__MascotRemoteQuery_t {
    QByteArrayData data[14];
    char stringdata0[161];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_OpenMS__MascotRemoteQuery_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_OpenMS__MascotRemoteQuery_t qt_meta_stringdata_OpenMS__MascotRemoteQuery = {
    {
QT_MOC_LITERAL(0, 0, 25), // "OpenMS::MascotRemoteQuery"
QT_MOC_LITERAL(1, 26, 11), // "gotRedirect"
QT_MOC_LITERAL(2, 38, 0), // ""
QT_MOC_LITERAL(3, 39, 14), // "QNetworkReply*"
QT_MOC_LITERAL(4, 54, 5), // "reply"
QT_MOC_LITERAL(5, 60, 4), // "done"
QT_MOC_LITERAL(6, 65, 3), // "run"
QT_MOC_LITERAL(7, 69, 8), // "timedOut"
QT_MOC_LITERAL(8, 78, 12), // "readResponse"
QT_MOC_LITERAL(9, 91, 16), // "downloadProgress"
QT_MOC_LITERAL(10, 108, 10), // "bytes_read"
QT_MOC_LITERAL(11, 119, 11), // "bytes_total"
QT_MOC_LITERAL(12, 131, 14), // "uploadProgress"
QT_MOC_LITERAL(13, 146, 14) // "followRedirect"

    },
    "OpenMS::MascotRemoteQuery\0gotRedirect\0"
    "\0QNetworkReply*\0reply\0done\0run\0timedOut\0"
    "readResponse\0downloadProgress\0bytes_read\0"
    "bytes_total\0uploadProgress\0followRedirect"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_OpenMS__MascotRemoteQuery[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       8,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       2,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    1,   54,    2, 0x06 /* Public */,
       5,    0,   57,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
       6,    0,   58,    2, 0x0a /* Public */,
       7,    0,   59,    2, 0x08 /* Private */,
       8,    1,   60,    2, 0x08 /* Private */,
       9,    2,   63,    2, 0x08 /* Private */,
      12,    2,   68,    2, 0x08 /* Private */,
      13,    1,   73,    2, 0x08 /* Private */,

 // signals: parameters
    QMetaType::Void, 0x80000000 | 3,    4,
    QMetaType::Void,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 3,    4,
    QMetaType::Void, QMetaType::LongLong, QMetaType::LongLong,   10,   11,
    QMetaType::Void, QMetaType::LongLong, QMetaType::LongLong,   10,   11,
    QMetaType::Void, 0x80000000 | 3,    4,

       0        // eod
};

void OpenMS::MascotRemoteQuery::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        MascotRemoteQuery *_t = static_cast<MascotRemoteQuery *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->gotRedirect((*reinterpret_cast< QNetworkReply*(*)>(_a[1]))); break;
        case 1: _t->done(); break;
        case 2: _t->run(); break;
        case 3: _t->timedOut(); break;
        case 4: _t->readResponse((*reinterpret_cast< QNetworkReply*(*)>(_a[1]))); break;
        case 5: _t->downloadProgress((*reinterpret_cast< qint64(*)>(_a[1])),(*reinterpret_cast< qint64(*)>(_a[2]))); break;
        case 6: _t->uploadProgress((*reinterpret_cast< qint64(*)>(_a[1])),(*reinterpret_cast< qint64(*)>(_a[2]))); break;
        case 7: _t->followRedirect((*reinterpret_cast< QNetworkReply*(*)>(_a[1]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        switch (_id) {
        default: *reinterpret_cast<int*>(_a[0]) = -1; break;
        case 0:
            switch (*reinterpret_cast<int*>(_a[1])) {
            default: *reinterpret_cast<int*>(_a[0]) = -1; break;
            case 0:
                *reinterpret_cast<int*>(_a[0]) = qRegisterMetaType< QNetworkReply* >(); break;
            }
            break;
        case 4:
            switch (*reinterpret_cast<int*>(_a[1])) {
            default: *reinterpret_cast<int*>(_a[0]) = -1; break;
            case 0:
                *reinterpret_cast<int*>(_a[0]) = qRegisterMetaType< QNetworkReply* >(); break;
            }
            break;
        case 7:
            switch (*reinterpret_cast<int*>(_a[1])) {
            default: *reinterpret_cast<int*>(_a[0]) = -1; break;
            case 0:
                *reinterpret_cast<int*>(_a[0]) = qRegisterMetaType< QNetworkReply* >(); break;
            }
            break;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        {
            using _t = void (MascotRemoteQuery::*)(QNetworkReply * );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&MascotRemoteQuery::gotRedirect)) {
                *result = 0;
                return;
            }
        }
        {
            using _t = void (MascotRemoteQuery::*)();
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&MascotRemoteQuery::done)) {
                *result = 1;
                return;
            }
        }
    }
}

QT_INIT_METAOBJECT const QMetaObject OpenMS::MascotRemoteQuery::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_OpenMS__MascotRemoteQuery.data,
      qt_meta_data_OpenMS__MascotRemoteQuery,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *OpenMS::MascotRemoteQuery::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *OpenMS::MascotRemoteQuery::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_OpenMS__MascotRemoteQuery.stringdata0))
        return static_cast<void*>(this);
    if (!strcmp(_clname, "DefaultParamHandler"))
        return static_cast< DefaultParamHandler*>(this);
    return QObject::qt_metacast(_clname);
}

int OpenMS::MascotRemoteQuery::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 8)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 8;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 8)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 8;
    }
    return _id;
}

// SIGNAL 0
void OpenMS::MascotRemoteQuery::gotRedirect(QNetworkReply * _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void OpenMS::MascotRemoteQuery::done()
{
    QMetaObject::activate(this, &staticMetaObject, 1, nullptr);
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
