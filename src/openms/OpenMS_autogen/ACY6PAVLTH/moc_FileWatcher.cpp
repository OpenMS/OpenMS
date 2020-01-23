/****************************************************************************
** Meta object code from reading C++ file 'FileWatcher.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.11.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../../../../../openms/OpenMS/src/openms/include/OpenMS/SYSTEM/FileWatcher.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'FileWatcher.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.11.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_OpenMS__FileWatcher_t {
    QByteArrayData data[7];
    char stringdata0[81];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_OpenMS__FileWatcher_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_OpenMS__FileWatcher_t qt_meta_stringdata_OpenMS__FileWatcher = {
    {
QT_MOC_LITERAL(0, 0, 19), // "OpenMS::FileWatcher"
QT_MOC_LITERAL(1, 20, 11), // "fileChanged"
QT_MOC_LITERAL(2, 32, 0), // ""
QT_MOC_LITERAL(3, 33, 6), // "String"
QT_MOC_LITERAL(4, 40, 19), // "monitorFileChanged_"
QT_MOC_LITERAL(5, 60, 4), // "name"
QT_MOC_LITERAL(6, 65, 15) // "timerTriggered_"

    },
    "OpenMS::FileWatcher\0fileChanged\0\0"
    "String\0monitorFileChanged_\0name\0"
    "timerTriggered_"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_OpenMS__FileWatcher[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       3,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    1,   29,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
       4,    1,   32,    2, 0x09 /* Protected */,
       6,    0,   35,    2, 0x09 /* Protected */,

 // signals: parameters
    QMetaType::Void, 0x80000000 | 3,    2,

 // slots: parameters
    QMetaType::Void, QMetaType::QString,    5,
    QMetaType::Void,

       0        // eod
};

void OpenMS::FileWatcher::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        FileWatcher *_t = static_cast<FileWatcher *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->fileChanged((*reinterpret_cast< const String(*)>(_a[1]))); break;
        case 1: _t->monitorFileChanged_((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 2: _t->timerTriggered_(); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        {
            using _t = void (FileWatcher::*)(const String & );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&FileWatcher::fileChanged)) {
                *result = 0;
                return;
            }
        }
    }
}

QT_INIT_METAOBJECT const QMetaObject OpenMS::FileWatcher::staticMetaObject = {
    { &QFileSystemWatcher::staticMetaObject, qt_meta_stringdata_OpenMS__FileWatcher.data,
      qt_meta_data_OpenMS__FileWatcher,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *OpenMS::FileWatcher::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *OpenMS::FileWatcher::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_OpenMS__FileWatcher.stringdata0))
        return static_cast<void*>(this);
    return QFileSystemWatcher::qt_metacast(_clname);
}

int OpenMS::FileWatcher::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QFileSystemWatcher::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 3)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 3;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 3)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 3;
    }
    return _id;
}

// SIGNAL 0
void OpenMS::FileWatcher::fileChanged(const String & _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
