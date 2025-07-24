/****************************************************************************
** Meta object code from reading C++ file 'TSGDialog_test.h'
**
** Created by: The Qt Meta Object Compiler version 68 (Qt 6.6.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../include/GUI/TSGDialog_test.h"
#include <QtNetwork/QSslError>
#include <QtGui/qtextcursor.h>
#include <QtGui/qscreen.h>
#include <QtCore/qmetatype.h>

#if __has_include(<QtCore/qtmochelpers.h>)
#include <QtCore/qtmochelpers.h>
#else
QT_BEGIN_MOC_NAMESPACE
#endif


#include <memory>

#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'TSGDialog_test.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 68
#error "This file was generated using the moc from 6.6.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

#ifndef Q_CONSTINIT
#define Q_CONSTINIT
#endif

QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
QT_WARNING_DISABLE_GCC("-Wuseless-cast")
namespace {

#ifdef QT_MOC_HAS_STRINGDATA
struct qt_meta_stringdata_CLASSOpenMSSCOPETestTSGDialogENDCLASS_t {};
constexpr auto qt_meta_stringdata_CLASSOpenMSSCOPETestTSGDialogENDCLASS = QtMocHelpers::stringData(
    "OpenMS::TestTSGDialog",
    "testConstruction",
    "",
    "testGui",
    "testParameterImport",
    "testSpectrumCalculation",
    "testErrors"
);
#else  // !QT_MOC_HAS_STRING_DATA
struct qt_meta_stringdata_CLASSOpenMSSCOPETestTSGDialogENDCLASS_t {
    uint offsetsAndSizes[14];
    char stringdata0[22];
    char stringdata1[17];
    char stringdata2[1];
    char stringdata3[8];
    char stringdata4[20];
    char stringdata5[24];
    char stringdata6[11];
};
#define QT_MOC_LITERAL(ofs, len) \
    uint(sizeof(qt_meta_stringdata_CLASSOpenMSSCOPETestTSGDialogENDCLASS_t::offsetsAndSizes) + ofs), len 
Q_CONSTINIT static const qt_meta_stringdata_CLASSOpenMSSCOPETestTSGDialogENDCLASS_t qt_meta_stringdata_CLASSOpenMSSCOPETestTSGDialogENDCLASS = {
    {
        QT_MOC_LITERAL(0, 21),  // "OpenMS::TestTSGDialog"
        QT_MOC_LITERAL(22, 16),  // "testConstruction"
        QT_MOC_LITERAL(39, 0),  // ""
        QT_MOC_LITERAL(40, 7),  // "testGui"
        QT_MOC_LITERAL(48, 19),  // "testParameterImport"
        QT_MOC_LITERAL(68, 23),  // "testSpectrumCalculation"
        QT_MOC_LITERAL(92, 10)   // "testErrors"
    },
    "OpenMS::TestTSGDialog",
    "testConstruction",
    "",
    "testGui",
    "testParameterImport",
    "testSpectrumCalculation",
    "testErrors"
};
#undef QT_MOC_LITERAL
#endif // !QT_MOC_HAS_STRING_DATA
} // unnamed namespace

Q_CONSTINIT static const uint qt_meta_data_CLASSOpenMSSCOPETestTSGDialogENDCLASS[] = {

 // content:
      12,       // revision
       0,       // classname
       0,    0, // classinfo
       5,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags, initial metatype offsets
       1,    0,   44,    2, 0x08,    1 /* Private */,
       3,    0,   45,    2, 0x08,    2 /* Private */,
       4,    0,   46,    2, 0x08,    3 /* Private */,
       5,    0,   47,    2, 0x08,    4 /* Private */,
       6,    0,   48,    2, 0x08,    5 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

Q_CONSTINIT const QMetaObject OpenMS::TestTSGDialog::staticMetaObject = { {
    QMetaObject::SuperData::link<QObject::staticMetaObject>(),
    qt_meta_stringdata_CLASSOpenMSSCOPETestTSGDialogENDCLASS.offsetsAndSizes,
    qt_meta_data_CLASSOpenMSSCOPETestTSGDialogENDCLASS,
    qt_static_metacall,
    nullptr,
    qt_incomplete_metaTypeArray<qt_meta_stringdata_CLASSOpenMSSCOPETestTSGDialogENDCLASS_t,
        // Q_OBJECT / Q_GADGET
        QtPrivate::TypeAndForceComplete<TestTSGDialog, std::true_type>,
        // method 'testConstruction'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        // method 'testGui'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        // method 'testParameterImport'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        // method 'testSpectrumCalculation'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        // method 'testErrors'
        QtPrivate::TypeAndForceComplete<void, std::false_type>
    >,
    nullptr
} };

void OpenMS::TestTSGDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<TestTSGDialog *>(_o);
        (void)_t;
        switch (_id) {
        case 0: _t->testConstruction(); break;
        case 1: _t->testGui(); break;
        case 2: _t->testParameterImport(); break;
        case 3: _t->testSpectrumCalculation(); break;
        case 4: _t->testErrors(); break;
        default: ;
        }
    }
    (void)_a;
}

const QMetaObject *OpenMS::TestTSGDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *OpenMS::TestTSGDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_CLASSOpenMSSCOPETestTSGDialogENDCLASS.stringdata0))
        return static_cast<void*>(this);
    return QObject::qt_metacast(_clname);
}

int OpenMS::TestTSGDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 5)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 5;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 5)
            *reinterpret_cast<QMetaType *>(_a[0]) = QMetaType();
        _id -= 5;
    }
    return _id;
}
QT_WARNING_POP
