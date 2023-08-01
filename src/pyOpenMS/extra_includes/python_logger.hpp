#ifndef __PYTHON_LOGSTREAM_HPP__
#define __PYTHON_LOGSTREAM_HPP__

#include <OpenMS/CONCEPT/LogStream.h>

class PyLogStream : public std::stringbuf
{
public:
    PyLogStream()
    {
        // import logging module on demand
        if (logging == NULL){
            logging = PyImport_ImportModuleNoBlock("logging");
            if (logging == NULL)
                PyErr_SetString(PyExc_ImportError,
                    "Could not import module 'logging'");
        }
    }

private:
    static PyObject *logging = NULL;
};

class PyLogWarn : public PyLogStream
{
public:
    virtual int sync() {
        static PyObject *string = NULL;

        // build msg-string
        string = Py_BuildValue("s", this->str());
        PyObject_CallMethod(logging, "warn", "O", string);
        this->clear();
    }
}

inline PyLogWarn Global_PyLogWarn_Buf();
inline std::ostream Global_PyLogWarn_Stream = std::ostream(*Global_PyLogWarn);

static void redirectLoggingToPython(){
    OpenMS::OpenMS_Log_warn.removeAllStreams();
    OpenMS::OpenMS_Log_warn.insert(*Global_PyLogWarn_Stream);

}
#endif
