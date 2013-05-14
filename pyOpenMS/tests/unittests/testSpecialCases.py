import pyopenms

def testMSSpectrumIter():
    # did raise segfault because of missing index check in __getitem__
    # should be resolved by implementing a proper iterator:
    s = pyopenms.MSSpectrum()
    assert list(s) == []

def testParamEntry():
    # as ParamEntry::isValid takes "String &" as input argument, which
    # can not be implemened by a Python string, here no automatic
    # conversion from a basestring should happen:
    p = pyopenms.ParamEntry()
    message = pyopenms.String()
    assert p.isValid(message)
    assert message.c_str() == ""


