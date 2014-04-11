import pyopenms

def testCVTermList():

    term = pyopenms.CVTerm()

    term.setAccession(b"ACC")
    assert term.getAccession() == b"ACC"

    term.setName(b"NAME")
    assert term.getName() == b"NAME"

    term.setCVIdentifierRef(b"CVREF")
    assert term.getCVIdentifierRef() == b"CVREF"

    term.setValue(123)
    assert term.getValue() == 123

    li = pyopenms.CVTermList()

    li.setCVTerms([term])

    (name, (term,)), =  li.getCVTerms().items()

    assert name == b"ACC"
    assert term.getName() == b"NAME"
    assert term.getCVIdentifierRef() == b"CVREF"
    assert term.getValue() == 123

    li.replaceCVTerm(term)
    (name, (term,)), =  li.getCVTerms().items()

    assert name == b"ACC"
    assert term.getName() == b"NAME"
    assert term.getCVIdentifierRef() == b"CVREF"
    assert term.getValue() == 123

    li.replaceCVTerms([term], b"ACC2")

    dd =  li.getCVTerms()
    assert len(dd) == 2
    assert b"ACC" in dd and b"ACC2" in dd
    term, = dd[b"ACC"]
    assert term.getName() == b"NAME"
    assert term.getCVIdentifierRef() == b"CVREF"
    assert term.getValue() == 123
    term, = dd[b"ACC2"]
    assert term.getName() == b"NAME"
    assert term.getCVIdentifierRef() == b"CVREF"
    assert term.getValue() == 123

    li.replaceCVTerms(li.getCVTerms())
    dd =  li.getCVTerms()
    assert len(dd) == 2
    assert b"ACC" in dd and b"ACC2" in dd
    term, = dd[b"ACC"]
    assert term.getName() == b"NAME"
    assert term.getCVIdentifierRef() == b"CVREF"
    assert term.getValue() == 123
    term, = dd[b"ACC2"]
    assert term.getName() == b"NAME"
    assert term.getCVIdentifierRef() == b"CVREF"
    assert term.getValue() == 123

    li.addCVTerm(term)
    dd = li.getCVTerms()

    assert len(dd) == 2
    assert b"ACC" in dd and b"ACC2" in dd

    term1, term2, = dd[b"ACC"]
    assert term1.getName() == b"NAME"
    assert term1.getCVIdentifierRef() == b"CVREF"
    assert term1.getValue() == 123

    assert term2.getName() == b"NAME"
    assert term2.getCVIdentifierRef() == b"CVREF"
    assert term2.getValue() == 123

    term, = dd[b"ACC2"]
    assert term.getName() == b"NAME"
    assert term.getCVIdentifierRef() == b"CVREF"
    assert term.getValue() == 123

    assert li.hasCVTerm(b"ACC")
    assert not li.hasCVTerm(b"ABC")
