import pyopenms

def testCVTermList():

    term = pyopenms.CVTerm()

    term.setAccession("ACC")
    assert term.getAccession() == "ACC"

    term.setName("NAME")
    assert term.getName() == "NAME"

    term.setCVIdentifierRef("CVREF")
    assert term.getCVIdentifierRef() == "CVREF"

    term.setValue(pyopenms.DataValue(123))
    assert term.getValue().toInt() == 123

    li = pyopenms.CVTermList()

    li.setCVTerms([term])

    (name, (term,)), =  li.getCVTerms().items()

    assert name == "ACC"
    assert term.getName() == "NAME"
    assert term.getCVIdentifierRef() == "CVREF"
    assert term.getValue().toInt() == 123

    li.replaceCVTerm(term)
    (name, (term,)), =  li.getCVTerms().items()

    assert name == "ACC"
    assert term.getName() == "NAME"
    assert term.getCVIdentifierRef() == "CVREF"
    assert term.getValue().toInt() == 123

    li.replaceCVTerms([term], "ACC2")

    dd =  li.getCVTerms()
    assert len(dd) == 2
    assert "ACC" in dd and "ACC2" in dd
    term, = dd["ACC"]
    assert term.getName() == "NAME"
    assert term.getCVIdentifierRef() == "CVREF"
    assert term.getValue().toInt() == 123
    term, = dd["ACC2"]
    assert term.getName() == "NAME"
    assert term.getCVIdentifierRef() == "CVREF"
    assert term.getValue().toInt() == 123

    li.replaceCVTerms(li.getCVTerms())
    dd =  li.getCVTerms()
    assert len(dd) == 2
    assert "ACC" in dd and "ACC2" in dd
    term, = dd["ACC"]
    assert term.getName() == "NAME"
    assert term.getCVIdentifierRef() == "CVREF"
    assert term.getValue().toInt() == 123
    term, = dd["ACC2"]
    assert term.getName() == "NAME"
    assert term.getCVIdentifierRef() == "CVREF"
    assert term.getValue().toInt() == 123

    li.addCVTerm(term)
    dd = li.getCVTerms()

    assert len(dd) == 2
    assert "ACC" in dd and "ACC2" in dd

    term1, term2, = dd["ACC"]
    assert term1.getName() == "NAME"
    assert term1.getCVIdentifierRef() == "CVREF"
    assert term1.getValue().toInt() == 123

    assert term2.getName() == "NAME"
    assert term2.getCVIdentifierRef() == "CVREF"
    assert term2.getValue().toInt() == 123

    term, = dd["ACC2"]
    assert term.getName() == "NAME"
    assert term.getCVIdentifierRef() == "CVREF"
    assert term.getValue().toInt() == 123

    assert li.hasCVTerm("ACC")
    assert not li.hasCVTerm("ABC")
