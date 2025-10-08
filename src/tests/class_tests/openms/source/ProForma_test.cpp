#include <iostream>
#include <OpenMS/CHEMISTRY/ProForma.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

void test_parseSimpleModification()
{
  OpenMS::AASequence seq = OpenMS::AASequence::fromString("ACDEFGHIK");
  OpenMS::ProForma proforma(seq);

  std::string proforma_str = "A[Phospho]CDEFGHIK";
  proforma.fromProFormaString(proforma_str);
  std::string output = proforma.toProFormaString();
  if (output == proforma_str)
  {
    std::cout << "test_parseSimpleModification PASSED" << std::endl;
  }
  else
  {
    std::cerr << "test_parseSimpleModification FAILED: Expected " << proforma_str << " but got " << output << std::endl;
  }
}

void test_parseMassShiftModification()
{
  OpenMS::AASequence seq = OpenMS::AASequence::fromString("ACDEFGHIK");
  OpenMS::ProForma proforma(seq);

  std::string proforma_str = "A[+15.99]CDEFGHIK";
  proforma.fromProFormaString(proforma_str);
  std::string output = proforma.toProFormaString();
  if (output == proforma_str)
  {
    std::cout << "test_parseMassShiftModification PASSED" << std::endl;
  }
  else
  {
    std::cerr << "test_parseMassShiftModification FAILED: Expected " << proforma_str << " but got " << output << std::endl;
  }
}

void test_parseNTerminalModification()
{
  OpenMS::AASequence seq = OpenMS::AASequence::fromString("ACDEFGHIK");
  OpenMS::ProForma proforma(seq);

  std::string proforma_str = "[Acetyl]-ACDEFGHIK";
  proforma.fromProFormaString(proforma_str);
  std::string output = proforma.toProFormaString();
  if (output == proforma_str)
  {
    std::cout << "test_parseNTerminalModification PASSED" << std::endl;
  }
  else
  {
    std::cerr << "test_parseNTerminalModification FAILED: Expected " << proforma_str << " but got " << output << std::endl;
  }
}

void test_parseCTerminalModification()
{
  OpenMS::AASequence seq = OpenMS::AASequence::fromString("ACDEFGHIK");
  OpenMS::ProForma proforma(seq);

  std::string proforma_str = "ACDEFGHIK-[Amidation]";
  proforma.fromProFormaString(proforma_str);
  std::string output = proforma.toProFormaString();
  if (output == proforma_str)
  {
    std::cout << "test_parseCTerminalModification PASSED" << std::endl;
  }
  else
  {
    std::cerr << "test_parseCTerminalModification FAILED: Expected " << proforma_str << " but got " << output << std::endl;
  }
}

void test_parseAmbiguousModification()
{
  OpenMS::AASequence seq = OpenMS::AASequence::fromString("ACDEFGHIK");
  OpenMS::ProForma proforma(seq);

  std::string proforma_str = "A[+15.99]?CDEFGHIK";
  proforma.fromProFormaString(proforma_str);
  std::string output = proforma.toProFormaString();
  if (output == proforma_str)
  {
    std::cout << "test_parseAmbiguousModification PASSED" << std::endl;
  }
  else
  {
    std::cerr << "test_parseAmbiguousModification FAILED: Expected " << proforma_str << " but got " << output << std::endl;
  }
}

void test_invalidModificationFormat()
{
  OpenMS::AASequence seq = OpenMS::AASequence::fromString("ACDEFGHIK");
  OpenMS::ProForma proforma(seq);

  std::string invalid_proforma_str = "A[PhosphoCDEFGHIK"; // Missing closing bracket
  try
  {
    proforma.fromProFormaString(invalid_proforma_str);
    std::cerr << "test_invalidModificationFormat FAILED: No exception thrown" << std::endl;
  }
  catch (const std::runtime_error&)
  {
    std::cout << "test_invalidModificationFormat PASSED" << std::endl;
  }
  catch (const std::exception& e)
  {
    std::cerr << "test_invalidModificationFormat FAILED with unexpected exception: " << e.what() << std::endl;
  }
}

void test_unsupportedCV()
{
  OpenMS::AASequence seq = OpenMS::AASequence::fromString("ACDEFGHIK");
  OpenMS::ProForma proforma(seq);

  std::string invalid_proforma_str = "A[XYZ:12345]CDEFGHIK"; // Unsupported CV
  try
  {
    proforma.fromProFormaString(invalid_proforma_str);
    std::cerr << "test_unsupportedCV FAILED: No exception thrown" << std::endl;
  }
  catch (const std::invalid_argument&)
  {
    std::cout << "test_unsupportedCV PASSED" << std::endl;
  }
  catch (const std::exception& e)
  {
    std::cerr << "test_unsupportedCV FAILED with unexpected exception: " << e.what() << std::endl;
  }
}

void test_removeModification()
{
  OpenMS::AASequence seq = OpenMS::AASequence::fromString("ACDEFGHIK");
  OpenMS::ProForma proforma(seq);

  std::string proforma_str = "A[Phospho]CDEFGHIK";
  proforma.fromProFormaString(proforma_str);
  proforma.removeModification(1);
  std::string output = proforma.toProFormaString();
  if (output == "ACDEFGHIK")
  {
    std::cout << "test_removeModification PASSED" << std::endl;
  }
  else
  {
    std::cerr << "test_removeModification FAILED: Expected ACDEFGHIK but got " << output << std::endl;
  }
}

void test_validProFormaPSIMOD()
{
  OpenMS::AASequence seq = OpenMS::AASequence::fromString("EMEVEESPEK");
  OpenMS::ProForma proforma(seq);

  std::string proforma_str = "EM[MOD:00719]EVEES[MOD:00046]PEK";
  proforma.fromProFormaString(proforma_str);
  std::string output = proforma.toProFormaString();
  if (output == proforma_str)
  {
    std::cout << "test_validProFormaPSIMOD PASSED" << std::endl;
  }
  else
  {
    std::cerr << "test_validProFormaPSIMOD FAILED: Expected " << proforma_str << " but got " << output << std::endl;
  }
}

void test_validProFormaUNIMOD()
{
  OpenMS::AASequence seq = OpenMS::AASequence::fromString("EMEVEESPEK");
  OpenMS::ProForma proforma(seq);

  std::string proforma_str = "EM[UNIMOD:35]EVEES[UNIMOD:56]PEK";
  proforma.fromProFormaString(proforma_str);
  std::string output = proforma.toProFormaString();
  if (output == proforma_str)
  {
    std::cout << "test_validProFormaUNIMOD PASSED" << std::endl;
  }
  else
  {
    std::cerr << "test_validProFormaUNIMOD FAILED: Expected " << proforma_str << " but got " << output << std::endl;
  }
}

void test_validProFormaRESID()
{
  OpenMS::AASequence seq = OpenMS::AASequence::fromString("EMEVEESPEK");
  OpenMS::ProForma proforma(seq);

  std::string proforma_str = "EM[RESID:AA0581]EVEES[RESID:AA0037]PEK";
  proforma.fromProFormaString(proforma_str);
  std::string output = proforma.toProFormaString();
  if (output == proforma_str)
  {
    std::cout << "test_validProFormaRESID PASSED" << std::endl;
  }
  else
  {
    std::cerr << "test_validProFormaRESID FAILED: Expected " << proforma_str << " but got " << output << std::endl;
  }
}

void test_invalidProFormaPSIMODSyntax()
{
  OpenMS::AASequence seq = OpenMS::AASequence::fromString("EMEVEESPEK");
  OpenMS::ProForma proforma(seq);

  std::string invalid_proforma_str = "EM[M:00719]EVEES[M:00046]PEK";
  try
  {
    proforma.fromProFormaString(invalid_proforma_str);
    std::cerr << "test_invalidProFormaPSIMODSyntax FAILED: No exception thrown" << std::endl;
  }
  catch (const std::invalid_argument&)
  {
    std::cout << "test_invalidProFormaPSIMODSyntax PASSED" << std::endl;
  }
  catch (const std::exception& e)
  {
    std::cerr << "test_invalidProFormaPSIMODSyntax FAILED with unexpected exception: " << e.what() << std::endl;
  }
}

void test_invalidProFormaUNIMODSyntax()
{
  OpenMS::AASequence seq = OpenMS::AASequence::fromString("EMEVEESPEK");
  OpenMS::ProForma proforma(seq);

  std::string invalid_proforma_str = "EM[U:35]EVEES[U:56]PEK";
  try
  {
    proforma.fromProFormaString(invalid_proforma_str);
    std::cerr << "test_invalidProFormaUNIMODSyntax FAILED: No exception thrown" << std::endl;
  }
  catch (const std::invalid_argument&)
  {
    std::cout << "test_invalidProFormaUNIMODSyntax PASSED" << std::endl;
  }
  catch (const std::exception& e)
  {
    std::cerr << "test_invalidProFormaUNIMODSyntax FAILED with unexpected exception: " << e.what() << std::endl;
  }
}

void test_invalidProFormaRESIDSyntax()
{
  OpenMS::AASequence seq = OpenMS::AASequence::fromString("EMEVEESPEK");
  OpenMS::ProForma proforma(seq);

  std::string invalid_proforma_str = "EM[R:AA0581]EVEES[R:AA0037]PEK";
  try
  {
    proforma.fromProFormaString(invalid_proforma_str);
    std::cerr << "test_invalidProFormaRESIDSyntax FAILED: No exception thrown" << std::endl;
  }
  catch (const std::invalid_argument&)
  {
    std::cout << "test_invalidProFormaRESIDSyntax PASSED" << std::endl;
  }
  catch (const std::exception& e)
  {
    std::cerr << "test_invalidProFormaRESIDSyntax FAILED with unexpected exception: " << e.what() << std::endl;
  }
}

void test_parseRangeModification()
{
  OpenMS::AASequence seq = OpenMS::AASequence::fromString("ACDEFGHIK");
  OpenMS::ProForma proforma(seq);

  std::string proforma_str = "A(CDE)[+12.5]FGHIK";
  proforma.fromProFormaString(proforma_str);
  std::string output = proforma.toProFormaString();
  if (output == proforma_str)
  {
    std::cout << "test_parseRangeModification PASSED" << std::endl;
  }
  else
  {
    std::cerr << "test_parseRangeModification FAILED: Expected " << proforma_str << " but got " << output << std::endl;
  }
}

int main()
{
  test_parseSimpleModification();
  test_parseMassShiftModification();
  test_parseNTerminalModification();
  test_parseCTerminalModification();
  test_parseAmbiguousModification();
  test_invalidModificationFormat();
  test_unsupportedCV();
  test_removeModification();

  test_validProFormaPSIMOD();
  test_validProFormaUNIMOD();
  test_validProFormaRESID();
  test_invalidProFormaPSIMODSyntax();
  test_invalidProFormaUNIMODSyntax();
  test_invalidProFormaRESIDSyntax();
  test_parseRangeModification();

  return 0;
}

