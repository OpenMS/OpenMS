
START_SECTION(BuiltInProteaseDataProvider - Percolator Specific Enzymes)
{
  BuiltInProteaseDataProvider provider;
  auto enzymes = provider.loadEnzymes();
  
  bool found_pepsin = false;
  bool found_elastase = false;
  bool found_gluc = false;

  for (const auto& enz : enzymes)
  {
    if (enz->getName() == "Percolator Pepsin")
    {
      found_pepsin = true;
      TEST_EQUAL(enz->getRegEx(), "(?<=[FLWY])(?!R)|(?<!R)(?=[FLWY])")
    }
    if (enz->getName() == "Percolator Elastase")
    {
      found_elastase = true;
      TEST_EQUAL(enz->getRegEx(), "(?<=[LVAG])(?!P)")
    }
    if (enz->getName() == "Percolator Glu-C")
    {
      found_gluc = true;
      TEST_EQUAL(enz->getRegEx(), "(?<=E)(?!P)")
    }
  }
  TEST_EQUAL(found_pepsin, true)
  TEST_EQUAL(found_elastase, true)
  TEST_EQUAL(found_gluc, true)
}
END_SECTION


START_SECTION(BuiltInProteaseDataProvider - Percolator Specific Enzymes)
{
  BuiltInProteaseDataProvider provider;
  auto enzymes = provider.loadEnzymes();
  
  bool found_pepsin = false;
  bool found_elastase = false;
  bool found_gluc = false;

  for (const auto& enz : enzymes)
  {
    if (enz->getName() == "Percolator Pepsin")
    {
      found_pepsin = true;
      TEST_EQUAL(enz->getRegEx(), "(?<=[FLWY])(?!R)|(?<!R)(?=[FLWY])")
    }
    if (enz->getName() == "Percolator Elastase")
    {
      found_elastase = true;
      TEST_EQUAL(enz->getRegEx(), "(?<=[LVAG])(?!P)")
    }
    if (enz->getName() == "Percolator Glu-C")
    {
      found_gluc = true;
      TEST_EQUAL(enz->getRegEx(), "(?<=E)(?!P)")
    }
  }
  TEST_EQUAL(found_pepsin, true)
  TEST_EQUAL(found_elastase, true)
  TEST_EQUAL(found_gluc, true)
}
END_SECTION

END_TEST
