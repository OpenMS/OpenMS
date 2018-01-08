//! [Register]

void registerOptionsAndFlags_()
{
  registerInputFile_("in", "<file>", "","Input FASTA file, containing a database.");
  setValidFormats_("in", ListUtils::create<String>("fasta"));

  registerInputFile_("accession", "<file>", "","Input IdXML file, containing the identified peptides.", true);
  setValidFormats_("accession", ListUtils::create<String>("idXML"));
 
  registerStringOption_("method", "<choice>", "whitelist", "Switch between white/blacklisting", false);
  setValidStrings_("method", ListUtils::create<String>("whitelist, blacklist"));

  registerOutputFile_("out", "<file>", "", "Output FASTA file where the reduced database will be written to.");
  setValidFormats_("out", ListUtils::create<String>("fasta"));
}

//! [Register]
