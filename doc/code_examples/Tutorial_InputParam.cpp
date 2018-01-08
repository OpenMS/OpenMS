//! [InputParam]

ExitCodes main_(int, const char **)
{
  //-------------------------------------------------------------
  // parsing parameters
  //-------------------------------------------------------------
  String in(getStringOption_("in"));  // read the database filename
  String ids(getStringOption_("accession")); // read the identification filename 
  String method(getStringOption_("method")); // check if user wants black- or whitelisting
  bool whitelist = (method == "whitelist"); // and store this in a boolean variable
  String out(getStringOption_("out")); // read the output filename

//! [InputParam]
