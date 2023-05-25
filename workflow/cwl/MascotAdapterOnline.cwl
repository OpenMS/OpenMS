inputs:
  in:
    doc: "input file in mzML format.\n"
    type: File
  out:
    doc: "output file in idXML format.\n"
    type: string
  log:
    doc: Name of log file (created only when specified)
    type: string?
  debug:
    doc: Sets the debug level
    type: long?
  threads:
    doc: Sets the number of threads allowed to be used by the TOPP tool
    type: long?
  no_progress:
    doc: Disables progress logging to command line
    type: boolean?
  force:
    doc: Overrides tool-specific checks
    type: boolean?
  Mascot_parameters__database:
    doc: Name of the sequence database
    type: string?
  Mascot_parameters__search_type:
    doc: Name of the search type for the query
    type: string?
  Mascot_parameters__enzyme:
    doc: The enzyme descriptor to the enzyme used for digestion. (Trypsin is default, None would be best for peptide input or unspecific digestion, for more please refer to your mascot server).
    type: string?
  Mascot_parameters__instrument:
    doc: Instrument definition which specifies the fragmentation rules
    type: string?
  Mascot_parameters__missed_cleavages:
    doc: Number of missed cleavages allowed for the enzyme
    type: long?
  Mascot_parameters__precursor_mass_tolerance:
    doc: Tolerance of the precursor peaks
    type: double?
  Mascot_parameters__precursor_error_units:
    doc: Units of the precursor mass tolerance
    type: string?
  Mascot_parameters__fragment_mass_tolerance:
    doc: Tolerance of the peaks in the fragment spectrum
    type: double?
  Mascot_parameters__fragment_error_units:
    doc: Units of the fragment peaks tolerance
    type: string?
  Mascot_parameters__charges:
    doc: Charge states to consider, given as a comma separated list of integers (only used for spectra without precursor charge information)
    type: string?
  Mascot_parameters__taxonomy:
    doc: Taxonomy specification of the sequences
    type: string?
  Mascot_parameters__fixed_modifications:
    doc: List of fixed modifications, according to UniMod definitions.
    type: string[]?
  Mascot_parameters__variable_modifications:
    doc: Variable modifications given as UniMod definitions.
    type: string[]?
  Mascot_parameters__special_modifications:
    doc: Modifications with specificity groups that are used by Mascot and have to be treated specially
    type: string?
  Mascot_parameters__mass_type:
    doc: Defines the mass type, either monoisotopic or average
    type: string?
  Mascot_parameters__number_of_hits:
    doc: Number of hits which should be returned, if 0 AUTO mode is enabled.
    type: long?
  Mascot_parameters__skip_spectrum_charges:
    doc: Sometimes precursor charges are given for each spectrum but are wrong, setting this to 'true' does not write any charge information to the spectrum, the general charge information is however kept.
    type: boolean?
  Mascot_parameters__decoy:
    doc: Set to true if mascot should generate the decoy database.
    type: boolean?
  Mascot_parameters__search_title:
    doc: Sets the title of the search.
    type: string?
  Mascot_parameters__username:
    doc: Sets the username which is mentioned in the results file.
    type: string?
  Mascot_parameters__email:
    doc: "Sets the email which is mentioned in the results file. Note: Some server require that a proper email is provided."
    type: string?
  Mascot_server__hostname:
    doc: Address of the host where Mascot listens, e.g. 'mascot-server' or '127.0.0.1'
    type: string?
  Mascot_server__host_port:
    doc: Port where the Mascot server listens, 80 should be a good guess
    type: long?
  Mascot_server__server_path:
    doc: Path on the server where Mascot server listens, 'mascot' should be a good guess
    type: string?
  Mascot_server__timeout:
    doc: Timeout in seconds, after which the query is declared as failed.This is NOT the whole time the search takes, but the time in between two progress steps. Some Mascot servers freeze during this (unstable network etc) and idle forever, the connection is killed. Set this to 0 to disable timeout!
    type: long?
  Mascot_server__boundary:
    doc: Boundary for the MIME section
    type: string?
  Mascot_server__use_proxy:
    doc: Flag which enables the proxy usage for the HTTP requests, please specify at least 'proxy_host' and 'proxy_port'
    type: boolean?
  Mascot_server__proxy_host:
    doc: Host where the proxy server runs on
    type: string?
  Mascot_server__proxy_port:
    doc: Port where the proxy server listens
    type: long?
  Mascot_server__proxy_username:
    doc: Login name for the proxy server, if needed
    type: string?
  Mascot_server__proxy_password:
    doc: Login password for the proxy server, if needed
    type: string?
  Mascot_server__login:
    doc: Flag which should be set 'true' if Mascot security is enabled; also set 'username' and 'password' then.
    type: boolean?
  Mascot_server__username:
    doc: Name of the user if login is used (Mascot security must be enabled!)
    type: string?
  Mascot_server__password:
    doc: Password of the user if login is used (Mascot security must be enabled!)
    type: string?
  Mascot_server__use_ssl:
    doc: Flag indicating whether you want to send requests to an HTTPS server or not (HTTP). Requires OpenSSL to be installed (see openssl.org)
    type: boolean?
  Mascot_server__export_params:
    doc: Adjustable export parameters (passed to Mascot's 'export_dat_2.pl' script). Generally only parameters that control which hits to export are safe to adjust/add. Many settings that govern what types of information to include are required by OpenMS and cannot be changed. Note that setting 'query_master' to 1 may lead to incorrect protein references for peptides.
    type: string?
  Mascot_server__skip_export:
    doc: "For use with an external Mascot Percolator (via GenericWrapper): Run the Mascot search, but do not export the results. The output file produced by MascotAdapterOnline will contain only the Mascot search number."
    type: boolean?
  Mascot_server__batch_size:
    doc: Number of spectra processed in one batch by Mascot (default 50000)
    type: long?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: MascotAdapterOnline
doc: Annotates MS/MS spectra using Mascot.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - MascotAdapterOnline
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json