inputs:
  in:
    doc: Input FASTA file(s), each containing a database. It is recommended to include a contaminant database as well.
    type: File[]
  out:
    doc: Output FASTA file where the decoy database will be written to.
    type: string
  decoy_string:
    doc: String that is combined with the accession of the protein identifier to indicate a decoy protein.
    type: string?
  decoy_string_position:
    doc: Should the 'decoy_string' be prepended (prefix) or appended (suffix) to the protein accession?
    type: string?
  only_decoy:
    doc: Write only decoy proteins to the output database instead of a combined database.
    type: boolean?
  type:
    doc: Type of sequence. RNA sequences may contain modification codes, which will be handled correctly if this is set to 'RNA'.
    type: string?
  method:
    doc: Method by which decoy sequences are generated from target sequences. Note that all sequences are shuffled using the same random seed, ensuring that identical sequences produce the same shuffled decoy sequences. Shuffled sequences that produce highly similar output sequences are shuffled again (see shuffle_sequence_identity_threshold).
    type: string?
  shuffle_max_attempts:
    doc: "shuffle: maximum attempts to lower the amino acid sequence identity between target and decoy for the shuffle algorithm"
    type: long?
  shuffle_sequence_identity_threshold:
    doc: "shuffle: target-decoy amino acid sequence identity threshold for the shuffle algorithm. If the sequence identity is above this threshold, shuffling is repeated. In case of repeated failure, individual amino acids are 'mutated' to produce a different amino acid sequence."
    type: double?
  seed:
    doc: Random number seed (use 'time' for system time)
    type: string?
  enzyme:
    doc: Enzyme used for the digestion of the sample. Only applicable if parameter 'type' is 'protein'.
    type: string?
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
  Decoy__non_shuffle_pattern:
    doc: Residues to not shuffle (keep at a constant position when shuffling). Separate by comma, e.g. use 'K,P,R' here.
    type: string?
  Decoy__keepPeptideNTerm:
    doc: Whether to keep peptide N terminus constant when shuffling / reversing.
    type: string?
  Decoy__keepPeptideCTerm:
    doc: Whether to keep peptide C terminus constant when shuffling / reversing.
    type: string?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: DecoyDatabase
doc: Create decoy sequence database from forward sequence database.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - DecoyDatabase
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json