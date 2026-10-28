if (NOT DEFINED INPUT_FILE)
  message(FATAL_ERROR "INPUT_FILE is required")
endif()

file(READ "${INPUT_FILE}" contents)

# With the RNA-UV (UCGA) preset, several precursor adducts share a sum formula (e.g., C-H3N1 and U-H2O1)
# and dinucleotide adducts have more than one cross-linkable nucleotide (e.g., A and U for AU).
# Every adduct and nucleotide must be scored, not only the first of each.

# Returns the PeptideHit element with the given sequence attribute, or fails if there is none.
function(get_hit sequence out_var)
  string(FIND "${contents}" "sequence=\"${sequence}\"" hit_start)
  if (hit_start EQUAL -1)
    message(FATAL_ERROR "No PeptideHit with sequence ${sequence}")
  endif()
  string(SUBSTRING "${contents}" ${hit_start} -1 hit)
  string(FIND "${hit}" "</PeptideHit>" hit_end)
  string(SUBSTRING "${hit}" 0 ${hit_end} hit)
  set(${out_var} "${hit}" PARENT_SCOPE)
endfunction()

function(expect_user_param hit sequence name value)
  string(FIND "${hit}" "name=\"${name}\" value=\"${value}\"" position)
  if (position EQUAL -1)
    message(FATAL_ERROR "PeptideHit ${sequence} does not have ${name} = ${value}")
  endif()
endfunction()

# The second adduct of the sum formula (U-H2O1 after C-H3N1) explains this spectrum best.
get_hit("AEADNLDDK(NuXL:U-H2O1)K" hit)
expect_user_param("${hit}" "AEADNLDDK(NuXL:U-H2O1)K" "NuXL:NA" "U-H2O1")
expect_user_param("${hit}" "AEADNLDDK(NuXL:U-H2O1)K" "NuXL:NT" "U")
# Candidates are counted once per peptide, sum formula and spectrum, not once per adduct and nucleotide.
expect_user_param("${hit}" "AEADNLDDK(NuXL:U-H2O1)K" "nr_candidates" "2")

# The second cross-linkable nucleotide of AU (U after A) explains this spectrum best.
get_hit("YTPY(NuXL:AU)LGR" hit)
expect_user_param("${hit}" "YTPY(NuXL:AU)LGR" "NuXL:NT" "U")

# C-H3N1 and U-H2O1 score equally on this spectrum. The tie goes to the first adduct of the sum formula,
# and the equally scoring copy with the other adduct does not take the second of the two reported hits.
get_hit("DDMNMQK(NuXL:C-H3N1)DAIFGGK" hit)
string(FIND "${contents}" "sequence=\"DDMNMQK(NuXL:U-H2O1)DAIFGGK\"" copy_position)
if (NOT copy_position EQUAL -1)
  message(FATAL_ERROR "DDMNMQK(NuXL:U-H2O1)DAIFGGK is reported next to its equally scoring C-H3N1 copy")
endif()
