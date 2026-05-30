#!/bin/bash
for wf in basic_peptide_*.knwf
do
 "${1}" -consoleLog --launcher.suppressErrors -nosave -nosplash -data "${2}" -application org.knime.product.KNIME_BATCH_APPLICATION -workflowFile="$(pwd)/${wf}" -reset -vmargs -Djava.io.tmpdir="${2}" > log_${wf}.txt
done
