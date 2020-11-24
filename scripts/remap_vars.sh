#!/bin/bash
for i in /home/pelmo/sscrofa_variant_remapping/data/var_coords/*.vars ; do 
	test -f $i.annot || perl /home/pelmo/sscrofa_variant_remapping/scripts/remap_api.pl --mode asm-asm --annotation $i --from GCF_000003025.5 --dest GCF_000003025.6 --annot_out $i.annot --report_out $i.report

done
