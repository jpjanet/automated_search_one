#! /bin/bash
for files in `/bin/ls jobs/*.in`; do
	qsub sol_gib_wrap.sh $files 
done
