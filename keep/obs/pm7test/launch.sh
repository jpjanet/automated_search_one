#! /bin/bash
for files in `/bin/ls jobs/*.in`; do
	qsub gib_wrap.sh $files 
done
