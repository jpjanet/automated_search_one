#! /bin/bash
for files in `/bin/ls molpac/*.mop`; do
	namebase=`echo $files | sed "s/[.]xyz" | sed "sed:.*/::"`
	echo " $namebase is running "
	mopac $files > $namebase.out
done

