for files in `/bin/ls completejobs/*.in`; do	
	echo "**********\n"
	namebase=`echo $files | sed "s/[.]in//"| sed "s:.*/::"`
	echo "this name is  = $namebase" 
	generalpathn="/home/jp/redox_search"
	generalpath="/home/jp/redox_search/"
	sourcepath=$generalpathn/jobs/$namebase.in
	inpath=$generalpathn/infiles/$namebase.in
	opt_geo_path=$generalpathn/optimized_geo/$namebase.xyz
	prog_geo_path=$generalpathn/prog_geo/$namebase.xyz
	opt_geo_path=`echo "$opt_geo_path" | awk '{print tolower($0)}'`
	echo "optgeo is  $opt_geo_path"

	initial_geo_path=$generalpathn/initial_geo/$namebase.xyz
	outpath=$generalpathn/geo_outfiles/$namebase.out
	completepath=$generalpathin/completejobs/$namebase.in
	scrpath=scr/geo/$namebase
	geoextr=$generalpathn/optgeo_extract.py
	stringtotest="$scrpath/optim.xyz"
	if [ -e $stringtotest ]; then	
		python $geoextr $stringtotest $opt_geo_path			
		echo "extracting opt geo"
	else
		echo " error! $stringtoest not found - caution needed"
	fi	

done

