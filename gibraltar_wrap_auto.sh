#$ -S bin/bash		#shell
#$ -cwd			#return job from cwd
#$ -l h_rt=100:00:00	#runtime max
#$ -l h_rss=8G		#memory req
#$ -q gpus	        #gpus
#$ -l gpus=1            #
#$ -j y
#$ -pe smp 1 		#number parrallel jobs
#  -fin optimized_geo 
#  -fin optgeo_extract.py
#  -fin *.sh
#  -fout scr    	# copy back to current dir
#  -fout infiles
#  -fout outfiles
#  -fout optimized_geo
#  -fout completejobs
##################CHOOSE A METAL ######################
fullpath="$1"

echo $fullpath
namebase=`echo $fullpath | sed "s/[.]in//"| sed "s:.*/::"`
echo "namebase = $namebase" 
generalpathn="/home/jp/redox_search"
generalpath="/home/jp/redox_search/"

echo "Begining Gibraltar geometry optimization run"
export OMP_NUM_THREADS=1
module load cuda
module load terachem/alpha-6.5
echo "gen path = $generalpath" 
echo "namebase = $namebase" 
sourcepath=$generalpathn/jobs/$namebase.in
inpath=$generalpathn/infiles/$namebase.in
opt_geo_path=$generalpathn/optimized_geo/$namebase.xyz
prog_geo_path=$generalpathn/prog_geo/$namebase.xyz

opt_geo_path=`echo "$opt_geo_path" | awk '{print tolower($0)}'`
echo $opt_geo_path

initial_geo_path=$generalpathn/initial_geo/$namebase.xyz
outpath=$generalpathn/geo_outfiles/$namebase.out
completepath=$generalpathin/completejobs/$namebase.in
scrpath=scr/geo/$namebase
geoextr=$generalpathn/optgeo_extract.py
echo "inpath is $inpath"
echo "Initializing local run, finding input files..."
mkdir -p scr
mkdir - scr/geo
spacer='_'
echo "begining"
echo "file is  $namebase"
echo "this current home: $HOME"
echo "outpath is $outpath"
echo "optGpath is $opt_geo_path"

#echo "this current env: $SGE_JOB_SPOOL_DIR"
wf_guess_flag=0
##begin geo-optimization

coordfile=$initial_geo_path
echo $coordfile
if [ -e $prog_geo_path ]; then
	echo "restarting from previously optimized geo"
	coordfile=$prog_geo_path #write for continutation in alpha
	guess_opt="scr/geo/$namebase/ca0 scr/geo/$namebase/cb0"
	wf_guess_flag=1
	echo "Since there is no hand-on guess, and optgeo exists \n"
	echo "I will try to use the scr value. \n"
fi

if [ $wf_guess_flag -eq 0 ]; then ## see if we load in a guess file
	guess_opt="generate"
	echo "wf from scratch"
fi
if [ -e $inpath ]; then
	rm $inpath
fi
echo "copying from $sourcepath to $inpath"
cp $sourcepath $inpath 
cat >> $inpath <<-EOF
	coordinates $coordfile
	scrdir $scrpath/
	guess $guess_opt
	end
EOF
echo "inpath is $inpath"
echo "Launching geo calc: $namebase"
terachem $inpath >  $outpath
msg=$?
if [ -e $scrpath ]; then
	echo "we found scr: $scrpath"
fi
stringtotest="$scrpath/optim.xyz"
echo "TC return status = $msg"
if [ $msg -eq 0 ]; then	
	python $geoextr $stringtotest $opt_geo_path			
	echo "the file is $namebase , the target is $completepath"
	echo "removing job $namebase"
	echo " the job was SUCCESSFUL! TC = 0 "
#	cp $sourcepath $completepath 
	printf "\n"
fi
if [ -e $stringtotest ]; then	
	python $geoextr $stringtotest $prog_geo_path			
	printf "\n"
	echo "extracting prog geo"
fi


echo "Complete"

