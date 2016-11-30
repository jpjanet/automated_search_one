#$ -S bin/bash		#shell
#$ -cwd			#return job from cwd
#$ -l h_rt=100:00:00	#runtime max
#$ -l h_rss=8G		#memory req
#$ -q gpus	        #gpus
#$ -l gpus=1            #
#$ -pe smp 2 		#number parrallel jobs
#$ -N solt_tests
#  -fin optimized_geo 
#  -fin *.py
#  -fin pcm_radii
#  -fout scr    	# copy back to current dir
#  -fout infiles
#  -fout outfiles
#  -fout optimized_geo
#  -fout completejobs
##################CHOOSE A METAL ######################
fullpath="$1"

echo $fullpath
generalpath=`echo $(dirname $fullpath) | sed "s,/*[^/]\*$,,"`
#gennumpath=$(basename $generalpath)
#generalpath=`echo $(dirname $generalpath) | sed "s,/*[^/]\*$,,"`
#generalpath=`echo $(dirname $generalpath) | sed "s,/*[^/]\*$,,"`

namebase=`echo $fullpath | sed "s/[.]in//"| sed "s:.*/::"`
echo "gen path = $generalpath" 
echo "namebase = $namebase" 
pcm_rad_path=$(readlink -e pcm_radii) # comso radii

generalpath="/home/jp/solpm7test"
echo "Begining Gibraltar geometry optimization run"
export OMP_NUM_THREADS=2
module load cuda
module load terachem/alpha-6.5
#outpath="$SGE_O_WORKDIR/outfiles/"
echo "gen path = $generalpath" 
echo "gen num = $gennumpath"
echo "namebase = $namebase" 
sourcepath=$generalpath/jobs/$namebase.in
inpath=$generalpath/infiles/$namebase.in
opt_geo_path=$generalpath/optimized_geo/$namebase.xyz
initial_geo_path=$generalpath/initial_geo/$namebase.xyz
outpath=$generalpath/outfiles/$namebase.out
completepath=$generalpath/completejobs/$namebase.in
echo "inpath is $inpath"
echo "Initializing local run, finding input files..."
mkdir -p scr
spacer='_'
echo "begining"
echo "file is  $namebase"
echo "this current home: $HOME"
echo "outpath is $outpath"

#echo "this current env: $SGE_JOB_SPOOL_DIR"
wf_guess_flag=0
##begin geo-optimization

coordfile=$initial_geo_path
echo $coordfile
if [ -e $opt_geo_path]; then
	echo "restarting from previously optimized geo"
	coordfile=$opt_geo_path #write for continutation in alpha
	guess_opt="scr/$namebase/ca0 scr/$namebase/cb0"
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
	pcm cosmo
	pcm_grid iswig
	epsilon 1.40
	print_ms yes
	pcm_radii read
	pcm_radii_file $pcm_rad_path
	new_minimizer yes
	guess $guess_opt
	end
EOF

echo "inpath is $inpath"
echo "Launching geo calc: $namebase"
terachem $inpath >  $outpath
msg=$?
scrpath=$(readlink -e scr)
echo "we found scr: $scrpath"
stringtotest="$scrpath/$namebase/optim.xyz"
echo "TC return status = $msg"
if [ $msg -eq 0 ]; then	
	python optgeo_extract.py $stringtotest $opt_geo_path			
	echo "the file is $namebase , the target is $completepath"
	echo "removing job $namebase"
	mv $sourcepath $completepath 
	printf "\n"
fi

echo "Complete"

