#!/bin/bash

# settings from input

width=${1:-14}
boundary=${2:-1}
liebdef=${3:-31}
flux0=${4:-1.0}  # This is a disorder loop with fix energy(in flux loop below) 
flux1=${5:-1.0}
dflux=${6:-0.1}

# settings for files

binary=tmseLMxD.GF

# settings for directories

currdir=`pwd`


binarydir=$currdir/../EXE
[ -d ${binarydir} ] || mkdir ${binarydir}

cp $currdir/../src/$binary $binarydir 


jobdir="TM-L$liebdef-M$width"
[ -d $jobdir ] || mkdir $jobdir

cd $jobdir
echo "present workdir " ${currdir}/${jobdir}

#if [ $(echo "$fluxFlag==1" | bc) =1 ]; then

#    echo "It is energy loop with fix disorder"
#else
#    echo "It is disorder loop with fix energy"
#fi


for flux in 15.5 16.0 16.5 17.0 17.5 # this is energy value 
do
    
jobname="L$liebdef-M$width-E$flux"
echo $jobname

jobfile=`printf "$jobname.sh"`
logfile=`printf "$jobname.log"`

inpfile=tmseLMxD-E$flux.inp
echo "binarydir=" $binarydir " jobdir=" $jobdir


# settings for parallel submission

cat > ${jobfile} << EOD
#!/bin/bash
#PBS -l nodes=${nodes}:ppn=16
#PBS -l pmem=${memory}
#PBS -l walltime=04:00:00

#       The jobname
#PBS -N ${jobname}

##############################
## comment out when INTERACTIVE
## SLURM_NTASKS=8
## SLURM_JOBID=1
##############################

# construct the input file

touch $inpfile

echo "ISeed         = 1277      "> $inpfile
echo "NOfIter       = 100000000.">> $inpfile
echo "NOfOrtho      = 10        ">> $inpfile
echo "NOfPrint      = 100000    ">> $inpfile
echo "NOfGamma      = 1         ">> $inpfile
echo "IDimenFlag    = $liebdef  ">> $inpfile
echo "IBCFlag       = $boundary ">> $inpfile
echo "IRNGFlag      = 10      ">> $inpfile
echo "IKeepFlag     = 0      ">> $inpfile
echo "IWriteFlag    = 2      ">> $inpfile
echo "ISortFlag     = 0      ">> $inpfile
echo "IFluxFlag     = 0      ">> $inpfile
echo "Width0        = $width ">> $inpfile
echo "Width1        = $width ">> $inpfile
echo "dWidth        = 2      ">> $inpfile
echo "DiagDis0      = $flux0 ">> $inpfile
echo "DiagDis1      = $flux1 ">> $inpfile
echo "dDiagDis      = $dflux ">> $inpfile
echo "Energy0       = $flux  ">> $inpfile
echo "Energy1       = $flux  ">> $inpfile
echo "dEnergy       = 0.02   ">> $inpfile
echo "Kappa         = 1.0    ">> $inpfile
echo "Epsilon       = 1.0E-2 ">> $inpfile

cat $inpfile

${binarydir}/${binary} <${inpfile} >& ${logfile}
#${binarydir}/${binary} <$inpfile >& ${logfile}

wait
#exit 0
EOD

chmod 755 ${currdir}/${jobdir}/${jobfile}
#(cd ${jobdir} ; msub -q devel ./${jobfile}) # for development/testing only
(source ${currdir}/${jobdir}/${jobfile} ) & # for PRODUCTION

sleep 1

done

wait
#done

