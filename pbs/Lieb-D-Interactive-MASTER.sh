#!/bin/bash

# settings from input

width=${1:-4}
#flux=${2:-10.0}
boundary=${2:-0}

echo "LM: boundary=" $boundary

# defining the number of nodes/cores to use

let "cores= 1"
let "nodes = ( $cores/16 )" # tinis 16 cores per node
#echo "nodes=" $nodes ", cores= " $cores
nodes=$(($nodes>1?$nodes:1))
nodes=$(($nodes<4?$nodes:4))
echo "nodes=" $nodes ", cores= " $cores #make sure at least 1 node is specified

# setting the memory to use

#memory="2048mb" # L=10
memory="4096mb" # L=20
echo "memory= " $memory

# settings for files

binary=tmseLMxD.GF

##############################
# comment out when INTERACTIVE
# binary=tmseLMxD.GF
##############################

inpfile=tmseLMxD.inp

# settings for directories

currdir=`pwd`

binarydir=~/Projects/LiebModelTMM/EXE

submitdir=$currdir
[ -d ${submitdir} ] || mkdir ${submitdir}

jobdir=${submitdir}

[ -d ${jobdir} ] || mkdir ${jobdir}

# compute jobs for gate sweep
# this will use individual jobs along each x-line

#for width in \
#10 #20
#do
echo "LM: loop-width=" $width

for flux in 5.0 2.0 1.0
do
echo "LM: loop-flux=" $flux

for liebdef in 21 22 23 24 31 32 33 34
do
echo "LM: loop-lieb=" $liebdef

jobname="L$liebdef-M$width-B$boundary-D$flux"
echo $jobname

jobfile=`printf "$jobname.sh"`

targetdir=${submitdir}/../DATA/${jobname}
[ -d ${targetdir} ] || mkdir ${targetdir}

tmpdir=~/RUNS/${jobname}

echo "binarydir=" $binarydir " submitdir=" $submitdir 
echo "tmpdir=" $tmpdir " targetdir=" $targetdir

# settings for parallel submission

cat > ${jobdir}/${jobfile} << EOD
#!/bin/bash
#PBS -l nodes=${nodes}:ppn=16
#PBS -l pmem=${memory}
#PBS -l walltime=04:00:00

#       The jobname
#PBS -N ${jobname}

##############################
# comment out when INTERACTIVE
SLURM_NTASKS=8
SLURM_JOBID=1
##############################

# The base directory is the directory that the job was submitted from.
basedir=\$PBS_O_WORKDIR
echo "basedir=" \${basedir}

# The binary directory is the directory where the scripts are
echo "binarydir=" ${binarydir}

# do the run in a TMPdir for safekeeping

[ -d ${tmpdir} ] || mkdir ${tmpdir}

cp ${binarydir}/${binary} ${tmpdir}/

# construct the input file

cd ${tmpdir}
echo "hostname=" $HOSTNAME
pwd

# doing to interactive runs

rm -rf $inpfile
touch $inpfile

echo "ISeed         = 1277">> $inpfile
echo "NOfIter       = 100000000.">> $inpfile
echo "NOfOrtho      = 10     ">> $inpfile
echo "NOfPrint      = 100000">> $inpfile
echo "NOfGamma      = 1">> $inpfile
echo "IDimenFlag    = $liebdef">> $inpfile
echo "IBCFlag       = $boundary"    >> $inpfile
echo "IRNGFlag      = 0"      >> $inpfile
echo "IKeepFlag     = 0">> $inpfile
echo "IWriteFlag    = 2">> $inpfile
echo "ISortFlag     = 0      ">> $inpfile
echo "IFluxFlag     = 1">> $inpfile
echo "Width0        = $width">> $inpfile
echo "Width1        = $width">> $inpfile
echo "dWidth        = 2      ">> $inpfile
echo "DiagDis0      = $flux  ">> $inpfile
echo "DiagDis1      = $flux ">> $inpfile
echo "dDiagDis      = 0.1">> $inpfile
echo "Energy0       = -10.0">> $inpfile
echo "Energy1       = 10.0">> $inpfile
echo "dEnergy       = 0.02">> $inpfile
echo "Kappa         = 1.0">> $inpfile
echo "Epsilon       = 1.0E-2">> $inpfile

cat $inpfile
ls -al \${tmpdir}

# sed 's/NAV/\\$iseed/g' ../$inpfile.NAV >>$inpfile
# cat $inpfile

${binarydir}/${binary} >& ${jobname}.log

# copy the result of the run in the destination for safekeeping

[ -d ${targetdir} ] || mkdir ${targetdir}

mv *.raw ${targetdir}
mv *.inp ${targetdir}

cd ${targetdir}
pwd
#gzip -fv9 *.dat 

cd ..
#rsync -vr LM-L$width-ND$flux-Ed\$liebdef $USER@godzilla.csc.warwick.ac.uk:/storage/disqs/TwoBandModel/DATA/

wait
#exit 0
EOD

chmod 755 ${jobdir}/${jobfile}
#(cd ${jobdir} ; msub -q devel ./${jobfile}) # for development/testing only
(cd ${jobdir} ; source ./${jobfile} ) & # for PRODUCTION

sleep 1

done
done
#done

