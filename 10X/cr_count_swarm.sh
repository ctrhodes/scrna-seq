#! /bin/bash
# this file is cr_count_swarm.sh
# CMD: sbatch cr_count_swarm.sh

set -o pipefail
set -e

function fail {
    echo "$@" >&2
    exit 1
}

#rm /home/rhodesct/bin/cr_counts.swarm

# Main directories
export BASE_DIR="/home/rhodesct/Data"
export REF_SUB="indices/cellranger/GRCm38_ai9"
export TRANSCRIPTOME=$BASE_DIR/$REF_SUB

# Project specific directories
export TYPE="Projects/scRNA-Seq-10x"
export PROJECT="Petros_May2018_SingleCell"
export FLOW_CELL="180508"
export FQ_SUB="outs/fastq_path/PETROS"
export FASTQ=$BASE_DIR/$TYPE/$PROJECT/$FLOW_CELL/$FQ_SUB

export RUNS=$(echo {A..C})

echo "base dir: "$BASE_DIR
echo "reference subdir: "$REF_SUB
echo "transcriptome path: "$TRANSCRIPTOME
echo "type: "$TYPE
echo "project: "$PROJECT
echo "flow cell: "$FLOW_CELL
echo "fastq subdir: "$FQ_SUB
echo "fastq path: "$FASTQ
echo "runs: "$RUNS


SCRIPT_DIR="/home/rhodesct/bin"
SCRIPT_NAME=`basename "$0"`

SWARM_NAME=$(echo ${SCRIPT_NAME} | sed 's/.sh//').swarm

SWARM_FILE=$SCRIPT_DIR/$SWARM_NAME
echo "swarmfile: "$SWARM_FILE

if [ -e $SWARM_FILE ]
then
    echo "rebuilding swarm file: "$SWARM_FILE
    rm $SWARM_FILE
else
    echo "creating swarmfile..."
fi

for run in $RUNS
do

echo "cellranger count --id=$LETTER \
--transcriptome=${TRANSCRIPTOME} \
--fastqs=${FASTQ}/$run \
--chemistry=SC3Pv2 \
--localcores="'$SLURM_CPUS_PER_TASK'" \
--localmem=70" >> $SWARM_FILE
done

chmod 755 $SWARM_FILE

ACTIVE=$BASE_DIR/$TYPE/$PROJECT
echo "going to: "$ACTIVE
cd $ACTIVE || fail "no such directory"

swarm -f $SWARM_FILE -g 73 -t 16 --verbose 1 --time=48:00:00 --module cellranger

