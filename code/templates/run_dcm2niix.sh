#!/usr/bin/bash

# load required modules...
module use /projects/ics/modules
module load fsl/6.0.3
module load dcmtk 
export DCMDICTPATH=/projects/ics/software/dcmtk/share/dcmtk/dicom.dic
module load freesurfer/7.1.0
module load dcm2niix/v1.0.20190410


# Organize variables
inputdicomdir=DICOM_PLACEHOLDER
bidsfile=BIDSFILE_PLACEHOLDER
tempfile=NIFTI_PLACEHOLDER

tmpdir=`dirname $tempfile`
mkdir -p $tmpdir
log=$tempfile.log
cmd="dcm2niix -z y -b y -f %f y -o $tmpdir -v y $inputdicomdir"
echo $cmd >> $log
$cmd >> $log 2>&1


# after converting to nifti move to bids format
outdir=`dirname $bidsfile`
mkdir -p $outdir

mv ${tempfile}.json $bidsfile.json
mv ${tempfile}.nii.gz $bidsfile.nii.gz

ln -s $bidsfile.nii.gz ${tempfile}.nii.gz
ln -s $bidsfile.json ${tempfile}.json

if [ -f ${tempfile}.bvec ]; then mv ${tempfile}.bvec ${bidsfile}.bvec ; fi
if [ -f ${tempfile}.bval ]; then mv ${tempfile}.bval ${bidsfile}.bval ; fi

