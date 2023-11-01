#! /bin/bash

## Generate random points uniformly distributed in a disk:
## pts_disk.sh PREC PTS RADIUS

## One generates two independent uniform variables A (on 0-R^2) and B (on 0-2*Pi)
## The points sqrt(A)*e^{i*B} are uniformly distributed on the disk of radius R

## No sanity checks: use responsively and at your own risks

PREC=$1
PTS=$2
RADIUS=$3
FastPolyEval="UNDEFINED"
#FastPolyEval="$HOME/FPE_Demo/bin/FastPolyEval"

if [[ "$FastPolyEval" =~ "UNDEFINED" ]]; then
    echo "Please update the FastPolyEval variable in this script."
    echo "The variable should contain the path of the executable file."
    echo "For example: <local_gitclone_path>/bin/FastPolyEval"
    exit 1
fi

function square(){ ## bash workaround to compute square with floats
    echo $1 | awk '{print $1*$1}'
}

TMP1=tmp1.csv
TMP2=tmp2.csv
TMP3=tmp3.csv
FILE="disk_${PREC}_${PTS}.csv"

echo "Generating $PTS random evaluation points on the unit disk (precision $PREC bits)"
echo ""


$FastPolyEval -rand $PREC $PTS $TMP1 0 $(square $RADIUS)
## square root is not implemented in FastPolyEval; using shell magic (precision may be lost)
cat $TMP1 | cut -d ',' -f1 | awk '{print sqrt($1)}' | sed -e 's/$/, 0/' > $TMP3
$FastPolyEval -rand $PREC $PTS $TMP2 0 6.283185307179586
$FastPolyEval -join $PREC $TMP3 $TMP2 $TMP1
$FastPolyEval -rot $PREC $TMP1 $FILE

echo ""
echo "Cleaning up $TMP1, $TMP2, $TMP3"
rm -f $TMP1 $TMP2 $TMP3
