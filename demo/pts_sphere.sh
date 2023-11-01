#! /bin/bash

## Generate random points on the Riemann sphere:
## pts_sphere.sh PREC SEEDS
## There will be about SEEDS^2/pi points on the sphere.
## 4*pi*R^2 = SEEDS^2/pi if R = SEEDS/(2*pi)

PREC=$1
SEEDS=$2
FastPolyEval="UNDEFINED"

if [[ "$FastPolyEval" =~ "UNDEFINED" ]]; then
    echo "Please update the FastPolyEval variable in this script."
    echo "The variable should contain the path of the executable file."
    echo "For example: <local_gitclone_path>/bin/FastPolyEval"
    exit 1
fi

TMP=tmp.csv
FILE="sph_${PREC}_s${SEEDS}.csv"

echo "Generating random evaluation points on the sphere (precision $PREC bits, $SEEDS seeds)"
echo "Expect $((SEEDS*SEEDS*100/314)) points."
echo ""

$FastPolyEval -sphere $PREC $SEEDS $TMP
$FastPolyEval -polar $PREC $TMP $FILE

echo ""
echo "Cleaning up $TMP"
rm -f $TMP
