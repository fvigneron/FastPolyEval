#! /bin/bash

## Benchmark polynomials with normal coefficients by increasing degree

## Update this variable to match your system
FastPolyEval="UNDEFINED"
#FastPolyEval="$HOME/FPE_Demo/bin/FastPolyEval"

if [[ "$FastPolyEval" =~ "UNDEFINED" ]]; then
    echo "Please update the FastPolyEval variable in this script."
    echo "The variable should contain the path of the executable file."
    echo "For example: <local_gitclone_path>/bin/FastPolyEval"
    exit 1
fi

echo -e "\nEasy benchmark of FastPolyEval\n"
echo "The coefficients of normal polynomials are randomly distributed with law N(0,1)."
echo "For this family the gain of PFE over Horner is about halfway between the best and worst case."
echo "In this benchmark, the degree will be doubled at each increment."

echo -e "\n--------------------------------------------------------\n"

## Precision of the computation in this benchmark run

PREC=54

## Generate evaluation points on the sphere

NBPTS=10000

TMP="tmp.csv"
SPHERE="sph_${PREC}_${NBPTS}.csv"
SEEDS=$(echo $NBPTS | awk '{print 1+int(sqrt(314*$1/100))}');
## The seed is an integer such that NBPTS = SEEDS*SEEDS*100/314
echo "Generating random evaluation points on the sphere (precision $PREC bits, $SEEDS seeds)"
echo "Expect about $NBPTS evaluation points."
echo ""

$FastPolyEval -sphere $PREC $SEEDS $TMP
$FastPolyEval -polar $PREC $TMP $SPHERE
ls -lh "$SPHERE"
echo -e "\nRemoving $TMP"
rm -f $TMP
echo ""

## Actual benchmark

N=64  ## benchmark starts with 2*N
OK=true
while $OK; do

    N=$((2*N))
    POLY="hyp_${N}.csv"
    echo "********* Benchmark of normal poly of deg $N with a requested precision of $PREC bits *********"
    echo ""
    read -n1 -p "Proceed ? (y/n) : "
    echo -e "\n"
    if [[ "$REPLY" =~ "y" || "$REPLY" = "" ]]; then
    
        date
        ## Generate polynomial
        TMP1="tmp1.csv"
        TMP2="tmp2.csv"
        $FastPolyEval -normal $PREC $N "$TMP1"
        $FastPolyEval -normal $PREC $N "$TMP2"
        $FastPolyEval -join $PREC "$TMP1" "$TMP2" "$POLY"
        rm -f "$TMP1" "$TMP2"
        ls -lh "$POLY"
        echo ""
        
        date
        ## Benchmark evaluation
        ## Modify this command if you want to keep the output files
        $FastPolyEval -eval $PREC "$POLY" "$SPHERE" /dev/null /dev/null 1 /dev/null 1

        ## If you want to only run the FPE algorith and not compare to Horner,
        ## comment the previous instruction and uncomment the following one instead.
        #$FastPolyEval -eval $PREC "$POLY" "$SPHERE" /dev/null

        date
        ## Cleanup
        rm -f "$POLY"
        echo ""
        
    else
        OK=false
        echo -e "Thank you for your interest in FastPolyEval. Please, come back soon!\n"
    fi
done
