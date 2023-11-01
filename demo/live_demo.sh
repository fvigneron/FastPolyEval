#! /bin/bash
#
## ----------------------------------------------------------------
## Live demo of the FastPolyEval app on linux/mac system
## ----------------------------------------------------------------

RED='\033[0;31m'
YLW='\033[1;33m'
CYN='\033[1;36m'
GRN='\033[1;32m'
NC='\033[0m' # no color

CLR=$CYN  ## choose best prompt color for your screen settings

function shell(){ ## display a fancy command prompt
   echo -e "${CLR}> $* ${NC}"
}

function prompt(){ ## display command then execute upon answer
   shell "$*"
   read -n1 -p ""
   $*
}

function sprompt(){ ## display command then execute upon answer
   shell "$*"
   read -n1 -p "Use 's' skip this step or any other key to proceed: "
   if [ "x$REPLY" != "xs" ]
   then
   	$* 2> /dev/null ## warning: error messages are hidden
   else
	echo
   fi
}

function comment(){ ## display a comment
	echo -e "${GRN}# $* ${NC}"
}

## ----------------------------------------------------------------

TGT="$HOME/FPE_Demo"
DEMODIR="$TGT/demo_tmp"

## Clone the repository into TGT
echo "FPE live demo will take place in $TGT."
echo -e "At each prompt ${CLR}> command${NC} press any key to execute the shell command\n"
read -n1 -p "Overwrite and download a fresh clone (n/y) ? "
if [ "x$REPLY" == "xy" ] ## prompt is not used here to streamline the flow of the demo
then
	rm -rf $TGT
	shell "\n\ngit clone http://github.com/fvigneron/FastPolyEval $TGT"
	git clone http://github.com/fvigneron/FastPolyEval "$TGT"
else
	echo -e "\n\nAutonomous mode: assuming that you have already run"
	shell "git clone http://github.com/fvigneron/FastPolyEval $TGT"
fi
echo

## Sanity check
prompt cd $TGT
if [ $(basename `pwd`) != $(basename "$TGT") ]
then
	echo "Something appears to be wrong with $TGT"
	exit 1
fi

## Explore the repository
comment 'Explore the repository'
which tree > /dev/null
if [ $? == 0 ]
then
    prompt tree -C
else
    prompt ls
fi
echo

## About the license
prompt cat LICENSE.txt
echo
prompt cat THIRDPARTY.md
echo

## Compile from the sources
comment 'Compile from scratch'
prompt cd code
shell make
shell make fpe
## sed -i -e 's/^GCC_OPTIONS=-Wall$/GCC_OPTIONS=-Wall -Wl,-no_pie/' Makefile ## MacOS harmless anoyance on some installs
sprompt make hardware 
cp "$TGT/bin/FastPolyEval_FP64" "$TGT/bin/FastPolyEval" ## streamlined to emulate running both targets fpe+hardware
echo
prompt ls -lh "$TGT/bin"
echo
shell PATH="\$PATH:$TGT/bin" ## add to the PATH
PATH="$PATH:$TGT/bin"
echo

## First contact
shell mkdir "$DEMODIR"
shell cd "$DEMODIR"
mkdir -p "$DEMODIR" && cd "$DEMODIR" ## streamline shell commands
echo
comment 'First contact'
prompt FastPolyEval
echo

prompt FastPolyEval -eval -help
echo

shell FastPolyEval_FP32
prompt FastPolyEval_FP32 | tail -2
echo

shell FastPolyEval_FP80
prompt FastPolyEval_FP80 | tail -2
echo

## Generate polynomials
comment 'Generate polynomials'
prompt FastPolyEval -Laguerre 53 15 laguerre_FP64_deg15.csv
echo
comment 'File format is CSV (list of complex numbers)'
prompt cat laguerre_FP64_deg15.csv  ## file format is CSV (list of complex numbers)
echo
prompt FastPolyEval -Legendre 53 15 legendre_FP64_deg15.csv
echo
prompt cat legendre_FP64_deg15.csv
echo
comment 'File manipulation tools'
prompt FastPolyEval -join 53 laguerre_FP64_deg15.csv legendre_FP64_deg15.csv poly.csv  ## extensive file manipulation tools
echo
prompt cat poly.csv
echo
comment 'File safe operations'
prompt FastPolyEval -prod 128 poly.csv poly.csv poly.csv  ## operations are file safe + precision change on the fly
echo
comment 'Easy precision changes (be self-conscious: no magic gain!)'
prompt cat poly.csv
echo

## Generate evaluation points
SEED=178
comment 'Generate evaluation points'
shell FastPolyEval -sphere 128 $SEED tmp_file.csv
shell FastPolyEval -polar 128 tmp_file.csv sphere.csv
read -n1 -p ""
FastPolyEval -sphere 128 $SEED tmp_file.csv
FastPolyEval -polar 128 tmp_file.csv sphere.csv
rm -f tmp_file.csv
echo

## Evaluation benchmarks
comment 'Evaluation'
prompt FastPolyEval -eval 128 poly.csv sphere.csv val_poly_sphere.csv
comment 'Benchmark call'
prompt FastPolyEval -eval 128 poly.csv sphere.csv val_poly_sphere.csv err.csv 1 horner.csv 1

## poly.csv is deg 30 ; raise the degree
comment 'Raise the degree of poly.csv to 240'
shell FastPolyEval -prod 128 poly.csv poly.csv poly.csv  ## deg 60
FastPolyEval -prod 128 poly.csv poly.csv poly.csv > /dev/null
shell FastPolyEval -prod 128 poly.csv poly.csv poly.csv  ## deg 120
FastPolyEval -prod 128 poly.csv poly.csv poly.csv > /dev/null
prompt FastPolyEval -prod 128 poly.csv poly.csv poly.csv  ## deg 240
echo

comment 'Balance with Horner'
prompt FastPolyEval -eval 128 poly.csv sphere.csv val_poly_sphere.csv err.csv 3 /dev/null 3 ## should barely beat Horner
comment 'Raise degree to 480'
prompt FastPolyEval -prod 128 poly.csv poly.csv poly.csv  ## deg 480
comment 'Beat Horner'
prompt FastPolyEval -eval 128 poly.csv sphere.csv val_poly_sphere.csv err.csv 1 /dev/null 1
comment 'Please, wait to the end to REALLY beat Horner :-)'
echo
comment 'FP64 exponents are limited to 10^308 (error message is normal)'
prompt FastPolyEval -eval 53 poly.csv sphere.csv val_poly_sphere.csv ## FP64 exponents are limited to 10^308
comment 'MPFR exponents are virtually unlimited'
prompt FastPolyEval -eval 54 poly.csv sphere.csv val_poly_sphere.csv err.csv 1 /dev/null 1 ## MPFR numbers have unlimited exponents
comment 'Four reasons to use FPE for your large polynomial evaluations:'
echo
comment 'Benefit 1: Generate individual reports for bit confidence and absolute error'
prompt head err.csv  ## we can generate individual reports for absolute error and confidence
tail err.csv
echo
comment 'Benefit 2: Let_s beat Horner by a factor 20... (FPE instantaneous, Horner 15s)'
prompt FastPolyEval -eval 54 sphere.csv sphere.csv val.csv err_report.csv 1 horner.csv 1  ## 10k evaluations of a poly of degree 10k : beats Horner by factor 20
rm -f err_report.csv horner.csv
comment 'Benefit 3: Polynomial analysis (e.g. for onboard electronics)'
prompt FastPolyEval -analyse 53 legendre_FP64_deg15.csv concave_cover.csv analysis_53bits.csv  ## bonus: polynomial analysis
prompt cat analysis_53bits.csv
echo
comment 'Benefit 4: CSV format + Library'
comment 'Easy to integrate the FastPolyEval app in your production pipeline (system call).'
comment 'For advanced applications, include the FastPolyEval library (standard C code).'
