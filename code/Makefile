# =============================================
#
# Makefile for FPE: A Fast Polynomial Evaluator
#
# =============================================

# C compiler
CC=gcc
# Automatically find all source files *.c (assuming Makefile is at the root of the 'code' directory).
# The content of TMP_DIR should be ignored.
SOURCE_DIR=.
SOURCE_FILES=tools/*.c apps/*.c eval/*.c numbers/*.c poly/*.c
DIR_INCLUDE=-I ./tools -I ./apps -I ./tmp -I ./eval -I ./numbers -I ./poly

# A disposable temporary directory, used for the build (created/deleted on demand).
TMP_DIR=./hw_tmp

# Default target for binaries (will create a new a directory alongside 'code').
TARGET_DIR=../bin
TARGET_NAME=FastPolyEval
TARGET=$(TARGET_DIR)/$(TARGET_NAME)
TARGET_FP32=$(TARGET)_FP32
TARGET_FP64=$(TARGET)_FP64
TARGET_FP80=$(TARGET)_FP80

# FastPolyEval uses machine numbers for low precision computations.
# The state of the #define variables MACHINE_EXTRA_PREC and MACHINE_LOW_PREC in numbers/ntypes.h
# defines which machine numbers FastPolyEval defaults to (either FP32, FP64 or FP80).
# The 'hardware' target in this Makefile generates all 3 possible binaries.
# The relevant code is duplicated in TMP_DIR and tweaked automatically with sed (toggle comment on the lines).
# The temporary directory is then included in priority by the compiler, which overrides the default code.
TMP_FILE=$(TMP_DIR)/ntypes.h

# Compiler's options
# GCC_INCLUDE contains the list of all sub-directories in SOURCE_DIR, each prefixed with -I.
# The -Wall flag enforces good programing practice by tracking ambiugous/questionable code.
# The -lm -lgmp -lmpfr are the 3 math libraries used by FPE (install GMP/MPFR first).
# The -O3 optimization flag is best for production code.
GCC_LIBRARY_PATH=#-L./<pkg>/lib
GCC_INCLUDE_PATH=#-I./<pkg>/include
GCC_INCLUDE=${DIR_INCLUDE}
GCC_OPTIONS=-Wall
GCC_LIBRARIES=-lm -lgmp -lmpfr
GCC_OPTIMIZATION=-O3

.PHONY: greetings fpe hardware hw32 hw64 hw80 hw_clean clean

greetings:
	@echo "\n"
	@echo '             * ********************************** *'
	@echo '             *  FPE: A Fast Polynomial Evaluator  *'
	@echo "             * ********************************** *\n"
	@echo 'Installation instructions:'
	@echo '   1. Install GMP (https://gmplib.org)'
	@echo '   2. Install MPFR (https://www.mpfr.org)'
	@echo "   3. Run 'make fpe' to compile from the source code"
	@echo "   4. Use '$(TARGET_NAME)' regularly and enjoy"
	@echo '   5. Recommend it to at least two friends ;-)'
	@echo "\nBinaries will be installed in:\n   $(shell pwd)/$(TARGET_DIR)"
	@echo 'You can move them to a new location at your convenience.'
	@echo 'Alternatively, consider adding this directory to your path.'
	@echo "\nTo start fresh, run 'make clean'\n"
	@echo "In the default build (step 3), FastPolyEval uses FP64 for low-precision computations."
	@echo "To generate binaries with FP32, FP64 or FP80 support, run 'make hardware'"
	@echo "Alternatively, (un)comment the relevant #define in numbers/ntypes.h and return to step 3.\n"
	@echo 'If your libraries are not installed in the default location, or if you want to adjust'
	@echo 'the compile options, or for good practice, check & tweak this Makefile before step 3.'
	@echo "\nConsider reading the documentation. Congratulations if you already did !\n"


fpe: greetings
	@echo "----------------------------------------------------------------------------\n"
	@echo "\nCompiling FastPolyEval from the current source code"
	@mkdir -p $(TARGET_DIR)
	@echo "The following binary will be generated: $(TARGET)\n"
	${CC} $(SOURCE_FILES) -o $(TARGET) ${GCC_INCLUDE_PATH} $(GCC_INCLUDE) $(GCC_OPTIONS) ${GCC_LIBRARY_PATH} $(GCC_LIBRARIES) $(GCC_OPTIMIZATION)

hardware: greetings hw32 hw64 hw80

hw32:
	@echo "----------------------------------------------------------------------------\n"
	@echo "\nCompiling a version of $(TARGET_NAME) with FP32 machine numbers"
	@mkdir -p $(TARGET_DIR)
	@echo "The following binary will be generated: $(TARGET_FP32)\n"
	@mkdir -p $(TMP_DIR)
	@sed -e 's|^.*#define MACHINE_EXTRA_PREC|//#define MACHINE_EXTRA_PREC|; s|^.*#define MACHINE_LOW_PREC|#define MACHINE_LOW_PREC|' $(SOURCE_DIR)/numbers/ntypes.h > $(TMP_FILE) 2>&1
	${CC} $(SOURCE_FILES) -o $(TARGET_FP32) -I $(TMP_DIR) ${GCC_INCLUDE_PATH} $(GCC_INCLUDE) $(GCC_OPTIONS) ${GCC_LIBRARY_PATH} $(GCC_LIBRARIES) $(GCC_OPTIMIZATION)
	@make hw_clean

hw64:
	@echo "----------------------------------------------------------------------------\n"
	@echo "\nCompiling a version of $(TARGET_NAME) with FP64 machine numbers"
	@mkdir -p $(TARGET_DIR)
	@echo "The following binary will be generated: $(TARGET_FP64)\n"
	@mkdir -p $(TMP_DIR)
	@sed -e 's|^.*#define MACHINE_EXTRA_PREC|//#define MACHINE_EXTRA_PREC|; s|^.*#define MACHINE_LOW_PREC|//#define MACHINE_LOW_PREC|' $(SOURCE_DIR)/numbers/ntypes.h > $(TMP_FILE) 2>&1
	${CC} $(SOURCE_FILES) -o $(TARGET_FP64) -I $(TMP_DIR) ${GCC_INCLUDE_PATH} $(GCC_INCLUDE) $(GCC_OPTIONS) ${GCC_LIBRARY_PATH} $(GCC_LIBRARIES) $(GCC_OPTIMIZATION)
	@make hw_clean

hw80:
	@echo "----------------------------------------------------------------------------\n"
	@echo "\nCompiling a version of $(TARGET_NAME) with FP80 machine numbers"
	@mkdir -p $(TARGET_DIR)
	@echo "The following binary will be generated: $(TARGET_FP80)\n"
	@mkdir -p $(TMP_DIR)
	@sed -e 's|^.*#define MACHINE_EXTRA_PREC|#define MACHINE_EXTRA_PREC|; s|^.*#define MACHINE_LOW_PREC|//#define MACHINE_LOW_PREC|' $(SOURCE_DIR)/numbers/ntypes.h > $(TMP_FILE) 2>&1
	${CC} $(SOURCE_FILES) -o $(TARGET_FP80) -I $(TMP_DIR) ${GCC_INCLUDE_PATH} $(GCC_INCLUDE) $(GCC_OPTIONS) ${GCC_LIBRARY_PATH} $(GCC_LIBRARIES) $(GCC_OPTIMIZATION)
	@make hw_clean

hw_clean:
	@rm -rf $(TMP_DIR)

clean: 
	@echo 'Removing executables from $(TARGET_DIR) subdirectory'
	rm -f $(TARGET)
	rm -f $(TARGET_FP32)
	rm -f $(TARGET_FP64)
	rm -f $(TARGET_FP80)
	@make hw_clean

