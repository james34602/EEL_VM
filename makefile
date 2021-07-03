#*******************************************************************************
# Name of application
#*******************************************************************************

TARGET = eel_CLI

#*******************************************************************************

.PHONY : all

TOOLPATH=/Applications/Xcode.app/Contents/Developer/usr/bin
PREFIX=
TOOLCHAIN=$(TOOLPATH)/$(PREFIX)
SDK=/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.15.sdk

ENV_MACRO = 

ARCH = -arch x86_64 -isysroot $(SDK)

F_OPTS = -F$(SDK)System/Library/Frameworks -F$(SDK)System/Library/PrivateFrameworks
W_OPTS = -Wimplicit -Wno-unused-value -Wno-incompatible-pointer-types \
    -Wno-implicit-function-declaration -Wno-multichar -Wno-logical-op-parentheses \
    -Wno-shift-count-overflow -Wno-missing-declarations -Wno-switch \
    -Wno-int-to-void-pointer-cast \
 
ARCH_DIR = $(platform)

#-------------------------------------------------------------------------------

CPP       = $(TOOLCHAIN)g++
CC        = $(TOOLCHAIN)gcc
LD        = $(TOOLCHAIN)ld
AR        = ar
STRIP     = strip
OBJCOPY   = $(TOOLCHAIN)objcopy
OBJDUMP   = $(TOOLCHAIN)objdump

#*******************************************************************************
# config
#*******************************************************************************

ROOT   = .
INCDIR = $(ROOT)/include

#*******************************************************************************
# Compiler and linker options
#*******************************************************************************

INCLUDE   = $(addprefix -I,$(INCDIR)) -I$(TOOLCHAIN)/include

CC_OPTS_NO_OPT = -g $(INCLUDE) $(F_OPTS) $(W_OPTS) $(ARCH) -pipe $(ENV_MACRO)

CC_OPTS   = -Ofast -ftree-vectorize $(CC_OPTS_NO_OPT)
CPP_OPTS  = $(CC_OPTS)
CC_OPTS_A = $(CC_OPTS) -D_ASSEMBLER_

#*******************************************************************************
# Directories & Files to be compiled & Rules
#*******************************************************************************

SRCS = $(ROOT)
OBJS = \
	numericSys/FilterDesign/generalFdesign.o \
	numericSys/FilterDesign/eqnerror.o \
	numericSys/FilterDesign/firls.o \
	numericSys/FilterDesign/cos_fib_paraunitary.o \
	numericSys/SolveLinearSystem/inv.o \
	numericSys/SolveLinearSystem/pinv.o \
	numericSys/SolveLinearSystem/mldivide.o \
	numericSys/SolveLinearSystem/mrdivide.o \
	numericSys/solvopt.o \
	numericSys/cpoly.o \
	numericSys/MersenneTwister.o \
	numericSys/quadprog.o \
	numericSys/codelet.o \
	numericSys/FFTConvolver.o \
	libsamplerate/src_sinc.o \
	libsamplerate/samplerate.o \
	cpthread.o \
	s_str.o \
	fft.o \
	nseel-compiler.o \
	nseel-ram.o \
	y.tab.o \
	loose_eel.o

#*******************************************************************************
# target
#*******************************************************************************

$(TARGET) : $(OBJS)
	$(CC) $(CC_OPTS) -o $(TARGET) -framework Accelerate $^

#*******************************************************************************

$(OBJS) : %.o : $(SRCS)/%.c
	$(CC) $(CC_OPTS) -c -o $@ $<

#*******************************************************************************

.PHONY : clean

clean :
	rm $(TARGET) $(OBJS) *.d

#*******************************************************************************
#                          E N D  O F  F I L E
#*******************************************************************************
