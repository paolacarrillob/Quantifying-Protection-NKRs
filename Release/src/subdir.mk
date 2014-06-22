################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Bitstring.cpp \
../src/Coevolution.cpp \
../src/Genes.cpp \
../src/Host.cpp \
../src/MathFunctions.cpp \
../src/World.cpp 

OBJS += \
./src/Bitstring.o \
./src/Coevolution.o \
./src/Genes.o \
./src/Host.o \
./src/MathFunctions.o \
./src/World.o 

CPP_DEPS += \
./src/Bitstring.d \
./src/Coevolution.d \
./src/Genes.d \
./src/Host.d \
./src/MathFunctions.d \
./src/World.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


