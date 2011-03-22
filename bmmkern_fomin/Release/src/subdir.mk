################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/__moldefs.cpp \
../src/bmmkern.cpp 

OBJS += \
./src/__moldefs.o \
./src/bmmkern.o 

CPP_DEPS += \
./src/__moldefs.d \
./src/bmmkern.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/home/fomin/projects/bmmkern/include -O3 -fno-strict-aliasing -finline-limit=15000 -Wall -c -fmessage-length=0 -Wno-deprecated -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


