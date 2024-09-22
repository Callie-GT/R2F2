AUTOPILOT_ROOT :=/tools/software/xilinx/Vitis_HLS/2023.1

ASSEMBLE_SRC_ROOT := .
IFLAG += -I "${AUTOPILOT_ROOT}/include"
# IFLAG += -I "${ASSEMBLE_SRC_ROOT}"
# IFLAG += -I "/usr/include/x86_64-linux-gnu"
IFLAG += -D__SIM_FPO__ -D__SIM_OPENCV__ -D__SIM_FFT__ -D__SIM_FIR__ -D__SIM_DDS__ -D__DSP48E1__ -DHLS_NO_XIL_FPO_LIB


# IFLAG += -DDEBUG_FILE_PRINT=1
IFLAG += -g 
CFLAG += -fPIC -O0 #-fsanitize=address
CFLAG += -lm
CFLAG += -std=c++11 -Wno-unused-result 


all:
	g++ flexible_floating_point.cpp templated_floating_point.cpp helper.cpp tb.cpp -o result $(CFLAG) $(IFLAG)
	# g++ heat_equation.cpp -o result $(CFLAG) $(IFLAG)
	# g++ shallow_water_equation.cpp -o result $(CFLAG) $(IFLAG)
	
clean:
	rm -f *.o result
