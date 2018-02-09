#include <iostream>
#include "iFPGA.h"

using namespace std;

int main(int argc, char* argv[]) {

	uint32_t page_size_in_cache_lines = 65536; // 65536 x 64B = 4 MB
	uint32_t pages_to_allocate = 8;

	RuntimeClient runtimeClient;
	iFPGA* interfaceFPGA = new iFPGA(&runtimeClient, pages_to_allocate, page_size_in_cache_lines);
	if(!runtimeClient.isOK()) {
		cout << "FPGA runtime failed to start" << endl;
		exit(1);
	}

	// Example communation with the FPGA

	uint32_t num_cache_lines_to_send = 100;
	// 1. Prepare memory to send
	for (uint32_t i = 0; i < num_cache_lines_to_send*8; i++) {
		interfaceFPGA->writeToMemory64('i', i+1, i);
	}
	// 2. Configure registers on the FPGA
	CSR_WRITE32(interfaceFPGA, CSR_READ_OFFSET, 0);
	CSR_WRITE32(interfaceFPGA, CSR_WRITE_OFFSET, 0);
	CSR_WRITE32(interfaceFPGA, CSR_NUM_LINES, num_cache_lines_to_send);
	// 3. Trigger the FPGA
	interfaceFPGA->doTransaction();
	
	delete interfaceFPGA;

	return 0;
}