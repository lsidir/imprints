CPPFLAGS = -O3 -std=c++11 -Wall
CXX	= g++
LDFLAGS	= -lpthread -lOSAL -lAAS

all: HARPv2


HARPv2: main.cpp iFPGA.cpp RuntimeClient.cpp
	$(CXX) -D HARPv2 -I$(AALSDK)/include $(CPPFLAGS) -o t main.cpp iFPGA.cpp RuntimeClient.cpp -L$(AALSDK)/lib $(LDFLAGS) -laalrt