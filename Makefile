CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall -pthread -march=native
TARGET = brkga_solver
OBJS = brkga.o scp_cs_data.o decodificador.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

brkga.o: brkga.cpp scp_cs_data.hpp decodificador.hpp
	$(CXX) $(CXXFLAGS) -c brkga.cpp

scp_cs_data.o: scp_cs_data.cpp scp_cs_data.hpp
	$(CXX) $(CXXFLAGS) -c scp_cs_data.cpp

decodificador.o: decodificador.cpp decodificador.hpp scp_cs_data.hpp
	$(CXX) $(CXXFLAGS) -c decodificador.cpp

clean:
	rm -f $(OBJS) $(TARGET) 2>NUL || del $(OBJS) $(TARGET) 2>NUL

run: $(TARGET)
	./$(TARGET)