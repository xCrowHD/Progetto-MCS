# Nome del compilatore
CXX = g++

# File sorgenti
SRC = main.cpp DCT.hpp

# Nome dell'eseguibile
TARGET = programma.o

# Directory di Eigen
EIGEN_DIR = ../eigen-3.4.0

# Opzioni di compilazione
CXXFLAGS = -I$(EIGEN_DIR) -std=c++17 -Wall -O3

# Regole Phony
.PHONY: build clean run all

# Regola predefinita
all: build

# Compilazione
build: $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)
	@echo "✅ Compilazione completata: $(TARGET)"

# Esecuzione
run: build
	./$(TARGET)

# Pulizia
clean:
	rm -f $(TARGET)
	@echo "🧹 Pulizia completata"
