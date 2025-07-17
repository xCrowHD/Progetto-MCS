# Nome del compilatore
CXX = g++

# File sorgenti
SRC = main.cpp linear_resolver.hpp lr_utils.hpp lr_test.hpp

# Nome dell'eseguibile
TARGET = programma.o

CXXFLAGS  = -std=c++17 -Wall -O3             \
            -I.                              \
            -I../eigen-3.4.0                 

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
	@./$(TARGET)

# Pulizia
clean:
	rm -f $(TARGET)
	@echo "🧹 Pulizia completata"
