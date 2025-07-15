# -------------------
#  Impostazioni
# -------------------
CXX       := g++
CXXSTD    := c++17
TARGET    := programma

# Sorgenti progetto
SRC := main.cpp DCT.hpp DCT_FDCT.hpp DCT_basic_test.hpp

# Eigen
EIGEN_DIR := ../eigen-3.4.0

# Tabulate (single header in ../)
TABULATE_HEADER := ../tabulate.hpp

# OpenCV
OPENCV_CXXFLAGS := $(shell pkg-config --cflags opencv4)
OPENCV_LIBS     := $(shell pkg-config --libs opencv4)

# Flags compilazione
CXXFLAGS := -std=$(CXXSTD) -Wall -O3 \
            -I$(EIGEN_DIR) \
            -I.. \
            $(OPENCV_CXXFLAGS)

LIBS := $(OPENCV_LIBS) -lpthread -lz

# -------------------
#  Phony targets
# -------------------
.PHONY: all check build run clean

all: check build

check:
	@echo "🔍 Controllo dipendenze..."
	@test -d "$(EIGEN_DIR)" || (echo "❌ Eigen non trovato in $(EIGEN_DIR)"; exit 1)
	@test -f "$(TABULATE_HEADER)" || (echo "❌ tabulate.hpp non trovato in ../"; exit 1)
	@pkg-config --exists opencv4 || (echo "❌ OpenCV non trovato"; exit 1)
	@echo "✅ Tutto ok."

# -------------------
#  Compilazione finale
# -------------------
build:
#@echo "🛠 Compilo con:"
#@echo "   CXXFLAGS   = $(CXXFLAGS)"
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) $(LIBS)

run: build
	./$(TARGET)

clean:
	-rm -f $(TARGET)
	@echo "🧹 Pulizia completata."
