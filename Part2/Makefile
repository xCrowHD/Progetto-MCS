# Compilatore
CXX = g++
CXXSTD = c++17

# Eseguibile e sorgenti
TARGET = programma.o
SRC = main.cpp DCT.hpp

# Percorsi librerie
EIGEN_DIR = ../eigen-3.4.0
MATPLOTLIBCPP_DIR = ../

# Flag di compilazione comuni
CXXFLAGS = -std=$(CXXSTD) -Wall -O3 -I$(EIGEN_DIR) -I$(MATPLOTLIBCPP_DIR)

# ===== Python config (per matplotlib-cpp) =====
PYTHON_INCLUDE := $(shell python3 -c "from sysconfig import get_paths as gp; print(gp()['include'])")
PYTHON_LIB := $(shell python3-config --ldflags)

CXXFLAGS += -I$(PYTHON_INCLUDE)
LDFLAGS += $(PYTHON_LIB)

# ===== OpenCV config =====
OPENCV_CFLAGS := $(shell pkg-config --cflags opencv4 2>/dev/null)
OPENCV_LIBS := $(shell pkg-config --libs opencv4 2>/dev/null)

# Phony targets
.PHONY: all build run clean check_all check_eigen check_opencv check_python check_matplotlib_cpp

# Target di default
all: check_all build

# Controllo completo
check_all: check_eigen check_opencv check_python check_matplotlib_cpp

# Controllo Eigen (header-only)
check_eigen:
	@if [ ! -d "$(EIGEN_DIR)" ]; then \
		echo "❌ Directory Eigen NON trovata in $(EIGEN_DIR)."; \
		echo "   Scaricala da: https://gitlab.com/libeigen/eigen/-/releases"; \
		exit 1; \
	else \
		echo "✅ Eigen trovato"; \
	fi

# Controllo OpenCV
check_opencv:
	@if [ -z "$(OPENCV_CFLAGS)" ]; then \
		echo "❌ OpenCV NON trovato tramite pkg-config."; \
		echo "   Installa con: sudo apt install libopencv-dev"; \
		exit 1; \
	else \
		echo "✅ OpenCV trovato"; \
	fi

# Controllo Python + matplotlib
check_python:
	@command -v python3 >/dev/null 2>&1 || (echo "❌ Python3 non trovato. Installa con: sudo apt install python3-dev" && exit 1)
	@python3 -c "import matplotlib" 2>/dev/null || (echo "❌ matplotlib non trovato. Installa con: pip3 install matplotlib" && exit 1)
	@echo "✅ Python3 e matplotlib trovati"

# Controllo matplotlib-cpp header
check_matplotlib_cpp:
	@if [ ! -f "$(MATPLOTLIBCPP_DIR)/matplotlibcpp.h" ]; then \
		echo "❌ matplotlibcpp.h non trovato in $(MATPLOTLIBCPP_DIR)"; \
		exit 1; \
	else \
		echo "✅ matplotlib-cpp header trovato"; \
	fi

# Compilazione
build: $(SRC)
	$(CXX) $(CXXFLAGS) $(OPENCV_CFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS) $(OPENCV_LIBS)
	@echo "✅ Compilazione completata: $(TARGET)"

# Esecuzione
run: build
	./$(TARGET)

# Pulizia
clean:
	rm -f $(TARGET)
	@echo "🧹 Pulizia completata"
