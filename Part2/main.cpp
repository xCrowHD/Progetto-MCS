#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <random>

#include "DCT_basic_test.hpp"
#include "DCT_FDCT.hpp"
#include "BMP_FDCT.hpp"
#include "tinyfiledialogs.h"

// Colori ANSI (puoi rimuoverli se non vuoi)
#define COLOR_RED     "\033[1;31m"
#define COLOR_GREEN   "\033[1;32m"
#define COLOR_YELLOW  "\033[1;33m"
#define COLOR_CYAN    "\033[1;36m"
#define COLOR_RESET   "\033[0m"

void Basic_DCT_testing()
{
    DCTBasicTest test;
    test.test_vector(test.vector);
    test.test_matrix(test.matrix);
}

void DCT_FDCT_plots_testing()
{
    DCTFDCT dct_fdct;
    dct_fdct.test_matrix();
}

void Bmp_FDCT_testing()
{
    BMPFDCT fdct;
    fdct.test_it();
}

void show_menu()
{
    std::cout << COLOR_CYAN << "\n===== MENU PRINCIPALE =====\n" << COLOR_RESET;
    std::cout << COLOR_GREEN << "1. Test su vettori e matrici con DCT base\n";
    std::cout << "2. Plot DCT/FDCT su CSV\n";
    std::cout << "3. DCT su immagine BMP\n";
    std::cout << "q. Esci\n" << COLOR_RESET;
    std::cout << COLOR_YELLOW << "Scegli un'opzione: " << COLOR_RESET;
}

int main(int argc, char const *argv[])
{

    // Modalità interattiva
    char choice;
    do
    {
        show_menu();
        std::cin >> choice;

        switch (choice)
        {
        case '1':
            Basic_DCT_testing();
            break;
        case '2':
            DCT_FDCT_plots_testing();
            break;
        case '3':
            Bmp_FDCT_testing();
            break;
        case 'q':
            std::cout << COLOR_GREEN << "Chiusura in corso...\n"
                      << COLOR_RESET;
            break;
        default:
            std::cout << COLOR_RED << "Scelta non valida. Riprova.\n"
                      << COLOR_RESET;
        }

    } while (choice != 'q');

    return 0;
}
