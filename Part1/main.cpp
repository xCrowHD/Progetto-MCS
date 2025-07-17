#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include "linear_resolver.hpp"
#include "tabulate.hpp"
#include "lr_test.hpp"

using namespace tabulate;
using Row_t = Table::Row_t;

// Colori ANSI (puoi rimuoverli se non vuoi)
#define COLOR_RED "\033[1;31m"
#define COLOR_GREEN "\033[1;32m"
#define COLOR_YELLOW "\033[1;33m"
#define COLOR_CYAN "\033[1;36m"
#define COLOR_RESET "\033[0m"

void show_menu()
{
    std::cout << COLOR_CYAN << "\n===== MENU PRINCIPALE =====\n"
              << COLOR_RESET;
    std::cout << COLOR_GREEN << "1. Test su matrici .mtx \n";
    std::cout << "2. Test su matrice a scelta\n";
    std::cout << "q. Esci\n"
              << COLOR_RESET;
    std::cout << COLOR_YELLOW << "Scegli un'opzione: " << COLOR_RESET;
}

int main()
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
            lr_test<Eigen::SparseMatrix<double>> lr_t;
            lr_t.test_lr();
            break;
        case '2':
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
