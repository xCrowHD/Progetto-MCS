#include <iostream>
#include <string>
#include <filesystem>
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

double ask_tolerance_from_user(const std::string &prompt = "Inserisci la tolleranza (es. 1e-4): ")
{
    double tol;
    std::string input;

    while (true)
    {
        std::cout << prompt;
        std::getline(std::cin, input);

        std::stringstream ss(input);
        if ((ss >> tol) && ss.eof() && tol > 0.0)
        {
            return tol;
        }
        else
        {
            std::cout << "Input non valido. Riprova con un numero positivo in notazione scientifica (es. 1e-4).\n";
        }
    }
}

inline bool ends_with(const std::string &str, const std::string &suffix)
{
    return str.size() >= suffix.size() &&
           str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

bool load_csv_to_dense(const std::string &filename, Eigen::MatrixXd &matrix, const char &sep)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Impossibile aprire il file: " << filename << std::endl;
        return false;
    }

    std::vector<std::vector<double>> values;
    std::string line;

    while (std::getline(file, line))
    {
        std::stringstream ss(line);
        std::string cell;
        std::vector<double> row;

        while (std::getline(ss, cell, sep))
        {
            try
            {
                row.push_back(std::stod(cell));
            }
            catch (...)
            {
                std::cerr << "Errore di conversione nel file CSV." << std::endl;
                return false;
            }
        }

        if (!row.empty())
        {
            values.push_back(row);
        }
    }

    if (values.empty())
    {
        std::cerr << "Il file CSV è vuoto o mal formattato." << std::endl;
        return false;
    }

    size_t rows = values.size();
    size_t cols = values[0].size();

    for (const auto &row : values)
    {
        if (row.size() != cols)
        {
            std::cerr << "Numero di colonne non coerente tra le righe." << std::endl;
            return false;
        }
    }

    matrix = Eigen::MatrixXd(rows, cols);
    for (size_t i = 0; i < rows; ++i)
        for (size_t j = 0; j < cols; ++j)
            matrix(i, j) = values[i][j];

    return true;
}

int main()
{
    // Modalità interattiva
    char choice;
    do
    {
        show_menu();
        std::cin >> choice;
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');  // <-- aggiungi questa riga

        switch (choice)
        {
        case '1':
        {
            lr_test<Eigen::SparseMatrix<double>> lr_t_dm;
            lr_t_dm.test_lr_data_matrix();

            break;
        }
        case '2':
        {

            std::string path;
            std::cout << "Inserisci il path al file della matrice (.csv o .mtx): ";
            std::getline(std::cin, path);

            double tol = ask_tolerance_from_user();

            std::cout << path << " " << tol << std::endl;

            std::string matrix_name = fs::path(path).filename().string();

            if (ends_with(path, ".csv"))
            {
                lr_test<Eigen::MatrixXd> lr_t_m_dense;
                Eigen::MatrixXd d_matrix;
                if (!load_csv_to_dense(path, d_matrix, ','))
                {
                    std::cerr << "Errore durante il caricamento." << std::endl;
                }

                lr_t_m_dense.create_tables(d_matrix, tol, matrix_name);
            }
            else if (ends_with(path, ".mtx"))
            {
                lr_test<Eigen::SparseMatrix<double>> lr_t_m_sparse;
                Eigen::SparseMatrix<double> s_matrix;
                if (!Eigen::loadMarket(s_matrix, path))
                {
                    std::cerr << "Errore nel caricamento della matrice sparsa." << std::endl;
                }

                lr_t_m_sparse.create_tables(s_matrix, tol, matrix_name);
            }
            else
            {
                std::cerr << "Formato file non supportato: " << path << "\n";
            }

            break;
        }
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
