#ifndef LR_TEST_c
#define LR_TEST_c
#include "linear_resolver.hpp"
#include <filesystem>
#include <vector>
#include <string>
#include <limits>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include "tabulate.hpp"

using namespace tabulate;
using Row_t = Table::Row_t;
namespace fs = std::filesystem;

std::map<std::string, std::vector<int>> metodo_min_count{
    {"Jacobi", std::vector<int>(5, 0)},
    {"Gauss", std::vector<int>(5, 0)},
    {"Gradiente", std::vector<int>(5, 0)},
    {"Gradiente Coniugato", std::vector<int>(5, 0)}};

template <typename MatrixType>
class lr_test
{

public:
    lr_test() = default;
    ~lr_test() = default;

    void test_lr()
    {

        std::string directory = "./dati";
        std::vector<std::string> mtx_files;
        std::vector<double> tols = {1e-4, 1e-6, 1e-8, 1e-10};
        // Scansione directory per file .mtx
        for (const auto &entry : fs::directory_iterator(directory))
        {
            if (entry.is_regular_file() && entry.path().extension() == ".mtx")
            {
                mtx_files.push_back(entry.path().string());
            }
        }

        if (mtx_files.empty())
        {
            std::cerr << "Nessun file .mtx trovato in " << directory << std::endl;
            return;
        }

        // std::cout << mtx_files.size() << std::endl;

        // Per ogni file mtx trovato, esegui dtc_fdct e raccogli risultati
        for (const auto &file_path : mtx_files)
        {
            MatrixType matrix;
            if (!Eigen::loadMarket(matrix, file_path))
            {
                std::cerr << "Errore nel caricamento della matrice sparsa." << std::endl;
            }

            std::string matrix_name = fs::path(file_path).filename().string();

            for (const auto &t : tols)
            {
                // TITLE
                Table title;
                title.format().font_align(FontAlign::center);
                std::ostringstream oss;
                size_t rows = matrix.rows();
                size_t cols = matrix.cols();
                size_t nnz = 0;
                if constexpr (std::is_same<MatrixType, Eigen::SparseMatrix<double>>::value)
                {
                    nnz = matrix.nonZeros();
                }

                std::string matrix_type = std::is_same<MatrixType, Eigen::MatrixXd>::value ? "Densa" : "Sparsa";
                oss << matrix_name << " (" << matrix_type << ", "
                    << rows << "x" << cols;
                if (nnz > 0)
                    oss << ", NNZ = " << nnz;
                oss << ", tol = " << std::scientific << std::setprecision(1) << t << ")";
                title.add_row(Row_t{oss.str()}); // title

                // CONTENT
                Table table;
                table.format().border_color(Color::red);

                std::string nome_metodo = "Nome Metodo";
                std::string n_iter = "Numero iterazioni";
                std::string res_rel = "Residuo relativo ||b - Ax|| / ||b||";
                std::string err_rel = "Errore relativo ||x - x_exact|| / ||x_exact||";
                std::string err_ass = "Errore assoluto ||x - x_exact||";
                std::string time = "Tempo in secondi";

                table.add_row({nome_metodo, n_iter, res_rel, err_rel, err_ass, time});
                linear_resolver<MatrixType> lr(matrix, t);

                std::tuple<Eigen::VectorXd, int, double, double, double, double> jac_result = lr.jacobi_resolver(lr.getA(), lr.getb());
                add_to_table(table, jac_result, std::string("Jacobi"));

                std::tuple<Eigen::VectorXd, int, double, double, double, double> gas_result = lr.gauss_resolver(lr.getA(), lr.getb());
                add_to_table(table, gas_result, std::string("Gauss"));

                std::tuple<Eigen::VectorXd, int, double, double, double, double> grad_result = lr.gradient_resolver(lr.getA(), lr.getb());
                add_to_table(table, grad_result, std::string("Gradiente"));

                std::tuple<Eigen::VectorXd, int, double, double, double, double> cgrad_result = lr.conjugate_gradient_resolver(lr.getA(), lr.getb());
                add_to_table(table, cgrad_result, std::string("Gradiente Coniugato"));

                title.format().font_color(Color::cyan).font_style({FontStyle::bold});
                table[0].format().font_color(Color::yellow).font_style({FontStyle::italic});

                color_min_values_per_column(table);

                std::cout << title << "\n";
                std::cout << table << "\n\n";
            }
        }

        stats_table();
    }

private:
    void add_to_table(Table &t, const std::tuple<Eigen::VectorXd, int, double, double, double, double> &result, std::string nome)
    {
        // Eigen::VectorXd aprox_sol = std::get<0>(result);
        int iter_n = std::get<1>(result);
        double residuo_rel = std::get<2>(result);
        double time = std::get<3>(result);
        double err_abs = std::get<4>(result);
        double err_rel = std::get<5>(result);

        std::ostringstream oss_res, oss_err_rel, oss_err_abs, oss_time;

        oss_res << std::scientific << std::setprecision(2) << residuo_rel;
        oss_err_rel << std::scientific << std::setprecision(2) << err_rel;
        oss_err_abs << std::scientific << std::setprecision(2) << err_abs;
        oss_time << std::scientific << std::setprecision(2) << time;

        t.add_row({nome,
                   std::to_string(iter_n),
                   oss_res.str(),
                   oss_err_rel.str(),
                   oss_err_abs.str(),
                   oss_time.str()});
    }

    void color_min_values_per_column(Table &table, bool update_stats = true, bool find_min = true)
    {
        if (table.size() <= 1)
            return; // Nessuna riga dati

        size_t num_rows = table.size();
        size_t num_cols = table[0].size();

        for (size_t col = 1; col < num_cols; ++col)
        {
            double extreme_val = find_min ? std::numeric_limits<double>::max()
                                          : std::numeric_limits<double>::lowest(); // lowest = valore più piccolo rappresentabile
            std::vector<size_t> extreme_row_indices;

            // Cerca minimo o massimo
            for (size_t row = 1; row < num_rows; ++row)
            {
                try
                {
                    double val = std::stod(table[row][col].get_text());

                    if ((find_min && val < extreme_val) || (!find_min && val > extreme_val))
                    {
                        extreme_val = val;
                        extreme_row_indices.clear();
                        extreme_row_indices.push_back(row);
                    }
                    else if (val == extreme_val)
                    {
                        extreme_row_indices.push_back(row);
                    }
                }
                catch (...)
                {
                    // ignora celle non numeriche
                }
            }

            // Colora di verde il valore minimo/max
            for (size_t row_index : extreme_row_indices)
            {
                table[row_index][col].format().font_color(Color::green).font_style({FontStyle::bold});

                // 3. Aggiorna statistiche solo per minimi
                if (update_stats && find_min)
                {
                    std::string metodo = table[row_index][0].get_text();
                    switch (col)
                    {
                    case 1:
                        metodo_min_count[metodo][0]++;
                        break;
                    case 2:
                        metodo_min_count[metodo][1]++;
                        break;
                    case 3:
                        metodo_min_count[metodo][2]++;
                        break;
                    case 4:
                        metodo_min_count[metodo][3]++;
                        break;
                    case 5:
                        metodo_min_count[metodo][4]++;
                        break;
                    default:
                        break;
                    }
                }
            }
        }
    }

    void stats_table()
    {
        Table title;
        title.add_row({"Statistiche Finali (Conteggi Minimi)"});

        Table final_table;
        final_table.format().border_color(Color::blue);

        final_table.add_row({"Metodo", "Iter min", "Residuo min", "Err Rel min", "Err Ass min", "Tempo min"});

        for (const auto &[metodo, counts] : metodo_min_count)
        {
            final_table.add_row({metodo,
                                 std::to_string(counts[0]),
                                 std::to_string(counts[1]),
                                 std::to_string(counts[2]),
                                 std::to_string(counts[3]),
                                 std::to_string(counts[4])});
        }

        final_table[0].format().font_color(Color::yellow).font_style({FontStyle::bold});
        color_min_values_per_column(final_table, false, false);
        std::cout << "\n"
                  << title << "\n";
        std::cout << final_table << "\n";
    }
};
#endif