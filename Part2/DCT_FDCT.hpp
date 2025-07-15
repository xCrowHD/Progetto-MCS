#ifndef DCT_FDCT_c
#define DCT_FDCT_c

#include <fstream>
#include <sstream>
#include <tuple>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include <opencv2/core.hpp>
#include"DCT.hpp"
#include "json.hpp" 
#include <cstdlib>

using json = nlohmann::json;

class DCTFDCT
{

public:

    DCTFDCT() = default;
    ~DCTFDCT() = default;


    void test_matrix()
    {
        int count = 5;

        std::string matrix[] {"./matrix/matrice_16x16.csv", "./matrix/matrice_32x32.csv",
             "./matrix/matrice_64x64.csv", "./matrix/matrice_128x128.csv", "./matrix/matrice_256x256.csv"};
        
        json j;
        j["results"] = json::array();

        for (int i = 0; i < count; i++)
        {
            std::tuple<std::string, double, double> result = dtc_fdct(matrix[i]);
            std::string label = std::get<0>(result);
            double hm_time = std::get<1>(result);
            double cv_time = std::get<2>(result);


            j["results"].push_back({
                {"label", label},
                {"homemade_time", hm_time},
                {"opencv_time", cv_time}
            });

        }
        
        // Scrivo su file JSON
        std::ofstream file("./data/results.json");
        file << j.dump(4);  // indentazione 4 spazi
        file.close();
        int ret = system("python3 data/dct_fdct_plot.py &");
        if(ret != 0) {
            std::cerr << "Errore nell'esecuzione dello script python, codice: " << ret << "\n";
        }
        std::cout << "test" << std::endl;
        // std::cout << label << std::endl;
        // std::cout << "Homemade DCT time: " << hm_time << std::endl;
        // std::cout << "OpenCV DCT time: " << cv_time << std::endl;
        
    }

    std::tuple<std::string, double, double> dtc_fdct(const std::string& matrix_csv_name)
    {
        Eigen::MatrixXd mat = read_csv_to_eigen(matrix_csv_name);
        if (mat.size() == 0) {
            throw std::runtime_error("Input matrix is empty!");
        }
        cv::Mat cv_mat_input = eigen_to_cv_mat(mat);
        cv::Mat cv_mat_output;

        DCT dct;
        Eigen::MatrixXd a = dct.run_DCT2(mat);

        auto start = std::chrono::high_resolution_clock::now();
        cv::dct(cv_mat_input, cv_mat_output);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = end - start;
        double cv_time = elapsed.count();

        int r = cv_mat_input.rows;
        int c = cv_mat_input.cols;
        //std::cout << a << std::endl;
        //std::cout << cv_mat_output << std::endl;
        std::string s = std::to_string(r) + "x" + std::to_string(c);
        return {s, dct.time, cv_time};
    }

private:

    Eigen::MatrixXd read_csv_to_eigen(const std::string& filename) {
        std::ifstream file(filename);
        std::string line;
        std::vector<std::vector<double>> data;
    
        while (std::getline(file, line)) {
            std::stringstream line_stream(line);
            std::string cell;
            std::vector<double> row;
    
            while (std::getline(line_stream, cell, ',')) {
                row.push_back(std::stod(cell));
            }
            data.push_back(row);
        }
    
        if (data.empty()) return Eigen::MatrixXd();
    
        size_t rows = data.size();
        size_t cols = data[0].size();
    
        Eigen::MatrixXd mat(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            if (data[i].size() != cols) {
                throw std::runtime_error("CSV row size mismatch");
            }
            for (size_t j = 0; j < cols; ++j) {
                mat(i, j) = data[i][j];
            }
        }
        return mat;
    }

    cv::Mat eigen_to_cv_mat(const Eigen::MatrixXd& eigen_mat) {
        cv::Mat cv_mat(eigen_mat.rows(), eigen_mat.cols(), CV_64F);
        for (int i = 0; i < eigen_mat.rows(); ++i) {
            for (int j = 0; j < eigen_mat.cols(); ++j) {
                cv_mat.at<double>(i, j) = eigen_mat(i, j);
            }
        }
        return cv_mat;
    }

};

#endif