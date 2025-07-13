#ifndef DCT_FDCT_c
#define DCT_FDCT_c

#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include <opencv2/core.hpp>

class DCT_FDCT
{

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