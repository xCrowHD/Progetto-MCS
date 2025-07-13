#ifndef DCT_c
#define DCT_c

#include<iostream>
#include<chrono>
#include<math.h>
#include<Eigen/Dense>
#include<cmath>

class DCT
{
    
public:
    /* data */
    double time;
    
public:
    DCT() = default;
    ~DCT() = default;

    Eigen::VectorXd run_DCT(const Eigen::VectorXd& v)
    {
        auto start = std::chrono::high_resolution_clock::now();
        const int n = v.size();

        Eigen::VectorXd a(n);

        for (int k = 0; k < n; k++)
        {
            double d = (k == 0) ? std::sqrt(1.0 / n) : std::sqrt(2.0 / n);
            double sum = 0;
            
            for (int i = 0; i < n; i++)
            {
                double angle = k * M_PI * (2.0 * i + 1) / (2.0 * n);
                sum = sum + (v(i) * cos(angle));
            }

            a(k) = sum * d;
        }
        auto end = std::chrono::high_resolution_clock::now();
        time = (end - start).count();
        return a;
    }

    Eigen::MatrixXd run_DCT2(const Eigen::MatrixXd& mat)
    {
        auto start = std::chrono::high_resolution_clock::now();
        int n = mat.rows(); // assume quadrato n x n
        Eigen::MatrixXd tmp(n, n);
        Eigen::MatrixXd res(n, n);

        // DCT sulle righe
        for(int r = 0; r < n; ++r)
        {
            tmp.row(r) = run_DCT(mat.row(r).transpose()).transpose();
        }
        
        // DCT sulle colonne
        for(int c = 0; c < n; ++c)
        {
            res.col(c) = run_DCT(tmp.col(c));
        }
        auto end = std::chrono::high_resolution_clock::now();
        time = (end - start).count();
        return res;
    }

    
    
};

#endif
