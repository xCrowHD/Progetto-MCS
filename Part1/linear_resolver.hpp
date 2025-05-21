#ifndef LINEAR_RESOLVER
#define LINEAR_RESOLVER

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include<iostream>
#include <chrono>
#include"lr_utils.hpp"

template<typename MatrixType>
class linear_resolver {

private:
    MatrixType _A;
    Eigen::VectorXd _x, _b;
    double _tol;

    int maxIter = 20000;

public:
    linear_resolver(const MatrixType& A, const double& tol)
    {
        _A = A;
        _tol = tol;
        //Create vector X
        int n = _A.rows();  
        _x = Eigen::VectorXd::Ones(n);
        _b = _A * _x;  // genera b coerente con x

    };    // Costruttore

    ~linear_resolver() = default;   // Distruttore

    void run_resolvers()
    {
        std::cout << "\n=== Test JACOBI ===\n";
        jacobi_resolver(_A, _b);
        std::cout << "\n=== Test GAUSS ===\n";
        gauss_resolver(_A, _b);
        std::cout << "\n=== Test GRADIENT ===\n";
        gradient_resolver(_A, _b);
        std::cout << "\n=== Test CONJUGATE GRADIENT ===\n";
        conjugate_gradient_resolver(_A, _b);
    }

private:
    void jacobi_resolver(const MatrixType& A, const Eigen::VectorXd& b)
    {
        //Start Timer
        auto start = std::chrono::high_resolution_clock::now();

        //Calcolo D
        Eigen::VectorXd D = lr_utils::getD(A);

        //Situazione Pre Jacobi
        Eigen::VectorXd x = Eigen::VectorXd::Zero(b.size());
        Eigen::VectorXd r = b - A * x;
        int n_it = 0;
        Eigen::MatrixXd D_inv = D.cwiseInverse().asDiagonal();

        //Jacobi
        while(n_it < maxIter && r.norm() / b.norm() > _tol)
        {
            x = x + D_inv * r;
            r = b - A * x;
            n_it++;
        }

        //End Timer
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration_sec = end - start;

        //std::cout << "Souluzione: " << x << std::endl;
        std::cout << "Iterazioni: " << n_it << std::endl;
        std::cout << "Tolleranza: " << r.norm() / b.norm() << std::endl;
        std::cout << "Tempo in sec: " << duration_sec.count() << std::endl;
        std::cout << "Errore assoluto: " << (x - _x).norm() << std::endl;
        std::cout << "Errore relativo: " << (x - _x).norm() / _x.norm() << std::endl;

    };

    void gauss_resolver(const MatrixType& A, const Eigen::VectorXd& b)
    {
        //Start Timer
        auto start = std::chrono::high_resolution_clock::now();

        //Situazione Pre Gauss
        Eigen::VectorXd x = Eigen::VectorXd::Zero(b.size());
        Eigen::VectorXd r = b - A * x;
        int n_it = 0;

        //Build P
        MatrixType P = lr_utils::getLowerTriangular(A);
        //Build first y
        Eigen::VectorXd y = lr_utils::forward_substitution(P, r);

        while(n_it < maxIter && r.norm() / b.norm() > _tol)
        {
            x = x + y;
            r = b - A * x;
            y = lr_utils::forward_substitution(P, r);
            n_it++;
        }

        //End Timer
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration_sec = end - start;

        //std::cout << "Souluzione: " << x << std::endl;
        std::cout << "Iterazioni: " << n_it << std::endl;
        std::cout << "Tolleranza: " << r.norm() / b.norm() << std::endl;
        std::cout << "Tempo in sec: " << duration_sec.count() << std::endl;
        std::cout << "Errore assoluto: " << (x - _x).norm() << std::endl;
        std::cout << "Errore relativo: " << (x - _x).norm() / _x.norm() << std::endl;

    };

    void gradient_resolver(const MatrixType& A, const Eigen::VectorXd& b)
    {
        // Start timer
        auto start = std::chrono::high_resolution_clock::now();

        Eigen::VectorXd x = Eigen::VectorXd::Zero(b.size());
        Eigen::VectorXd r = b - A * x;  // <-- Calcolo iniziale del residuo

        Eigen::VectorXd y = Eigen::VectorXd::Zero(b.size());

        int n_it = 0;

        while(n_it < maxIter && r.norm() / b.norm() > _tol)
        {
            y = A * r;
            double a = r.dot(r);
            double b_val = r.dot(y);  // non facciamo la trasposta ma direttamente .dot() che è piu veloce e semplice
            double alpha = a / b_val;

            x = x + alpha * r;
            r = b - A * x; 

            n_it++;
        }

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration_sec = end - start;

        //std::cout << "Souluzione: " << x << std::endl;
        std::cout << "Iterazioni: " << n_it << std::endl;
        std::cout << "Tolleranza: " << r.norm() / b.norm() << std::endl;
        std::cout << "Tempo in sec: " << duration_sec.count() << std::endl;
        std::cout << "Errore assoluto: " << (x - _x).norm() << std::endl;
        std::cout << "Errore relativo: " << (x - _x).norm() / _x.norm() << std::endl;

    };

    void conjugate_gradient_resolver(const MatrixType& A, const Eigen::VectorXd& b)
    {
        // Start timer
        auto start = std::chrono::high_resolution_clock::now();

        Eigen::VectorXd x = Eigen::VectorXd::Zero(b.size());
        Eigen::VectorXd r = b - A * x;  // <-- Calcolo iniziale del residuo
        Eigen::VectorXd d = r;

        Eigen::VectorXd y = Eigen::VectorXd::Zero(b.size());
        Eigen::VectorXd w = Eigen::VectorXd::Zero(b.size());

        int n_it = 0;

        while(n_it < maxIter && r.norm() / b.norm() > _tol)
        {
            y = A * d;
            double alpha = r.dot(d) / d.dot(y);

            x = x + (alpha * d);

            r = b - A * x;

            w = A * r;
            double beta = d.dot(w) / d.dot(y);

            d = r - beta * d;

            n_it++;
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration_sec = end - start;

        //std::cout << "Souluzione: " << x << std::endl;
        std::cout << "Iterazioni: " << n_it << std::endl;
        std::cout << "Tolleranza: " << r.norm() / b.norm() << std::endl;
        std::cout << "Tempo in sec: " << duration_sec.count() << std::endl;
        std::cout << "Errore assoluto: " << (x - _x).norm() << std::endl;
        std::cout << "Errore relativo: " << (x - _x).norm() / _x.norm() << std::endl;

    };
   
};

#endif // JACOBI_ITER_RES
