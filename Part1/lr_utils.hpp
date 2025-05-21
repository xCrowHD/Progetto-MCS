#ifndef LR_UTILS
#define LR_UTILS

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace lr_utils
{
    Eigen::VectorXd getD(const Eigen::SparseMatrix<double>& A)
    {
        Eigen::VectorXd D(A.rows());
        for (int k = 0; k < A.rows(); ++k) 
        {
            D[k] = A.coeff(k, k);
        }
        return D;
    }

    Eigen::VectorXd getD(const Eigen::MatrixXd& A)
    {
        return A.diagonal();
    }

    //Custom solver perchè quelli gia presenti faticano sui metodi iterativi e non sono pensati per essere usati con essi
    Eigen::VectorXd forward_substitution(const Eigen::MatrixXd& P, const Eigen::VectorXd& r) 
    {
        const int n = P.rows();
        Eigen::VectorXd y(n);
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < i; ++j) {
                sum += P(i, j) * y(j);
            }
            y(i) = (r(i) - sum) / P(i, i);
        }
        return y;
    }

    Eigen::VectorXd forward_substitution(const Eigen::SparseMatrix<double>& P, const Eigen::VectorXd& r) 
    {
        const int n = P.rows();
        Eigen::VectorXd y(n);
        y.setZero();
    
        // Convertiamo temporaneamente P in formato RowMajor per iterare per riga
        Eigen::SparseMatrix<double, Eigen::RowMajor> P_row = P;
    
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            double diag = 0.0;
    
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(P_row, i); it; ++it) {
                int j = it.col();
                if (j < i) {
                    sum += it.value() * y(j);
                } else if (j == i) {
                    diag = it.value();
                }
            }
            y(i) = (r(i) - sum) / diag;
        }
    
        return y;
    }
    
    
    //Funzioni per fare P per il metodo di gauss
    Eigen::MatrixXd getLowerTriangular(const Eigen::MatrixXd& A) 
    {
        return A.triangularView<Eigen::Lower>();
    }

    Eigen::SparseMatrix<double> getLowerTriangular(const Eigen::SparseMatrix<double>& A) 
    {
        using Triplet = Eigen::Triplet<double>;
        std::vector<Triplet> triplets;
    
        for (int k = 0; k < A.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
                if (it.row() >= it.col()) {  // parte inferiore (inclusa diagonale)
                    triplets.emplace_back(it.row(), it.col(), it.value());
                }
            }
        }
    
        Eigen::SparseMatrix<double> P(A.rows(), A.cols());
        P.setFromTriplets(triplets.begin(), triplets.end());
        P.makeCompressed(); // opzionale ma consigliato
        return P;
    }
    
} // namespace name

#endif