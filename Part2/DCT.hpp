#include<iostream>
#include<chrono>
#include<math.h>
#include<Eigen/Dense>
#include <cmath>

class DCT
{
    
private:
    /* data */
public:
    DCT() = default;
    ~DCT() = default;

    Eigen::VectorXd run_DCT(const Eigen::VectorXd& v)
    {
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
        
        return a;
    }

    
};
