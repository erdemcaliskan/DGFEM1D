#include "util/Configuration.h"
#include "util/globalMacros.h"

#include <Eigen/Dense>

using namespace std;

int main(int, char **)
{
    std::cout << "Hello, world!\n";
    setGlobalMembers();

    string configName = "";
    Configuration conf;
    conf.Main_SolveDomain(configName);
    return 0;

    // Example
    Eigen::MatrixXcd m(2, 2); // MatrixXcd typedef Matrix< std::complex< double >,
                              // Dynamic, Dynamic >
    m << 1.0f + 2.0if, 2.0f + 1.0if, 3.0f - 1.0if, 4.0f - 2.0if;
    cout << "Matrix m:\n" << m << endl;

    Eigen::MatrixXcd n(2, 2);
    n(0, 0) = 3.0f + 4.0if;
    n(1, 0) = 2.5f + 10.0if;
    n(0, 1) = -1.0f + 0.35if;
    n(1, 1) = n(1, 0) + n(0, 1);
    cout << "Matrix n:\n" << n << endl;

    // Finding Eigen Values of m*(n)T
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
    ces.compute(m * n.transpose());
    cout << "Eigen Values:\n" << ces.eigenvalues() << endl;

    return 0;
}
