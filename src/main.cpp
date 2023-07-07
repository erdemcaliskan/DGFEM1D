#include "util/Configuration.h"
#include "util/globalMacros.h"

#if USE_COMPLEX
#if !VCPP
#include <Eigen/Dense>
#endif
#endif

using namespace std;

// usage windows: add -s (serialNumber) [e.g. -s 1] to command arguments
//			from cmd: use the compiled Target.exe in Release folder for this, may need to move it one folder up and run
//			Target.exe 	-s (serialNumber) 
// other: ./solver -s (serialNumber) [e.g. ./solver -s 1] 
int main(int argc, char *argv[])
{
	serialNumber = 0;
	string configName = "configFile.txt";
	// For Ali: uncomment the line below
	// configName = "configFileMM.txt";
	if (argc > 0)
	{
		for (int i = 1; i < argc; ++i)
		{
			if (strcmp(argv[i], "-s") == 0)
			{
				serialNumber = (unsigned int)atoi(argv[++i]);
			}
			else if (strcmp(argv[i], "-c") == 0)
			{
				configName = (string)argv[++i];
			}
			else
			{
				cout << "i\t" << i << '\n';
				cout << "argv[i]\t" << argv[i] << '\n';
				THROW("ERROR: Invalid flag--value option");
			}
		}
	}

//    std::cout << "Hello, world!\n";
    setGlobalMembers();

    Configuration conf;
    conf.Main_SolveDomain(configName, serialNumber);


    return 0;

    // Example
#if !VCPP
#if USE_COMPLEX
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
#endif
#endif
#if VCPP
	cout << "successful solution\n";
	getchar();
#endif
    return 0;
}
