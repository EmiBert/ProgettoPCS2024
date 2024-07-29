//#include <iostream>
//#include<iomanip>
#include "Eigen/Eigen"
//#include "Utils.hpp"
//#include "DFNlibrary.hpp"
#include "TestDFN.hpp"


using namespace std;
using namespace Eigen;
using namespace DFNlibrary;

int main(int argc, char ** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
