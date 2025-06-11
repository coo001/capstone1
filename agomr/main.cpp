#include "include/sophomr.h"
#include <string>

std::string OMR = "OMR";
std::string OMD = "OMD";

int main(int argc, char* argv[]) 
{
    if (argc == 4) {
        if (OMR.compare(argv[1]) == 0 && atoi(argv[2]) == 65536 && atoi(argv[3]) == 50) {
            param_OMR_65536_50();
        } else if (OMD.compare(argv[1]) == 0 && atoi(argv[2]) == 65536 && atoi(argv[3]) == 50) {
            param_OMD_65536_50();
        } else if (OMR.compare(argv[1]) == 0 && atoi(argv[2]) == 524288 && atoi(argv[3]) == 50) {
            param_OMR_524288_50();
        } else if (OMD.compare(argv[1]) == 0 && atoi(argv[2]) == 524288 && atoi(argv[3]) == 50) {
            param_OMD_524288_50();
        } else {
            std::cout << "No predefined parameters detected. Using the default parameters set in global.h" << std::endl;
        }
    }
    sophomr();
}