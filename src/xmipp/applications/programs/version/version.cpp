// Auto-generated code to get compilation Info
#include <sys/utsname.h>
#include <iostream>
#include <string.h>
using namespace std;

int main(int argc, char** argv){

    if (argc>2)
    {
        std::cout << "Incorrect parameter" << std::endl;
        return 1;
    }
    int shrt = 0;
    if (argc>1)
    {
        if((strcmp(argv[1], "--short") == 0))
        {
            shrt = 1;
        }else{
            std::cout << "Incorrect parameter: " << argv[1] << std::endl;
            return 2;
        }
    }

    if (shrt==1)
    {
        std::cout << "Xmipp version: devel" << std::endl;
    }else{
        struct utsname utsname; // stores the data returned by uname()
        struct utsname *utsname_ptr = &utsname; // pointer to the struct holding the data returned by uname()
    
        int ret;
        ret = uname(utsname_ptr);
    
        std::cout << std::endl;
        std::cout << "  [4mXmipp version[24m: [1mdevel[0m" << std::endl;
        std::cout << std::endl;
        std::cout << "  Release date:     not released yet" << std::endl;
        std::cout << "  Xmipp branch:     (no git repo detected)" << std::endl;
        std::cout << "  Core branch :     devel (417728f)" << std::endl;
        std::cout << "  Viz branch:       devel (c733d39)" << std::endl;
        std::cout << "  Compilation date: 27/04/2020" << std::endl;
        std::cout << "  Compiler:         g++ " << __VERSION__ << std::endl;
        std::cout << "  Compiling system: " << utsname.machine << " " << utsname.sysname
                  << " " << utsname.release << std::endl 
                  << "                    " << utsname.version << std::endl;
        std::cout << std::endl;
    }
    return 0;
}
