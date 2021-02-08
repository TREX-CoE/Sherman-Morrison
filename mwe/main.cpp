#include <iostream>


int main()
{

    typedef int (*to2D)[3]; //pint2 is a pointer to an array of 2 ints
    int linArray[9] = {0,1,2,3,4,5,6,7,8};
    
    to2D dArray = (to2D)linArray;
    
    std::cout << dArray[0][0] << std::endl;
    std::cout << dArray[0][1] << std::endl;
    std::cout << dArray[0][2] << std::endl;
    std::cout << dArray[1][0] << std::endl;
    std::cout << dArray[1][1] << std::endl;
    std::cout << dArray[1][2] << std::endl;
    std::cout << dArray[2][0] << std::endl;
    std::cout << dArray[2][1] << std::endl;
    std::cout << dArray[2][2] << std::endl;

    return 0;
}
