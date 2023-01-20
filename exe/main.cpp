#include "interpolate.h"
#include <iostream>     // std::cout

int main()
{
    std::vector<std::pair<float, float>> coordinates{ { {0.5f, 2.0f}, {1.0f, 2.3f}, {3.0f, 3.44f}, {5.0f, 3.657f} } };
    std::pair<float, float> out1;
    
    out1.first = 0.5f;
    
    out1.second = interpolate::hermite(coordinates, out1.first);
    
    std::cout << "x1 = " << out1.first << ", y1 = " << out1.second << '\n';

    return 0;
}
