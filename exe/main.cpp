#include "interpolate.h"
#include <iostream>     // std::cout

int main()
{
    std::vector<std::pair<float, float>> coordinates{ { {1.0f, 5.0f}, {2.0f, 2.0f}, {3.0f, 2.5f}, {4.0f, 5.0f} } };
    std::pair<float, float> out1;
    std::pair<float, float> out2;
    std::pair<float, float> out3;
    
    out1.first = 2.25f;
    out2.first = 2.5f;
    out3.first = 2.75f;
    
    out1.second = interpolate::hermite(coordinates, out1.first);
    out2.second = interpolate::hermite(coordinates, out2.first);
    out3.second = interpolate::hermite(coordinates, out3.first);
    
    std::cout << "Hermite: \n";
    std::cout << "x1 = " << out1.first << ", y1 = " << out1.second << '\n';
    std::cout << "x2 = " << out2.first << ", y2 = " << out2.second << '\n';
    std::cout << "x3 = " << out3.first << ", y3 = " << out3.second << '\n';

    out1.second = interpolate::hyman(coordinates, out1.first);
    out2.second = interpolate::hyman(coordinates, out2.first);
    out3.second = interpolate::hyman(coordinates, out3.first);

    std::cout << "Hyman: \n";
    std::cout << "x1 = " << out1.first << ", y1 = " << out1.second << '\n';
    std::cout << "x2 = " << out2.first << ", y2 = " << out2.second << '\n';
    std::cout << "x3 = " << out3.first << ", y3 = " << out3.second << '\n';

    return 0;
}
