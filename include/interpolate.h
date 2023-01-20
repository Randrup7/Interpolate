
#include <vector>       // std::vector
#include <algorithm>    // std::upper_bound
#include <iterator>     // std::prev
#include <utility>      // std::pair
#include <cmath>        // std::pow

namespace interpolate 
{
    float hermite(std::vector<std::pair<float, float>>& coord, float outX);
}