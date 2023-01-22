#include "interpolate.h"

// This function returns the general formula of a cubicSpline. tau is within unit interval
float cubicSpline(float a, float b, float c, float d, float tau)
{
        return a + b * tau + c * std::pow(tau, 2) + d * std::pow(tau, 3);
}

// Hermite interpolation with a vector of pairs
float interpolate::hermite(std::vector<std::pair<float, float>>& coord, float outX)
{
        using iteratorType = std::vector<std::pair<float, float>>::iterator;

        // uBound is first iterator with first-value above outX using lambda function
        iteratorType uBound{ std::upper_bound(coord.begin(), coord.end(), outX, 
                [](float value,  std::pair<float, float>& coord) { return value < coord.first; } ) };

        // lBound is iterator before uBound
        iteratorType lBound = std::prev(uBound, 1);

        // dist is distance from first pair to uBound
        std::iterator_traits<iteratorType>::difference_type dist{ std::distance(coord.begin(), uBound) };
                
        float bi;
        float bk;

        if (dist == coord.size()) { return coord.back().second; }   // Flat extrapolation
        else if (dist == 0) { return coord.front().second; }        // Flat extrapolation

        else if (dist == 1) // when lbound - 1 doesn't exist
        {   
                bi = ((std::next(uBound, 1)->first + uBound->first - 2 * lBound->first) * 
                        (uBound->second - lBound->second) /
                        (uBound->first - lBound->first) -
                        (uBound->first - lBound->first) *
                        (std::next(uBound, 1)->second - uBound->second) /
                        (std::next(uBound, 1)->first - uBound->first)) /
                        (std::next(uBound, 1)->first - lBound->first);  // Hagan West 2006 equation 18

                bk = ((std::next(uBound, 1)->first - uBound->first) *
                        (uBound->second - lBound->second) /
                        (uBound->first - lBound->first) +
                        (uBound->first - lBound->first) *
                        (std::next(uBound, 1)->second - uBound->second) /
                        (std::next(uBound, 1)->first - uBound->first)) /
                        (std::next(uBound, 1)->first - lBound->first);  // Hagan West 2006 equation 19
        }

        else if (dist == coord.size() - 1) // when ubound + 1 doesn't exist
        {        
                bi = ((uBound->first - lBound->first) * 
                        (lBound->second - std::prev(lBound, 1)->second) /
                        (lBound->first - std::prev(lBound, 1)->first) +
                        (lBound->first - std::prev(lBound, 1)->first) *
                        (uBound->second - lBound->second) /
                        (uBound->first - lBound->first)) /
                        (uBound->first - std::prev(lBound, 1)->first);  // Hagan West 2006 equation 19

                bk = -((uBound->first - lBound->first) * 
                        (lBound->second - std::prev(lBound, 1)->second) /
                        (lBound->first - std::prev(lBound, 1)->first) -
                        (2 * uBound->first - lBound->first - std::prev(lBound, 1)->first) * 
                        (uBound->second - lBound->second) /
                        (uBound->first - lBound->first)) /
                        (uBound->first - std::prev(lBound, 1)->first);  // Hagan West 2006 equation 20
        }

        else // uBound = k, lBound = i
        {
                bi = ((uBound->first - lBound->first) * 
                        (lBound->second - std::prev(lBound, 1)->second) /
                        (lBound->first - std::prev(lBound, 1)->first) +
                        (lBound->first - std::prev(lBound, 1)->first) *
                        (uBound->second - lBound->second) /
                        (uBound->first - lBound->first)) /
                        (uBound->first - std::prev(lBound, 1)->first);  // Hagan West 2006 equation 19

                bk = ((std::next(uBound, 1)->first - uBound->first) *
                        (uBound->second - lBound->second) /
                        (uBound->first - lBound->first) +
                        (uBound->first - lBound->first) *
                        (std::next(uBound, 1)->second - uBound->second) /
                        (std::next(uBound, 1)->first - uBound->first)) /
                        (std::next(uBound, 1)->first - lBound->first);  // Hagan West 2006 equation 19
        }

        float hi = uBound->first - lBound->first;           // lenght between x-coordinates lBound and uBound
        float mi = (uBound->second - lBound->second) / hi;  // Hagan West 2006 equation 15
        float ci = (3 * mi - bk - 2 * bi) / hi;             // Hagan West 2006 equation 16
        float di = (bk + bi - 2 * mi) * std::pow(hi, -2);   // Hagan West 2006 equation 17

        
        //float outY = lBound->second + bi * (outX - lBound->first) + ci * std::pow(outX - lBound->first, 2) + di * std::pow(outX - lBound->first, 3);
        // this Bessel Hermite Cubic Spline may not preserve monotonicity between knot points. See Hyman filter for this.

        return cubicSpline(lBound->second, bi, ci, di, outX - lBound->first);
}

// generally results seem too low, check formulas
float interpolate::hyman(std::vector<std::pair<float, float>>& coord, float outX)
{
        using iteratorType = std::vector<std::pair<float, float>>::iterator;

        // uBound is first iterator with first-value above outX using lambda function
        iteratorType uBound{ std::upper_bound(coord.begin(), coord.end(), outX, 
        [](float value,  std::pair<float, float>& coord) { return value < coord.first; } ) };

        std::iterator_traits<iteratorType>::difference_type dist{ std::distance(coord.begin(), uBound) };

        if (dist == coord.size()) { return coord.back().second; }   // Flat extrapolation
        else if (dist == 0) { return coord.front().second; }        // Flat extrapolation

        // lBound is iterator before uBound
        iteratorType lBound = std::prev(uBound, 1);

        float hh = lBound->first - std::prev(lBound, 1)->first;
        float hi = uBound->first - lBound->first;
        float hk = std::next(uBound, 1)->first - uBound->first;

        float mh = (lBound->second - std::prev(lBound, 1)->second) / hh;
        float mi = (uBound->second - lBound->second) / hi;
        float mk = (std::next(uBound, 1)->second - uBound->second) / hk;

        float bi;
        float bk;
        
        if (dist == 1)
        {
                bi = 0.0f; // Hagen West 2006 equation 25

                if (mi * mk < 0.0f)
                        bk = 0.0f;
                else
                        bk = (3.0f * mi * mk) / (std::max(mi, mk) + 2.0f * std::min(mi, mk)); // Hagen West 2006 equation 26
        }

        else if (dist == coord.size() - 1)
        {
                if (mh * mi < 0.0f)
                        bi = 0.0f;
                else
                        bi = (3.0f * mh * mi) / (std::max(mh, mi) + 2.0f * std::min(mh, mi)); // Hagen West 2006 equation 26

                bk = 0.0f; // Hagen West 2006 equation 25
        }

        else
        {
                if (mh * mi < 0.0f)
                        bi = 0.0f;
                else
                        bi = (3.0f * mh * mi) / (std::max(mh, mi) + 2.0f * std::min(mh, mi)); // Hagen West 2006 equation 26

                if (mi * mk < 0.0f)
                        bk = 0.0f;
                else
                        bk = (3.0f * mi * mk) / (std::max(mi, mk) + 2.0f * std::min(mi, mk)); // Hagen West 2006 equation 26
        }

        // Adjustment to b, Hagen West 2006 equation 27
        // Make sure that the first below is sufficient to conlcude local decreasing monotonicity (second vice versa)
        if (mi < 0)
        {
                bi = std::max(std::min(0.0f, bi), 3.0f * std::max(mh, mi));
                bk = std::max(std::min(0.0f, bk), 3.0f * std::max(mi, mk));
        }
        if (mi > 0)
        {
                bi = std::min(std::max(0.0f, bi), 3.0f * std::min(mh, mi));
                bk = std::min(std::max(0.0f, bk), 3.0f * std::min(mi, mk));
        }

        float ci = (3 * mi - bk - 2 * bi) / hi;             // Hagan West 2006 equation 16
        float di = (bk + bi - 2 * mi) * std::pow(hi, -2);   // Hagan West 2006 equation 17

        return cubicSpline(lBound->second, bi, ci, di, outX - lBound->first);
}

