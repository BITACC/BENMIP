#pragma once
#include <algorithm>    // std::min_element, std::max_element
#include <functional>   // std::negate
#include <math.h>       /* fabs */

inline std::string trim(const std::string& str,
	const std::string& whitespace = " \t")
{
	const auto strBegin = str.find_first_not_of(whitespace);
	if (strBegin == std::string::npos)
		return ""; // no content

	const auto strEnd = str.find_last_not_of(whitespace);
	const auto strRange = strEnd - strBegin + 1;

	return str.substr(strBegin, strRange);
}
enum ConstraintType { UNIQUELY_INTEGER=10, UNIQUELY_CONTINUOUS, MIXED };
enum VariableType		{ INTEGER=20, CONTINUOUS };
enum cutType     { OPTIMALITY, FEASIBILITY, COMBINATORIAL };
enum Status     { OPTIMAL, INFEASIBILE, UNKNOWN};
//enum RowType		{ MASTER_ROW = 13, SUBPROBLEM_ROW };
//enum ColType		{ MASTER_COL = 13, SUBPROBLEM_COL };


/**********************************************************
Let's say you want to scale a range [min,max] to [a,b]. You're looking for a(continuous) function that satisfies
f(min) = a
f(max) = b
In your case, a would be 1 and b would be 30, but let's start with something simpler and try to map [min,max] into the range [0,1].
Putting min into a function and getting out 0 could be accomplished with
f(x) = x - min == = > f(min) = min - min = 0
So that's almost what we want. But putting in max would give us max - min when we actually want 1. So we'll have to scale it :
		x - min                                  max - min
f(x) = -------- - == = > f(min) = 0;  f(max) = -------- - = 1
		max - min                                 max - min
which is what we want.So we need to do a translation and a scaling.Now if instead we want to get arbitrary values of a and b, we need 
something a little more complicated :
	(b - a)(x - min)
f(x) = -------------- + a
		max - min
You can verify that putting in min for x now gives a, and putting in max gives b.
You might also notice that(b - a) / (max - min) is a scaling factor between the size of the new range and the size of the original range.
So really we are first translating x by - min, scaling it to the correct factor, and then translating it back up to the new minimum value of a.
***********************************************************/

inline double scale(double valueIn, double baseMin, double baseMax, double limitMin, double limitMax) {
	return ((limitMax - limitMin) * (valueIn - baseMin) / (baseMax - baseMin)) + limitMin;
}
template<typename T>
bool isIn(std::vector<T> vec, T element)
{
	return std::find(vec.begin(), vec.end(), element) == vec.end();
}