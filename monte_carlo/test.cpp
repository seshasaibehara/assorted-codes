#include <iostream> 
#include <cmath>
#include "headers.h"
#include "../eigen-git-mirror/Eigen/Core"
#include <vector>
#include <numeric>
int main(){
int i=0;
int mp = 10;
std::vector<double> test;
test.push_back(0.33);
test.push_back(0.33);
test.push_back(0.33);
test.push_back(0.33);
test.push_back(0.33);
std::cout << std::accumulate(test.begin(),test.end(),0) << std::endl;
}


