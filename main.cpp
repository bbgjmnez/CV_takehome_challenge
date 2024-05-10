#include <iostream>

#include "Eigen/Dense"

int main(int, char**){
    std::cout << "Hello, from cv_interview_challenge!\n";

    Eigen::Vector3f vec(0,1,3);
    std::cout << "vec: " << vec.transpose() << "\n";

}
