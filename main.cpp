#include <iostream>

#include "Eigen/Dense"
#include "nlohmann/json.hpp"

int main(int, char**){
    std::cout << "Hello, from cv_interview_challenge!\n";

    Eigen::Vector3f vec(0,1,3);
    std::cout << "vec: " << vec.transpose() << "\n";

    nlohmann::json json;
    json["mood"] = "happy";
    json["number"] = 6.3;
    std::cout << "json:\n" << json.dump(4) << "\n";

}
