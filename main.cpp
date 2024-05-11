#include <iostream>
#include <fstream>
#include <vector>

#include "Eigen/Dense"
#include "nlohmann/json.hpp"

// this is a pinhole camera model, with zero distortion or skew
struct Intrinsic {
    Eigen::Matrix3f K; // intrinsic matrix K
    int rows; // height of the image
    int cols; // width of the image
};

Intrinsic ExtractIntrinsicFromJson(const nlohmann::json &intrinsic_json) {
    
    Intrinsic intrinsic;

    intrinsic.K << 
        intrinsic_json["value"]["focal_length"], 0, intrinsic_json["value"]["principal_point"][0],
        0, intrinsic_json["value"]["focal_length"], intrinsic_json["value"]["principal_point"][1],
        0, 0, 1;

    intrinsic.rows = intrinsic_json["value"]["height"];
    intrinsic.cols = intrinsic_json["value"]["width"];

    return intrinsic;
}

struct Extrinsic {
    Eigen::Matrix3f R; // rotation matrix
    Eigen::Vector3f c; // camera center
};

Extrinsic ExtractExtrinsicFromJson(const nlohmann::json &extrinsic_json) {
    Extrinsic extrinsic;    
    for (int row = 0; row < 3; row++) {
        for (int col = 0; col < 3; col++) {
            extrinsic.R(row,col) = extrinsic_json["value"]["rotation"][row][col];
        }
    }
    for (int i = 0; i < 3; i++) {
        extrinsic.c(i) = extrinsic_json["value"]["center"][i];
    }
    return extrinsic;
}

Eigen::Vector3f ExtractPointFromJson(const nlohmann::json &point_json) {
    Eigen::Vector3f point;
    for (int i = 0; i < 3; i++) {
        point(i) = point_json["value"]["X"][i];
    }
    return point;
}

void ParseSceneDataJson(
    Intrinsic &intrinsic, // output
    std::vector<Extrinsic> &extrinsics, // output
    std::vector<Eigen::Vector3f> &points, // output
    const std::string &json_file_path // input
) {
    std::ifstream f(json_file_path);
    nlohmann::json scene_data = nlohmann::json::parse(f);

    // extract camera intrinsics.
    intrinsic = ExtractIntrinsicFromJson(scene_data["intrinsic"]);
    std::cout << "intrinsic K:\n" << intrinsic.K << "\n";
    std::cout << "image size: " << intrinsic.rows << ", " << intrinsic.cols << "\n";

    // extract the camera extrinsics for each view:
    for (const auto &extrinsic_json : scene_data["extrinsics"]) {
        extrinsics.emplace_back(ExtractExtrinsicFromJson(extrinsic_json));
    }
    std::cout << "extracted " << extrinsics.size() << " extrinsics\n";

    // extract the point cloud:
    for (const auto &point_json : scene_data["structure"]) {
        points.emplace_back(ExtractPointFromJson(point_json));
    }
    std::cout << "extracted " << points.size() << " points\n";
}

int main(int, char**){

    // read in the scene data:
    const std::string json_file_path = "/workspaces/CV_takehome_challenge/scene_data.json";
    Intrinsic intrinsic;
    std::vector<Extrinsic> extrinsics;
    std::vector<Eigen::Vector3f> points;
    ParseSceneDataJson(intrinsic, extrinsics, points, json_file_path);

    // TODO: do stuff here


}
