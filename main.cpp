#include <iostream>
#include <fstream>
#include <vector>

#include "Eigen/Dense"
#include "nlohmann/json.hpp"

#include "ply_stream.hpp"

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

Eigen::Vector4f EstimateGroundPlane(
    const std::vector<Eigen::Vector3f> &points
){
    // centroid 
    Eigen::Vector3f centroid(0.0f, 0.0f, 0.0f);
    for (const auto &pt : points) {
        centroid += pt; // running sum
    }
    centroid /= points.size(); //avergae by num samples

    // covar matrix
    Eigen::Matrix3f covariance_matrix = Eigen::Matrix3f::Zero();
    for (const auto &pt : points) {
        Eigen::Vector3f distance = pt - centroid;
        covariance_matrix += distance * distance.transpose(); // gete cov matrix of distance from cent
    }
    covariance_matrix /= points.size(); // avg

    // PCA
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigen_solver(covariance_matrix, Eigen::ComputeEigenvectors);
    Eigen::Vector3f normal = eigen_solver.eigenvectors().col(2).normalized();

    // plane eq coeff (ax + by + cz + d = 0)
    float d = -normal.dot(centroid);

    return Eigen::Vector4f(normal.x(), normal.y(), normal.z(), d);

    /*
        Ultimately I went with using PCA to create a flat "ground plane", if I had more time I might of used another approach that would have allowed for 
        a less uniform ground plane. In my reading I found out about RANSAC, which seemed promising since it allowed for more customization of the parameters used to 
        create a plane. Meaning the ground could have represented a less "flat" surface. I think this would have been my approach if there were more variance along the "z"
        direction. For example, the UAS covers a small area, but enclosed in the area maybe there are skyrises or other tall structures, which could throw PCA off. 
    */
}

// backproject "this" image corners onto the ground plane
std::vector<Eigen::Vector2f> BackprojectCorners(const Intrinsic &intrinsic, const Extrinsic &extrinsic, const Eigen::Vector4f &ground_plane) {
    // corners
    std::vector<Eigen::Vector2f> corners = {Eigen::Vector2f(0, 0), Eigen::Vector2f(intrinsic.cols - 1, 0),
                                            Eigen::Vector2f(0, intrinsic.rows - 1), Eigen::Vector2f(intrinsic.cols - 1, intrinsic.rows - 1)};
    
    // matrix to go from 2D (pixel) to 3D (camera) domains
    Eigen::Matrix3f K_inv = intrinsic.K.inverse();
    Eigen::Matrix<float, 3, 4> P;
    P.block<3, 3>(0, 0) = extrinsic.R.transpose();
    P.col(3) = -extrinsic.R.transpose() * extrinsic.c;

    // translate image corners to the world's ground plane
    std::vector<Eigen::Vector2f> projected_corners;
    for (const auto &corner : corners) {
        // norm pixel coor
        Eigen::Vector3f pixel_homogeneous(corner.x(), corner.y(), 1.0f);
        Eigen::Vector3f ray_direction = K_inv * pixel_homogeneous;
        ray_direction.normalize();

        // intersection w/ ground 
        float t = -(ground_plane.head<3>().dot(Eigen::Vector3f(extrinsic.c)) + ground_plane[3]) / ground_plane.head<3>().dot(ray_direction);
        Eigen::Vector3f intersection = extrinsic.c + t * ray_direction;

        // project intersection to world coordinates
        Eigen::Vector3f intersection_homogeneous = P * Eigen::Vector4f(intersection.x(), intersection.y(), intersection.z(), 1.0f);
        Eigen::Vector2f projected_point = intersection_homogeneous.head<2>() / intersection_homogeneous[2];

        projected_corners.push_back(projected_point);
    }

    return projected_corners;
}

// write view coverage polygons to a ply file
void WriteViewCoveragePolygons(const std::vector<Extrinsic> &extrinsics, const Intrinsic &intrinsic, const std::vector<Eigen::Vector3f> &points, const Eigen::Vector4f &ground_plane, const std::string &output_folder) {
    PointPlyStream ply_stream(output_folder + "view_coverage_polygons.ply");
    ply_stream.WriteHeader({"float x", "float y", "float z", "uchar red", "uchar green", "uchar blue"});

    // og point cloud
    for (const auto &pt : points) {
        ply_stream << pt.x() << pt.y() << pt.z() << (unsigned char)255 << (unsigned char)0 << (unsigned char)0; // red color for points
    }

    // camera positions (based on ext parameters)
    for (const auto &extrinsic : extrinsics) {
        ply_stream << extrinsic.c.x() << extrinsic.c.y() << extrinsic.c.z() << (unsigned char)0 << (unsigned char)0 << (unsigned char)255; // blue color for cameras
    }
	float count = 0.0f;
    // view coverage polygons based on backprojection
    for (const auto &extrinsic : extrinsics) {
        std::vector<Eigen::Vector2f> projected_corners = BackprojectCorners(intrinsic, extrinsic, ground_plane);
		
        // 4 coreners of each "view" polygon
        for (int i = 0; i < 4; ++i) {
            ply_stream << projected_corners[i].x() << projected_corners[i].y() << count << (unsigned char)0 << (unsigned char)255 << (unsigned char)0; // green color for polygon edges
            ply_stream << projected_corners[(i + 1) % 4].x() << projected_corners[(i + 1) % 4].y() << count << (unsigned char)0 << (unsigned char)255 << (unsigned char)0;
        }
		count -= 10.0f;
    }

    /*
        I i would have had more time I would have liked to acutally line up the polygons representing the FOV with where the camera was located at the time. Currently I have them stacked 
        along the "z" direction (with a spacing of -10f) as that is how I interepreted the sentence:
        "You can write the view coverage polygons by constructing each edge of points evenly spaced along the line."
        Although this method does not seem super useful for visualizing what the corners actually correspond to in the "world" domain.

        Furthermore, the scaling of the polygons "feels" wrong. 
        It seems to me that the "FOV cone" from the cameras to the ground should be a lot larger, so there may be and error in my code -_-
    */
}

int main(int, char**){

    // read in the scene data:
    const std::string json_file_path = "/workspaces/CV_takehome_challenge/scene_data.json";
    Intrinsic intrinsic;
    std::vector<Extrinsic> extrinsics;
    std::vector<Eigen::Vector3f> points;
    Eigen::Vector4f ground_plane;

    ParseSceneDataJson(intrinsic, extrinsics, points, json_file_path);


    // write results to this folder:
    const std::string output_folder = "/workspaces/CV_takehome_challenge/";

    // // here's an example of how to use the PointPlyStream class (ply_stream.hpp) to write the point cloud to a .ply file:
    // PointPlyStream ply_stream(output_folder + "pt_cloud.ply");
    // ply_stream.WriteHeader({"float x", "float y", "float z", "uchar red", "uchar green", "uchar blue"});
    // for (const auto &pt : points) {
    //     ply_stream << pt.x() << pt.y() << pt.z() << (unsigned char)0 << (unsigned char)255 << (unsigned char)255;
    // }

    // Brian's Code
    // Estimate Ground plane
    ground_plane = EstimateGroundPlane(points);
    // Write view coverage polygons to ply file
    WriteViewCoveragePolygons(extrinsics, intrinsic, points, ground_plane, output_folder);
    /*
        If I had the opportunity to work on this more I would probably have also thrown the functions in another .cpp/.hpp file, rather than jamming them above. 
    */

}
