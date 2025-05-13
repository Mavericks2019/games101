#include "Triangle.hpp"
#include "rasterizer.hpp"
#include <eigen3/Eigen/Eigen>
#include <iostream>
#include <opencv2/opencv.hpp>

constexpr double MY_PI = 3.1415926;

Eigen::Matrix4f get_view_matrix(Eigen::Vector3f eye_pos)
{
    Eigen::Matrix4f view = Eigen::Matrix4f::Identity();

    Eigen::Matrix4f translate;
    translate << 1, 0, 0, -eye_pos[0],
                 0, 1, 0, -eye_pos[1],
                 0, 0, 1,-eye_pos[2], 
                 0, 0, 0, 1;

    view = translate * view;

    return view;
}

Eigen::Matrix4f get_model_matrix(float rotation_angle)
{
    Eigen::Matrix4f model = Eigen::Matrix4f::Identity();


    float a=rotation_angle / 180.0 * MY_PI;
    Eigen::Matrix4f rotate;
    rotate << cos(a), -sin(a), 0.0, 0.0, 
            sin(a), cos(a), 0.0, 0.0, 
            0.0, 0.0, 1.0, 0.0, 
            0.0, 0.0, 0.0, 1.0;
    model = rotate * model;

    return model;
}

Eigen::Matrix4f get_projection_matrix(float eye_fov, float aspect_ratio,
                                      float zNear, float zFar)
{

    Eigen::Matrix4f projection = Eigen::Matrix4f::Identity();
    Eigen::Matrix4f persp_to_ortho;
    persp_to_ortho<< -zNear, 0.0, 0.0, 0.0, 
                     0.0, -zNear, 0.0, 0.0, 
                     0.0, 0.0, -zNear + -zFar, -zNear * zFar, 
                     0.0, 0.0, 1.0, 0.0;
    Eigen::Matrix4f ortho = Eigen::Matrix4f::Identity() , a, b;
    float z = zFar - zNear;
    float y_half = zNear * tan(eye_fov / 2.0 / 180.0 * MY_PI);
    float x_half = y_half * aspect_ratio;
    //标准化
    a << 2.0/(x_half*2), 0.0, 0.0, 0.0,
        0.0, 2.0/(y_half*2), 0.0, 0.0, 
        0.0, 0.0, 2.0/z, 0.0, 
        0.0, 0.0, 0.0, 1.0;
    //平移方块中心到原点
    b << 1.0, 0.0, 0.0, 0.0,
         0.0, 1.0, 0.0, 0.0,
         0.0, 0.0, 1.0, -(-zFar+-zNear)/2.0,
         0.0, 0.0, 0.0, 1.0;

    ortho = a * b * ortho;

    projection = ortho * persp_to_ortho * projection;

    return projection;
}

Eigen::Matrix4f get_rotation(Vector3f axis, float angle)
{
    Eigen::Matrix4f rotation=Eigen::Matrix4f::Identity();
    Eigen::Matrix3f n = Eigen::Matrix3f::Identity();
    Eigen::Matrix3f tmp = Eigen::Matrix3f::Identity();
    n << 0, -axis[2], axis[1],
         axis[2], 0, -axis[0],
         -axis[1], axis[0], 0;
    float angle_pi = angle / 180.0 * MY_PI;
    tmp = cos(angle_pi) * tmp + (1 - cos(angle_pi)) * axis * axis.transpose() + sin(angle_pi) * n;
    rotation<<  tmp(0,0), tmp(0,1), tmp(0,2), 0,
                tmp(1,0), tmp(1,1), tmp(1,2), 0,
                tmp(2,0), tmp(2,1), tmp(2,2), 0,
                0, 0, 0, 1;
    return rotation;
}

int main(int argc, const char** argv)
{
    float angle = 0;
    bool command_line = false;
    std::string filename = "output.png";

    if (argc >= 3) {
        command_line = true;
        angle = std::stof(argv[2]);
        if (argc == 4) {
            filename = std::string(argv[3]);
        }
    }

    rst::rasterizer r(700, 700);

    Eigen::Vector3f eye_pos = {0, 0, 5};
    
    Eigen::Vector3f axis={0,0,1};

    std::vector<Eigen::Vector3f> pos{{2, 0, -2}, {0, 2, -2}, {-2, 0, -2}};
    std::vector<Eigen::Vector3i> ind{{0, 1, 2}};

    auto pos_id = r.load_positions(pos);
    auto ind_id = r.load_indices(ind);

    int key = 0;
    int frame_count = 0;

    if (command_line) {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_rotation(axis,angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));

        r.draw(pos_id, ind_id, rst::Primitive::Triangle);
        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);

        cv::imwrite(filename, image);

        return 0;
    }

    while (key != 27) {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);
                
        r.set_model(get_rotation(axis,angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(90, 1, 0.1, 50));
        
        r.draw(pos_id, ind_id, rst::Primitive::Triangle);

        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::imshow("image", image);
        key = cv::waitKey(10);

        std::cout << "frame count: " << frame_count++ << '\n';
        std::cout << key << '\n';

        if (key == 'a') {
            std::cout << "aa" << '\n';
            angle += 10;
        }
        else if (key == 'd') {
            angle -= 10;
        }
        else if(key == 'w') {
            axis = {1, 1, 1};
            angle += 10;
        }
        else if(key == 's') {
            axis = {0, 0, 1};
        }
    }

    return 0;
}
