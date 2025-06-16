#include "Triangle.hpp"
#include "rasterizer.hpp"
#include <eigen3/Eigen/Eigen>
#include <iostream>
#include <opencv2/opencv.hpp>

constexpr double MY_PI = 3.1415926;

//视图变换——放好相机
Eigen::Matrix4f get_view_matrix(Eigen::Vector3f eye_pos)
{
    Eigen::Matrix4f view = Eigen::Matrix4f::Identity();

    Eigen::Matrix4f translate;
    translate << 1, 0, 0, -eye_pos[0], 0, 1, 0, -eye_pos[1], 0, 0, 1,
        -eye_pos[2], 0, 0, 0, 1;

    view = translate * view;

    return view;
}

//模型变换——绕Z轴旋转
Eigen::Matrix4f get_model_matrix(float rotation_angle)
{
    Eigen::Matrix4f model = Eigen::Matrix4f::Identity();

    // TODO: Implement this function
    // Create the model matrix for rotating the triangle around the Z axis.
    // Then return it.

    float angle_pi=rotation_angle/180.0*MY_PI;
    Eigen::Matrix4f rotate;
    rotate<< cos(angle_pi),-sin(angle_pi),0.0,0.0, 
            sin(angle_pi),cos(angle_pi),0.0,0.0, 
            0.0,0.0,1.0,0.0, 
            0.0,0.0,0.0,1.0;
    model=rotate*model;

    return model;
}

//投影变换——透视投影
Eigen::Matrix4f get_projection_matrix(float eye_fov, float aspect_ratio,
                                      float zNear, float zFar)
{
    //传入的zNear和zFar为正值 但公式实际用的nf都是负值（实际坐标），所以套公式时其实是-zNear -zFar!!!
    // Students will implement this function

    Eigen::Matrix4f projection = Eigen::Matrix4f::Identity();

    // TODO: Implement this function
    // Create the projection matrix for the given parameters.
    // Then return it.

    //“挤压” 成为zNear平面长为zFar-Znear(输入的是正的)的方体
    Eigen::Matrix4f persp_to_ortho;
    persp_to_ortho<< -zNear,0.0,0.0,0.0, 0.0,-zNear,0.0,0.0, 0.0,0.0,-zNear+-zFar,-zNear*zFar, 0.0,0.0,1.0,0.0;

    Eigen::Matrix4f ortho=Eigen::Matrix4f::Identity() ,a,b;
    float z=zFar-zNear,y_half=zNear*tan(eye_fov/2.0/180.0*MY_PI),x_half=y_half*aspect_ratio;
    //将方体变为标准方体
    a<< 2.0/(x_half*2),0.0,0.0,0.0, 0.0,2.0/(y_half*2),0.0,0.0, 0.0,0.0,2.0/z,0.0, 0.0,0.0,0.0,1.0;
    //平移到以原点为中心
    b<< 1.0,0.0,0.0,0.0, 0.0,1.0,0.0,0.0, 0.0,0.0,1.0,-(-zFar+-zNear)/2.0, 0.0,0.0,0.0,1.0;
    ortho=a*b*ortho;

    projection=ortho*persp_to_ortho*projection;

    return projection;
}

//提高部分——按任意过原点轴旋转  Rodrigues’ Rotation Formula
Eigen::Matrix4f get_rotation(Vector3f axis, float angle)
{
    Eigen::Matrix4f rotation=Eigen::Matrix4f::Identity();
    Eigen::Matrix3f n,tmp=Eigen::Matrix3f::Identity();
    n<< 0,-axis[2],axis[1],
        axis[2],0,-axis[0],
        -axis[1],axis[0],0;
    float angle_pi=angle/180.0*MY_PI;
    tmp=cos(angle_pi)*tmp+(1-cos(angle_pi))*axis*axis.transpose()+sin(angle_pi)*n;
    rotation<<  tmp(0,0),tmp(0,1),tmp(0,2),0,
                tmp(1,0),tmp(1,1),tmp(1,2),0,
                tmp(2,0),tmp(2,1),tmp(2,2),0,
                0,0,0,1;
    return rotation;
}

int main(int argc, const char** argv)
{
    float angle = 0;
    bool command_line = false;
    std::string filename = "output.png";

    if (argc >= 3) {//参数多于三个，则是静态模式保存图片 1 ./rasterizer 2 -r 3 角度 4 名称
        command_line = true;
        angle = std::stof(argv[2]); // -r by default  stof()将字符串转换为float型
        if (argc == 4) {
            filename = std::string(argv[3]);
        }
        // else
            // return 0;
    }

    rst::rasterizer r(700, 700);//确定宽高参数的屏幕(定义一个光栅器)

    Eigen::Vector3f eye_pos = {0, 0, 5};//人眼（相机）坐标位置
    
    //提高——绕axis选择
    Eigen::Vector3f axis={0,0,1};

    std::vector<Eigen::Vector3f> pos{{2, 0, -2}, {0, 2, -2}, {-2, 0, -2}};//三角形三个顶点坐标

    std::vector<Eigen::Vector3i> ind{{0, 1, 2}};//三个点的索引值

    auto pos_id = r.load_positions(pos);
    auto ind_id = r.load_indices(ind);

    int key = 0;
    int frame_count = 0;

    if (command_line) {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);//清除整个显示屏

        // r.set_model(get_model_matrix(angle));
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

        //分别将main.cpp中用三个函数返回的m,v,p三个矩阵加载到光栅化器（rasterizer）中的model,view,projection这三个矩阵中去
        //三个变换矩阵也就是模型变换矩阵，视图变换矩阵，投影变换矩阵
        
        // r.set_model(get_model_matrix(angle));
        
        r.set_model(get_rotation(axis,angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));
        
        //绘制三角形
        r.draw(pos_id, ind_id, rst::Primitive::Triangle);

        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());//mat是opencv中存储图像的容器
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::imshow("image", image);
        key = cv::waitKey(10);//获取键盘内容

        std::cout << "frame count: " << frame_count++ << '\n';

        if (key == 'a') {
            angle += 10;
        }
        else if (key == 'd') {
            angle -= 10;
        }
    }

    return 0;
}
