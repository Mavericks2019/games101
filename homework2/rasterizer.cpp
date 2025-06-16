// clang-format off
//
// Created by goksu on 4/6/19.
//

#include <algorithm>
#include <vector>
#include "rasterizer.hpp"
#include <opencv2/opencv.hpp>
#include <math.h>
using namespace std;

rst::pos_buf_id rst::rasterizer::load_positions(const std::vector<Eigen::Vector3f> &positions)
{
    auto id = get_next_id();
    pos_buf.emplace(id, positions);

    return {id};
}

rst::ind_buf_id rst::rasterizer::load_indices(const std::vector<Eigen::Vector3i> &indices)
{
    auto id = get_next_id();
    ind_buf.emplace(id, indices);

    return {id};
}

rst::col_buf_id rst::rasterizer::load_colors(const std::vector<Eigen::Vector3f> &cols)
{
    auto id = get_next_id();
    col_buf.emplace(id, cols);

    return {id};
}

auto to_vec4(const Eigen::Vector3f& v3, float w = 1.0f)
{
    return Vector4f(v3.x(), v3.y(), v3.z(), w);
}


static bool insideTriangle(float x, float y, const Vector3f* _v)
{   
    Eigen::Vector3f point(x,y,1.0f);
    
    Eigen::Vector3f a,b,c;
    a=_v[1]-_v[0];b=_v[2]-_v[1];c=_v[0]-_v[2];

    Eigen::Vector3f p0,p1,p2;
    p0=point-_v[0];p1=point-_v[1];p2=point-_v[2];
    if((a.cross(p0)).z()>0&&(b.cross(p1)).z()>0&&(c.cross(p2)).z()>0 || (a.cross(p0)).z()<0&&(b.cross(p1)).z()<0&&(c.cross(p2)).z()<0)return true;
    return false; 
}

static std::tuple<float, float, float> computeBarycentric2D(float x, float y, const Vector3f* v)
{
    float c1 = (x*(v[1].y() - v[2].y()) + (v[2].x() - v[1].x())*y + v[1].x()*v[2].y() - v[2].x()*v[1].y()) / (v[0].x()*(v[1].y() - v[2].y()) + (v[2].x() - v[1].x())*v[0].y() + v[1].x()*v[2].y() - v[2].x()*v[1].y());
    float c2 = (x*(v[2].y() - v[0].y()) + (v[0].x() - v[2].x())*y + v[2].x()*v[0].y() - v[0].x()*v[2].y()) / (v[1].x()*(v[2].y() - v[0].y()) + (v[0].x() - v[2].x())*v[1].y() + v[2].x()*v[0].y() - v[0].x()*v[2].y());
    float c3 = (x*(v[0].y() - v[1].y()) + (v[1].x() - v[0].x())*y + v[0].x()*v[1].y() - v[1].x()*v[0].y()) / (v[2].x()*(v[0].y() - v[1].y()) + (v[1].x() - v[0].x())*v[2].y() + v[0].x()*v[1].y() - v[1].x()*v[0].y());
    return {c1, c2, c3};
}

void rst::rasterizer::draw(pos_buf_id pos_buffer, ind_buf_id ind_buffer, col_buf_id col_buffer, Primitive type)
{
    auto& buf = pos_buf[pos_buffer.pos_id];
    auto& ind = ind_buf[ind_buffer.ind_id];
    auto& col = col_buf[col_buffer.col_id];

    float f1 = (50 - 0.1) / 2.0;
    float f2 = (50 + 0.1) / 2.0;

    Eigen::Matrix4f mvp = projection * view * model;
    for (auto& i : ind)
    {
        Triangle t;
        Eigen::Vector4f v[] = {
                mvp * to_vec4(buf[i[0]], 1.0f),
                mvp * to_vec4(buf[i[1]], 1.0f),
                mvp * to_vec4(buf[i[2]], 1.0f)
        };
        //Homogeneous division
        for (auto& vec : v) {
            vec /= vec.w();
        }
        //Viewport transformation
        for (auto & vert : v)
        {
            vert.x() = 0.5*width*(vert.x()+1.0);
            vert.y() = 0.5*height*(vert.y()+1.0);
            // vert.z() = vert.z() * f1 + f2; 
            vert.z() = -vert.z() * f1 + f2;
        }

        for (int i = 0; i < 3; ++i)
        {
            t.setVertex(i, v[i].head<3>());
            t.setVertex(i, v[i].head<3>());
            t.setVertex(i, v[i].head<3>());
        }

        auto col_x = col[i[0]];
        auto col_y = col[i[1]];
        auto col_z = col[i[2]];

        t.setColor(0, col_x[0], col_x[1], col_x[2]);
        t.setColor(1, col_y[0], col_y[1], col_y[2]);
        t.setColor(2, col_z[0], col_z[1], col_z[2]);

        rasterize_triangle(t);
    }
}

//Screen space rasterization
void rst::rasterizer::rasterize_triangle(const Triangle& t) {
    auto v = t.toVector4();
    std::array<Eigen::Vector2d, 2> bounding_box;
    float tmp_min[2], tmp_max[2];
    tmp_min[0] = std::min(v[0].x(), std::min(v[1].x(), v[2].x()));
    tmp_min[1] = std::min(v[0].y(), std::min(v[1].y(), v[2].y()));
    tmp_max[0] = std::max(v[0].x(), std::max(v[1].x(), v[2].x()));
    tmp_max[1] = std::max(v[0].y(), std::max(v[1].y(), v[2].y()));
    bounding_box[0] << floor(tmp_min[0]), floor(tmp_min[1]);
    bounding_box[1] << ceil(tmp_max[0]), ceil(tmp_max[1]);

    for(int x = bounding_box[0].x(); x < bounding_box[1].x(); x++)
    {
        for(int y = bounding_box[0].y(); y < bounding_box[1].y(); y++)
        {
            vector<vector<float>> tmp = {{0.25, 0.25}, {0.25, 0.75}, {0.75, 0.25}, {0.75, 0.75}};
            float z_pixel = 0x3f3f3f3f;
            int pix_index = get_index(x, y) * 4;
            int flag = 4;
            for(int i = 0; i < 4; i++) {
                float x_t = x + tmp[i][0];
                float y_t = y + tmp[i][1];
                if(insideTriangle(x_t, y_t, t.v)) {
                    flag--;
                    auto[alpha, beta, gamma] = computeBarycentric2D(x_t, y_t, t.v);
                    float w_reciprocal = 1.0/(alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
                    float z_interpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
                    z_interpolated *= w_reciprocal;
                    if(z_interpolated < depth_samp_buf[i + pix_index])
                    {
                        depth_samp_buf[i + pix_index] = z_interpolated;
                        color_samp_buf[i + pix_index] = t.getColor() / 4;
                    }
                    z_pixel = min(z_pixel, depth_samp_buf[i + pix_index]);
                }
            }
            if(flag != 4 && flag != 0) {
                std::cout << flag << std::endl;
            }
            auto color = color_samp_buf[pix_index] + color_samp_buf[pix_index + 1] + color_samp_buf[pix_index + 2] +  color_samp_buf[pix_index + 3];
            depth_buf[get_index(x, y)] = z_pixel;
            set_pixel(Eigen::Vector3f(x, y, z_pixel), color);
        }
    }

}

void rst::rasterizer::set_model(const Eigen::Matrix4f& m)
{
    model = m;
}

void rst::rasterizer::set_view(const Eigen::Matrix4f& v)
{
    view = v;
}

void rst::rasterizer::set_projection(const Eigen::Matrix4f& p)
{
    projection = p;
}

void rst::rasterizer::clear(rst::Buffers buff)
{
    if ((buff & rst::Buffers::Color) == rst::Buffers::Color)
    {
        std::fill(frame_buf.begin(), frame_buf.end(), Eigen::Vector3f{0, 0, 0});
    }
    if ((buff & rst::Buffers::Color) == rst::Buffers::Color)
    {
        std::fill(color_samp_buf.begin(), color_samp_buf.end(), Eigen::Vector3f{0, 0, 0});
    }
    if ((buff & rst::Buffers::Depth) == rst::Buffers::Depth)
    {
        std::fill(depth_buf.begin(), depth_buf.end(), std::numeric_limits<float>::infinity());
    }
    if ((buff & rst::Buffers::Depth) == rst::Buffers::Depth)
    {
        std::fill(depth_samp_buf.begin(), depth_samp_buf.end(), std::numeric_limits<float>::infinity());
    }
}

rst::rasterizer::rasterizer(int w, int h) : width(w), height(h)
{
    frame_buf.resize(w * h);
    depth_buf.resize(w * h);
    color_samp_buf.resize(w * h * 4);
    depth_samp_buf.resize(w * h * 4);
}

int rst::rasterizer::get_index(int x, int y)
{
    return (height-1-y)*width + x;
}

void rst::rasterizer::set_pixel(const Eigen::Vector3f& point, const Eigen::Vector3f& color)
{
    //old index: auto ind = point.y() + point.x() * width;
    auto ind = (height-1-point.y())*width + point.x();
    frame_buf[ind] = color;

}

// clang-format on