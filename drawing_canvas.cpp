#include "drawing_canvas.h"
#include <FL/fl_draw.H>
#include <cmath>
#include <iostream>

drawing_canvas::drawing_canvas(int x, int y, int w, int h) : 
    Fl_Box(x, y, w, h),
    current_color(FL_RED),
    current_shape(RECTANGLE), // 默认圆形
    is_dragging(false),
    drag_index(-1),
    last_x(0),
    last_y(0)
{
    color(FL_BLACK);
    box(FL_FLAT_BOX);
}

void drawing_canvas::clear_points() {
    points.clear();
    redraw();
}

void drawing_canvas::print_points() const {
    std::cout << "Total Points: " << points.size() << "\n";
    for (const auto &p : points) {
        std::cout << "Point(" << p.x << ", " << p.y << ") - Color: ";
        
        if (p.color == FL_WHITE) std::cout << "WHITE";
        else if (p.color == FL_RED) std::cout << "RED";
        else if (p.color == FL_GREEN) std::cout << "GREEN";
        else if (p.color == FL_BLUE) std::cout << "BLUE";
        else if (p.color == FL_YELLOW) std::cout << "YELLOW";
        else if (p.color == FL_CYAN) std::cout << "CYAN";
        else if (p.color == FL_MAGENTA) std::cout << "MAGENTA";
        else {
            unsigned char r, g, b;
            Fl::get_color(p.color, r, g, b);
            std::cout << "RGB(" << static_cast<int>(r) << "," 
                      << static_cast<int>(g) << "," 
                      << static_cast<int>(b) << ")";
        }
        
        // 添加形状信息输出
        std::cout << " - Shape: ";
        switch (p.shape) {
            case CIRCLE: std::cout << "CIRCLE"; break;
            case RECTANGLE: std::cout << "RECTANGLE"; break;
            case TRIANGLE: std::cout << "TRIANGLE"; break;
        }
        
        std::cout << "\n";
    }
}

void drawing_canvas::set_color(Fl_Color c) {
    current_color = c;
}

// 设置形状方法
void drawing_canvas::set_shape(ShapeType s) {
    current_shape = s;
}

int drawing_canvas::find_point_at(int x, int y) {
    for (int i = 0; i < points.size(); i++) {
        const point& p = points[i];
        // 计算鼠标位置与点的距离
        int dx = x - p.x;
        int dy = y - p.y;
        int distance_squared = dx * dx + dy * dy;
        
        // 如果距离在阈值内，则认为选中该点
        if (distance_squared <= selection_radius * selection_radius) {
            return i;
        }
    }
    return -1; // 没有找到点
}

void drawing_canvas::draw() {
    // 绘制黑色背景
    fl_color(FL_BLACK);
    fl_rectf(x(), y(), w(), h());
    
    // 绘制所有点
    for (const auto &p : points) {
        fl_color(p.color);
        
        // 根据形状类型绘制不同的图形
        switch (p.shape) {
            case CIRCLE: {
                fl_begin_polygon();
                for (int i = 0; i < 360; i += 10) {
                    double angle = i * M_PI / 180.0;
                    fl_vertex(
                        p.x + point_size * cos(angle),
                        p.y + point_size * sin(angle)
                    );
                }
                fl_end_polygon();
                break;
            }
            
            case RECTANGLE: {
                fl_rectf(p.x - point_size, p.y - point_size, 
                         point_size * 2, point_size * 2);
                break;
            }
            
            case TRIANGLE: {
                fl_begin_polygon();
                fl_vertex(p.x, p.y - point_size); // 上顶点
                fl_vertex(p.x - point_size, p.y + point_size); // 左下
                fl_vertex(p.x + point_size, p.y + point_size); // 右下
                fl_end_polygon();
                break;
            }
        }
    }
}

int drawing_canvas::handle(int event) { 
    int canvas_x = Fl::event_x(); 
    int canvas_y = Fl::event_y(); 

    bool in_canvas = (canvas_x >= this->x() && canvas_x < this->x() + w() && 
                      canvas_y >= this->y() && canvas_y < this->y() + h()); 

    switch (event) { 
        case FL_PUSH:
            if (Fl::event_button() == FL_LEFT_MOUSE && in_canvas) {
                // 首先检查是否点击了已有的点
                int index = find_point_at(canvas_x, canvas_y);
                if (index != -1) {
                    // 开始拖动点
                    is_dragging = true;
                    drag_index = index;
                    last_x = canvas_x;
                    last_y = canvas_y;
                    return 1;
                } else {
                    // 没有点到已有图形，记录位置准备添加新点
                    last_x = canvas_x;
                    last_y = canvas_y;
                    return 1;
                }
            }
            else if (Fl::event_button() == FL_RIGHT_MOUSE && in_canvas && !points.empty()) {
                // 删除点：首先尝试删除鼠标下的点
                int index = find_point_at(canvas_x, canvas_y);
                if (index != -1) {
                    points.erase(points.begin() + index);
                    redraw();
                    return 1;
                }
                // 如果没点到点，删除最后一个点
                points.pop_back();
                redraw();
                return 1;
            }
            break;

        case FL_DRAG:
            if (is_dragging && drag_index >= 0 && drag_index < points.size()) {
                // 更新被拖动的点的位置
                points[drag_index].x += Fl::event_x() - last_x;
                points[drag_index].y += Fl::event_y() - last_y;
                last_x = Fl::event_x();
                last_y = Fl::event_y();
                redraw();
                return 1;
            }
            break;

        case FL_RELEASE:
            if (is_dragging) {
                // 结束拖动
                is_dragging = false;
                drag_index = -1;
                return 1;
            } else if (Fl::event_button() == FL_LEFT_MOUSE && in_canvas) {
                // 添加新点
                points.push_back({last_x, last_y, current_color, current_shape});
                redraw();
                return 1;
            }
            break;
    }

    return Fl_Box::handle(event); 
}