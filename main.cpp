#include "drawing_canvas.h"
#include <FL/Fl_Window.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Round_Button.H>
#include <FL/Fl_Group.H>
#include <iostream>
#include <chrono>
#include "solver.h"
using namespace Eigen;

// 安全传递形状类型的辅助结构
struct ShapeData {
    drawing_canvas* canvas;
    ShapeType shape;
};

// 前向声明回调函数
static void clear_callback(Fl_Widget *widget, void *data);
static void print_callback(Fl_Widget *widget, void *data);
static void color_callback(Fl_Widget *widget, void *data);
static void shape_callback(Fl_Widget *widget, void *data);
static void new_button_callback(Fl_Widget *widget, void *data); 
static void poly_button_callback(Fl_Widget *widget, void *data); 

int main() {
    // 创建更大的窗口以容纳新按钮
    Fl_Window *window = new Fl_Window(1100, 810, "Enhanced Drawing Canvas");
    
    // 设置窗口不可调整大小和最大化
    window->resizable(nullptr); // 禁止调整大小
    window->size_range(1100, 810, 1100, 810); // 固定窗口大小
    
    // 创建画布（占据窗口左侧）
    drawing_canvas *canvas = new drawing_canvas(0, 0, 810, 810);
    
    // 工具栏（占据窗口右侧）
    Fl_Box *toolbar = new Fl_Box(810, 0, 300, 810);
    toolbar->box(FL_ENGRAVED_BOX);
    toolbar->color(fl_rgb_color(240, 240, 240));
    
    // ================== 基本功能按钮 ==================
    int y_pos = 30;
    const int btn_height = 40;
    const int btn_width = 260;
    const int btn_spacing = 10;
    
    Fl_Button *clear_btn = new Fl_Button(820, y_pos, btn_width, btn_height, "Clear Canvas");
    clear_btn->callback(clear_callback, canvas);
    y_pos += btn_height + btn_spacing;
    
    Fl_Button *print_btn = new Fl_Button(820, y_pos, btn_width, btn_height, "Print Points");
    print_btn->callback(print_callback, canvas);
    y_pos += btn_height + btn_spacing;
    
    // ================== 新增功能按钮 ==================
    Fl_Button *poly_btn = new Fl_Button(820, y_pos, btn_width, btn_height, "Ploynomial Fitting");
    poly_btn->callback(poly_button_callback, canvas);
    y_pos += btn_height + btn_spacing;
    
    Fl_Button *load_btn = new Fl_Button(820, y_pos, btn_width, btn_height, "Gaussian Fitting");
    load_btn->callback(new_button_callback, (void*)"Load");
    y_pos += btn_height + btn_spacing;
    
    Fl_Button *undo_btn = new Fl_Button(820, y_pos, btn_width, btn_height, "Ploynomial Regression");
    undo_btn->callback(new_button_callback, (void*)"Undo");
    y_pos += btn_height + btn_spacing;
    
    Fl_Button *size_btn = new Fl_Button(820, y_pos, btn_width, btn_height, "Ridge Regression");
    size_btn->callback(new_button_callback, (void*)"Size");
    y_pos += btn_height + btn_spacing;
    
    Fl_Button *bg_btn = new Fl_Button(820, y_pos, btn_width, btn_height, "All");
    bg_btn->callback(new_button_callback, (void*)"Background");
    y_pos += btn_height + btn_spacing * 2;
    
    // ================== 点大小信息 ==================
    Fl_Box *size_info = new Fl_Box(820, y_pos, btn_width, 30, "Point Size: 5px");
    size_info->labelsize(16);
    y_pos += 40;
    
    // ================== 颜色选择组 ==================
    Fl_Group *color_group = new Fl_Group(820, y_pos, btn_width, 180, "Drawing Color");
    color_group->box(FL_ENGRAVED_FRAME);
    color_group->align(FL_ALIGN_TOP | FL_ALIGN_INSIDE);
    
    int color_y = y_pos + 30;
    const int color_btn_height = 30;
    
    Fl_Round_Button *white_btn = new Fl_Round_Button(830, color_y, btn_width - 20, color_btn_height, "White");
    white_btn->type(FL_RADIO_BUTTON);
    white_btn->value(1);
    white_btn->color(FL_WHITE);
    white_btn->selection_color(FL_WHITE);
    white_btn->callback(color_callback, canvas);
    color_y += color_btn_height;
    
    Fl_Round_Button *red_btn = new Fl_Round_Button(830, color_y, btn_width - 20, color_btn_height, "Red");
    red_btn->type(FL_RADIO_BUTTON);
    red_btn->color(FL_RED);
    red_btn->selection_color(FL_RED);
    red_btn->callback(color_callback, canvas);
    color_y += color_btn_height;
    
    Fl_Round_Button *green_btn = new Fl_Round_Button(830, color_y, btn_width - 20, color_btn_height, "Green");
    green_btn->type(FL_RADIO_BUTTON);
    green_btn->color(FL_GREEN);
    green_btn->selection_color(FL_GREEN);
    green_btn->callback(color_callback, canvas);
    color_y += color_btn_height;
    
    Fl_Round_Button *blue_btn = new Fl_Round_Button(830, color_y, btn_width - 20, color_btn_height, "Blue");
    blue_btn->type(FL_RADIO_BUTTON);
    blue_btn->color(FL_BLUE);
    blue_btn->selection_color(FL_BLUE);
    blue_btn->callback(color_callback, canvas);
    color_y += color_btn_height;
    
    Fl_Round_Button *yellow_btn = new Fl_Round_Button(830, color_y, btn_width - 20, color_btn_height, "Yellow");
    yellow_btn->type(FL_RADIO_BUTTON);
    yellow_btn->color(FL_YELLOW);
    yellow_btn->selection_color(FL_YELLOW);
    yellow_btn->callback(color_callback, canvas);
    
    color_group->end();
    y_pos += 200;
    
    // ================== 形状选择组 ==================
    Fl_Group *shape_group = new Fl_Group(820, y_pos, btn_width, 120, "Drawing Shape");
    shape_group->box(FL_ENGRAVED_FRAME);
    shape_group->align(FL_ALIGN_TOP | FL_ALIGN_INSIDE);
    
    int shape_y = y_pos + 30;
    
    // 为每个形状创建数据对象
    ShapeData* circle_data = new ShapeData{canvas, CIRCLE};
    ShapeData* rect_data = new ShapeData{canvas, RECTANGLE};
    ShapeData* triangle_data = new ShapeData{canvas, TRIANGLE};
    
    Fl_Round_Button *circle_btn = new Fl_Round_Button(830, shape_y, btn_width - 20, color_btn_height, "Circle");
    circle_btn->type(FL_RADIO_BUTTON);
    circle_btn->value(1); // 默认选中
    circle_btn->callback(shape_callback, circle_data);
    shape_y += color_btn_height;
    
    Fl_Round_Button *rect_btn = new Fl_Round_Button(830, shape_y, btn_width - 20, color_btn_height, "Rectangle");
    rect_btn->type(FL_RADIO_BUTTON);
    rect_btn->callback(shape_callback, rect_data);
    shape_y += color_btn_height;
    
    Fl_Round_Button *triangle_btn = new Fl_Round_Button(830, shape_y, btn_width - 20, color_btn_height, "Triangle");
    triangle_btn->type(FL_RADIO_BUTTON);
    triangle_btn->callback(shape_callback, triangle_data);
    
    shape_group->end();
    y_pos += 130;
    
    // ================== 退出按钮 ==================
    Fl_Button *quit_btn = new Fl_Button(820, y_pos, btn_width, btn_height, "Quit");
    quit_btn->callback([](Fl_Widget *, void *w) { 
        static_cast<Fl_Window *>(w)->hide(); 
    }, window);
    
    window->end();
    window->show();
    
    int ret = Fl::run();
    
    // 清理形状数据对象
    delete circle_data;
    delete rect_data;
    delete triangle_data;
    
    return ret;
}

// 回调函数实现
static void clear_callback(Fl_Widget *widget, void *data) {
    drawing_canvas *canvas = static_cast<drawing_canvas *>(data);
    canvas->clear_points();
}

static void print_callback(Fl_Widget *widget, void *data) {
    drawing_canvas *canvas = static_cast<drawing_canvas *>(data);
    canvas->print_points();
}

static void color_callback(Fl_Widget *widget, void *data) {
    Fl_Round_Button *btn = static_cast<Fl_Round_Button *>(widget);
    if (btn->value()) {
        drawing_canvas *canvas = static_cast<drawing_canvas *>(data);
        canvas->set_color(btn->color());
    }
}

static void poly_button_callback(Fl_Widget *widget, void *data) {
    Fl_Round_Button *btn = static_cast<Fl_Round_Button *>(widget);
    drawing_canvas *canvas = static_cast<drawing_canvas *>(data);
    canvas->set_color(btn->color());
    int n = canvas->points.size();
    MatrixXd A(n, n);
    VectorXd b(n);
    int minx = 810;
    int maxx = 0;
    for(int j = 0; j < canvas->points.size(); j++) {
        auto point = canvas->points[j];
        b(j) = point.y;
        VectorXd tmp(n);
        minx = min(minx, point.x);
        maxx = max(minx, point.x);
        for(int i = 0; i < n; i++) {
            A(j, i) = pow(point.x, i);
        }
    }

    auto res = solve_linear_system(A, b);
    for(float m = minx; m < maxx; m += 0.1){
        double tmpy = 0;
        for(int n = 0; n < res.size(); n++) {
            tmpy += res(n) * pow(m, n);
        }
        std::cout << tmpy << std::endl;
        canvas->add_point({(int)m, (int)tmpy, FL_RED, CIRCLE});
    }
}

// 形状选择回调函数
static void shape_callback(Fl_Widget *widget, void *data) {
    Fl_Round_Button *btn = static_cast<Fl_Round_Button *>(widget);
    if (btn->value()) {
        ShapeData* shape_data = static_cast<ShapeData*>(data);
        if (shape_data && shape_data->canvas) {
            shape_data->canvas->set_shape(shape_data->shape);
        }
    }
}

// 新增按钮回调函数
static void new_button_callback(Fl_Widget *widget, void *data) {
    const char* action = static_cast<const char*>(data);
    std::cout << "Button clicked: " << action << std::endl;
    
    // 这里可以添加实际功能代码
    if (strcmp(action, "Save") == 0) {
        // 实现保存功能
    }
    else if (strcmp(action, "Load") == 0) {
        // 实现加载功能
    }
    else if (strcmp(action, "Undo") == 0) {
        // 实现撤销功能
    }
    else if (strcmp(action, "Size") == 0) {
        // 实现改变点大小功能
    }
    else if (strcmp(action, "Background") == 0) {
        // 实现切换背景功能
    }
}