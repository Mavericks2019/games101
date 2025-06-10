#include "drawing_canvas.h"
#include <FL/Fl_Window.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Round_Button.H>
#include <FL/Fl_Group.H>
#include <iostream> // 添加用于调试

// 前向声明回调函数
static void clear_callback(Fl_Widget *widget, void *data);
static void print_callback(Fl_Widget *widget, void *data);
static void color_callback(Fl_Widget *widget, void *data);
static void shape_callback(Fl_Widget *widget, void *data); // 形状回调

// 安全传递形状类型的辅助结构
struct ShapeData {
    drawing_canvas* canvas;
    ShapeType shape;
};

int main() {
    Fl_Window *window = new Fl_Window(1000, 600, "FLTK Drawing Canvas");
    
    drawing_canvas *canvas = new drawing_canvas(0, 0, 800, 600);
    
    Fl_Box *toolbar = new Fl_Box(800, 0, 200, 600);
    toolbar->box(FL_ENGRAVED_BOX);
    toolbar->color(fl_rgb_color(240, 240, 240));
    
    Fl_Button *clear_btn = new Fl_Button(820, 30, 160, 40, "Clear Canvas");
    clear_btn->callback(clear_callback, canvas);
    
    Fl_Button *print_btn = new Fl_Button(820, 80, 160, 40, "Print Points");
    print_btn->callback(print_callback, canvas);
    
    Fl_Box *size_info = new Fl_Box(820, 130, 160, 30, "Point Size: 5px");
    size_info->labelsize(16);
    
    // ================== 颜色选择组 ==================
    Fl_Group *color_group = new Fl_Group(820, 170, 160, 180, "Drawing Color");
    color_group->box(FL_ENGRAVED_FRAME);
    color_group->align(FL_ALIGN_TOP | FL_ALIGN_INSIDE);
    
    int y_pos = 200;
    const int btn_height = 30;
    
    Fl_Round_Button *white_btn = new Fl_Round_Button(830, y_pos, 140, btn_height, "White");
    white_btn->type(FL_RADIO_BUTTON);
    white_btn->value(1);
    white_btn->color(FL_WHITE);
    white_btn->selection_color(FL_WHITE);
    white_btn->callback(color_callback, canvas);
    y_pos += btn_height;
    
    Fl_Round_Button *red_btn = new Fl_Round_Button(830, y_pos, 140, btn_height, "Red");
    red_btn->type(FL_RADIO_BUTTON);
    red_btn->color(FL_RED);
    red_btn->selection_color(FL_RED);
    red_btn->callback(color_callback, canvas);
    y_pos += btn_height;
    
    Fl_Round_Button *green_btn = new Fl_Round_Button(830, y_pos, 140, btn_height, "Green");
    green_btn->type(FL_RADIO_BUTTON);
    green_btn->color(FL_GREEN);
    green_btn->selection_color(FL_GREEN);
    green_btn->callback(color_callback, canvas);
    y_pos += btn_height;
    
    Fl_Round_Button *blue_btn = new Fl_Round_Button(830, y_pos, 140, btn_height, "Blue");
    blue_btn->type(FL_RADIO_BUTTON);
    blue_btn->color(FL_BLUE);
    blue_btn->selection_color(FL_BLUE);
    blue_btn->callback(color_callback, canvas);
    y_pos += btn_height;
    
    Fl_Round_Button *yellow_btn = new Fl_Round_Button(830, y_pos, 140, btn_height, "Yellow");
    yellow_btn->type(FL_RADIO_BUTTON);
    yellow_btn->color(FL_YELLOW);
    yellow_btn->selection_color(FL_YELLOW);
    yellow_btn->callback(color_callback, canvas);
    y_pos += btn_height;
    
    color_group->end();
    
    // ================== 形状选择组 ==================
    Fl_Group *shape_group = new Fl_Group(820, 360, 160, 120, "Drawing Shape");
    shape_group->box(FL_ENGRAVED_FRAME);
    shape_group->align(FL_ALIGN_TOP | FL_ALIGN_INSIDE);
    
    y_pos = 390;
    
    // 为每个形状创建数据对象
    ShapeData* circle_data = new ShapeData{canvas, CIRCLE};
    ShapeData* rect_data = new ShapeData{canvas, RECTANGLE};
    ShapeData* triangle_data = new ShapeData{canvas, TRIANGLE};
    
    Fl_Round_Button *circle_btn = new Fl_Round_Button(830, y_pos, 140, btn_height, "Circle");
    circle_btn->type(FL_RADIO_BUTTON);
    circle_btn->value(1); // 默认选中
    circle_btn->callback(shape_callback, circle_data);
    y_pos += btn_height;
    
    Fl_Round_Button *rect_btn = new Fl_Round_Button(830, y_pos, 140, btn_height, "Rectangle");
    rect_btn->type(FL_RADIO_BUTTON);
    rect_btn->callback(shape_callback, rect_data);
    y_pos += btn_height;
    
    Fl_Round_Button *triangle_btn = new Fl_Round_Button(830, y_pos, 140, btn_height, "Triangle");
    triangle_btn->type(FL_RADIO_BUTTON);
    triangle_btn->callback(shape_callback, triangle_data);
    
    shape_group->end();
    
    // ================== 退出按钮 ==================
    Fl_Button *quit_btn = new Fl_Button(820, 500, 160, 40, "Quit");
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

// 修复的形状选择回调函数
static void shape_callback(Fl_Widget *widget, void *data) {
    Fl_Round_Button *btn = static_cast<Fl_Round_Button *>(widget);
    if (btn->value()) {
        ShapeData* shape_data = static_cast<ShapeData*>(data);
        if (shape_data && shape_data->canvas) {
            shape_data->canvas->set_shape(shape_data->shape);
        }
    }
}