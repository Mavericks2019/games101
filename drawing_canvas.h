#ifndef DRAWING_CANVAS_H
#define DRAWING_CANVAS_H

#include <FL/Fl.H>
#include <FL/Fl_Box.H>
#include <vector>

// 增加形状枚举类型
enum ShapeType {
    CIRCLE,
    RECTANGLE,
    TRIANGLE
};

struct point {
    int x;
    int y;
    Fl_Color color;
    ShapeType shape; // 增加形状类型
};

class drawing_canvas : public Fl_Box {
private:
    std::vector<point> points;
    Fl_Color current_color;
    ShapeType current_shape; // 当前选择的形状
    bool is_drawing;
    int last_x;
    int last_y;
    const int point_size = 5;

public:
    drawing_canvas(int x, int y, int w, int h);
    void clear_points();
    void print_points() const;
    void set_color(Fl_Color c);
    void set_shape(ShapeType s); // 设置形状方法
    void draw() override;
    int handle(int event) override;
};

#endif // DRAWING_CANVAS_H