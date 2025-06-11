cmake_minimum_required(VERSION 3.10)
project(FLTK_Drawing_Canvas)

# 查找FLTK库
find_package(FLTK REQUIRED)

# 添加可执行文件
add_executable(drawing_app 
    main.cpp
    drawing_canvas.cpp
    drawing_canvas.h
    solver.h
)

# 链接FLTK库
target_link_libraries(drawing_app ${FLTK_LIBRARIES})