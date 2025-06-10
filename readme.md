# FLTK Drawing Application

This is a drawing application built with the FLTK (Fast Light Toolkit) library. It features a black canvas where users can draw with different colors, clear the canvas, and print point coordinates to the console.

## Features

- Black drawing canvas with adjustable point size
- 7-color palette selection via radio buttons
- Clear canvas functionality
- Print all drawn points to console
- Responsive drawing with mouse

## Requirements

### Operating System
- Linux (Ubuntu/Debian recommended)
- Windows (with MinGW)
- macOS (with Homebrew)

### Dependencies
- FLTK 1.3.x or newer
- CMake 3.10+
- C++17 compatible compiler
- Math library (`-lm`)

## Installation

### Linux (Ubuntu/Debian)
```bash
sudo apt update
sudo apt install build-essential cmake libfltk1.3-dev

## Building with CMake
cmake_minimum_required(VERSION 3.10)
project(FLTK_Drawing_App)

set(CMAKE_CXX_STANDARD 17)

find_package(FLTK REQUIRED)
include_directories(${FLTK_INCLUDE_DIRS})

add_executable(drawing_app kernel_style_drawing_app.cpp)
target_link_libraries(drawing_app ${FLTK_LIBRARIES} m)

## Build and run the application:
mkdir build
cd build
cmake ..
make
./drawing_app