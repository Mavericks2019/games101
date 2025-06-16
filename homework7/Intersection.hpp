//
// Created by LEI XU on 5/16/19.
//

#ifndef RAYTRACING_INTERSECTION_H
#define RAYTRACING_INTERSECTION_H
#include "Vector.hpp"
#include "Material.hpp"
class Object;
class Sphere;

struct Intersection
{
    Intersection(){
        happened=false;
        coords=Vector3f();
        normal=Vector3f();
        distance= std::numeric_limits<double>::max();
        obj =nullptr;
        m=nullptr;
    }
    bool happened;
    Vector3f coords;
    Vector3f tcoords;
    Vector3f normal;
    Vector3f emit;
    double distance;
    Object* obj;
    Material* m;
    //happened 是否相交
    //obj 相交的物体
    //distance 相交点离相机的距离
    //normal 三角形法向量
    //coords 相交点的坐标
    //m 相交点表面材质
};
#endif //RAYTRACING_INTERSECTION_H
