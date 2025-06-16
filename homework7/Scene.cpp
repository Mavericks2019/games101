//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"


void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const
{
    float emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
        }
    }
    float p = get_random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
            if (p <= emit_area_sum){
                objects[k]->Sample(pos, pdf);
                break;
            }
        }
    }
}

bool Scene::trace(
        const Ray &ray,
        const std::vector<Object*> &objects,
        float &tNear, uint32_t &index, Object **hitObject)
{
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
        }
    }


    return (*hitObject != nullptr);
}

// Implementation of Path Tracing
Vector3f Scene::castRay(const Ray &ray, int depth) const
{
    Intersection inter = intersect(ray);// 获取相交信息
    Vector3f l_dir(0, 0, 0);
    Vector3f l_indir(0, 0, 0);
    if(!inter.happened)//光线无相交
    {
        return {};
    }
    if (inter.m->hasEmission()) {
        return inter.m->getEmission();// 打到光源 直接返回光源颜色
    }

    Intersection x;
    float pdf_light;
    sampleLight(x, pdf_light);
    auto objtolight = x.coords - inter.coords;
    auto objtolightDir = objtolight.normalized();
    float dist = objtolight.x * objtolight.x + objtolight.y * objtolight.y + objtolight.z * objtolight.z;

    Ray objtolightRay(inter.coords, objtolightDir);
    Intersection check = intersect(objtolightRay);
    if (abs(check.distance - objtolight.norm()) < EPSILON)// 浮点数的比较
    {
        l_dir = x.emit * inter.m->eval(ray.direction, objtolightDir, inter.normal)* 
        dotProduct(objtolightDir, inter.normal) * 
        dotProduct(-objtolightDir, x.normal) / dist / pdf_light;
    }

    float ksi = get_random_float();
    
    if(ksi <= RussianRoulette)//俄罗斯轮盘赌 需要继续发射光线
    {
        Vector3f wi = inter.m->sample(ray.direction, inter.normal).normalized();
        Ray ray_intertoobj = Ray(inter.coords, wi);

        Intersection tmp = intersect(ray_intertoobj);
        if(tmp.happened==true && !tmp.m->hasEmission())
        {
            float pdf = inter.m->pdf(ray.direction, wi, inter.normal);
            if(pdf > EPSILON){
                l_indir = castRay(ray_intertoobj ,depth + 1) * inter.m->eval(ray.direction, wi, inter.normal)
                        * dotProduct(wi, inter.normal)/ pdf / RussianRoulette;
            }

        }
    }
    return l_dir + l_indir;


    // TO DO Implement Path Tracing Algorithm here
}