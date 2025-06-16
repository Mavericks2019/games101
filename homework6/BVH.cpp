#include <algorithm>
#include <cassert>
#include "BVH.hpp"

BVHAccel::BVHAccel(std::vector<Object*> p, int maxPrimsInNode,
                   SplitMethod splitMethod)
    : maxPrimsInNode(std::min(255, maxPrimsInNode)), splitMethod(splitMethod),
      primitives(std::move(p))
{
    time_t start, stop;
    time(&start);
    if (primitives.empty())
        return;

    //root = recursiveBuild(primitives);
    root = recursiveBuildSvh(primitives);

    time(&stop);
    double diff = difftime(stop, start);
    int hrs = (int)diff / 3600;
    int mins = ((int)diff / 60) - (hrs * 60);
    int secs = (int)diff - (hrs * 3600) - (mins * 60);

    printf(
        "\rBVH Generation complete: \nTime Taken: %i hrs, %i mins, %i secs\n\n",
        hrs, mins, secs);
}

BVHBuildNode* BVHAccel::recursiveBuild(std::vector<Object*> objects)
{
    BVHBuildNode* node = new BVHBuildNode();

    // Compute bounds of all primitives in BVH node
    Bounds3 bounds;
    for (int i = 0; i < objects.size(); ++i)
        bounds = Union(bounds, objects[i]->getBounds());
    if (objects.size() == 1) {
        // Create leaf _BVHBuildNode_
        node->bounds = objects[0]->getBounds();
        node->object = objects[0];
        node->left = nullptr;
        node->right = nullptr;
        return node;
    }
    else if (objects.size() == 2) {
        node->left = recursiveBuild(std::vector{objects[0]});
        node->right = recursiveBuild(std::vector{objects[1]});

        node->bounds = Union(node->left->bounds, node->right->bounds);
        return node;
    }
    else {
        Bounds3 centroidBounds;
        for (int i = 0; i < objects.size(); ++i)
            centroidBounds =
                Union(centroidBounds, objects[i]->getBounds().Centroid());
        int dim = centroidBounds.maxExtent();
        switch (dim) {
        case 0:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                return f1->getBounds().Centroid().x <
                       f2->getBounds().Centroid().x;
            });
            break;
        case 1:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                return f1->getBounds().Centroid().y <
                       f2->getBounds().Centroid().y;
            });
            break;
        case 2:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                return f1->getBounds().Centroid().z <
                       f2->getBounds().Centroid().z;
            });
            break;
        }

        auto beginning = objects.begin();
        auto middling = objects.begin() + (objects.size() / 2);
        auto ending = objects.end();

        auto leftshapes = std::vector<Object*>(beginning, middling);
        auto rightshapes = std::vector<Object*>(middling, ending);

        assert(objects.size() == (leftshapes.size() + rightshapes.size()));

        node->left = recursiveBuild(leftshapes);
        node->right = recursiveBuild(rightshapes);

        node->bounds = Union(node->left->bounds, node->right->bounds);
    }

    return node;
}

BVHBuildNode* BVHAccel::recursiveBuildSvh(std::vector<Object*> objects)
{
    BVHBuildNode* node = new BVHBuildNode();
    int size = objects.size();

    // Compute bounds of all primitives in BVH node
    Bounds3 bounds;
    for (int i = 0; i < objects.size(); ++i)
        bounds = Union(bounds, objects[i]->getBounds());
    if (objects.size() == 1) {
        // Create leaf _BVHBuildNode_
        node->bounds = objects[0]->getBounds();
        node->object = objects[0];
        node->left = nullptr;
        node->right = nullptr;
        return node;
    }
    else if (objects.size() == 2) {
        node->left = recursiveBuildSvh(std::vector{objects[0]});
        node->right = recursiveBuildSvh(std::vector{objects[1]});

        node->bounds = Union(node->left->bounds, node->right->bounds);
        return node;
    }
    else 
    {
        int bestDiv = 0;
        int bestDim = 0;
        float minCost = std::numeric_limits<float>::max();
        int div_num = 5;
        float div[] = { 1.0 / 6, 2.0 / 6, 3.0 / 6, 4.0 / 6, 5.0 / 6 };

        for (int i = 0; i < 5; i++){
             div[i] *= size;
        }

        for (int dim = 0; dim < 3; dim++) {
            std::sort(objects.begin(), objects.end(), cmp[dim]);
            for (int i = 0; i < 5; i++){
                auto l = objects.begin();
                auto r = objects.end();
                auto mid = l + (int)div[i];
    
                auto leftshapes = std::vector<Object*>(l, mid);
                auto rightshapes = std::vector<Object*>(mid, r);

                Bounds3 uBounds, vBounds;
                for (int i = 0; i < leftshapes.size(); ++i){
                    uBounds = Union(uBounds, leftshapes[i]->getBounds().Centroid());
                }
                
                for (int i = 0; i < rightshapes.size(); ++i){
                    vBounds = Union(vBounds, rightshapes[i]->getBounds().Centroid());
                }
                auto leftBoxSize = uBounds.SurfaceArea();
                auto rightBoxSize = vBounds.SurfaceArea();

                double cost = 100.0 + (leftBoxSize * leftshapes.size() + rightBoxSize * rightshapes.size()) 
                                            / bounds.SurfaceArea();

                auto u = std::vector<Object*>(l, mid), v = std::vector<Object*>(mid, r);
                if (cost < minCost)
                {
                    minCost = cost;
                    bestDiv = (int)div[i];
                    bestDim = dim;
                }
            }
        }
                    // 根据选择更新并建子树
        if(bestDim != 2){
            std::sort(objects.begin(), objects.end(), cmp[bestDim]);
        }
        auto l = objects.begin();
        auto r = objects.end();
        auto mid = l + bestDiv;

        auto u = std::vector<Object*>(l, mid);
        auto v= std::vector<Object*>(mid, r);

        node->left = recursiveBuildSvh(u);
        node->right = recursiveBuildSvh(v);

        node->bounds = Union(node->left->bounds, node->right->bounds);
    }
    return node;
}


Intersection BVHAccel::Intersect(const Ray& ray) const
{
    Intersection isect;
    if (!root)
        return isect;
    isect = BVHAccel::getIntersection(root, ray);
    return isect;
}

Intersection BVHAccel::getIntersection(BVHBuildNode* node, const Ray& ray) const
{
    // TODO Traverse the BVH to find intersection
    std::array<int,3> dirIsNeg;
    dirIsNeg[0] = ray.direction.x > 0 ? 1:0;
    dirIsNeg[1] = ray.direction.y > 0 ? 1:0;
    dirIsNeg[2] = ray.direction.z > 0 ? 1:0;

    if (!node->bounds.IntersectP(ray, ray.direction_inv, dirIsNeg)){
        return Intersection();
    }
    if (node->left == nullptr && node->right == nullptr){
         return  node->object->getIntersection(ray);
    }
    
    auto u = getIntersection(node->left, ray);
    auto v = getIntersection(node->right, ray);
    return u.distance < v.distance ? u : v;

}