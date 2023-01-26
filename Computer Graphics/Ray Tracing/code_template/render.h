#ifndef CENG477_HW1_RENDER_H
#define CENG477_HW1_RENDER_H

#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <cmath>

struct RGB{
    unsigned char r;
    unsigned char g;
    unsigned char b;
};

#define EPSILON 1e-7

using namespace parser;
using namespace std;

class Renderer{
public:
    Vec3i background_color;
    float shadow_ray_epsilon;
    int max_recursion_depth;
    std::vector<Camera> cameras;
    Vec3f ambient_light;
    std::vector<PointLight> point_lights;
    std::vector<Material> materials;
    std::vector<Vec3f> vertex_data;
    std::vector<Mesh> meshes;
    std::vector<Triangle> triangles;
    std::vector<Sphere> spheres;



public:
    Renderer() = default;
    ~Renderer();
    Renderer(Scene* scene_ptr);
    void Render(Camera& cam);
    Vec3f ray_value(int x, int y, Camera *cam, Vec3f u);
    void write_to_ppm(const char *filename, unsigned char *image, int width, int height);
    Vec3f cross(Vec3f v1, Vec3f v2);
    float dot(Vec3f v1, Vec3f v2);
    Vec3f sub(Vec3f v1, Vec3f v2);
    Vec3f add(Vec3f v1, Vec3f v2);
    Vec3f normalize(Vec3f v);
    Vec3f mult(float constant, Vec3f v1);
    float intersect_sphere(const Vec3f &ray_direction, const Vec3f &ray_origin, const Sphere &sphere, Vec3f *normal);
    float intersect_triangle(const Vec3f &ray_direction, const Vec3f &ray_origin, const Vec3f &v0, const Vec3f &v1,
                             const Vec3f &v2, Vec3f *normal);
    bool shadow_intersect(const Vec3f &shadow_origin, const Vec3f &light_origin, const Vec3f &normal);
    Vec3f pixel_value(const Vec3f &ray_direction, const Vec3f &ray_origin, int material_id, float t_val,
                      const Vec3f &normal);
    float intersection_check(const Vec3f& direction, const Vec3f& origin, Vec3f *normal, int *material_id);
};




#endif //CENG477_HW1_RENDER_H
