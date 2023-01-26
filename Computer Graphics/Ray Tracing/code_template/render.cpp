#include "render.h"
#include <limits>


Renderer::Renderer(Scene *scene_ptr) {
    background_color = scene_ptr->background_color;
    shadow_ray_epsilon = scene_ptr->shadow_ray_epsilon;
    max_recursion_depth = scene_ptr->max_recursion_depth;
    cameras = scene_ptr->cameras;
    ambient_light = scene_ptr->ambient_light;
    point_lights = scene_ptr->point_lights;
    materials = scene_ptr->materials;
    vertex_data = scene_ptr->vertex_data;
    meshes = scene_ptr->meshes;
    triangles = scene_ptr->triangles;
    spheres = scene_ptr->spheres;
}


Renderer::~Renderer() {

}


void Renderer::Render(Camera& cam) {
    int width = cam.image_width;
    int height = cam.image_height;
    Vec3f origin = cam.position;

    /* CALCULATE U VECTOR*/
    Vec3f gaze_op = (Vec3f) {.x = cam.gaze.x, .y = cam.gaze.y, .z = cam.gaze.z};
    Vec3f u = cross(gaze_op, cam.up);

    auto* image = new unsigned char [width * height * 3];
    auto* function_color = new unsigned char[3];
    Vec3f* normal = new Vec3f;
    int *material_id = new int;
    Vec3f min_normal;
    for(int y = 0 ; y < height ; y++){
        for(int x = 0 ; x < width ; x++){
            bool is_intersect = false;
            Vec3f direction = ray_value(x, y, &cam, u);

            float t = intersection_check(direction, origin,normal,material_id);
            if(t != -1){is_intersect = true;}
            if (is_intersect){
                Vec3f val = pixel_value(direction,origin,*material_id,t,*normal);
                image[(y * width + x)*3] = (int) val.x;
                image[(y * width + x)*3 + 1] = (int) val.y;
                image[(y * width + x)*3 + 2] = (int) val.z;
            }
            // if no intersection than assign background color
            else {
                image[(y * width + x) * 3] = background_color.x;
                image[(y * width + x) * 3 + 1] = background_color.y;
                image[(y * width + x) * 3 + 2] = background_color.z;
            }
        }

    }
    write_to_ppm(cam.image_name.c_str(), image, width, height);
    delete[] function_color;
    delete[] image; // free image
    delete normal;
    delete material_id;
}


Vec3f Renderer::ray_value(int x, int y, Camera *cam, Vec3f u) {
    float l = cam->near_plane.x;
    float r = cam->near_plane.y;
    float b = cam->near_plane.z;
    float t = cam->near_plane.w;
    Vec3f v = cam->up;
    Vec3f position = cam->position;

    Vec3f m = (Vec3f) {.x = position.x + (cam->gaze.x * cam->near_distance),
            .y = position.y + (cam->gaze.y * cam->near_distance),
            .z = position.z + (cam->gaze.z * cam->near_distance)};

    //left,right,bottom,top   x, y, z, w
    Vec3f q = (Vec3f) { .x = m.x + l * u.x + t * cam->up.x,
            .y = m.y + l * u.y + t * cam->up.y,
            .z = m.z + l * u.z+ t * cam->up.z};

    float su = (float)(x + 0.5) * ((r-l) / cam->image_width);
    float sv = (float)(y + 0.5) * ((t-b) / cam->image_height);

    Vec3f s = (Vec3f) { .x = q.x + su * u.x - sv * v.x,
            .y = q.y + su * u.y - sv * v.y,
            .z = q.z + su * u.z - sv * v.z};

    Vec3f d = (Vec3f) { .x = s.x - position.x,
            .y = s.y - position.y,
            .z = s.z - position.z};
    return d;
}

void Renderer::write_to_ppm(const char* filename, unsigned char *image, int width, int height) {
    write_ppm(filename, image, width, height);
}

Vec3f Renderer::add(Vec3f v1, Vec3f v2) {
    Vec3f ret;
    ret.x = v1.x + v2.x;
    ret.y = v1.y + v2.y;
    ret.z = v1.z + v2.z;
    return ret;
}

Vec3f Renderer::sub(Vec3f v1, Vec3f v2) {
    Vec3f ret;
    ret.x = v1.x - v2.x;
    ret.y = v1.y - v2.y;
    ret.z = v1.z - v2.z;
    return ret;
}

Vec3f Renderer::cross(Vec3f v1, Vec3f v2) {
    Vec3f cross;
    cross.x = v1.y * v2.z - v1.z * v2.y;
    cross.y = v1.z * v2.x - v1.x * v2.z;
    cross.z = v1.x * v2.y - v1.y * v2.x;

    return cross;
}

float Renderer::dot(Vec3f v1, Vec3f v2) {
    float total = 0;
    total += v1.x * v2.x;
    total += v1.y * v2.y;
    total += v1.z * v2.z;
    return total;
}

Vec3f Renderer::mult(float constant, Vec3f v1) {
    Vec3f ret_vec = {
            .x = v1.x * constant,
            .y = v1.y * constant,
            .z = v1.z * constant,
    };
    return ret_vec;
}

float Renderer::intersect_sphere(const Vec3f &ray_direction, const Vec3f &ray_origin, const Sphere &sphere, Vec3f *normal) {

    Vec3f sphere_center = vertex_data.at(sphere.center_vertex_id-1);
    float radius = sphere.radius;

    // discriminant = sqrt((d.(o-c))^2 - (d.d)((o-c).(o-c)-R^2))

    //o-c
    Vec3f origin_difference = this->sub(ray_origin, sphere_center);
    float disc = pow(dot(ray_direction,origin_difference),2) - (dot(ray_direction, ray_direction) * (dot(origin_difference,origin_difference) - pow(radius, 2)));

    if(disc < EPSILON){
        return -1;
    }
    float t1 = (-1*dot(ray_direction, sub(ray_origin, sphere_center)) + sqrt(disc)) / dot(ray_direction, ray_direction);
    float t2 = (-1*dot(ray_direction, sub(ray_origin, sphere_center)) - sqrt(disc)) / dot(ray_direction, ray_direction);
    if(t2 < t1){
        t1 = t2;
    }

    //calculate normal
    Vec3f hit_point = add(ray_origin, mult(t1, ray_direction));
    *normal = normalize(mult((1/(radius*radius)), sub(hit_point, sphere_center)));

    return t1;
}

float
Renderer::intersect_triangle(const Vec3f &ray_direction, const Vec3f &ray_origin, const Vec3f &v0, const Vec3f &v1,
                             const Vec3f &v2, Vec3f *normal) {

    const Vec3f v0v1 = sub(v1, v0);
    const Vec3f v0v2 = sub(v2, v0);
    const Vec3f pvec = cross(ray_direction, v0v2);
    float det = dot(v0v1, pvec);

    // parallel
    if(fabs(det) < EPSILON){return -1;}
    float inv_det = 1.0/det;

    // v, to origin distance
    const Vec3f tvec = sub(ray_origin, v0);
    float u = dot(tvec, pvec) * inv_det;
    if(u < 0.0 || u > 1.0){return -1;}

    const Vec3f qvec = cross(tvec, v0v1);
    float v = dot(ray_direction, qvec) * inv_det;
    if(v < 0.0 || (u + v) > 1.0){return -1;}

    float t = dot(v0v2, qvec) * inv_det;
    if(t < 0.0){return -1;}

    /*t *= inv_det;
    u *= inv_det;
    v *= inv_det;*/

    // normal calculation
    Vec3f hit_point = add(ray_origin, mult(t, ray_direction));
    *normal = normalize(cross(v0v1, v0v2));


    return t;
}

Vec3f Renderer::normalize(Vec3f v) {
    float len = sqrt((v.x * v.x) + (v.y * v.y) + (v.z * v.z));
    Vec3f normalized;
    normalized = {    .x = v.x / len,
            .y = v.y / len,
            .z = v.z / len};
    return normalized;
}

bool Renderer::shadow_intersect(const Vec3f &shadow_origin, const Vec3f &light_origin, const Vec3f &normal) {

    Vec3f epsilon_origin = add(shadow_origin, mult(shadow_ray_epsilon, normal));
    Vec3f shadow_direction = sub(light_origin, epsilon_origin);
    for(auto sphere: spheres){

        Vec3f sphere_center = vertex_data.at(sphere.center_vertex_id-1);
        float radius = sphere.radius;
        // discriminant = sqrt((d.(o-c))^2 - (d.d)((o-c).(o-c)-R^2))

        //o-c
        float dist1 = sqrt(pow(light_origin.x-sphere_center.x,2) + pow(light_origin.y-sphere_center.y,2) + pow(light_origin.z-sphere_center.z,2));
        float dist2 = sqrt(pow(light_origin.x-epsilon_origin.x,2) + pow(light_origin.y-epsilon_origin.y,2) + pow(light_origin.z-epsilon_origin.z,2));
        if(dist1 > dist2){
            continue;
        }
        Vec3f origin_difference = this->sub(epsilon_origin, sphere_center);
        float disc = pow(dot(shadow_direction,origin_difference),2) -
                (dot(shadow_direction, shadow_direction) * (dot(origin_difference,origin_difference)
                - pow(radius, 2)));

        if(disc >= EPSILON){
            return true;
        }
    }
    for(auto mesh : meshes){

        for(auto face: mesh.faces){
            Vec3f v0 = vertex_data.at(face.v0_id-1);
            Vec3f v1 = vertex_data.at(face.v1_id-1);
            Vec3f v2 = vertex_data.at(face.v2_id-1);

            const Vec3f v0v1 = sub(v1, v0);
            const Vec3f v0v2 = sub(v2, v0);
            const Vec3f pvec = cross(shadow_direction, v0v2);
            float det = dot(v0v1, pvec);

            // parallel
            if(fabs(det) < EPSILON){ continue;}
            float inv_det = 1.0/det;

            // v, to origin distance
            const Vec3f tvec = sub(epsilon_origin, v0);
            float u = dot(tvec, pvec) * inv_det;
            if(u < 0.0 || u > 1.0){ continue;}

            const Vec3f qvec = cross(tvec, v0v1);
            float v = dot(shadow_direction, qvec) * inv_det;
            if(v < 0.0 || (u + v) > 1.0){ continue;}

            float t = dot(v0v2, qvec) * inv_det;
            if(t < 0.0){ continue;}
            return true;
        }
    }
    for(auto triangle: triangles){
        Face face = triangle.indices;
        Vec3f v0 = vertex_data.at(face.v0_id-1);
        Vec3f v1 = vertex_data.at(face.v1_id-1);
        Vec3f v2 = vertex_data.at(face.v2_id-1);

        const Vec3f v0v1 = sub(v1, v0);
        const Vec3f v0v2 = sub(v2, v0);
        const Vec3f pvec = cross(shadow_direction, v0v2);
        float det = dot(v0v1, pvec);

        // parallel
        if(fabs(det) < EPSILON){ continue;}
        float inv_det = 1.0/det;

        // v, to origin distance
        const Vec3f tvec = sub(epsilon_origin, v0);
        float u = dot(tvec, pvec) * inv_det;
        if(u < 0.0 || u > 1.0){ continue;}

        const Vec3f qvec = cross(tvec, v0v1);
        float v = dot(shadow_direction, qvec) * inv_det;
        if(v < 0.0 || (u + v) > 1.0){ continue;}

        float t = dot(v0v2, qvec) * inv_det;
        if(t < 0.0){ continue;}
        return true;
    }
    return false;
}

Vec3f Renderer::pixel_value(const Vec3f &ray_direction, const Vec3f &ray_origin, int material_id, float t_val,
                            const Vec3f &normal) {

    Material material = materials.at(material_id-1);
    // DIFFUSE SHADING
    Vec3f hit_point = add(ray_origin, mult(t_val, ray_direction));
    Vec3f value = {.x=0.0,.y=0.0,.z=0.0};
    for (auto light: point_lights){

        // SHADOW CHECK
        bool is_shadow = shadow_intersect(hit_point, light.position, normal);
        if(is_shadow){ continue;}

        Vec3f light_dir = normalize(sub(light.position, hit_point));
        float angle = dot(light_dir, normal);
        float cos_alpha = (angle < 0) ? 0.0 : angle;
        if (cos_alpha == 0.0){
            continue;
        }
        Vec3f v = sub(light.position, hit_point);
        float light_distance = sqrt((v.x*v.x) + (v.y*v.y) + (v.z*v.z));
        Vec3f received_irradiance = mult(1/(light_distance * light_distance), light.intensity);
        Vec3f this_light = mult(cos_alpha, received_irradiance);
        Vec3f diffuse = materials.at(material_id-1).diffuse;
        this_light.x = diffuse.x * this_light.x;
        this_light.y = diffuse.y * this_light.y;
        this_light.z = diffuse.z * this_light.z;

        // SPECULAR SHADING
        Vec3f specular = materials.at(material_id-1).specular;
        Vec3f h = normalize(add(normalize(sub(light.position, hit_point)), normalize(sub(ray_origin, hit_point))));
        float angle2 = dot(normal, h);
        float cos_beta = (angle2 < 0) ? 0.0 : angle2;
        // phong exponent
        if(cos_beta != 0){
            cos_beta = pow(cos_beta, materials.at(material_id-1).phong_exponent);
            this_light.x += specular.x * cos_beta * received_irradiance.x;
            this_light.y += specular.y * cos_beta * received_irradiance.y;
            this_light.z += specular.z * cos_beta * received_irradiance.z;
        }
        value = add(value, this_light);
    }

    // material reflecting light
    /*Vec3f *reflect_normal = new Vec3f;
    int *reflect_material_id = new int;
    if(material.is_mirror){
        float cos_theta = dot(sub(ray_origin, hit_point),normal); //w0.n
        Vec3f reflect_direction = add(mult(-1,normal), mult(2*cos_theta,normal)); // -w0 + 2cos_theta*n
        Vec3f reflect_origin = add(hit_point, mult(shadow_ray_epsilon, reflect_direction));

        float t = intersection_check(reflect_direction, reflect_origin, reflect_normal, reflect_material_id);

        if(t != -1){
            Vec3f pixel_value = this->pixel_value(reflect_direction, reflect_origin, material_id, t, *reflect_normal);
            value.x += material.mirror.x * pixel_value.x;
            value.y += material.mirror.y * pixel_value.y;
            value.z += material.mirror.z * pixel_value.z;
        }

    }*/
    //delete reflect_normal;
    // AMBIENT SHADING
    Vec3f ambient = materials.at(material_id-1).ambient;
    value.x += ambient_light.x * ambient.x;
    value.y += ambient_light.y * ambient.y;
    value.z += ambient_light.z * ambient.z;
    // END OF SHADING

    /*if (value.x>=255){value.x=255.0;}
    if (value.y>=255){value.y=255.0;}
    if (value.z>=255){value.z=255.0;}*/

    value.x = min(max(0.0f, value.x), 255.0f);
    value.y = min(max(0.0f, value.y), 255.0f);
    value.z = min(max(0.0f, value.z), 255.0f);

    return value;
}

float Renderer::intersection_check(const Vec3f& direction, const Vec3f& origin, Vec3f *normal, int *material_id) {

    float t = numeric_limits<float>::max();
    Vec3f min_normal;
    bool is_intersect;
    // SPHERE INTERSECTION
    for(auto sphere : this->spheres){
        float intersect = intersect_sphere(direction, origin, sphere, normal);
        if(intersect != -1){
            if(intersect < t){
                t = intersect;
                is_intersect = true;
                *material_id = sphere.material_id;
                min_normal = *normal;
            }
        }
    }
    // MESH INTERSECTION
    for(const auto& mesh : this->meshes){
        for(const auto& face: mesh.faces){
            Vec3f v0 = vertex_data.at(face.v0_id-1);
            Vec3f v1 = vertex_data.at(face.v1_id-1);
            Vec3f v2 = vertex_data.at(face.v2_id-1);
            float intersect = intersect_triangle(direction, origin, v0, v1, v2,
                                                 normal);
            if(intersect != -1){
                if(intersect < t){
                    t = intersect;
                    is_intersect = true;
                    *material_id = mesh.material_id;
                    min_normal = *normal;
                }
            }
        }
    }
    // TRIANGLE INTERSECTION
    for(auto triangle: this->triangles){
        Face face = triangle.indices;
        Vec3f v0 = vertex_data.at(face.v0_id-1);
        Vec3f v1 = vertex_data.at(face.v1_id-1);
        Vec3f v2 = vertex_data.at(face.v2_id-1);
        float intersect = intersect_triangle(direction, origin, v0, v1, v2,
                                             normal);
        if(intersect != -1){
            if(intersect < t){
                t = intersect;
                is_intersect = true;
                *material_id = triangle.material_id;
                min_normal = *normal;
            }
        }
    }
    *normal = min_normal;
    if (is_intersect){
        return t;
    }
    else {
        return -1;
    }
}

