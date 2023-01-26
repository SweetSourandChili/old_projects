#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>

#include "Scene.h"
#include "Camera.h"
#include "Color.h"
#include "Mesh.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Triangle.h"
#include "Vec3.h"
#include "tinyxml2.h"
#include "Helpers.h"

using namespace tinyxml2;
using namespace std;

/*
	Transformations, clipping, culling, rasterization are done here.
	You may define helper functions.
*/

void Scene::forwardRenderingPipeline(Camera *camera)
{
	// TODO: Implement this function.

    /************************ VERTEX PROCESSING****************/

    // get eye transformation matrix for the Param:camera
    Matrix4 M_cam = camera->eyeTransformationMatrix();
    // get perspective transformation matrix for the Param:camera
    Matrix4 M_per = camera->projectionMatrix();
    // get viewport transformation matrix for the Param:camera
    Matrix4 M_vp = camera->viewportMatrix(0,0);

    Matrix4 M_per_cam = multiplyMatrixWithMatrix(M_per, M_cam);

    for(auto *mesh: meshes){
        // get model transformation matrix(s,t,r operations) for param:mesh
        Matrix4 M_model = meshTranslationMatrix(*mesh);

        // bind camera and model transformations
        Matrix4 M_per_cam_model = multiplyMatrixWithMatrix(M_per_cam,M_model);

        for(auto triangle: mesh->triangles){
            Vec3 * v0 = this->vertices[triangle.getFirstVertexId()-1];
            Vec3 * v1 = this->vertices[triangle.getSecondVertexId()-1];
            Vec3 * v2 = this->vertices[triangle.getThirdVertexId()-1];
            Color * c0 = this->colorsOfVertices[v0->colorId-1];
            Color * c1 = this->colorsOfVertices[v1->colorId-1];
            Color * c2 = this->colorsOfVertices[v2->colorId-1];

            Vec4 homogeneous0 = multiplyMatrixWithVec4(M_per_cam_model, Vec4(v0->x, v0->y, v0->z, 1, -1));
            Vec4 homogeneous1 = multiplyMatrixWithVec4(M_per_cam_model, Vec4(v1->x, v1->y, v1->z, 1, -1));
            Vec4 homogeneous2 = multiplyMatrixWithVec4(M_per_cam_model, Vec4(v2->x, v2->y, v2->z, 1, -1));

            // perspective division
            homogeneous0 = Vec4(homogeneous0.x / homogeneous0.t, homogeneous0.y / homogeneous0.t, homogeneous0.z / homogeneous0.t, homogeneous0.t / homogeneous0.t, -1);
            homogeneous1 = Vec4(homogeneous1.x / homogeneous1.t, homogeneous1.y / homogeneous1.t, homogeneous1.z / homogeneous1.t, homogeneous1.t / homogeneous1.t, -1);
            homogeneous2 = Vec4(homogeneous2.x / homogeneous2.t, homogeneous2.y / homogeneous2.t, homogeneous2.z / homogeneous2.t, homogeneous2.t / homogeneous2.t, -1);

            /******* CULLING *************/
            if(cullingEnabled){
                if(backFaceCulling(homogeneous0, homogeneous1, homogeneous2, camera->pos)){
                    continue;
                }
            }

            // solid mode - rasTriangles
            if(mesh->type){

                // no clipping
                /******* VIEWPORT TRANSFORM *************/
                homogeneous0 = multiplyMatrixWithVec4(M_vp, homogeneous0);
                homogeneous1 = multiplyMatrixWithVec4(M_vp, homogeneous1);
                homogeneous2 = multiplyMatrixWithVec4(M_vp, homogeneous2);

                /******* RASTERIZATION *************/
                rasTriangle(homogeneous0, homogeneous1, homogeneous2, *c0, *c1, *c2, camera->horRes, camera->verRes);
            }
            // wireframe mode - rasLine
            else{

                /******* CLIPPING *************/
                // create copy of the vectors since they are passed liangBarsky with reference
                Vec4 line1_0(homogeneous0);
                Vec4 line1_1(homogeneous1);
                Vec4 line2_0(homogeneous0);
                Vec4 line2_1(homogeneous2);
                Vec4 line3_0(homogeneous1);
                Vec4 line3_1(homogeneous2);

                Color color1_0(*c0);
                Color color1_1(*c1);
                Color color2_0(*c0);
                Color color2_1(*c2);
                Color color3_0(*c1);
                Color color3_1(*c2);

                /******* VIEWPORT TRANSFORM AND RASTERIZATION*************/
                if(liangBarsky(line1_0, line1_1, &color1_0, &color1_1)){
                    line1_0 = multiplyMatrixWithVec4(M_vp, line1_0);
                    line1_1 = multiplyMatrixWithVec4(M_vp, line1_1);

                    rasLine(line1_0, line1_1, color1_0, color1_1, camera->horRes, camera->verRes);
                }

                if(liangBarsky(line2_0, line2_1, &color2_0, &color2_1)){
                    line2_0 = multiplyMatrixWithVec4(M_vp, line2_0);
                    line2_1 = multiplyMatrixWithVec4(M_vp, line2_1);

                    rasLine(line2_0, line2_1, color2_0, color2_1, camera->horRes, camera->verRes);
                }

                if(liangBarsky(line3_0, line3_1, &color3_0, &color3_1)){
                    line3_0 = multiplyMatrixWithVec4(M_vp, line3_0);
                    line3_1 = multiplyMatrixWithVec4(M_vp, line3_1);

                    rasLine(line3_0, line3_1, color3_0, color3_0, camera->horRes, camera->verRes);
                }
            }
        }
    }
}


Matrix4 Scene::meshTranslationMatrix(Mesh const& mesh) {
    Matrix4 retMatrix = getIdentityMatrix();
    for( int i = 0 ; i < mesh.numberOfTransformations ; i++){
        int id = mesh.transformationIds[i];
        char type = mesh.transformationTypes[i];
        switch (type) {
            case 't': {
                Translation *t_values = translations[id - 1];
                double T_arr[4][4] = {{1, 0, 0, t_values->tx},
                                      {0, 1, 0, t_values->ty},
                                      {0, 0, 1, t_values->tz},
                                      {0, 0, 0, 1}};
                Matrix4 T(T_arr);
                retMatrix = multiplyMatrixWithMatrix(T, retMatrix);
                break;
            }
            case 's': {
                Scaling *s_values = scalings[id - 1];
                double S_arr[4][4] = {{s_values->sx, 0,            0,            0},
                                      {0,            s_values->sy, 0,            0},
                                      {0,            0,            s_values->sz, 0},
                                      {0,            0,            0,            1}};
                Matrix4 S(S_arr);
                retMatrix = multiplyMatrixWithMatrix(S, retMatrix);
                break;
            }
            case 'r':{
                // alternate method for rotation around an axis
                Rotation *r = rotations[id-1];
                Vec3 u(r->ux, r->uy, r->uz, -1);
                //find |the smallest| then construct v
                Vec3 v;
                double absolute_min = min(abs(u.x), min(abs(u.y), abs(u.z)));
                if(absolute_min == abs(u.x)){
                    v = Vec3(0, -u.z, u.y,0);

                }
                else if(absolute_min == abs(u.y)){
                    v = Vec3(-u.z, 0, u.x,-1);
                }
                else{
                    v = Vec3(-u.y, u.x, 0,-1);
                }
                Vec3 w = crossProductVec3(u, v);
                v = normalizeVec3(v);
                w = normalizeVec3(w);
                double M_normal_arr[4][4] = {{u.x,u.y,u.z,0},
                                             {v.x,v.y,v.z,0},
                                             {w.x,w.y,w.z,0},
                                             {0,0,0,1}};
                Matrix4 M(M_normal_arr);
                double M_transpose_arr[4][4] = {{u.x,v.x,w.x,0},
                                                {u.y,v.y,w.y,0},
                                                {u.z,v.z,w.z,0},
                                                {0,0,0,1}};
                Matrix4 M_T(M_transpose_arr);
                double sin_angle = std::sin(r->angle * 3.14159/180);      // convert radian
                double cos_angle = std::cos(r->angle * 3.14159/180);      // convert radian
                double R_arr[4][4] = {{1,0,0,0},
                                      {0,cos_angle,-sin_angle,0},
                                      {0,sin_angle,cos_angle,0},
                                      {0,0,0,1}};
                Matrix4 R(R_arr);
                // Mt*R*M
                Matrix4 R_axis = multiplyMatrixWithMatrix(M_T, multiplyMatrixWithMatrix(R, M));
                retMatrix = multiplyMatrixWithMatrix(R_axis, retMatrix);
                break;
            }
            default: {
                std::cout << "Wrong translation type" << std::endl;
                break;
            }
        }

    }
    return retMatrix;
}

void Scene::rasLine(Vec4 v1, Vec4 v2, Color color1, Color color2, int nx, int ny) {

    double dx = v2.x - v1.x;
    double dy = v2.y - v1.y;


    // m between 0-1
    if(abs(dy) < abs(dx)){
        // should start with smaller x
        if(v2.x < v1.x){
            Vec4 temp_vec = v1;
            v1 = v2;
            v2 = temp_vec;
            Color temp_color = color1;
            color1 = color2;
            color2 = temp_color;
            dx = -dx;
            dy = -dy;
        }
        double dr = (color2.r - color1.r) / dx, dg = (color2.g - color1.g) / dx, db = (color2.b - color1.b) / dx;

        int loop_val = v2.y >= (int)v1.y ? 1:-1;
        double error = -2 * dy + loop_val * dx;
        int y = (int) v1.y;
        for(int x = (int) v1.x ; x < v2.x ; x++){
            error -= dy;

            // NE chosen
            if(error * loop_val < 0){
                error += loop_val * dx ;
                y += loop_val;
            }
            color1.r+=dr, color1.g+=dg, color1.b+=db;

            // control y value
            if(y < 0)
                y = 0;
            if(y >= ny)
                y = ny;

            image[x][y] = color1;
        }
    }
    // m bigger than 1 swap x - y
    else{
        if(v2.y < v1.y){
            Vec4 temp_vec = v1;
            v1 = v2;
            v2 = temp_vec;
            Color temp_color = color1;
            color1 = color2;
            color2 = temp_color;
            dx = -dx;
            dy = -dy;
        }
        double dr = (color2.r - color1.r) / dy, dg = (color2.g - color1.g) / dy, db = (color2.b - color1.b) / dy;

        int loop_val = v2.x >= v1.x ? 1:-1;
        double error = loop_val * dy - 2 * dx;
        int x = (int) v1.x;

        for(int y = (int) v1.y ; y < (int)v2.y ; y++){
            error -= dx;

            // NE chosen
            if(error * loop_val < 0){
                error += loop_val * dy;
                x += loop_val;
            }
            color1.r+=dr, color1.g+=dg, color1.b+=db;

            // control x value
            if(x < 0)
                x = 0;
            if(x >= nx)
                x = nx;

            image[x][y] = color1;
        }
    }

}

auto line_eq = []  (const Vec4& v0, const Vec4& v1, double x, double y)
        -> double {
    double x0 = v0.x;
    double y0 = v0.y;
    double x1 = v1.x;
    double y1 = v1.y;
    return x * (y0 - y1) + y * (x1 - x0) + (x0 * y1) - y0 * x1;
};

void Scene::rasTriangle(Vec4 &v0, Vec4 &v1, Vec4 &v2, const Color& c0, const Color& c1, const Color& c2, int nx, int ny) {

    double alpha, beta, gamma, r, g, b;
    double f01, f02, f12, f01v2, f02v1, f12v0;
    // calculate box coordinates
    int box_x_min = min(max((int) min(min(v0.x,v1.x),v2.x),0), nx-1);
    int box_y_min = min(max((int) min(min(v0.y,v1.y),v2.y),0), ny-1);
    int box_x_max = max(min((int) max(max(v0.x, v1.x), v2.x), nx-1),0);
    int box_y_max = max(min((int) max(max(v0.y, v1.y), v2.y), ny-1),0);

    for(int y = box_y_min ; y <= box_y_max; y++){
        for(int x = box_x_min ; x <= box_x_max ; x++){
            // barycentric coordinate calculation
            f01 = line_eq(v0,v1,x,y);            f01v2 = line_eq(v0,v1,v2.x,v2.y);
            f02 = line_eq(v0,v2,x,y);        f02v1 = line_eq(v0,v2,v1.x,v1.y);
            f12 = line_eq(v1,v2,x,y);    f12v0 = line_eq(v1,v2,v0.x,v0.y);

            alpha = f12 / f12v0;
            beta = f02 / f02v1;
            gamma = f01 / f01v2;
            // check if inside the triangle using barycentric coordinates
            if((alpha >= 0 && beta >= 0) && gamma >= 0){
                r = c0.r * alpha + c1.r * beta + c2.r * gamma;
                g = c0.g * alpha + c1.g * beta + c2.g * gamma;
                b = c0.b * alpha + c1.b * beta + c2.b * gamma;

                image[x][y] = Color(r, g, b);
            }
        }
    }

}
/*
bool visible(den, num, tE, tL):
if (den > 0): // potentially entering
t = num / den;
if (t > tL):
return false;
if (t > tE)
tE = t;
else if (den < 0): // potentially leaving
t = num / den;
if (t < tE):
return false;
if (t < tL)
tL = t;
else if num > 0: // line parallel to edge
return false;
return true;
 * */
bool Scene::visible(double den, double num, double &tE, double &tL) {
    double t;
    if (den > 0) { // potentially entering
        t = num / den;
        if (t > tL)
            return false;
        if (t > tE)
            tE = t;
    }
    else if (den < 0) { // potentially leaving
        t = num / den;
        if (t < tE)
            return false;
        if (t < tL)
            tL = t;
    }
    else if (num > 0){ // line parallel to edge
        return false;
    }
    return true;
}

/*
tE = 0; tL = 1;
visible = false;
if visible(dx, xmin – x0,tE, tL): // left
if visible (-dx, x0 – xmax,tE,tL): // right
if visible (dy, ymin – y0,tE,tL): // bottom
if visible (-dy, y0 – ymax,tE,tL): // top
if visible (dz, zmin – z0,tE,tL): // front
if visible (-dz, z0 – zmax,tE,tL): // back
visible = true;
if (tL < 1):
x1 = x0 + dxtL; y1 = y0 + dytL; z1 = z0 + dztL;
if (tE > 0):
x0 = x0 + dxtE; y0 = y0 + dytE; z0 = z0 + dztE;
 * */

bool Scene::liangBarsky(Vec4 &v0, Vec4 &v1, Color *c0, Color *c1) {

    double dx,dy,dz,d_r,d_g,d_b,x_min,y_min,z_min,x_max,y_max,z_max,tE,tL;
    dx = v1.x - v0.x; dy = v1.y - v0.y; dz = v1.z - v0.z;
    d_r = c1->r - c0->r; d_g = c1->g - c0->g; d_b = c1->b - c0->b;
    x_min = -1; y_min = -1; z_min = -1; x_max = 1; y_max = 1; z_max = 1;

    tE = 0; tL = 1;
    //check visibility on edges
    if(!visible(dx, x_min-v0.x,tE,tL)){return false;}
    if(!visible(-dx, v0.x-x_max,tE,tL)){return false;}
    if(!visible(dy, y_min-v0.y,tE,tL)){return false;}
    if(!visible(-dy, v0.y-y_max,tE,tL)){return false;}
    if(!visible(dz, z_min-v0.z,tE,tL)){return false;}
    if(!visible(-dz, v0.z-z_max,tE,tL)){return false;}

    if (tL < 1) {
        v1.x = v0.x + dx * tL;
        v1.y = v0.y + dy * tL;
        v1.z = v0.z + dz * tL;
        //color
        c1->r = c0->r + d_r * tL; c1->g = c0->g + d_g * tL; c1->b = c0->b + d_b * tL;
    }
    if (tE > 0) {
        v0.x = v0.x + dx * tE;
        v0.y = v0.y + dy * tE;
        v0.z = v0.z + dz * tE;
        //color
        c0->r = c0->r + d_r * tE; c0->g = c0->g + d_g * tE; c0->b = c0->b + d_b * tE;
    }
    return true;
}


bool Scene::backFaceCulling(const Vec4& v0, const Vec4& v1, const Vec4& v2, const Vec3& pos) {
    Vec3 edge1(v1.x-v0.x,v1.y-v0.y,v1.z-v0.z,-1);
    Vec3 edge2(v2.x-v0.x,v2.y-v0.y,v2.z-v0.z,-1);
    Vec3 normal = crossProductVec3(edge1, edge2);
    Vec3 eye_to_view(subtractVec3(Vec3(v0.x,v0.y,v0.z,-1), pos));
    return dotProductVec3(normal, eye_to_view) < 0;
}



/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
    const char *str;
    XMLDocument xmlDoc;
    XMLElement *pElement;

    xmlDoc.LoadFile(xmlPath);

    XMLNode *pRoot = xmlDoc.FirstChild();

    // read background color
    pElement = pRoot->FirstChildElement("BackgroundColor");
    str = pElement->GetText();
    sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

    // read culling
    pElement = pRoot->FirstChildElement("Culling");
    if (pElement != NULL) {
        str = pElement->GetText();

        if (strcmp(str, "enabled") == 0) {
            cullingEnabled = true;
        }
        else {
            cullingEnabled = false;
        }
    }

    // read cameras
    pElement = pRoot->FirstChildElement("Cameras");
    XMLElement *pCamera = pElement->FirstChildElement("Camera");
    XMLElement *camElement;
    while (pCamera != NULL)
    {
        Camera *cam = new Camera();

        pCamera->QueryIntAttribute("id", &cam->cameraId);

        // read projection type
        str = pCamera->Attribute("type");

        if (strcmp(str, "orthographic") == 0) {
            cam->projectionType = 0;
        }
        else {
            cam->projectionType = 1;
        }

        camElement = pCamera->FirstChildElement("Position");
        str = camElement->GetText();
        sscanf(str, "%lf %lf %lf", &cam->pos.x, &cam->pos.y, &cam->pos.z);

        camElement = pCamera->FirstChildElement("Gaze");
        str = camElement->GetText();
        sscanf(str, "%lf %lf %lf", &cam->gaze.x, &cam->gaze.y, &cam->gaze.z);

        camElement = pCamera->FirstChildElement("Up");
        str = camElement->GetText();
        sscanf(str, "%lf %lf %lf", &cam->v.x, &cam->v.y, &cam->v.z);

        cam->gaze = normalizeVec3(cam->gaze);
        cam->u = crossProductVec3(cam->gaze, cam->v);
        cam->u = normalizeVec3(cam->u);

        cam->w = inverseVec3(cam->gaze);
        cam->v = crossProductVec3(cam->u, cam->gaze);
        cam->v = normalizeVec3(cam->v);

        camElement = pCamera->FirstChildElement("ImagePlane");
        str = camElement->GetText();
        sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
               &cam->left, &cam->right, &cam->bottom, &cam->top,
               &cam->near, &cam->far, &cam->horRes, &cam->verRes);

        camElement = pCamera->FirstChildElement("OutputName");
        str = camElement->GetText();
        cam->outputFileName = string(str);

        cameras.push_back(cam);

        pCamera = pCamera->NextSiblingElement("Camera");
    }

    // read vertices
    pElement = pRoot->FirstChildElement("Vertices");
    XMLElement *pVertex = pElement->FirstChildElement("Vertex");
    int vertexId = 1;

    while (pVertex != NULL)
    {
        Vec3 *vertex = new Vec3();
        Color *color = new Color();

        vertex->colorId = vertexId;

        str = pVertex->Attribute("position");
        sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

        str = pVertex->Attribute("color");
        sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

        vertices.push_back(vertex);
        colorsOfVertices.push_back(color);

        pVertex = pVertex->NextSiblingElement("Vertex");

        vertexId++;
    }

    // read translations
    pElement = pRoot->FirstChildElement("Translations");
    XMLElement *pTranslation = pElement->FirstChildElement("Translation");
    while (pTranslation != NULL)
    {
        Translation *translation = new Translation();

        pTranslation->QueryIntAttribute("id", &translation->translationId);

        str = pTranslation->Attribute("value");
        sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

        translations.push_back(translation);

        pTranslation = pTranslation->NextSiblingElement("Translation");
    }

    // read scalings
    pElement = pRoot->FirstChildElement("Scalings");
    XMLElement *pScaling = pElement->FirstChildElement("Scaling");
    while (pScaling != NULL)
    {
        Scaling *scaling = new Scaling();

        pScaling->QueryIntAttribute("id", &scaling->scalingId);
        str = pScaling->Attribute("value");
        sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

        scalings.push_back(scaling);

        pScaling = pScaling->NextSiblingElement("Scaling");
    }

    // read rotations
    pElement = pRoot->FirstChildElement("Rotations");
    XMLElement *pRotation = pElement->FirstChildElement("Rotation");
    while (pRotation != NULL)
    {
        Rotation *rotation = new Rotation();

        pRotation->QueryIntAttribute("id", &rotation->rotationId);
        str = pRotation->Attribute("value");
        sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

        rotations.push_back(rotation);

        pRotation = pRotation->NextSiblingElement("Rotation");
    }

    // read meshes
    pElement = pRoot->FirstChildElement("Meshes");

    XMLElement *pMesh = pElement->FirstChildElement("Mesh");
    XMLElement *meshElement;
    while (pMesh != NULL)
    {
        Mesh *mesh = new Mesh();

        pMesh->QueryIntAttribute("id", &mesh->meshId);

        // read projection type
        str = pMesh->Attribute("type");

        if (strcmp(str, "wireframe") == 0) {
            mesh->type = 0;
        }
        else {
            mesh->type = 1;
        }

        // read mesh transformations
        XMLElement *pTransformations = pMesh->FirstChildElement("Transformations");
        XMLElement *pTransformation = pTransformations->FirstChildElement("Transformation");

        while (pTransformation != NULL)
        {
            char transformationType;
            int transformationId;

            str = pTransformation->GetText();
            sscanf(str, "%c %d", &transformationType, &transformationId);

            mesh->transformationTypes.push_back(transformationType);
            mesh->transformationIds.push_back(transformationId);

            pTransformation = pTransformation->NextSiblingElement("Transformation");
        }

        mesh->numberOfTransformations = mesh->transformationIds.size();

        // read mesh faces
        char *row;
        char *clone_str;
        int v1, v2, v3;
        XMLElement *pFaces = pMesh->FirstChildElement("Faces");
        str = pFaces->GetText();
        clone_str = strdup(str);

        row = strtok(clone_str, "\n");
        while (row != NULL)
        {
            int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

            if (result != EOF) {
                mesh->triangles.push_back(Triangle(v1, v2, v3));
            }
            row = strtok(NULL, "\n");
        }
        mesh->numberOfTriangles = mesh->triangles.size();
        meshes.push_back(mesh);

        pMesh = pMesh->NextSiblingElement("Mesh");
    }
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
    if (this->image.empty())
    {
        for (int i = 0; i < camera->horRes; i++)
        {
            vector<Color> rowOfColors;

            for (int j = 0; j < camera->verRes; j++)
            {
                rowOfColors.push_back(this->backgroundColor);
            }

            this->image.push_back(rowOfColors);
        }
    }
    else
    {
        for (int i = 0; i < camera->horRes; i++)
        {
            for (int j = 0; j < camera->verRes; j++)
            {
                this->image[i][j].r = this->backgroundColor.r;
                this->image[i][j].g = this->backgroundColor.g;
                this->image[i][j].b = this->backgroundColor.b;
            }
        }
    }
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
    if (value >= 255.0)
        return 255;
    if (value <= 0.0)
        return 0;
    return (int)(value);
}

/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
    ofstream fout;

    fout.open(camera->outputFileName.c_str());

    fout << "P3" << endl;
    fout << "# " << camera->outputFileName << endl;
    fout << camera->horRes << " " << camera->verRes << endl;
    fout << "255" << endl;

    for (int j = camera->verRes - 1; j >= 0; j--)
    {
        for (int i = 0; i < camera->horRes; i++)
        {
            fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
                 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
                 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
        }
        fout << endl;
    }
    fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
	os_type == 1 		-> Ubuntu
	os_type == 2 		-> Windows
	os_type == other	-> No conversion
*/
void Scene::convertPPMToPNG(string ppmFileName, int osType)
{
    string command;

    // call command on Ubuntu
    if (osType == 1)
    {
        command = "convert " + ppmFileName + " " + ppmFileName + ".png";
        system(command.c_str());
    }

        // call command on Windows
    else if (osType == 2)
    {
        command = "magick convert " + ppmFileName + " " + ppmFileName + ".png";
        system(command.c_str());
    }

        // default action - don't do conversion
    else
    {
    }
}