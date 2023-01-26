#include "Camera.h"
#include <string>
#include <iostream>
#include <iomanip>

using namespace std;

Camera::Camera() {}

Camera::Camera(int cameraId,
               int projectionType,
               Vec3 pos, Vec3 gaze,
               Vec3 u, Vec3 v, Vec3 w,
               double left, double right, double bottom, double top,
               double near, double far,
               int horRes, int verRes,
               string outputFileName)
{

    this->cameraId = cameraId;
    this->projectionType = projectionType;
    this->pos = pos;
    this->gaze = gaze;
    this->u = u;
    this->v = v;
    this->w = w;
    this->left = left;
    this->right = right;
    this->bottom = bottom;
    this->top = top;
    this->near = near;
    this->far = far;
    this->horRes = horRes;
    this->verRes = verRes;
    this->outputFileName = outputFileName;
}

Camera::Camera(const Camera &other)
{
    this->cameraId = other.cameraId;
    this->projectionType = other.projectionType;
    this->pos = other.pos;
    this->gaze = other.gaze;
    this->u = other.u;
    this->v = other.v;
    this->w = other.w;
    this->left = other.left;
    this->right = other.right;
    this->bottom = other.bottom;
    this->top = other.top;
    this->near = other.near;
    this->far = other.far;
    this->horRes = other.horRes;
    this->verRes = other.verRes;
    this->outputFileName = other.outputFileName;
}

ostream &operator<<(ostream &os, const Camera &c)
{
    const char *camType = c.projectionType ? "perspective" : "orthographic";

    os << fixed << setprecision(6) << "Camera " << c.cameraId << " (" << camType << ") => pos: " << c.pos << " gaze: " << c.gaze << endl
       << "\tu: " << c.u << " v: " << c.v << " w: " << c.w << endl
       << fixed << setprecision(3) << "\tleft: " << c.left << " right: " << c.right << " bottom: " << c.bottom << " top: " << c.top << endl
       << "\tnear: " << c.near << " far: " << c.far << " resolutions: " << c.horRes << "x" << c.verRes << " fileName: " << c.outputFileName;

    return os;
}

Matrix4 Camera::eyeTransformationMatrix() {
    // res = RT

    double R_arr[4][4] = {{this->u.x,this->u.y,this->u.z,0},
                          {this->v.x,this->v.y,this->v.z,0},
                          {this->w.x,this->w.y,this->w.z,0},
                          {0,0,0,1}};
    double T_arr[4][4] = {{1, 0, 0, -(this->pos.x)},
                          {0, 1, 0, -(this->pos.y)},
                          {0, 0, 1, -(this->pos.z)},
                          {0, 0, 0, 1}};
    Matrix4 R(R_arr);
    Matrix4 T(T_arr);
    return multiplyMatrixWithMatrix(R, T);
}

Matrix4 Camera::projectionMatrix() {

    double O_arr[4][4] = {{2/(right-left), 0, 0, -((right+left)/(right-left))},
                          {0, 2/(top-bottom), 0, -((top+bottom)/(top-bottom))},
                          {0, 0, -2/(far-near), -(far+near)/(far-near)},
                          {0, 0, 0, 1}};
    Matrix4 M_orth(O_arr);
    if(projectionType == 0){        // orthographic
        return M_orth;
    }
    else{   // perspective
        double P2O_arr[4][4] = {{near, 0, 0, 0},
                        {0, near, 0, 0},
                        {0, 0, far+near, far*near},
                        {0, 0, -1, 0}};
        Matrix4 M_p2o(P2O_arr);
        Matrix4 M_pers = multiplyMatrixWithMatrix(M_orth, M_p2o);
        return M_pers;
    }
}

Matrix4 Camera::viewportMatrix(int x_min, int y_min) {
    double V_arr[4][4] = {{(horRes-1)/2.0, 0, 0, (horRes-1)/2.0 + x_min},
                          {0, (verRes-1)/2.0, 0, (verRes-1)/2.0 + y_min},
                          {0, 0, 0.5, 0.5},
                          {0, 0, 0, 1}};

    return {V_arr};
}
