#ifndef _SCENE_H_
#define _SCENE_H_

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "Camera.h"
#include "Color.h"
#include "Mesh.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Triangle.h"
#include "Vec3.h"
#include "Vec4.h"

using namespace std;

class Scene
{
public:
	Color backgroundColor;
	bool cullingEnabled;

	vector< vector<Color> > image;
	vector< Camera* > cameras;
	vector< Vec3* > vertices;
	vector< Color* > colorsOfVertices;
	vector< Scaling* > scalings;
	vector< Rotation* > rotations;
	vector< Translation* > translations;
	vector< Mesh* > meshes;

	Scene(const char *xmlPath);

	void initializeImage(Camera* camera);
	void forwardRenderingPipeline(Camera* camera);
	int makeBetweenZeroAnd255(double value);
	void writeImageToPPMFile(Camera* camera);
	void convertPPMToPNG(string ppmFileName, int osType);

    Matrix4 meshTranslationMatrix(Mesh const& mesh);

    void rasLine(Vec4 vect1, Vec4 vect2, Color color1, Color color2, int nx, int ny);
    void rasTriangle(Vec4 &v0, Vec4 &v1, Vec4 &v2, const Color& c0, const Color& c1, const Color& c2, int nx, int ny);
    bool visible(double den, double num, double &tE, double &tL);
    bool liangBarsky(Vec4& v0, Vec4& v1, Color* c0, Color* c1);
    bool backFaceCulling(const Vec4& v0, const Vec4& v1, const Vec4& v2, const Vec3& pos);
};

#endif
