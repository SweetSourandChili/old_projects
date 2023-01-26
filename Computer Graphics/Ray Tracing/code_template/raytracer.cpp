#include <iostream>
#include "parser.h"
#include "render.h"
#include <ctime>

using namespace std;
using namespace parser;


int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;

    scene.loadFromXml(argv[1]);

    // The code below creates a test pattern and writes
    // it to a PPM file to demonstrate the usage of the
    // ppm_write function.
    //
    // Normally, you would be running your ray tracing
    // code here to produce the desired image.



    auto* r = new Renderer(&scene);
    int c = 1;
    for (auto cam : r->cameras){
        const clock_t begin_time = clock();
        r->Render(cam);
        printf("Image : %s, camNo: Cam%d -- Time : %f\n",argv[1],c,float( clock () - begin_time )/ CLOCKS_PER_SEC);
        c++;
    }
    delete r;
}
