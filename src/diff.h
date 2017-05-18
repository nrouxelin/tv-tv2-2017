#ifndef DIFF_H
#define DIFF_H 0

#include "CImg.h"
namespace cil = cimg_library;
typedef cil::CImg<float> FloatImg;
typedef cil::CImgList<float> FloatImgList;
#include <stdexcept>

class DiffImg: public FloatImg{
public:
    DiffImg(FloatImg img);
    DiffImg(const char *const filename);
    float fdx(int x, int y, int c);
    float fdy(int x, int y, int c);
    float bdx(int x, int y, int c);
    float bdy(int x, int y, int c);
    float dxx(int x, int y, int c);
    float dyy(int x, int y, int c);
    float dxy(int x, int y, int c);
    float bdxy(int x, int y, int c);
    float at(int x,int y, int c);



};
#endif
