#include "diff.h"
using namespace cimg_library;
using namespace std;


DiffImg::DiffImg(const char *const filename):FloatImg(filename){};
DiffImg::DiffImg(FloatImg I):FloatImg(I){};

float DiffImg::at(int x, int y, int c){
    if(spectrum() && c>=0){
        int w = width();
        int h = height();

        if(x>=w){
            x = x-w;
        }else if(x<0){
            x = w+x;
        }

        if(y>=h){
            y = y-h;
        }else if(y<0){
            y = h+y;
        }

        return (*this)(x,y,c);
    }else{
        throw domain_error("This channel does not exist.");
    }
}

float DiffImg::fdx(int x, int y, int c){
    return at(x+1,y,c)-at(x,y,c);
}

float DiffImg::fdy(int x, int y, int c){
    return at(x,y+1,c)-at(x,y,c);
}

float DiffImg::bdx(int x, int y, int c){
    return at(x,y,c)-at(x-1,y,c);
}

float DiffImg::bdy(int x, int y, int c){
    return at(x,y,c)-at(x,y-1,c);
}

float DiffImg::dxx(int x, int y, int c){
    return at(x-1,y,c)-2.0*at(x,y,c)+at(x+1,y,c);
}

float DiffImg::dyy(int x, int y, int c){
    return at(x,y-1,c)-2.0*at(x,y,c)+at(x,y+1,c);
}

float DiffImg::dxy(int x, int y, int c){
    return at(x,y,c)-at(x+1,y,c)-at(x,y+1,c)+at(x+1,y+1,c);
}
