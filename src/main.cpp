#include <iostream>
#include "CImg.h"
#include "diff.h"
#include "denoiser.h"
#include "inpainter.h"
using namespace cimg_library;
using namespace std;
#include <fenv.h>




int main(){
  try{
      feenableexcept(FE_INVALID | FE_OVERFLOW);
      //Original image
      char* filename = "img/input/lena_masque.bmp";
      DiffImg img(filename);
      FloatImg comp(filename);
      FloatImg mask(img);

      mask.fill(1.0);
      cimg_forXY(mask,x,y){
          if(img(x,y,0)==0 && img(x,y,1)==0 && img(x,y,2)==0){
              mask(x,y,0) = 0.0;
          }
      }
      img = img.get_RGBtoYCbCr().get_channel(0);
      mask = mask.get_channel(0);
      //cout << "MSE : " << img.MSE(comp) << endl;
      BregmanInpainter solver(img,mask,0.001f,0.001f,0.01f,0.01f,0.01f);
      solver.solve();
      solver.save("lena.bmp");


      cout << "MSE : " << img.MSE(solver.get_reconstructed_image()) << endl;


  }catch(exception const& e){
      cerr << "ERROR : " << e.what() << endl;
  }
  return 0;
}




/**int main(){
  try{


      feenableexcept(FE_INVALID | FE_OVERFLOW);
      //Original image
      char* filename = "img/input/lena.bmp";
      DiffImg img(filename);
      FloatImg comp(filename);
      img.noise(30);//.normalize(0,255);
      img.save("lena_bruit.bmp");
      //img = img.get_RGBtoYCbCr().get_channel(0);
      //cout << "MSE : " << img.MSE(comp) << endl;
      BregmanDenoiser solver(img,0.06f,0.05f,0.6f,0.5f);
      solver.solve();
      solver.save("lena.bmp");


      cout << "MSE : " << img.MSE(solver.get_reconstructed_image()) << endl;


  }catch(exception const& e){
      cerr << "ERROR : " << e.what() << endl;
  }
  return 0;
}**/
