#include <iostream>
#include "CImg.h"
#include "diff.h"
#include "denoiser.h"
using namespace cimg_library;
using namespace std;
#include <fenv.h>

int main(){
  try{


      feenableexcept(FE_INVALID | FE_OVERFLOW);
      //Original image
      char* filename = "img/input/lena.bmp";
      DiffImg img(filename);
      FloatImg comp(filename);
      img.noise(15,3);//.normalize(0,255);
      img = img.get_RGBtoYCbCr().get_channel(0);
      img.save("lena_bruit.bmp");
      //cout << "MSE : " << img.MSE(comp) << endl;
      BregmanDenoiser solver(img,0.06f,0.05f,0.6f,0.5f);
      solver.solve();
      solver.save("lena.bmp");


      cout << "MSE : " << img.MSE(solver.get_reconstructed_image()) << endl;


  }catch(exception const& e){
      cerr << "ERROR : " << e.what() << endl;
  }
  return 0;
}
