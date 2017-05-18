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
      string filename;
      cin >> filename;
      filename = "img/input/"+filename;
      DiffImg img(filename.c_str());
      FloatImg comp(filename.c_str());
      img.noise(30).normalize(0,255);
      img = img.get_RGBtoYCbCr().get_channel(0);
      img.save("lena_bruit.bmp");
      //cout << "MSE : " << img.MSE(comp) << endl;
      float l1, l2, a, b;
      cin >> l1 >> l2 >> a >> b;
      // 0.06 0.05 0.6 0.5
      BregmanDenoiser solver(img,l1,l2,a,b);
      solver.solve();
      solver.save("lena.bmp");


      cout << "MSE : " << img.MSE(solver.get_reconstructed_image()) << endl;


  }catch(exception const& e){
      cerr << "ERROR : " << e.what() << endl;
  }
  return 0;
}
