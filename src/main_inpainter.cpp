#include <iostream>
#include "CImg.h"
#include "diff.h"
using namespace cimg_library;
using namespace std;
#include "inpainter.h"
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
      FloatImg mask(img);

      mask.fill(1.0);
      cimg_forXYC(mask,x,y,c){
          if(img(x,y,0)<4 && img(x,y,1)<4 && img(x,y,2)<4){
              mask(x,y,c) = 0.0;
          }
      }
      //img = img.get_RGBtoYCbCr().get_channel(0);
      //mask = mask.get_channel(0);
      //cout << "MSE : " << img.MSE(comp) << endl;
      float l0, l1, l2, a, b;
      cin >> l0 >> l1 >> l2 >> a >> b;
      //0.001f,0.001f,0.0001f,0.01f,0.01f
      BregmanInpainter solver(img,mask,l0,l1,l2,a,b);
      solver.solve();
      solver.save("out.png");


      cout << "MSE : " << img.MSE(solver.get_reconstructed_image()) << endl;

  }catch(exception const& e){
      cerr << "ERROR : " << e.what() << endl;
  }
  return 0;
}
