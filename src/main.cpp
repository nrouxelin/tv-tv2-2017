#include <iostream>
#include "CImg.h"
#include "diff.h"
#include "bregman.h"
using namespace cimg_library;
using namespace std;



int main(){
  try{
      //Original image
      char* filename = "img/input/lena.bmp";
      DiffImg img(filename);
      FloatImg comp(filename);
      img.noise(10).normalize(0,255);
      cout << "MSE : " << img.MSE(comp) << endl;
      BregmanSolver solver(img,0.06f,0.05f,0.6f,0.5f);
      solver.solve();
      solver.save("lena.bmp");


      cout << "MSE : " << img.MSE(solver.get_reconstructed_image()) << endl;


  }catch(exception const& e){
      cerr << "ERROR : " << e.what() << endl;
  }
  return 0;
}
