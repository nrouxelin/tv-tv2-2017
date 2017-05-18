#include <iostream>
#include "CImg.h"
using namespace cimg_library;
using namespace std;

int main(){
  string filename("img/input/lena.bmp");//Nom du fichier
  CImg<float> ud(filename.c_str());//Création de l'image
  ud.noise(30);//Ajout de bruit
  ud = ud.get_RGBtoYCbCr().get_channel(0);//RGB vers niveaux de gris
  CImg<float> u(ud);//On copie  l'image observée
  ud.save("img/output/h1_start.png");
  CImgDisplay disp(u,"Regularisation H^1",0,false,false);

  //Parametres
  float a  = 5.0;
  float dt = 0.01;

  //Dimensions de l'Image
  int M = u.height();//Largeur
  int N = u.width();//Longueur

  const float white[] = { 255,255,255 };

  int t = 0;

  while(!disp.is_closed()){
      t++;

      cimg_forXYC(u,x,y,c){
          //Interieur de l'image
          if(x>1 && x<N-1 && y>1 && y<M-1){
              u(x,y,c) += dt*(ud(x,y,c)-u(x,y,c))
                +(a*dt)*(u(x+1,y,c)-4.0*u(x,y,c)+u(x-1,y,c)+u(x,y+1,c)+u(x,y-1,c));
          }else{
            //Conditions de Neumann
            if(x==0)
                u(x,y,c) = u(1,y,c);
            else if(x==N-1)
                u(x,y,c) = u(N-2,y,c);
            if(y==0)
                u(x,y,c) = u(x,1,c);
            else if(y==M-1)
                u(x,y,c) = u(x,M-2,c);
          }
      }

      cimg_library::CImg<>(u).draw_text(2,2,"iter = %d",white,0,1,13,t).display(disp.wait(100));
  }
  cimg_library::CImg<>(u).save("img/output/h1.png",t);

  return 0;
}
