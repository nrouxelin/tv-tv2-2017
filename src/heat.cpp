/**void heat() {
  CImg<> img("img/input/lena.bmp");  // Image originale
  CImg<> tmp(img);
  CImgDisplay disp(img,"Equation de la chaleur",0,false,false);
  //img.noise(50);
  int t = 0;//Temps
  float dt = 0.1;//Pas de temps
  float c = 2;
  const float white[] = { 255,255,255 };
  CImg_3x3(I,float); //Voisinnage 3x3
  while(!disp.is_closed()){
      cimg_forC(img,k){ //Boucle sur les couleurs
          cimg_for3x3(img,x,y,0,k,I,float){ //Boucle sur les voisinnages
              tmp(x,y,k) = Ipc+Inc+Icp+Icn-4*Icc;
              //tmp(x,y,k) = Inc-Ipc-Icn+Icp;
          }
      }
      img += dt*tmp;
      cimg_library::CImg<>(img).draw_text(2,2,"iter = %d",white,0,1,13,t).display(disp.wait(100));
      cimg_library::CImg<>(img).draw_text(2,2,"iter = %d",white,0,1,13,t).save("img/output/deconv.bmp",t);
      t++;
  }
}
**/
