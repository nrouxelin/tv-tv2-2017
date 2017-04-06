#include "bregman.h"
using namespace std;
using namespace cimg_library;

BregmanSolver::BregmanSolver(DiffImg f, float a, float b, float l1, float l2):
m_f(f),m_u(f),m_rfd(f),m_ifd(f),m_l1(l1),m_l2(l2),m_a(a),m_b(b){
    //Creating v:=\nabla u
      DiffImg vx(f), vy(f);
      //Creating w:=\nabla^2 u
      DiffImg wxx(f), wyy(f), wxy(f);

      //Loop over x,y and c (color channels)
      cimg_forXYC(f,x,y,c){
          //filling v
          vx(x,y,c) = f.fdx(x,y,c);
          vy(x,y,c) = f.fdy(x,y,c);

          //filling w
          wxx(x,y,c) = f.dxx(x,y,c);
          wxy(x,y,c) = f.dxy(x,y,c);
          wyy(x,y,c) = f.dyy(x,y,c);
      }

      //Allocating the members m_v and m_w
      m_v.push_back(vx);
      m_v.push_back(vy);

      m_w.push_back(wxx);
      m_w.push_back(wyy);
      m_w.push_back(wxy);

      //-----WIP-----
      DiffImg b1_1(f), b1_2(f), b2_1(f), b2_2(f), b2_3(f);
      b1_1.fill(0.01);
      b1_2.fill(0.01);
      m_b1.push_back(b1_1);
      m_b1.push_back(b1_2);

      b2_1.fill(0.01);
      b2_2.fill(0.01);
      b2_3.fill(0.01);
      m_b2.push_back(b2_1);
      m_b2.push_back(b2_2);
      m_b2.push_back(b2_3);

      m_rfd.fill(0.0);
      m_ifd.fill(0.0);


      m_iter = 0;

}


void BregmanSolver::solve_subproblem2(){
    FloatImg s1(m_u), s2(m_u);
    //Filling s1 and s2
    cimg_forXYC(m_u,x,y,c){
        s1(x,y,c) = m_b1[0](x,y,c)+m_u.fdx(x,y,c);
        s2(x,y,c) = m_b1[1](x,y,c)+m_u.fdy(x,y,c);
    }
    //Computing v^{n+1}
    float norm, m;
    cimg_forXYC(m_u,x,y,c){
        norm          = sqrt(s1(x,y,c)*s1(x,y,c)+s2(x,y,c)*s2(x,y,c));
        m             = max(norm-m_a/m_l1,0.0f);
        m_v[0](x,y,c) = m*(s1(x,y,c)/norm);
        m_v[1](x,y,c) = m*(s2(x,y,c)/norm);
    }
}

void BregmanSolver::solve_subproblem3(){
    FloatImg t1(m_u), t2(m_u), t3(m_u);
    //Filling t1, t2 and t3
    cimg_forXYC(m_u,x,y,c){
        t1(x,y,c) = m_b2[0](x,y,c)+m_u.dxx(x,y,c);
        t2(x,y,c) = m_b2[1](x,y,c)+m_u.dyy(x,y,c);
        t3(x,y,c) = m_b2[2](x,y,c)+m_u.dxy(x,y,c);
    }

    //Computing w^{n+1}
    float norm, m;
    cimg_forXYC(m_u,x,y,c){
        norm = t1(x,y,c)*t1(x,y,c)+t2(x,y,c)*t2(x,y,c)+2.0f*t3(x,y,c)*t3(x,y,c);
        norm = sqrt(norm);

        m             = max(norm-m_b/m_l2,0.0f);
        m_w[0](x,y,c) = m*(t1(x,y,c)/norm);
        m_w[1](x,y,c) = m*(t2(x,y,c)/norm);
        m_w[2](x,y,c) = m*(t3(x,y,c)/norm);
    }
}

void BregmanSolver::update_b(){
    cimg_forXYC(m_u,x,y,c){
        //Updating b_1
        m_b1[0](x,y,c) = m_b1[0](x,y,c)+m_u.fdx(x,y,c)-m_v[0](x,y,c);
        m_b1[1](x,y,c) = m_b1[1](x,y,c)+m_u.fdy(x,y,c)-m_v[1](x,y,c);

        //Updating b_2
        m_b2[0](x,y,c) = m_b2[0](x,y,c)+m_u.dxx(x,y,c)-m_w[0](x,y,c);
        m_b2[1](x,y,c) = m_b2[1](x,y,c)+m_u.dyy(x,y,c)-m_w[1](x,y,c);
        m_b2[2](x,y,c) = m_b2[2](x,y,c)+m_u.dxy(x,y,c)-m_w[2](x,y,c);
    }
}

void BregmanSolver::solve_subproblem1_GS(){
    //Compute right side of the equation
    FloatImg rs(m_u);
    cimg_forXYC(rs,x,y,c){
        rs(x,y,c) = m_f(x,y,c)+m_l1*(m_b1[0].bdx(x,y,c)-m_v[0].bdx(x,y,c)+m_b1[1].bdy(x,y,c)-m_v[1].bdy(x,y,c))
                    -m_l2*(m_b2[0].dxx(x,y,c)-m_w[0].dxx(x,y,c)+m_b2[1].dyy(x,y,c)
                    -m_w[1].dyy(x,y,c)+2.0f*(m_b2[2].dxy(x,y,c)-m_w[2].dxy(x,y,c)));
        if(rs(x,y,c)!=rs(x,y,c))
            throw runtime_error("Result of RS is NaN");
    }

    //Gauss-Seidel method
    float diviser   = (1.0f+4.0f*m_l1+20.0f*m_l2);
    float coeff     = -m_l1-8.0f*m_l2;
    float tmp;
    FloatImg uOld(m_u);
    int i=0;
    while(i<1){
        cimg_forXYC(m_u,x,y,c){
            tmp  = m_l2*m_u.at(x-2,y,c);
            tmp += 2.0f*m_l2*m_u.at(x-1,y-1,c)+coeff*m_u.at(x-1,y,c)+2.0f*m_l2*m_u.at(x-1,y+1,c);
            tmp += m_l2*m_u.at(x,y-2,c)+coeff*m_u.at(x,y-1,c)+coeff*m_u.at(x,y+1,c)+m_l2*m_u.at(x,y+2,c);
            tmp += 2.0f*m_l2*m_u.at(x+1,y-1,c)+coeff*m_u.at(x+1,y,c)+2.0f*m_l2*m_u.at(x+1,y+1,c);
            tmp += m_l2*m_u.at(x+2,y,c);

            uOld(x,y,c) = m_u(x,y,c);
            m_u(x,y,c)  = (rs(x,y,c)-tmp)/diviser;
            if(m_u(x,y,c)!=m_u(x,y,c))
                throw runtime_error("Result of GS is NaN");
        }
        i++;
    }

}


void BregmanSolver::solve(){
    CImgDisplay disp((m_u,m_f),"TV-TV2 denoising",0,false,false);
    int t = 0;//Temps
    cout << "Iter : " << t << ", noise variance = " << m_u.variance_noise() << endl;
    const float white[] = { 255,255,255 };
    while(!disp.is_closed()){
        solve_subproblem1_GS();
        solve_subproblem2();
        solve_subproblem3();
        update_b();

        //Display iteration number
        t++;
        cout << "Iter : " << t << ", noise variance = " << m_u.variance_noise() << endl;
        (cimg_library::CImg<>(m_u).draw_text(2,2,"iter = %d",white,0,1,13,t),m_f).display(disp.wait(100));
  }
  m_iter = t;
  cimg_library::CImg<>(m_u).draw_text(2,2,"iter = %d",white,0,1,13,t+1).save("img/output/test.bmp",t);
}


DiffImg BregmanSolver::get_reconstructed_image(){
    return m_u;
}

void BregmanSolver::save(const char* const filename){
    const float white[] = { 255,255,255 };
    string str = string("img/output/")+filename;
    cimg_library::CImg<>(m_u).draw_text(2,2,"iter = %d",white,0,1,13,m_iter+1).save(str.c_str());
}

void BregmanSolver::solve_subproblem1(){
    //Compute right side of the equation
    FloatImg rs(m_u), rsRealFFT(m_u), rsImagFFT(m_u);
    cimg_forXYC(rs,x,y,c){
        rs(x,y,c) = m_f(x,y,c)+m_l1*(m_b1[0].bdx(x,y,c)-m_v[0].bdx(x,y,c)+m_b1[1].bdy(x,y,c)-m_v[1].bdy(x,y,c))
                    -m_l2*(m_b2[0].dxx(x,y,c)-m_w[0].dxx(x,y,c)+m_b2[1].dyy(x,y,c)
                    -m_w[1].dyy(x,y,c)+2.0f*(m_b2[2].dxy(x,y,c)-m_w[2].dxy(x,y,c)));
        if(rs(x,y,c)!=rs(x,y,c))
            throw runtime_error("Result of RS is NaN");
    }
    //FFT
    rs.FFT(rsRealFFT,rsImagFFT);



}

void BregmanSolver::compute_fourier_denominator(){
    int width  = m_u.width();
    int height = m_u.height();

    DiffImg rtmp(m_u), itmp(m_u);
    //dxx matrix
    DiffImg mxx(m_u);
    mxx.fill(0);
    cimg_forC(mxx,c){
        mxx(0,0,c)       = -2.0;
        mxx(1,0,c)       = 1.0;
        mxx(width-1,0,c) = 1.0;
    }
    mxx.FFT(rtmp,itmp);
    m_rfd -= m_l1*rtmp;
    m_ifd -= m_l1*itmp;
    delete &mxx;

    //dyy matrix
    DiffImg myy(m_u);
    myy.fill(0.0);
    cimg_forC(myy,c){
        myy(0,0,c)        = -2.0;
        myy(0,1,c)        = 1.0;
        myy(0,height-1,c) = 1.0;
    }
    myy.FFT(rtmp,itmp);
    m_rfd -= m_l1*rtmp;
    m_ifd -= m_l1*itmp;
    delete &myy;

    //d4x matrix
    DiffImg m4x(m_u);
}
