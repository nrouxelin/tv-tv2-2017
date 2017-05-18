#include "inpainter.h"
using namespace std;
using namespace cimg_library;

BregmanInpainter::BregmanInpainter(DiffImg f, FloatImg mask, float a, float b, float l0, float l1, float l2):
m_mask(mask),m_f(f),m_u(f),m_tu(f),m_l0(l0),m_l1(l1),m_l2(l2),m_a(a),m_b(b),m_b0(f){
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
      DiffImg b1_1(f), b1_2(f), b2_1(f), b2_2(f), b2_3(f), b0(f);
      m_b0.fill(1.0);

      b1_1.fill(1.0);
      b1_2.fill(1.0);
      m_b1.push_back(b1_1);
      m_b1.push_back(b1_2);

      b2_1.fill(1.0);
      b2_2.fill(1.0);
      b2_3.fill(1.0);
      m_b2.push_back(b2_1);
      m_b2.push_back(b2_2);
      m_b2.push_back(b2_3);

      m_iter = 0;
      compute_fourier_denominator();

}

void BregmanInpainter::compute_fourier_denominator(){
    int width  = m_u.width();
    int height = m_u.height();
    FloatImgList tmp;

    //dirac
    FloatImg id(m_u);
    id.fill(0.0);
    cimg_forC(id,c){
        id(0,0,c) = 1;
    }
    tmp = id.get_FFT();
    m_fd    = tmp;
    m_fd[0] = m_l0*m_fd[0];
    m_fd[1] = m_l0*m_fd[1];



    //dxx matrix
    DiffImg mxx(m_u);
    mxx.fill(0);
    cimg_forC(mxx,c){
        mxx(0,0,c)       = -2.0;
        mxx(1,0,c)       = 1.0;
        mxx(width-1,0,c) = 1.0;
    }

    tmp      = mxx.get_FFT();
    m_fd[0] -= m_l1*tmp[0];
    m_fd[1] -= m_l1*tmp[1];


    //dyy matrix
    DiffImg myy(m_u);
    myy.fill(0.0);
    cimg_forC(myy,c){
        myy(0,0,c)        = -2.0;
        myy(0,1,c)        = 1.0;
        myy(0,height-1,c) = 1.0;
    }

    tmp      = myy.get_FFT();
    m_fd[0] -= m_l1*tmp[0];
    m_fd[1] -= m_l1*tmp[1];
    //delete &myy;

    //d4x matrix
    DiffImg m4x(m_u);
    m4x.fill(0.0);
    cimg_forC(m4x,c){
        m4x(0,0,c)       = 6;
        m4x(1,0,c)       = -4;
        m4x(2,0,c)       = 1;
        m4x(width-2,0,c) = 1;
        m4x(width-1,0,c) = -4;
    }
    tmp      = m4x.get_FFT();
    m_fd[0] += m_l2*tmp[0];
    m_fd[1] += m_l2*tmp[1];

    //d4y matri
    DiffImg m4y(m_u);
    m4y.fill(0.0);
    cimg_forC(m4y,c){
        m4y(0,0,c)        = 6;
        m4y(0,1,c)        = -4;
        m4y(0,2,c)        = 1;
        m4y(0,height-2,c) = 1;
        m4y(0,height-1,c) = -4;
    }
    tmp      = m4y.get_FFT();
    m_fd[0] += 2.0*m_l2*tmp[0];
    m_fd[1] += 2.0*m_l2*tmp[1];

    //d2x2y matrix
    DiffImg m2x2y(m_u);
    m2x2y.fill(0.0);
    cimg_forC(m2x2y,c){
        m2x2y(0,0,c)                = 4;
        m2x2y(0,1,c)                = -2;
        m2x2y(1,0,c)                = -2;
        m2x2y(1,1,c)                = 1;
        m2x2y(width-1,0,c)          = -2;
        m2x2y(width-1,1,c)          = 1;
        m2x2y(0,height-1,c)         = -2;
        m2x2y(1,height-1,c)         = 1;
        m2x2y(width-1,height-1,c)   = 1;
    }
    tmp      = m2x2y.get_FFT();
    m_fd[0] += m_l2*tmp[0];
    m_fd[1] += m_l2*tmp[1];

    if(m_fd[0].is_nan() || m_fd[1].is_nan())
        throw runtime_error("m_fd is NaN");
}


void BregmanInpainter::solve_subproblem1(){
    float tmp;
    float coeff = 1.0/(m_l0+2.0);
    cimg_forXYC(m_u,x,y,c){
        tmp        = (2.0*m_mask(x,y,c))*m_f(x,y,c)
                        +(m_l0+2.0*(1.0-m_mask(x,y,c)))*(m_b0(x,y,c)+m_tu(x,y,c));
        m_u(x,y,c) = tmp*coeff;
    }
}



void BregmanInpainter::solve_subproblem2(){
    //Computing the right side of the equation
    DiffImg rs(m_u);
    cimg_forXYC(rs,x,y,c){
        rs(x,y,c) = m_l0*(m_u(x,y,c)-m_b0(x,y,c))
            +m_l1*(m_b1[0].bdx(x,y,c)-m_v[0].bdx(x,y,c)+m_b1[1].bdy(x,y,c)-m_v[1].bdy(x,y,c))
            -m_l2*(m_b2[0].dxx(x,y,c)-m_w[0].dxx(x,y,c)+m_b2[1].dyy(x,y,c)-m_w[1].dyy(x,y,c)
                +2.0*(m_b2[2].dxy(x,y,c)-m_w[2].dxy(x,y,c)));
    }

    //FFT
    FloatImgList frs = rs.get_FFT();

    FloatImg denom(frs[0]);
    denom = m_fd[0].get_pow(2)+m_fd[1].get_pow(2);
    //denom = m_fd[0].get_mul(m_fd[0])+m_fd[1].get_mul(m_fd[1]);
    cimg_forXYC(denom,x,y,c){
        if(denom(x,y,c) == 0){
            denom(x,y,c) = 0.001;
            cout << x << " " << y << " " << c << endl;
        }
    }

    frs[0] = frs[0].get_mul(m_fd[0])+frs[1].get_mul(m_fd[1]);
    frs[0].div(denom);

    frs[1] = frs[0].get_mul(m_fd[1])-frs[1].get_mul(m_fd[0]);
    frs[1].div(denom);

    FloatImg::FFT(frs[0],frs[1],true);
    m_tu = frs[0].normalize(0,255).rotate(180);

}

void BregmanInpainter::solve_subproblem3(){
    //Filling s1 and s2
    FloatImg s1(m_u), s2(m_u);
    cimg_forXYC(s1,x,y,c){
        s1(x,y,c) = m_b1[0](x,y,c)+m_tu.fdx(x,y,c);
        s2(x,y,c) = m_b1[1](x,y,c)+m_tu.fdy(x,y,c);
    }

    //Computing v^{n+1}
    float norm, m;
    cimg_forXYC(m_v[0],x,y,c){
        norm  = sqrt(s1(x,y,c)*s1(x,y,c)+s2(x,y,c)*s2(x,y,c));
        if(norm!=0){
            m             = max(norm-(m_a/m_l1),0.0f);
            m_v[0](x,y,c) = m*(s1(x,y,c)/norm);
            m_v[1](x,y,c) = m*(s2(x,y,c)/norm);
        }else{
            m_v[0](x,y,c) = 0.0f;
            m_v[1](x,y,c) = 0.0f;
        }

    }
}


void BregmanInpainter::solve_subproblem4(){
    //Filling t1, t2 and t3
    FloatImg t1(m_u), t2(m_u), t3(m_u);
    cimg_forXYC(t1,x,y,c){
        t1(x,y,c) = m_b2[0](x,y,c)+m_tu.dxx(x,y,c);
        t2(x,y,c) = m_b2[1](x,y,c)+m_tu.dyy(x,y,c);
        t3(x,y,c) = m_b2[2](x,y,c)+m_tu.dxy(x,y,c);
    }

    //Computing w^{n+1}
    float norm, m;
    cimg_forXYC(m_w[0],x,y,c){
        norm = sqrt(t1(x,y,c)*t1(x,y,c)+t2(x,y,c)*t2(x,y,c)
                            +2*t3(x,y,c)*t3(x,y,c));
        if(norm!=0){
            m             = max(norm-(m_b/m_l2),0.0f);
            m_w[0](x,y,c) = m*(t1(x,y,c)/norm);
            m_w[1](x,y,c) = m*(t2(x,y,c)/norm);
            m_w[2](x,y,c) = m*(t3(x,y,c)/norm);
        }else{
            m_w[0](x,y,c) = 0.0f;
            m_w[1](x,y,c) = 0.0f;
            m_w[2](x,y,c) = 0.0f;
        }

    }
}

void BregmanInpainter::update_b(){
    cimg_forXYC(m_b0,x,y,c){
        m_b0(x,y,c) += m_tu(x,y,c)-m_u(x,y,c);

        m_b1[0](x,y,c) += m_tu.fdx(x,y,c)-m_v[0](x,y,c);
        m_b1[1](x,y,c) += m_tu.fdy(x,y,c)-m_v[1](x,y,c);

        m_b2[0](x,y,c) += m_tu.dxx(x,y,c)-m_w[0](x,y,c);
        m_b2[1](x,y,c) += m_tu.dyy(x,y,c)-m_w[1](x,y,c);
        m_b2[2](x,y,c) += m_tu.dxy(x,y,c)-m_w[2](x,y,c);
    }
    /**m_b0.normalize(0,255);
    m_b1[1].normalize(0,255);
    m_b1[0].normalize(0,255);
    m_b2[0].normalize(0,255);
    m_b2[1].normalize(0,255);
    m_b2[2].normalize(0,255);**/
}

void BregmanInpainter::solve(){
    CImgDisplay disp((m_u,m_f),"TV-TV2 inpainting",0,false,false);
    int t = 0;//Temps
    const float white[] = { 255,255,255 };
    DiffImg uOld(m_u);
    float diff,tmp;
    diff = 2000;
    while(diff>8e-5){ //!disp.is_closed()
        diff = 0;
        uOld = m_u;
        solve_subproblem1();
        solve_subproblem2();
        solve_subproblem3();
        solve_subproblem4();
        update_b();
        cimg_forXYC(m_u,x,y,c){
            tmp = abs(m_u(x,y,c)-uOld(x,y,c));
            diff = (tmp > diff) ? tmp : diff;
        }

        //Display iteration number
        t++;
        (cimg_library::CImg<>(m_u).draw_text(2,2,"iter = %d",white,0,1,13,t),m_f).display(disp.wait(100));
        cimg_library::CImg<>(m_u).save("img/output/test.png",t);
  }
  m_iter = t;
  cimg_library::CImg<>(m_u).save("img/output/test.png",t);
}

DiffImg BregmanInpainter::get_reconstructed_image(){
    return m_u;
}

void BregmanInpainter::save(const char* const filename){
    const float white[] = { 255,255,255 };
    string str = string("img/output/")+filename;
    cimg_library::CImg<>(m_u).save(str.c_str());
}
