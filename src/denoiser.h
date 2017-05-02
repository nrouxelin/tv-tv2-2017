#ifndef BREGMAN_DENOISER_H
#define BREGMAN_DENOISER_H 0

#include "CImg.h"
namespace cil = cimg_library;
#include <stdexcept>
#include "diff.h"
#include <vector>
#include <iostream>
#include <string>



typedef std::vector<DiffImg> ArrayDiffImg;

class BregmanDenoiser{
public:
    BregmanDenoiser(DiffImg f, float a, float b, float l1, float l2);
    void solve_subproblem1();
    void solve_subproblem1_GS();
    void solve_subproblem2();
    void solve_subproblem3();
    void update_b();
    void solve();
    void compute_fourier_denominator();
    DiffImg get_reconstructed_image();
    void save(const char* const filename);


private:
    DiffImg m_f;//Original image
    DiffImg m_u;//Reconstructed image
    FloatImgList m_fd;//Fourier denominator
    ArrayDiffImg m_v;//v:=\nabla u
    ArrayDiffImg m_w; //w:=\nabla^2u
    ArrayDiffImg m_b1;//b_1^n
    ArrayDiffImg m_b2;//b_2^n
    DiffImg m_b0;
    float m_l1, m_l2,m_a,m_b;
    int m_iter;

};

#endif
