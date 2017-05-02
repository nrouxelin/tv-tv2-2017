#ifndef BREGMAN_INPAINTER_H
#define BREGMAN_INPAINTER_H 0

#include "CImg.h"
namespace cil = cimg_library;
#include <stdexcept>
#include "diff.h"
#include <vector>
#include <iostream>
#include <string>



typedef std::vector<DiffImg> ArrayDiffImg;

class BregmanInpainter{
public:
    BregmanInpainter(DiffImg f, FloatImg mask, float a, float b, float l0, float l1, float l2);
    void solve_subproblem1();
    void solve_subproblem2();
    void solve_subproblem3();
    void solve_subproblem4();
    void update_b();
    void solve();
    void compute_fourier_denominator();
    DiffImg get_reconstructed_image();
    void save(const char* const filename);


private:
    FloatImg m_mask;//Mask for inpainting domain
    DiffImg m_f;//Original image
    DiffImg m_u;//Reconstructed image
    DiffImg m_tu;//tu:=u
    ArrayDiffImg m_v;//v:=\nabla u
    ArrayDiffImg m_w; //w:=\nabla^2u
    ArrayDiffImg m_b1;//b_1^n
    ArrayDiffImg m_b2;//b_2^n
    float m_l0, m_l1, m_l2, m_a, m_b;//Parameters for the model
    DiffImg m_b0;//b_0^n
    int m_iter;
    FloatImgList m_fd;//Fourier denominator

};

#endif
