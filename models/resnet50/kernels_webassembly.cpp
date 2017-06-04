
#include <stdlib.h>
#include <math.h>

float data_buffer[28691560];


void im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(const int * meta_buffer)
{
    const float *im = data_buffer + meta_buffer[0];
    float *col = data_buffer + meta_buffer[1];

    const int N = meta_buffer[2];
    const int C1 = meta_buffer[3];
    const int H1 = meta_buffer[4];
    const int W1 = meta_buffer[5];
    const int H2 = meta_buffer[6];
    const int W2 = meta_buffer[7];
    const int KH = meta_buffer[8];
    const int KW = meta_buffer[9];
    const int SH = meta_buffer[10];
    const int SW = meta_buffer[11];
    const int PH = meta_buffer[12];
    const int PW = meta_buffer[13];

    for (int gid = 0; gid < N*H2*W2*KH*KW*C1; gid += 1) {
        const int c1 = gid % C1;
        const int kw = gid / C1 % KW;
        const int kh = gid / C1 / KW % KH;
        const int w2 = gid / C1 / KW / KH % W2;
        const int h2 = gid / C1 / KW / KH / W2 % H2;
        const int  n = gid / C1 / KW / KH / W2 / H2;
        
        const int h1 = h2 * SH - PH + kh;
        const int w1 = w2 * SW - PW + kw;

        col[gid] = (h1 < 0 || h1 >= H1 || w1 < 0 || w1 >= W1) ? 0 : im[((n*H1+h1)*W1+w1)*C1+c1];
    }
}


#ifndef INCLUDE_EIGEN
#define INCLUDE_EIGEN
#include <Eigen/Dense>
#endif

void sgemm_fa64ce9b8c79c14b329984cb2174bfc5303bf88d4f4ac1c14844d102(const int * meta_buffer)
{
    float *A = data_buffer + meta_buffer[0];
    float *B = data_buffer + meta_buffer[1];
    float *C = data_buffer + meta_buffer[2];

    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > a_mat(A, 12544, 147);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > b_mat(B, 147, 64);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > c_mat(C, 12544, 64);
    //Eigen::Map<Eigen::Matrix<float, 12544, 147, Eigen::RowMajor> > a_mat(A, 12544, 147);
    //Eigen::Map<Eigen::Matrix<float, 147, 64, Eigen::RowMajor> > b_mat(B, 147, 64);
    //Eigen::Map<Eigen::Matrix<float, 12544, 64, Eigen::RowMajor> > c_mat(C, 12544, 64);

    c_mat.noalias() = a_mat * b_mat;
}


void axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(const int * meta_buffer)
{
    const float *X = data_buffer + meta_buffer[0];
    float *Y = data_buffer + meta_buffer[1];
    const float *B = data_buffer + meta_buffer[2];
    const int N = meta_buffer[3];
    const int C = meta_buffer[4];
  
    for (int gid = 0; gid < N * C; gid += 1) {
        int c = gid % C;
        int n = gid / C;
        float result = X[gid] + B[c];

        Y[n * C + c] = result;
    }
}


void relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(const int * meta_buffer)
{
    const float *X = data_buffer + meta_buffer[0];
    float *Y = data_buffer + meta_buffer[1];

    const int N = meta_buffer[2];
  
    for (int gid = 0; gid < N; gid += 1) {
        float result = X[gid];
        result = result < 0.0 ? 0.0 : result;
        
        Y[gid] = result;
    }
}


void maxpooling2d_dc4fb0430710858abcfbe8b612e9209838d9c5742cd57322d6ac9d6e(const int * meta_buffer)
{
    const float *X = data_buffer + meta_buffer[0];
    float *Y = data_buffer + meta_buffer[1];
    const int N = meta_buffer[2];
    const int H1 = meta_buffer[3];
    const int W1 = meta_buffer[4];
    const int C = meta_buffer[5];
    const int H2 = meta_buffer[6];
    const int W2 = meta_buffer[7];
    const int K = meta_buffer[8];
    const int S = meta_buffer[9];
    const int P = meta_buffer[10];

    for (int gid = 0; gid < N * H2 * W2 * C; gid += 1) {
        const int c = gid % C;
        const int w2 = gid / C % W2;
        const int h2 = gid / C / W2 % H2;
        const int n = gid / C / W2 / H2;

        float v = -1e7;
        for (int kh = 0; kh < K; kh++) {
            const int h1 = h2 * S - P + kh;
            if (h1 < 0 || h1 >= H1) continue;
            
            for (int kw = 0; kw < K; kw++) {
                const int w1 = w2 * S - P + kw;
                if (w1 < 0 || w1 >= W1) continue;

                v = v > X[((n * H1 + h1) * W1 + w1) * C + c] ? v : X[((n * H1 + h1) * W1 + w1) * C + c];
            }
        }

        Y[gid] = v;
    }
}


#ifndef INCLUDE_EIGEN
#define INCLUDE_EIGEN
#include <Eigen/Dense>
#endif

void sgemm_b4e0ea82c8a51819f004d0d142d0cd64a3efee7e10ad400c69c24b03(const int * meta_buffer)
{
    float *A = data_buffer + meta_buffer[0];
    float *B = data_buffer + meta_buffer[1];
    float *C = data_buffer + meta_buffer[2];

    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > a_mat(A, 3136, 64);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > b_mat(B, 64, 256);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > c_mat(C, 3136, 256);
    //Eigen::Map<Eigen::Matrix<float, 3136, 64, Eigen::RowMajor> > a_mat(A, 3136, 64);
    //Eigen::Map<Eigen::Matrix<float, 64, 256, Eigen::RowMajor> > b_mat(B, 64, 256);
    //Eigen::Map<Eigen::Matrix<float, 3136, 256, Eigen::RowMajor> > c_mat(C, 3136, 256);

    c_mat.noalias() = a_mat * b_mat;
}


#ifndef INCLUDE_EIGEN
#define INCLUDE_EIGEN
#include <Eigen/Dense>
#endif

void sgemm_76efe07eff429e52715530e8757d30fca197eb005af636f1660ce936(const int * meta_buffer)
{
    float *A = data_buffer + meta_buffer[0];
    float *B = data_buffer + meta_buffer[1];
    float *C = data_buffer + meta_buffer[2];

    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > a_mat(A, 3136, 64);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > b_mat(B, 64, 64);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > c_mat(C, 3136, 64);
    //Eigen::Map<Eigen::Matrix<float, 3136, 64, Eigen::RowMajor> > a_mat(A, 3136, 64);
    //Eigen::Map<Eigen::Matrix<float, 64, 64, Eigen::RowMajor> > b_mat(B, 64, 64);
    //Eigen::Map<Eigen::Matrix<float, 3136, 64, Eigen::RowMajor> > c_mat(C, 3136, 64);

    c_mat.noalias() = a_mat * b_mat;
}


#ifndef INCLUDE_EIGEN
#define INCLUDE_EIGEN
#include <Eigen/Dense>
#endif

void sgemm_854e5419af5e54efa8eb16243cb2d043e975e31e191a276cd7ea2b4a(const int * meta_buffer)
{
    float *A = data_buffer + meta_buffer[0];
    float *B = data_buffer + meta_buffer[1];
    float *C = data_buffer + meta_buffer[2];

    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > a_mat(A, 3136, 576);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > b_mat(B, 576, 64);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > c_mat(C, 3136, 64);
    //Eigen::Map<Eigen::Matrix<float, 3136, 576, Eigen::RowMajor> > a_mat(A, 3136, 576);
    //Eigen::Map<Eigen::Matrix<float, 576, 64, Eigen::RowMajor> > b_mat(B, 576, 64);
    //Eigen::Map<Eigen::Matrix<float, 3136, 64, Eigen::RowMajor> > c_mat(C, 3136, 64);

    c_mat.noalias() = a_mat * b_mat;
}


void elementwisesum_54db89b5b7579805e7f7abdf1dd69f582f55da0ff40728276f251c8f(const int * meta_buffer)
{
    const float *X0 = data_buffer + meta_buffer[0];
    const float *X1 = data_buffer + meta_buffer[1];
    float *Y = data_buffer + meta_buffer[2];
    const int N = meta_buffer[3];
  
    for (int gid = 0; gid < N; gid += 1) {
        float result = X0[gid] + X1[gid];

        Y[gid] = result;
    }
}


#ifndef INCLUDE_EIGEN
#define INCLUDE_EIGEN
#include <Eigen/Dense>
#endif

void sgemm_7e37df37d70e21c18906ab234c3c8000ebc717a84946295168cde9b5(const int * meta_buffer)
{
    float *A = data_buffer + meta_buffer[0];
    float *B = data_buffer + meta_buffer[1];
    float *C = data_buffer + meta_buffer[2];

    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > a_mat(A, 3136, 256);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > b_mat(B, 256, 64);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > c_mat(C, 3136, 64);
    //Eigen::Map<Eigen::Matrix<float, 3136, 256, Eigen::RowMajor> > a_mat(A, 3136, 256);
    //Eigen::Map<Eigen::Matrix<float, 256, 64, Eigen::RowMajor> > b_mat(B, 256, 64);
    //Eigen::Map<Eigen::Matrix<float, 3136, 64, Eigen::RowMajor> > c_mat(C, 3136, 64);

    c_mat.noalias() = a_mat * b_mat;
}


#ifndef INCLUDE_EIGEN
#define INCLUDE_EIGEN
#include <Eigen/Dense>
#endif

void sgemm_bbc9374155ac4b56b56eca8e080210ccbcad3e120301193aa501f171(const int * meta_buffer)
{
    float *A = data_buffer + meta_buffer[0];
    float *B = data_buffer + meta_buffer[1];
    float *C = data_buffer + meta_buffer[2];

    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > a_mat(A, 784, 256);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > b_mat(B, 256, 512);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > c_mat(C, 784, 512);
    //Eigen::Map<Eigen::Matrix<float, 784, 256, Eigen::RowMajor> > a_mat(A, 784, 256);
    //Eigen::Map<Eigen::Matrix<float, 256, 512, Eigen::RowMajor> > b_mat(B, 256, 512);
    //Eigen::Map<Eigen::Matrix<float, 784, 512, Eigen::RowMajor> > c_mat(C, 784, 512);

    c_mat.noalias() = a_mat * b_mat;
}


#ifndef INCLUDE_EIGEN
#define INCLUDE_EIGEN
#include <Eigen/Dense>
#endif

void sgemm_eb85e5a4b2c5104c45e087e3078a5d782b912829c77b07cd4b3a6385(const int * meta_buffer)
{
    float *A = data_buffer + meta_buffer[0];
    float *B = data_buffer + meta_buffer[1];
    float *C = data_buffer + meta_buffer[2];

    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > a_mat(A, 784, 256);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > b_mat(B, 256, 128);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > c_mat(C, 784, 128);
    //Eigen::Map<Eigen::Matrix<float, 784, 256, Eigen::RowMajor> > a_mat(A, 784, 256);
    //Eigen::Map<Eigen::Matrix<float, 256, 128, Eigen::RowMajor> > b_mat(B, 256, 128);
    //Eigen::Map<Eigen::Matrix<float, 784, 128, Eigen::RowMajor> > c_mat(C, 784, 128);

    c_mat.noalias() = a_mat * b_mat;
}


#ifndef INCLUDE_EIGEN
#define INCLUDE_EIGEN
#include <Eigen/Dense>
#endif

void sgemm_0eb94687b795f23b3fa8aabd9019c2c0d564a8fc9f72a5971ea7956e(const int * meta_buffer)
{
    float *A = data_buffer + meta_buffer[0];
    float *B = data_buffer + meta_buffer[1];
    float *C = data_buffer + meta_buffer[2];

    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > a_mat(A, 784, 1152);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > b_mat(B, 1152, 128);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > c_mat(C, 784, 128);
    //Eigen::Map<Eigen::Matrix<float, 784, 1152, Eigen::RowMajor> > a_mat(A, 784, 1152);
    //Eigen::Map<Eigen::Matrix<float, 1152, 128, Eigen::RowMajor> > b_mat(B, 1152, 128);
    //Eigen::Map<Eigen::Matrix<float, 784, 128, Eigen::RowMajor> > c_mat(C, 784, 128);

    c_mat.noalias() = a_mat * b_mat;
}


#ifndef INCLUDE_EIGEN
#define INCLUDE_EIGEN
#include <Eigen/Dense>
#endif

void sgemm_1ff0b62f802dcdbc5760e78562eff821a8c810d25b97fcad8be73a8f(const int * meta_buffer)
{
    float *A = data_buffer + meta_buffer[0];
    float *B = data_buffer + meta_buffer[1];
    float *C = data_buffer + meta_buffer[2];

    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > a_mat(A, 784, 128);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > b_mat(B, 128, 512);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > c_mat(C, 784, 512);
    //Eigen::Map<Eigen::Matrix<float, 784, 128, Eigen::RowMajor> > a_mat(A, 784, 128);
    //Eigen::Map<Eigen::Matrix<float, 128, 512, Eigen::RowMajor> > b_mat(B, 128, 512);
    //Eigen::Map<Eigen::Matrix<float, 784, 512, Eigen::RowMajor> > c_mat(C, 784, 512);

    c_mat.noalias() = a_mat * b_mat;
}


#ifndef INCLUDE_EIGEN
#define INCLUDE_EIGEN
#include <Eigen/Dense>
#endif

void sgemm_8312b38283c7d5d1fcdef78686cd27a7a523deaa3be8d59b751471c4(const int * meta_buffer)
{
    float *A = data_buffer + meta_buffer[0];
    float *B = data_buffer + meta_buffer[1];
    float *C = data_buffer + meta_buffer[2];

    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > a_mat(A, 784, 512);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > b_mat(B, 512, 128);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > c_mat(C, 784, 128);
    //Eigen::Map<Eigen::Matrix<float, 784, 512, Eigen::RowMajor> > a_mat(A, 784, 512);
    //Eigen::Map<Eigen::Matrix<float, 512, 128, Eigen::RowMajor> > b_mat(B, 512, 128);
    //Eigen::Map<Eigen::Matrix<float, 784, 128, Eigen::RowMajor> > c_mat(C, 784, 128);

    c_mat.noalias() = a_mat * b_mat;
}


#ifndef INCLUDE_EIGEN
#define INCLUDE_EIGEN
#include <Eigen/Dense>
#endif

void sgemm_5dff80f9525dc6f6d35c28060669accfa4d0aa39fc596947630600d5(const int * meta_buffer)
{
    float *A = data_buffer + meta_buffer[0];
    float *B = data_buffer + meta_buffer[1];
    float *C = data_buffer + meta_buffer[2];

    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > a_mat(A, 196, 512);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > b_mat(B, 512, 1024);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > c_mat(C, 196, 1024);
    //Eigen::Map<Eigen::Matrix<float, 196, 512, Eigen::RowMajor> > a_mat(A, 196, 512);
    //Eigen::Map<Eigen::Matrix<float, 512, 1024, Eigen::RowMajor> > b_mat(B, 512, 1024);
    //Eigen::Map<Eigen::Matrix<float, 196, 1024, Eigen::RowMajor> > c_mat(C, 196, 1024);

    c_mat.noalias() = a_mat * b_mat;
}


#ifndef INCLUDE_EIGEN
#define INCLUDE_EIGEN
#include <Eigen/Dense>
#endif

void sgemm_0d6862002b1c18ab90f4d5491388ca0df80f697fca56f05bd47a4616(const int * meta_buffer)
{
    float *A = data_buffer + meta_buffer[0];
    float *B = data_buffer + meta_buffer[1];
    float *C = data_buffer + meta_buffer[2];

    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > a_mat(A, 196, 512);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > b_mat(B, 512, 256);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > c_mat(C, 196, 256);
    //Eigen::Map<Eigen::Matrix<float, 196, 512, Eigen::RowMajor> > a_mat(A, 196, 512);
    //Eigen::Map<Eigen::Matrix<float, 512, 256, Eigen::RowMajor> > b_mat(B, 512, 256);
    //Eigen::Map<Eigen::Matrix<float, 196, 256, Eigen::RowMajor> > c_mat(C, 196, 256);

    c_mat.noalias() = a_mat * b_mat;
}


#ifndef INCLUDE_EIGEN
#define INCLUDE_EIGEN
#include <Eigen/Dense>
#endif

void sgemm_76b68f11ef9ec398cfc28347c2d754b41e7ca8b38c0a19dcc17db0ac(const int * meta_buffer)
{
    float *A = data_buffer + meta_buffer[0];
    float *B = data_buffer + meta_buffer[1];
    float *C = data_buffer + meta_buffer[2];

    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > a_mat(A, 196, 2304);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > b_mat(B, 2304, 256);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > c_mat(C, 196, 256);
    //Eigen::Map<Eigen::Matrix<float, 196, 2304, Eigen::RowMajor> > a_mat(A, 196, 2304);
    //Eigen::Map<Eigen::Matrix<float, 2304, 256, Eigen::RowMajor> > b_mat(B, 2304, 256);
    //Eigen::Map<Eigen::Matrix<float, 196, 256, Eigen::RowMajor> > c_mat(C, 196, 256);

    c_mat.noalias() = a_mat * b_mat;
}


#ifndef INCLUDE_EIGEN
#define INCLUDE_EIGEN
#include <Eigen/Dense>
#endif

void sgemm_308de8e5de116f30813582dedc9b735bb25b3358790c59e1afdcdc25(const int * meta_buffer)
{
    float *A = data_buffer + meta_buffer[0];
    float *B = data_buffer + meta_buffer[1];
    float *C = data_buffer + meta_buffer[2];

    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > a_mat(A, 196, 256);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > b_mat(B, 256, 1024);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > c_mat(C, 196, 1024);
    //Eigen::Map<Eigen::Matrix<float, 196, 256, Eigen::RowMajor> > a_mat(A, 196, 256);
    //Eigen::Map<Eigen::Matrix<float, 256, 1024, Eigen::RowMajor> > b_mat(B, 256, 1024);
    //Eigen::Map<Eigen::Matrix<float, 196, 1024, Eigen::RowMajor> > c_mat(C, 196, 1024);

    c_mat.noalias() = a_mat * b_mat;
}


#ifndef INCLUDE_EIGEN
#define INCLUDE_EIGEN
#include <Eigen/Dense>
#endif

void sgemm_b5cb420578d588f8bb861fbb408652619c0d633383fc83ab9f53082e(const int * meta_buffer)
{
    float *A = data_buffer + meta_buffer[0];
    float *B = data_buffer + meta_buffer[1];
    float *C = data_buffer + meta_buffer[2];

    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > a_mat(A, 196, 1024);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > b_mat(B, 1024, 256);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > c_mat(C, 196, 256);
    //Eigen::Map<Eigen::Matrix<float, 196, 1024, Eigen::RowMajor> > a_mat(A, 196, 1024);
    //Eigen::Map<Eigen::Matrix<float, 1024, 256, Eigen::RowMajor> > b_mat(B, 1024, 256);
    //Eigen::Map<Eigen::Matrix<float, 196, 256, Eigen::RowMajor> > c_mat(C, 196, 256);

    c_mat.noalias() = a_mat * b_mat;
}


#ifndef INCLUDE_EIGEN
#define INCLUDE_EIGEN
#include <Eigen/Dense>
#endif

void sgemm_201e41d345cb0eebe77bed8c1bffe4ef2348e01d6d321ae317404457(const int * meta_buffer)
{
    float *A = data_buffer + meta_buffer[0];
    float *B = data_buffer + meta_buffer[1];
    float *C = data_buffer + meta_buffer[2];

    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > a_mat(A, 49, 1024);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > b_mat(B, 1024, 512);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > c_mat(C, 49, 512);
    //Eigen::Map<Eigen::Matrix<float, 49, 1024, Eigen::RowMajor> > a_mat(A, 49, 1024);
    //Eigen::Map<Eigen::Matrix<float, 1024, 512, Eigen::RowMajor> > b_mat(B, 1024, 512);
    //Eigen::Map<Eigen::Matrix<float, 49, 512, Eigen::RowMajor> > c_mat(C, 49, 512);

    c_mat.noalias() = a_mat * b_mat;
}


#ifndef INCLUDE_EIGEN
#define INCLUDE_EIGEN
#include <Eigen/Dense>
#endif

void sgemm_675f43f38f11c7e8b176da63f208afead50e062474e037d7e5569bbf(const int * meta_buffer)
{
    float *A = data_buffer + meta_buffer[0];
    float *B = data_buffer + meta_buffer[1];
    float *C = data_buffer + meta_buffer[2];

    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > a_mat(A, 49, 4608);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > b_mat(B, 4608, 512);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > c_mat(C, 49, 512);
    //Eigen::Map<Eigen::Matrix<float, 49, 4608, Eigen::RowMajor> > a_mat(A, 49, 4608);
    //Eigen::Map<Eigen::Matrix<float, 4608, 512, Eigen::RowMajor> > b_mat(B, 4608, 512);
    //Eigen::Map<Eigen::Matrix<float, 49, 512, Eigen::RowMajor> > c_mat(C, 49, 512);

    c_mat.noalias() = a_mat * b_mat;
}


#ifndef INCLUDE_EIGEN
#define INCLUDE_EIGEN
#include <Eigen/Dense>
#endif

void sgemm_4a7b51005273e75f1cab8ba4b77dc86c5ce53b149a963daaa5d423e1(const int * meta_buffer)
{
    float *A = data_buffer + meta_buffer[0];
    float *B = data_buffer + meta_buffer[1];
    float *C = data_buffer + meta_buffer[2];

    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > a_mat(A, 49, 512);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > b_mat(B, 512, 2048);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > c_mat(C, 49, 2048);
    //Eigen::Map<Eigen::Matrix<float, 49, 512, Eigen::RowMajor> > a_mat(A, 49, 512);
    //Eigen::Map<Eigen::Matrix<float, 512, 2048, Eigen::RowMajor> > b_mat(B, 512, 2048);
    //Eigen::Map<Eigen::Matrix<float, 49, 2048, Eigen::RowMajor> > c_mat(C, 49, 2048);

    c_mat.noalias() = a_mat * b_mat;
}


#ifndef INCLUDE_EIGEN
#define INCLUDE_EIGEN
#include <Eigen/Dense>
#endif

void sgemm_bcfd30b9710b12cd25b0a3bd8ea938966f7ab10cc9262ab99d63696b(const int * meta_buffer)
{
    float *A = data_buffer + meta_buffer[0];
    float *B = data_buffer + meta_buffer[1];
    float *C = data_buffer + meta_buffer[2];

    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > a_mat(A, 49, 1024);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > b_mat(B, 1024, 2048);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > c_mat(C, 49, 2048);
    //Eigen::Map<Eigen::Matrix<float, 49, 1024, Eigen::RowMajor> > a_mat(A, 49, 1024);
    //Eigen::Map<Eigen::Matrix<float, 1024, 2048, Eigen::RowMajor> > b_mat(B, 1024, 2048);
    //Eigen::Map<Eigen::Matrix<float, 49, 2048, Eigen::RowMajor> > c_mat(C, 49, 2048);

    c_mat.noalias() = a_mat * b_mat;
}


#ifndef INCLUDE_EIGEN
#define INCLUDE_EIGEN
#include <Eigen/Dense>
#endif

void sgemm_187235ada597108e1eb448f6650b020d7557ddf4bb3b5910355d58ca(const int * meta_buffer)
{
    float *A = data_buffer + meta_buffer[0];
    float *B = data_buffer + meta_buffer[1];
    float *C = data_buffer + meta_buffer[2];

    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > a_mat(A, 49, 2048);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > b_mat(B, 2048, 512);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > c_mat(C, 49, 512);
    //Eigen::Map<Eigen::Matrix<float, 49, 2048, Eigen::RowMajor> > a_mat(A, 49, 2048);
    //Eigen::Map<Eigen::Matrix<float, 2048, 512, Eigen::RowMajor> > b_mat(B, 2048, 512);
    //Eigen::Map<Eigen::Matrix<float, 49, 512, Eigen::RowMajor> > c_mat(C, 49, 512);

    c_mat.noalias() = a_mat * b_mat;
}


void averagepooling2d_203c4e21c29b022faa72777a679127ac4aec4267b6caaf17b363a270(const int * meta_buffer)
{
    const float *X = data_buffer + meta_buffer[0];
    float *Y = data_buffer + meta_buffer[1];
    const int N = meta_buffer[2];
    const int H1 = meta_buffer[3];
    const int W1 = meta_buffer[4];
    const int C = meta_buffer[5];
    const int H2 = meta_buffer[6];
    const int W2 = meta_buffer[7];
    const int K = meta_buffer[8];
    const int S = meta_buffer[9];
    const int P = meta_buffer[10];
    
    for (int gid = 0; gid < N * H2 * W2 * C; gid += 1) {
        const int c = gid % C;
        const int w2 = gid / C % W2;
        const int h2 = gid / C / W2 % H2;
        const int n = gid / C / W2 / H2;

        float v = 0;
        for (int kh = 0; kh < K; kh++) {
            const int h1 = h2 * S - P + kh;
            if (h1 < 0 || h1 >= H1) continue;
            
            for (int kw = 0; kw < K; kw++) {
                const int w1 = w2 * S - P + kw;
                if (w1 < 0 || w1 >= W1) continue;

                v += X[((n * H1 + h1) * W1 + w1) * C + c];
            }
        }
        v /= K * K;

        Y[gid] = v;
    }
}


void flatten_143a2725e849669b2cf7cd13c75f75d2068fc136b9235e3aaf6cab38(const int * meta_buffer )
{
    const float *x = data_buffer + meta_buffer[0];
    float *y = data_buffer + meta_buffer[1];

    const int N = meta_buffer[2];

    for (int gid = 0; gid < N; gid += 1) {
        y[gid] = x[gid];
    }
}


#ifndef INCLUDE_EIGEN
#define INCLUDE_EIGEN
#include <Eigen/Dense>
#endif

void sgemm_814127e7734e5b079739c0a1d49f4231b555135cd2d42721305bb3f6(const int * meta_buffer)
{
    float *A = data_buffer + meta_buffer[0];
    float *B = data_buffer + meta_buffer[1];
    float *C = data_buffer + meta_buffer[2];

    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > a_mat(A, 1, 2048);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > b_mat(B, 2048, 1000);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > c_mat(C, 1, 1000);
    //Eigen::Map<Eigen::Matrix<float, 1, 2048, Eigen::RowMajor> > a_mat(A, 1, 2048);
    //Eigen::Map<Eigen::Matrix<float, 2048, 1000, Eigen::RowMajor> > b_mat(B, 2048, 1000);
    //Eigen::Map<Eigen::Matrix<float, 1, 1000, Eigen::RowMajor> > c_mat(C, 1, 1000);

    c_mat.noalias() = a_mat * b_mat;
}

extern "C" void init() {
    //data_buffer = (float*)malloc(28691560 * sizeof(float));
}

extern "C" float* get_data_buffer(void) {
    return data_buffer;
}



extern "C" void run() {
const int meta_buf_0[] = {25530472,26483816,1,3,224,224,112,112,7,7,2,2,3,3};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_0);
const int meta_buf_1[] = {26483816,25489408,25681000,12544,64,147};
sgemm_fa64ce9b8c79c14b329984cb2174bfc5303bf88d4f4ac1c14844d102(meta_buf_1);
const int meta_buf_2[] = {25681000,25681000,25530088,12544,64};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_2);
const int meta_buf_3[] = {25681000,25681000,802816};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_3);
const int meta_buf_4[] = {25681000,26483816,1,112,112,64,56,56,3,2,0};
maxpooling2d_dc4fb0430710858abcfbe8b612e9209838d9c5742cd57322d6ac9d6e(meta_buf_4);
const int meta_buf_5[] = {26483816,25391104,27688040,3136,256,64};
sgemm_b4e0ea82c8a51819f004d0d142d0cd64a3efee7e10ad400c69c24b03(meta_buf_5);
const int meta_buf_6[] = {27688040,27688040,25526440,3136,256};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_6);
const int meta_buf_7[] = {26483816,25498816,25681000,3136,64,64};
sgemm_76efe07eff429e52715530e8757d30fca197eb005af636f1660ce936(meta_buf_7);
const int meta_buf_8[] = {25681000,25681000,25530024,3136,64};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_8);
const int meta_buf_9[] = {25681000,25681000,200704};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_9);
const int meta_buf_10[] = {25681000,25881704,1,64,56,56,56,56,3,3,1,1,1,1};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_10);
const int meta_buf_11[] = {25881704,25321472,28490856,3136,64,576};
sgemm_854e5419af5e54efa8eb16243cb2d043e975e31e191a276cd7ea2b4a(meta_buf_11);
const int meta_buf_12[] = {28490856,28490856,25530280,3136,64};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_12);
const int meta_buf_13[] = {28490856,28490856,200704};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_13);
const int meta_buf_14[] = {28490856,25440256,25681000,3136,256,64};
sgemm_b4e0ea82c8a51819f004d0d142d0cd64a3efee7e10ad400c69c24b03(meta_buf_14);
const int meta_buf_15[] = {25681000,25681000,25525160,3136,256};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_15);
const int meta_buf_16[] = {25681000,27688040,25681000,802816};
elementwisesum_54db89b5b7579805e7f7abdf1dd69f582f55da0ff40728276f251c8f(meta_buf_16);
const int meta_buf_17[] = {25681000,25681000,802816};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_17);
const int meta_buf_18[] = {25681000,25407488,26483816,3136,64,256};
sgemm_7e37df37d70e21c18906ab234c3c8000ebc717a84946295168cde9b5(meta_buf_18);
const int meta_buf_19[] = {26483816,26483816,25530152,3136,64};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_19);
const int meta_buf_20[] = {26483816,26483816,200704};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_20);
const int meta_buf_21[] = {26483816,26684520,1,64,56,56,56,56,3,3,1,1,1,1};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_21);
const int meta_buf_22[] = {26684520,25284608,26483816,3136,64,576};
sgemm_854e5419af5e54efa8eb16243cb2d043e975e31e191a276cd7ea2b4a(meta_buf_22);
const int meta_buf_23[] = {26483816,26483816,25530344,3136,64};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_23);
const int meta_buf_24[] = {26483816,26483816,200704};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_24);
const int meta_buf_25[] = {26483816,25456640,27888744,3136,256,64};
sgemm_b4e0ea82c8a51819f004d0d142d0cd64a3efee7e10ad400c69c24b03(meta_buf_25);
const int meta_buf_26[] = {27888744,27888744,25527208,3136,256};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_26);
const int meta_buf_27[] = {27888744,25681000,27888744,802816};
elementwisesum_54db89b5b7579805e7f7abdf1dd69f582f55da0ff40728276f251c8f(meta_buf_27);
const int meta_buf_28[] = {27888744,27888744,802816};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_28);
const int meta_buf_29[] = {27888744,25423872,25681000,3136,64,256};
sgemm_7e37df37d70e21c18906ab234c3c8000ebc717a84946295168cde9b5(meta_buf_29);
const int meta_buf_30[] = {25681000,25681000,25530408,3136,64};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_30);
const int meta_buf_31[] = {25681000,25681000,200704};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_31);
const int meta_buf_32[] = {25681000,25881704,1,64,56,56,56,56,3,3,1,1,1,1};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_32);
const int meta_buf_33[] = {25881704,25247744,27688040,3136,64,576};
sgemm_854e5419af5e54efa8eb16243cb2d043e975e31e191a276cd7ea2b4a(meta_buf_33);
const int meta_buf_34[] = {27688040,27688040,25530216,3136,64};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_34);
const int meta_buf_35[] = {27688040,27688040,200704};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_35);
const int meta_buf_36[] = {27688040,25473024,25681000,3136,256,64};
sgemm_b4e0ea82c8a51819f004d0d142d0cd64a3efee7e10ad400c69c24b03(meta_buf_36);
const int meta_buf_37[] = {25681000,25681000,25527720,3136,256};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_37);
const int meta_buf_38[] = {25681000,27888744,25681000,802816};
elementwisesum_54db89b5b7579805e7f7abdf1dd69f582f55da0ff40728276f251c8f(meta_buf_38);
const int meta_buf_39[] = {25681000,25681000,802816};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_39);
const int meta_buf_40[] = {25681000,26483816,1,256,56,56,28,28,1,1,2,2,0,0};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_40);
const int meta_buf_41[] = {26483816,24657920,26784872,784,512,256};
sgemm_bbc9374155ac4b56b56eca8e080210ccbcad3e120301193aa501f171(meta_buf_41);
const int meta_buf_42[] = {26784872,26784872,25523880,784,512};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_42);
const int meta_buf_43[] = {25681000,26483816,1,256,56,56,28,28,1,1,2,2,0,0};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_43);
const int meta_buf_44[] = {26483816,25358336,25681000,784,128,256};
sgemm_eb85e5a4b2c5104c45e087e3078a5d782b912829c77b07cd4b3a6385(meta_buf_44);
const int meta_buf_45[] = {25681000,25681000,25529768,784,128};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_45);
const int meta_buf_46[] = {25681000,25681000,100352};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_46);
const int meta_buf_47[] = {25681000,25781352,1,128,28,28,28,28,3,3,1,1,1,1};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_47);
const int meta_buf_48[] = {25781352,24084480,26684520,784,128,1152};
sgemm_0eb94687b795f23b3fa8aabd9019c2c0d564a8fc9f72a5971ea7956e(meta_buf_48);
const int meta_buf_49[] = {26684520,26684520,25529384,784,128};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_49);
const int meta_buf_50[] = {26684520,26684520,100352};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_50);
const int meta_buf_51[] = {26684520,25182208,25681000,784,512,128};
sgemm_1ff0b62f802dcdbc5760e78562eff821a8c810d25b97fcad8be73a8f(meta_buf_51);
const int meta_buf_52[] = {25681000,25681000,25519784,784,512};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_52);
const int meta_buf_53[] = {25681000,26784872,25681000,401408};
elementwisesum_54db89b5b7579805e7f7abdf1dd69f582f55da0ff40728276f251c8f(meta_buf_53);
const int meta_buf_54[] = {25681000,25681000,401408};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_54);
const int meta_buf_55[] = {25681000,24854528,26082408,784,128,512};
sgemm_8312b38283c7d5d1fcdef78686cd27a7a523deaa3be8d59b751471c4(meta_buf_55);
const int meta_buf_56[] = {26082408,26082408,25529256,784,128};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_56);
const int meta_buf_57[] = {26082408,26082408,100352};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_57);
const int meta_buf_58[] = {26082408,26182760,1,128,28,28,28,28,3,3,1,1,1,1};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_58);
const int meta_buf_59[] = {26182760,24231936,26082408,784,128,1152};
sgemm_0eb94687b795f23b3fa8aabd9019c2c0d564a8fc9f72a5971ea7956e(meta_buf_59);
const int meta_buf_60[] = {26082408,26082408,25529896,784,128};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_60);
const int meta_buf_61[] = {26082408,26082408,100352};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_61);
const int meta_buf_62[] = {26082408,25051136,26784872,784,512,128};
sgemm_1ff0b62f802dcdbc5760e78562eff821a8c810d25b97fcad8be73a8f(meta_buf_62);
const int meta_buf_63[] = {26784872,26784872,25521832,784,512};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_63);
const int meta_buf_64[] = {26784872,25681000,26784872,401408};
elementwisesum_54db89b5b7579805e7f7abdf1dd69f582f55da0ff40728276f251c8f(meta_buf_64);
const int meta_buf_65[] = {26784872,26784872,401408};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_65);
const int meta_buf_66[] = {26784872,24985600,25681000,784,128,512};
sgemm_8312b38283c7d5d1fcdef78686cd27a7a523deaa3be8d59b751471c4(meta_buf_66);
const int meta_buf_67[] = {25681000,25681000,25529128,784,128};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_67);
const int meta_buf_68[] = {25681000,25681000,100352};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_68);
const int meta_buf_69[] = {25681000,25781352,1,128,28,28,28,28,3,3,1,1,1,1};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_69);
const int meta_buf_70[] = {25781352,23937024,26684520,784,128,1152};
sgemm_0eb94687b795f23b3fa8aabd9019c2c0d564a8fc9f72a5971ea7956e(meta_buf_70);
const int meta_buf_71[] = {26684520,26684520,25529000,784,128};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_71);
const int meta_buf_72[] = {26684520,26684520,100352};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_72);
const int meta_buf_73[] = {26684520,24920064,25681000,784,512,128};
sgemm_1ff0b62f802dcdbc5760e78562eff821a8c810d25b97fcad8be73a8f(meta_buf_73);
const int meta_buf_74[] = {25681000,25681000,25519272,784,512};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_74);
const int meta_buf_75[] = {25681000,26784872,25681000,401408};
elementwisesum_54db89b5b7579805e7f7abdf1dd69f582f55da0ff40728276f251c8f(meta_buf_75);
const int meta_buf_76[] = {25681000,25681000,401408};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_76);
const int meta_buf_77[] = {25681000,25116672,26082408,784,128,512};
sgemm_8312b38283c7d5d1fcdef78686cd27a7a523deaa3be8d59b751471c4(meta_buf_77);
const int meta_buf_78[] = {26082408,26082408,25529512,784,128};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_78);
const int meta_buf_79[] = {26082408,26082408,100352};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_79);
const int meta_buf_80[] = {26082408,26182760,1,128,28,28,28,28,3,3,1,1,1,1};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_80);
const int meta_buf_81[] = {26182760,24379392,27085928,784,128,1152};
sgemm_0eb94687b795f23b3fa8aabd9019c2c0d564a8fc9f72a5971ea7956e(meta_buf_81);
const int meta_buf_82[] = {27085928,27085928,25529640,784,128};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_82);
const int meta_buf_83[] = {27085928,27085928,100352};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_83);
const int meta_buf_84[] = {27085928,24788992,26082408,784,512,128};
sgemm_1ff0b62f802dcdbc5760e78562eff821a8c810d25b97fcad8be73a8f(meta_buf_84);
const int meta_buf_85[] = {26082408,26082408,25522344,784,512};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_85);
const int meta_buf_86[] = {26082408,25681000,26082408,401408};
elementwisesum_54db89b5b7579805e7f7abdf1dd69f582f55da0ff40728276f251c8f(meta_buf_86);
const int meta_buf_87[] = {26082408,26082408,401408};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_87);
const int meta_buf_88[] = {26082408,25681000,1,512,28,28,14,14,1,1,2,2,0,0};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_88);
const int meta_buf_89[] = {25681000,20529152,26483816,196,1024,512};
sgemm_5dff80f9525dc6f6d35c28060669accfa4d0aa39fc596947630600d5(meta_buf_89);
const int meta_buf_90[] = {26483816,26483816,25512128,196,1024};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_90);
const int meta_buf_91[] = {26082408,25731176,1,512,28,28,14,14,1,1,2,2,0,0};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_91);
const int meta_buf_92[] = {25731176,24526848,25681000,196,256,512};
sgemm_0d6862002b1c18ab90f4d5491388ca0df80f697fca56f05bd47a4616(meta_buf_92);
const int meta_buf_93[] = {25681000,25681000,25524904,196,256};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_93);
const int meta_buf_94[] = {25681000,25681000,50176};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_94);
const int meta_buf_95[] = {25681000,25731176,1,256,14,14,14,14,3,3,1,1,1,1};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_95);
const int meta_buf_96[] = {25731176,17055744,26182760,196,256,2304};
sgemm_76b68f11ef9ec398cfc28347c2d754b41e7ca8b38c0a19dcc17db0ac(meta_buf_96);
const int meta_buf_97[] = {26182760,26182760,25527464,196,256};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_97);
const int meta_buf_98[] = {26182760,26182760,50176};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_98);
const int meta_buf_99[] = {26182760,21839872,25681000,196,1024,256};
sgemm_308de8e5de116f30813582dedc9b735bb25b3358790c59e1afdcdc25(meta_buf_99);
const int meta_buf_100[] = {25681000,25681000,25513152,196,1024};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_100);
const int meta_buf_101[] = {25681000,26483816,25681000,200704};
elementwisesum_54db89b5b7579805e7f7abdf1dd69f582f55da0ff40728276f251c8f(meta_buf_101);
const int meta_buf_102[] = {25681000,25681000,200704};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_102);
const int meta_buf_103[] = {25681000,23150592,25881704,196,256,1024};
sgemm_b5cb420578d588f8bb861fbb408652619c0d633383fc83ab9f53082e(meta_buf_103);
const int meta_buf_104[] = {25881704,25881704,25528232,196,256};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_104);
const int meta_buf_105[] = {25881704,25881704,50176};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_105);
const int meta_buf_106[] = {25881704,25931880,1,256,14,14,14,14,3,3,1,1,1,1};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_106);
const int meta_buf_107[] = {25931880,17645568,25881704,196,256,2304};
sgemm_76b68f11ef9ec398cfc28347c2d754b41e7ca8b38c0a19dcc17db0ac(meta_buf_107);
const int meta_buf_108[] = {25881704,25881704,25528744,196,256};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_108);
const int meta_buf_109[] = {25881704,25881704,50176};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_109);
const int meta_buf_110[] = {25881704,23412736,26232936,196,1024,256};
sgemm_308de8e5de116f30813582dedc9b735bb25b3358790c59e1afdcdc25(meta_buf_110);
const int meta_buf_111[] = {26232936,26232936,25511104,196,1024};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_111);
const int meta_buf_112[] = {26232936,25681000,26232936,200704};
elementwisesum_54db89b5b7579805e7f7abdf1dd69f582f55da0ff40728276f251c8f(meta_buf_112);
const int meta_buf_113[] = {26232936,26232936,200704};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_113);
const int meta_buf_114[] = {26232936,21315584,25681000,196,256,1024};
sgemm_b5cb420578d588f8bb861fbb408652619c0d633383fc83ab9f53082e(meta_buf_114);
const int meta_buf_115[] = {25681000,25681000,25525416,196,256};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_115);
const int meta_buf_116[] = {25681000,25681000,50176};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_116);
const int meta_buf_117[] = {25681000,25731176,1,256,14,14,14,14,3,3,1,1,1,1};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_117);
const int meta_buf_118[] = {25731176,18825216,26182760,196,256,2304};
sgemm_76b68f11ef9ec398cfc28347c2d754b41e7ca8b38c0a19dcc17db0ac(meta_buf_118);
const int meta_buf_119[] = {26182760,26182760,25525928,196,256};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_119);
const int meta_buf_120[] = {26182760,26182760,50176};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_120);
const int meta_buf_121[] = {26182760,22888448,25681000,196,1024,256};
sgemm_308de8e5de116f30813582dedc9b735bb25b3358790c59e1afdcdc25(meta_buf_121);
const int meta_buf_122[] = {25681000,25681000,25516224,196,1024};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_122);
const int meta_buf_123[] = {25681000,26232936,25681000,200704};
elementwisesum_54db89b5b7579805e7f7abdf1dd69f582f55da0ff40728276f251c8f(meta_buf_123);
const int meta_buf_124[] = {25681000,25681000,200704};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_124);
const int meta_buf_125[] = {25681000,21577728,25881704,196,256,1024};
sgemm_b5cb420578d588f8bb861fbb408652619c0d633383fc83ab9f53082e(meta_buf_125);
const int meta_buf_126[] = {25881704,25881704,25528488,196,256};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_126);
const int meta_buf_127[] = {25881704,25881704,50176};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_127);
const int meta_buf_128[] = {25881704,25931880,1,256,14,14,14,14,3,3,1,1,1,1};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_128);
const int meta_buf_129[] = {25931880,19415040,25881704,196,256,2304};
sgemm_76b68f11ef9ec398cfc28347c2d754b41e7ca8b38c0a19dcc17db0ac(meta_buf_129);
const int meta_buf_130[] = {25881704,25881704,25526184,196,256};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_130);
const int meta_buf_131[] = {25881704,25881704,50176};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_131);
const int meta_buf_132[] = {25881704,22626304,26232936,196,1024,256};
sgemm_308de8e5de116f30813582dedc9b735bb25b3358790c59e1afdcdc25(meta_buf_132);
const int meta_buf_133[] = {26232936,26232936,25514176,196,1024};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_133);
const int meta_buf_134[] = {26232936,25681000,26232936,200704};
elementwisesum_54db89b5b7579805e7f7abdf1dd69f582f55da0ff40728276f251c8f(meta_buf_134);
const int meta_buf_135[] = {26232936,26232936,200704};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_135);
const int meta_buf_136[] = {26232936,22364160,25681000,196,256,1024};
sgemm_b5cb420578d588f8bb861fbb408652619c0d633383fc83ab9f53082e(meta_buf_136);
const int meta_buf_137[] = {25681000,25681000,25526696,196,256};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_137);
const int meta_buf_138[] = {25681000,25681000,50176};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_138);
const int meta_buf_139[] = {25681000,25731176,1,256,14,14,14,14,3,3,1,1,1,1};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_139);
const int meta_buf_140[] = {25731176,16465920,26182760,196,256,2304};
sgemm_76b68f11ef9ec398cfc28347c2d754b41e7ca8b38c0a19dcc17db0ac(meta_buf_140);
const int meta_buf_141[] = {26182760,26182760,25527976,196,256};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_141);
const int meta_buf_142[] = {26182760,26182760,50176};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_142);
const int meta_buf_143[] = {26182760,23674880,25681000,196,1024,256};
sgemm_308de8e5de116f30813582dedc9b735bb25b3358790c59e1afdcdc25(meta_buf_143);
const int meta_buf_144[] = {25681000,25681000,25517248,196,1024};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_144);
const int meta_buf_145[] = {25681000,26232936,25681000,200704};
elementwisesum_54db89b5b7579805e7f7abdf1dd69f582f55da0ff40728276f251c8f(meta_buf_145);
const int meta_buf_146[] = {25681000,25681000,200704};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_146);
const int meta_buf_147[] = {25681000,21053440,25881704,196,256,1024};
sgemm_b5cb420578d588f8bb861fbb408652619c0d633383fc83ab9f53082e(meta_buf_147);
const int meta_buf_148[] = {25881704,25881704,25526952,196,256};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_148);
const int meta_buf_149[] = {25881704,25881704,50176};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_149);
const int meta_buf_150[] = {25881704,25931880,1,256,14,14,14,14,3,3,1,1,1,1};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_150);
const int meta_buf_151[] = {25931880,18235392,25881704,196,256,2304};
sgemm_76b68f11ef9ec398cfc28347c2d754b41e7ca8b38c0a19dcc17db0ac(meta_buf_151);
const int meta_buf_152[] = {25881704,25881704,25525672,196,256};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_152);
const int meta_buf_153[] = {25881704,25881704,50176};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_153);
const int meta_buf_154[] = {25881704,22102016,25956968,196,1024,256};
sgemm_308de8e5de116f30813582dedc9b735bb25b3358790c59e1afdcdc25(meta_buf_154);
const int meta_buf_155[] = {25956968,25956968,25515200,196,1024};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_155);
const int meta_buf_156[] = {25956968,25681000,25956968,200704};
elementwisesum_54db89b5b7579805e7f7abdf1dd69f582f55da0ff40728276f251c8f(meta_buf_156);
const int meta_buf_157[] = {25956968,25956968,200704};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_157);
const int meta_buf_158[] = {25956968,25706088,1,1024,14,14,7,7,1,1,2,2,0,0};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_158);
const int meta_buf_159[] = {25706088,20004864,25681000,49,512,1024};
sgemm_201e41d345cb0eebe77bed8c1bffe4ef2348e01d6d321ae317404457(meta_buf_159);
const int meta_buf_160[] = {25681000,25681000,25520296,49,512};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_160);
const int meta_buf_161[] = {25681000,25681000,25088};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_161);
const int meta_buf_162[] = {25681000,25706088,1,512,7,7,7,7,3,3,1,1,1,1};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_162);
const int meta_buf_163[] = {25706088,2359296,25931880,49,512,4608};
sgemm_675f43f38f11c7e8b176da63f208afead50e062474e037d7e5569bbf(meta_buf_163);
const int meta_buf_164[] = {25931880,25931880,25520808,49,512};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_164);
const int meta_buf_165[] = {25931880,25931880,25088};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_165);
const int meta_buf_166[] = {25931880,12271616,25681000,49,2048,512};
sgemm_4a7b51005273e75f1cab8ba4b77dc86c5ce53b149a963daaa5d423e1(meta_buf_166);
const int meta_buf_167[] = {25681000,25681000,25504960,49,2048};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_167);
const int meta_buf_168[] = {25956968,25881704,1,1024,14,14,7,7,1,1,2,2,0,0};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_168);
const int meta_buf_169[] = {25881704,7077888,25781352,49,2048,1024};
sgemm_bcfd30b9710b12cd25b0a3bd8ea938966f7ab10cc9262ab99d63696b(meta_buf_169);
const int meta_buf_170[] = {25781352,25781352,25502912,49,2048};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_170);
const int meta_buf_171[] = {25681000,25781352,25681000,100352};
elementwisesum_54db89b5b7579805e7f7abdf1dd69f582f55da0ff40728276f251c8f(meta_buf_171);
const int meta_buf_172[] = {25681000,25681000,100352};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_172);
const int meta_buf_173[] = {25681000,15417344,25781352,49,512,2048};
sgemm_187235ada597108e1eb448f6650b020d7557ddf4bb3b5910355d58ca(meta_buf_173);
const int meta_buf_174[] = {25781352,25781352,25521320,49,512};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_174);
const int meta_buf_175[] = {25781352,25781352,25088};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_175);
const int meta_buf_176[] = {25781352,25806440,1,512,7,7,7,7,3,3,1,1,1,1};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_176);
const int meta_buf_177[] = {25806440,4718592,25781352,49,512,4608};
sgemm_675f43f38f11c7e8b176da63f208afead50e062474e037d7e5569bbf(meta_buf_177);
const int meta_buf_178[] = {25781352,25781352,25522856,49,512};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_178);
const int meta_buf_179[] = {25781352,25781352,25088};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_179);
const int meta_buf_180[] = {25781352,14368768,25956968,49,2048,512};
sgemm_4a7b51005273e75f1cab8ba4b77dc86c5ce53b149a963daaa5d423e1(meta_buf_180);
const int meta_buf_181[] = {25956968,25956968,25507008,49,2048};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_181);
const int meta_buf_182[] = {25956968,25681000,25956968,100352};
elementwisesum_54db89b5b7579805e7f7abdf1dd69f582f55da0ff40728276f251c8f(meta_buf_182);
const int meta_buf_183[] = {25956968,25956968,100352};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_183);
const int meta_buf_184[] = {25956968,11223040,25681000,49,512,2048};
sgemm_187235ada597108e1eb448f6650b020d7557ddf4bb3b5910355d58ca(meta_buf_184);
const int meta_buf_185[] = {25681000,25681000,25523368,49,512};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_185);
const int meta_buf_186[] = {25681000,25681000,25088};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_186);
const int meta_buf_187[] = {25681000,25706088,1,512,7,7,7,7,3,3,1,1,1,1};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_187);
const int meta_buf_188[] = {25706088,0,25931880,49,512,4608};
sgemm_675f43f38f11c7e8b176da63f208afead50e062474e037d7e5569bbf(meta_buf_188);
const int meta_buf_189[] = {25931880,25931880,25524392,49,512};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_189);
const int meta_buf_190[] = {25931880,25931880,25088};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_190);
const int meta_buf_191[] = {25931880,13320192,25681000,49,2048,512};
sgemm_4a7b51005273e75f1cab8ba4b77dc86c5ce53b149a963daaa5d423e1(meta_buf_191);
const int meta_buf_192[] = {25681000,25681000,25509056,49,2048};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_192);
const int meta_buf_193[] = {25681000,25956968,25681000,100352};
elementwisesum_54db89b5b7579805e7f7abdf1dd69f582f55da0ff40728276f251c8f(meta_buf_193);
const int meta_buf_194[] = {25681000,25681000,100352};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_194);
const int meta_buf_195[] = {25681000,25781352,1,7,7,2048,1,1,7,1,0};
averagepooling2d_203c4e21c29b022faa72777a679127ac4aec4267b6caaf17b363a270(meta_buf_195);
const int meta_buf_196[] = {25781352,25781352,2048};
flatten_143a2725e849669b2cf7cd13c75f75d2068fc136b9235e3aaf6cab38(meta_buf_196);
const int meta_buf_197[] = {25781352,9175040,25681000,1,1000,2048};
sgemm_814127e7734e5b079739c0a1d49f4231b555135cd2d42721305bb3f6(meta_buf_197);
const int meta_buf_198[] = {25681000,25681000,25518272,1,1000};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_198);

}

