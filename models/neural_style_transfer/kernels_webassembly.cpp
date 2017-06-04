
#include <stdlib.h>
#include <math.h>

float data_buffer[9506180];


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

void sgemm_4eaf1def8d81f03db9cf692d12bd6194c1a87e870b40959440607085(const int * meta_buffer)
{
    float *A = data_buffer + meta_buffer[0];
    float *B = data_buffer + meta_buffer[1];
    float *C = data_buffer + meta_buffer[2];

    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > a_mat(A, 27648, 243);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > b_mat(B, 243, 32);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > c_mat(C, 27648, 32);
    //Eigen::Map<Eigen::Matrix<float, 27648, 243, Eigen::RowMajor> > a_mat(A, 27648, 243);
    //Eigen::Map<Eigen::Matrix<float, 243, 32, Eigen::RowMajor> > b_mat(B, 243, 32);
    //Eigen::Map<Eigen::Matrix<float, 27648, 32, Eigen::RowMajor> > c_mat(C, 27648, 32);

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


void elu_ada847a5922e416d701697bb1608e40cdaa9da9fd5b7c6cfe1a12822(const int * meta_buffer)
{
    const float *X = data_buffer + meta_buffer[0];
    float *Y = data_buffer + meta_buffer[1];

    const int N = meta_buffer[2];
  
    for (int gid = 0; gid < N; gid += 1) {
        float result = X[gid];
        result = result < 0.0 ? (expf(result)-1) : result;      
        Y[gid] = result;
    }
}


void axiswisescale_888c8fbdc4e527683d1ad51de41d17069b1bef2e5b2751baf0a07767(const int * meta_buffer)
{
    const float *X = data_buffer + meta_buffer[0];
    float *Y = data_buffer + meta_buffer[1];
    const float *S = data_buffer + meta_buffer[2];
    const int N = meta_buffer[3];
    const int C = meta_buffer[4];
  
    for (int gid = 0; gid < N; gid += 1) {
        int c = gid % C;
        float result = X[gid] * S[c];

        Y[gid] = result;
    }
}


#ifndef INCLUDE_EIGEN
#define INCLUDE_EIGEN
#include <Eigen/Dense>
#endif

void sgemm_cdcf209b5504ee1dba4af4471a0d26801273c19d0f5fbd08fa814da8(const int * meta_buffer)
{
    float *A = data_buffer + meta_buffer[0];
    float *B = data_buffer + meta_buffer[1];
    float *C = data_buffer + meta_buffer[2];

    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > a_mat(A, 6912, 512);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > b_mat(B, 512, 64);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > c_mat(C, 6912, 64);
    //Eigen::Map<Eigen::Matrix<float, 6912, 512, Eigen::RowMajor> > a_mat(A, 6912, 512);
    //Eigen::Map<Eigen::Matrix<float, 512, 64, Eigen::RowMajor> > b_mat(B, 512, 64);
    //Eigen::Map<Eigen::Matrix<float, 6912, 64, Eigen::RowMajor> > c_mat(C, 6912, 64);

    c_mat.noalias() = a_mat * b_mat;
}


#ifndef INCLUDE_EIGEN
#define INCLUDE_EIGEN
#include <Eigen/Dense>
#endif

void sgemm_a3649a498ff69f8e59f76baec30f77ef3130069d36fd241d9e27296e(const int * meta_buffer)
{
    float *A = data_buffer + meta_buffer[0];
    float *B = data_buffer + meta_buffer[1];
    float *C = data_buffer + meta_buffer[2];

    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > a_mat(A, 1728, 1024);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > b_mat(B, 1024, 128);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > c_mat(C, 1728, 128);
    //Eigen::Map<Eigen::Matrix<float, 1728, 1024, Eigen::RowMajor> > a_mat(A, 1728, 1024);
    //Eigen::Map<Eigen::Matrix<float, 1024, 128, Eigen::RowMajor> > b_mat(B, 1024, 128);
    //Eigen::Map<Eigen::Matrix<float, 1728, 128, Eigen::RowMajor> > c_mat(C, 1728, 128);

    c_mat.noalias() = a_mat * b_mat;
}


#ifndef INCLUDE_EIGEN
#define INCLUDE_EIGEN
#include <Eigen/Dense>
#endif

void sgemm_98b4dd858fe746aad1c9191bf2cde6d0ed574219ed43ca663c6f93b0(const int * meta_buffer)
{
    float *A = data_buffer + meta_buffer[0];
    float *B = data_buffer + meta_buffer[1];
    float *C = data_buffer + meta_buffer[2];

    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > a_mat(A, 1728, 1152);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > b_mat(B, 1152, 128);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > c_mat(C, 1728, 128);
    //Eigen::Map<Eigen::Matrix<float, 1728, 1152, Eigen::RowMajor> > a_mat(A, 1728, 1152);
    //Eigen::Map<Eigen::Matrix<float, 1152, 128, Eigen::RowMajor> > b_mat(B, 1152, 128);
    //Eigen::Map<Eigen::Matrix<float, 1728, 128, Eigen::RowMajor> > c_mat(C, 1728, 128);

    c_mat.noalias() = a_mat * b_mat;
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

void sgemm_0b7d0836e58963ea16a2d82a9f6b4970180ab32ecf36fdbf088fe621(const int * meta_buffer)
{
    float *A = data_buffer + meta_buffer[0];
    float *B = data_buffer + meta_buffer[1];
    float *C = data_buffer + meta_buffer[2];

    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > a_mat(A, 1728, 128);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > b_mat(B, 128, 1024);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > c_mat(C, 1728, 1024);
    //Eigen::Map<Eigen::Matrix<float, 1728, 128, Eigen::RowMajor> > a_mat(A, 1728, 128);
    //Eigen::Map<Eigen::Matrix<float, 128, 1024, Eigen::RowMajor> > b_mat(B, 128, 1024);
    //Eigen::Map<Eigen::Matrix<float, 1728, 1024, Eigen::RowMajor> > c_mat(C, 1728, 1024);

    c_mat.noalias() = a_mat * b_mat;
}


void col2im_085c6157f908ef0c5ce49834baf86a62092a5fb7eade79e19107af95(const int * meta_buffer)
{
    const float *col = data_buffer + meta_buffer[1];
    float *im = data_buffer + meta_buffer[0];

    const int N = meta_buffer[2];
    const int C1 = meta_buffer[5];
    const int H1 = meta_buffer[6];
    const int W1 = meta_buffer[7];
    const int H2 = meta_buffer[3];
    const int W2 = meta_buffer[4];
    const int KH = meta_buffer[8];
    const int KW = meta_buffer[9];
    const int SH = meta_buffer[10];
    const int SW = meta_buffer[11];
    const int PH = meta_buffer[12];
    const int PW = meta_buffer[13];

    for (int gid = 0; gid < N*H1*W1*C1; gid += 1) {
        const int c1 = gid % C1;
        const int w1 = gid / C1 % W1;
        const int h1 = gid / C1 / W1 % H1;
        const int n = gid / C1 / W1 / H1;

        float sum = 0;
        for (int kh = 0; kh < KH; kh++) {
            const int h2 = (h1 + PH - kh) / SH;
            if ((h1 + PH - kh) % SH != 0 || h2 < 0 || h2 >= H2) continue;

            for (int kw = 0; kw < KW; kw++) {
                const int w2 = (w1 + PW - kw) / SW;
                if ((w1 + PW - kw) % SW != 0 || w2 < 0 || w2 >= W2) continue;
                
                sum += col[((((n * H2 + h2) * W2 + w2) * KH + kh) * KW + kw) * C1 + c1];
            }
        }
        
        im[gid] = sum; 
    }
}


#ifndef INCLUDE_EIGEN
#define INCLUDE_EIGEN
#include <Eigen/Dense>
#endif

void sgemm_8ff030cc50aadf1b4364f9709856b31e5d0cb9229552b9fb5db8a90e(const int * meta_buffer)
{
    float *A = data_buffer + meta_buffer[0];
    float *B = data_buffer + meta_buffer[1];
    float *C = data_buffer + meta_buffer[2];

    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > a_mat(A, 6912, 64);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > b_mat(B, 64, 512);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > c_mat(C, 6912, 512);
    //Eigen::Map<Eigen::Matrix<float, 6912, 64, Eigen::RowMajor> > a_mat(A, 6912, 64);
    //Eigen::Map<Eigen::Matrix<float, 64, 512, Eigen::RowMajor> > b_mat(B, 64, 512);
    //Eigen::Map<Eigen::Matrix<float, 6912, 512, Eigen::RowMajor> > c_mat(C, 6912, 512);

    c_mat.noalias() = a_mat * b_mat;
}


#ifndef INCLUDE_EIGEN
#define INCLUDE_EIGEN
#include <Eigen/Dense>
#endif

void sgemm_8bac6b06f06b81759c72f3595bb8554b912e7a586a2be63ace7cb9da(const int * meta_buffer)
{
    float *A = data_buffer + meta_buffer[0];
    float *B = data_buffer + meta_buffer[1];
    float *C = data_buffer + meta_buffer[2];

    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > a_mat(A, 27648, 32);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > b_mat(B, 32, 243);
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > c_mat(C, 27648, 243);
    //Eigen::Map<Eigen::Matrix<float, 27648, 32, Eigen::RowMajor> > a_mat(A, 27648, 32);
    //Eigen::Map<Eigen::Matrix<float, 32, 243, Eigen::RowMajor> > b_mat(B, 32, 243);
    //Eigen::Map<Eigen::Matrix<float, 27648, 243, Eigen::RowMajor> > c_mat(C, 27648, 243);

    c_mat.noalias() = a_mat * b_mat;
}


void tanh_a44cac43865bb3d61f261d70e71eac1e9efb1e219a8f6e72eec179c7(const int * meta_buffer)
{
    const float *X = data_buffer + meta_buffer[0];
    float *Y = data_buffer + meta_buffer[1];

    const int N = meta_buffer[2];
  
    for (int gid = 0; gid < N; gid += 1) {
        float result = X[gid];
        result = tanhf(result);      
        Y[gid] = result;
    }
}


void scalaraffine_dffbd2c541f7ae17936a12bebd7bd61ed5ea92c262d4fc4d57856c64(const int * meta_buffer)
{
    const float *X = data_buffer + meta_buffer[0];
    float *Y = data_buffer + meta_buffer[1];

    const float scale = *((const float *)(& meta_buffer[3]));
    const float bias = *((const float *)(& meta_buffer[4]));
    const int N = meta_buffer[2];

    for (int gid = 0; gid < N; gid += 1) {
        float result = X[gid];
        result = result * scale + bias;

        Y[gid] = result;
    }
}

extern "C" void init() {
    //data_buffer = (float*)malloc(9506180 * sizeof(float));
}

extern "C" float* get_data_buffer(void) {
    return data_buffer;
}



extern "C" void run() {
const int meta_buf_0[] = {1820036,2787716,1,3,144,192,144,192,9,9,1,1,4,4};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_0);
const int meta_buf_1[] = {2787716,1810016,1902980,27648,32,243};
sgemm_4eaf1def8d81f03db9cf692d12bd6194c1a87e870b40959440607085(meta_buf_1);
const int meta_buf_2[] = {1902980,1902980,1819840,27648,32};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_2);
const int meta_buf_3[] = {1902980,1902980,884736};
elu_ada847a5922e416d701697bb1608e40cdaa9da9fd5b7c6cfe1a12822(meta_buf_3);
const int meta_buf_4[] = {1902980,1902980,1819936,884736,32};
axiswisescale_888c8fbdc4e527683d1ad51de41d17069b1bef2e5b2751baf0a07767(meta_buf_4);
const int meta_buf_5[] = {1902980,1902980,1819968,27648,32};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_5);
const int meta_buf_6[] = {1902980,2787716,1,32,144,192,72,96,4,4,2,2,1,1};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_6);
const int meta_buf_7[] = {2787716,1769472,1902980,6912,64,512};
sgemm_cdcf209b5504ee1dba4af4471a0d26801273c19d0f5fbd08fa814da8(meta_buf_7);
const int meta_buf_8[] = {1902980,1902980,1819520,6912,64};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_8);
const int meta_buf_9[] = {1902980,1902980,442368};
elu_ada847a5922e416d701697bb1608e40cdaa9da9fd5b7c6cfe1a12822(meta_buf_9);
const int meta_buf_10[] = {1902980,1902980,1819712,442368,64};
axiswisescale_888c8fbdc4e527683d1ad51de41d17069b1bef2e5b2751baf0a07767(meta_buf_10);
const int meta_buf_11[] = {1902980,1902980,1819584,6912,64};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_11);
const int meta_buf_12[] = {1902980,2345348,1,64,72,96,36,48,4,4,2,2,1,1};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_12);
const int meta_buf_13[] = {2345348,1474560,1902980,1728,128,1024};
sgemm_a3649a498ff69f8e59f76baec30f77ef3130069d36fd241d9e27296e(meta_buf_13);
const int meta_buf_14[] = {1902980,1902980,1818560,1728,128};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_14);
const int meta_buf_15[] = {1902980,1902980,221184};
elu_ada847a5922e416d701697bb1608e40cdaa9da9fd5b7c6cfe1a12822(meta_buf_15);
const int meta_buf_16[] = {1902980,1902980,1818432,221184,128};
axiswisescale_888c8fbdc4e527683d1ad51de41d17069b1bef2e5b2751baf0a07767(meta_buf_16);
const int meta_buf_17[] = {1902980,1902980,1818816,1728,128};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_17);
const int meta_buf_18[] = {1902980,2345348,1,128,36,48,36,48,3,3,1,1,1,1};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_18);
const int meta_buf_19[] = {2345348,0,2124164,1728,128,1152};
sgemm_98b4dd858fe746aad1c9191bf2cde6d0ed574219ed43ca663c6f93b0(meta_buf_19);
const int meta_buf_20[] = {2124164,2124164,1819200,1728,128};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_20);
const int meta_buf_21[] = {2124164,2124164,221184};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_21);
const int meta_buf_22[] = {2124164,2345348,1,128,36,48,36,48,3,3,1,1,1,1};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_22);
const int meta_buf_23[] = {2345348,442368,4336004,1728,128,1152};
sgemm_98b4dd858fe746aad1c9191bf2cde6d0ed574219ed43ca663c6f93b0(meta_buf_23);
const int meta_buf_24[] = {4336004,4336004,1818048,1728,128};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_24);
const int meta_buf_25[] = {4336004,1902980,4336004,221184};
elementwisesum_54db89b5b7579805e7f7abdf1dd69f582f55da0ff40728276f251c8f(meta_buf_25);
const int meta_buf_26[] = {4336004,2124164,1,128,36,48,36,48,3,3,1,1,1,1};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_26);
const int meta_buf_27[] = {2124164,294912,1902980,1728,128,1152};
sgemm_98b4dd858fe746aad1c9191bf2cde6d0ed574219ed43ca663c6f93b0(meta_buf_27);
const int meta_buf_28[] = {1902980,1902980,1817920,1728,128};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_28);
const int meta_buf_29[] = {1902980,1902980,221184};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_29);
const int meta_buf_30[] = {1902980,2124164,1,128,36,48,36,48,3,3,1,1,1,1};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_30);
const int meta_buf_31[] = {2124164,589824,1902980,1728,128,1152};
sgemm_98b4dd858fe746aad1c9191bf2cde6d0ed574219ed43ca663c6f93b0(meta_buf_31);
const int meta_buf_32[] = {1902980,1902980,1819072,1728,128};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_32);
const int meta_buf_33[] = {1902980,4336004,1902980,221184};
elementwisesum_54db89b5b7579805e7f7abdf1dd69f582f55da0ff40728276f251c8f(meta_buf_33);
const int meta_buf_34[] = {1902980,2345348,1,128,36,48,36,48,3,3,1,1,1,1};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_34);
const int meta_buf_35[] = {2345348,1032192,2124164,1728,128,1152};
sgemm_98b4dd858fe746aad1c9191bf2cde6d0ed574219ed43ca663c6f93b0(meta_buf_35);
const int meta_buf_36[] = {2124164,2124164,1817792,1728,128};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_36);
const int meta_buf_37[] = {2124164,2124164,221184};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_37);
const int meta_buf_38[] = {2124164,2345348,1,128,36,48,36,48,3,3,1,1,1,1};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_38);
const int meta_buf_39[] = {2345348,884736,4336004,1728,128,1152};
sgemm_98b4dd858fe746aad1c9191bf2cde6d0ed574219ed43ca663c6f93b0(meta_buf_39);
const int meta_buf_40[] = {4336004,4336004,1818688,1728,128};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_40);
const int meta_buf_41[] = {4336004,1902980,4336004,221184};
elementwisesum_54db89b5b7579805e7f7abdf1dd69f582f55da0ff40728276f251c8f(meta_buf_41);
const int meta_buf_42[] = {4336004,2124164,1,128,36,48,36,48,3,3,1,1,1,1};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_42);
const int meta_buf_43[] = {2124164,737280,1902980,1728,128,1152};
sgemm_98b4dd858fe746aad1c9191bf2cde6d0ed574219ed43ca663c6f93b0(meta_buf_43);
const int meta_buf_44[] = {1902980,1902980,1819328,1728,128};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_44);
const int meta_buf_45[] = {1902980,1902980,221184};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_45);
const int meta_buf_46[] = {1902980,2124164,1,128,36,48,36,48,3,3,1,1,1,1};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_46);
const int meta_buf_47[] = {2124164,1179648,1902980,1728,128,1152};
sgemm_98b4dd858fe746aad1c9191bf2cde6d0ed574219ed43ca663c6f93b0(meta_buf_47);
const int meta_buf_48[] = {1902980,1902980,1818176,1728,128};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_48);
const int meta_buf_49[] = {1902980,4336004,1902980,221184};
elementwisesum_54db89b5b7579805e7f7abdf1dd69f582f55da0ff40728276f251c8f(meta_buf_49);
const int meta_buf_50[] = {1902980,2345348,1,128,36,48,36,48,3,3,1,1,1,1};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_50);
const int meta_buf_51[] = {2345348,147456,2124164,1728,128,1152};
sgemm_98b4dd858fe746aad1c9191bf2cde6d0ed574219ed43ca663c6f93b0(meta_buf_51);
const int meta_buf_52[] = {2124164,2124164,1818944,1728,128};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_52);
const int meta_buf_53[] = {2124164,2124164,221184};
relu_98667a12ddf6cd220e5632bfe6c4affc7f55eb5b3ad2138513486a8a(meta_buf_53);
const int meta_buf_54[] = {2124164,2345348,1,128,36,48,36,48,3,3,1,1,1,1};
im2col_bc00b53e478359c9136f25ec94004097a9f9a23b63fa471d7d94c024(meta_buf_54);
const int meta_buf_55[] = {2345348,1327104,2124164,1728,128,1152};
sgemm_98b4dd858fe746aad1c9191bf2cde6d0ed574219ed43ca663c6f93b0(meta_buf_55);
const int meta_buf_56[] = {2124164,2124164,1818304,1728,128};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_56);
const int meta_buf_57[] = {2124164,1902980,2124164,221184};
elementwisesum_54db89b5b7579805e7f7abdf1dd69f582f55da0ff40728276f251c8f(meta_buf_57);
const int meta_buf_58[] = {2124164,1605632,2345348,1728,1024,128};
sgemm_0b7d0836e58963ea16a2d82a9f6b4970180ab32ecf36fdbf088fe621(meta_buf_58);
const int meta_buf_59[] = {1902980,2345348,1,36,48,64,72,96,4,4,2,2,1,1};
col2im_085c6157f908ef0c5ce49834baf86a62092a5fb7eade79e19107af95(meta_buf_59);
const int meta_buf_60[] = {1902980,1902980,1819456,6912,64};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_60);
const int meta_buf_61[] = {1902980,1902980,442368};
elu_ada847a5922e416d701697bb1608e40cdaa9da9fd5b7c6cfe1a12822(meta_buf_61);
const int meta_buf_62[] = {1902980,1902980,1819648,442368,64};
axiswisescale_888c8fbdc4e527683d1ad51de41d17069b1bef2e5b2751baf0a07767(meta_buf_62);
const int meta_buf_63[] = {1902980,1902980,1819776,6912,64};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_63);
const int meta_buf_64[] = {1902980,1736704,2787716,6912,512,64};
sgemm_8ff030cc50aadf1b4364f9709856b31e5d0cb9229552b9fb5db8a90e(meta_buf_64);
const int meta_buf_65[] = {1902980,2787716,1,72,96,32,144,192,4,4,2,2,1,1};
col2im_085c6157f908ef0c5ce49834baf86a62092a5fb7eade79e19107af95(meta_buf_65);
const int meta_buf_66[] = {1902980,1902980,1820000,27648,32};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_66);
const int meta_buf_67[] = {1902980,1902980,884736};
elu_ada847a5922e416d701697bb1608e40cdaa9da9fd5b7c6cfe1a12822(meta_buf_67);
const int meta_buf_68[] = {1902980,1902980,1819872,884736,32};
axiswisescale_888c8fbdc4e527683d1ad51de41d17069b1bef2e5b2751baf0a07767(meta_buf_68);
const int meta_buf_69[] = {1902980,1902980,1819904,27648,32};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_69);
const int meta_buf_70[] = {1902980,1802240,2787716,27648,243,32};
sgemm_8bac6b06f06b81759c72f3595bb8554b912e7a586a2be63ace7cb9da(meta_buf_70);
const int meta_buf_71[] = {1902980,2787716,1,144,192,3,144,192,9,9,1,1,4,4};
col2im_085c6157f908ef0c5ce49834baf86a62092a5fb7eade79e19107af95(meta_buf_71);
const int meta_buf_72[] = {1902980,1902980,1820032,27648,3};
axiswisebias_17c5d4656f48efee4d8b12e4489079ce064a004f3cd635926548b2f6(meta_buf_72);
const int meta_buf_73[] = {1902980,1902980,82944};
tanh_a44cac43865bb3d61f261d70e71eac1e9efb1e219a8f6e72eec179c7(meta_buf_73);
const int meta_buf_74[] = {1902980,1902980,82944,1124007936,1124007936};
scalaraffine_dffbd2c541f7ae17936a12bebd7bd61ed5ea92c262d4fc4d57856c64(meta_buf_74);

}

