
#include <metal_stdlib>
using namespace metal;

#define OPTIMIZE 1


kernel void im2col_c11bc86fe77290c624886ebd52e1e88ff9cb3655ad3bfada8b15ad55(const device float *param_buffer[[buffer(0)]],
                          device float *data_buffer[[buffer(1)]],
                          const device int * meta_buffer [[buffer(2)]],
                          ushort index_thread[[thread_position_in_threadgroup]],
                          ushort index_group[[threadgroup_position_in_grid]])
{
#define SH_EQUAL_1 1
#define SW_EQUAL_1 1
#define C1_DIVIDABLE_BY_4 0


#if OPTIMIZE && C1_DIVIDABLE_BY_4
    const device float4 *im4 = (const device float4 *)(data_buffer + meta_buffer[0]);
    device float4 *col4 = (device float4 *)(data_buffer + meta_buffer[1]);
    const int C1_4 = (meta_buffer[3]) >> 2;
#else
    const device float *im = data_buffer + meta_buffer[0];
    device float *col = data_buffer + meta_buffer[1];
    const int C1 = meta_buffer[3];
#endif

    // const int N = meta_buffer[2];
    const int H1 = meta_buffer[4];
    const int W1 = meta_buffer[5];
    const int H2 = meta_buffer[6];
    const int W2 = meta_buffer[7];
    const int KH = meta_buffer[8];
    const int KW = meta_buffer[9];
    const int PH = meta_buffer[12];
    const int PW = meta_buffer[13];

#if !OPTIMIZE || !SH_EQUAL_1
    const int SH = meta_buffer[10];
#endif

#if !OPTIMIZE || !SW_EQUAL_1
    const int SW = meta_buffer[11];
#endif

    const int H1P = H1 + 2 * PH;
    const int W1P = W1 + 2 * PW;

    const int w1 = (index_group % W1P) - PW;
    const int h1 = (index_group / W1P % H1P) - PH;
    const int  n = index_group / W1P / H1P;

#if OPTIMIZE && C1_DIVIDABLE_BY_4
    for (int c1_4 = index_thread; c1_4 < C1_4; c1_4 += 64) {
        const float4 v4 = (h1 < 0 || h1 >= H1 || w1 < 0 || w1 >= W1) ? 0 : im4[((n * H1 + h1) * W1 + w1) * C1_4 + c1_4];
#else
    for (int c1 = index_thread; c1 < C1; c1 += 64) {
        const float v = (h1 < 0 || h1 >= H1 || w1 < 0 || w1 >= W1) ? 0 : im[((n * H1 + h1) * W1 + w1) * C1 + c1];
#endif

#if OPTIMIZE && SH_EQUAL_1
        for (int kh = 0; kh < KH; kh++) {
            const int h2 = h1 + PH - kh;
#else
        for (int kh = (h1 + PH) % SH; kh < KH; kh += SH) {
            const int h2 = (h1 + PH - kh) / SH;
#endif
            if (h2 < 0 || h2 >= H2) continue;

#if OPTIMIZE && SH_EQUAL_1
            for (int kw = 0; kw < KW; kw++) {
                const int w2 = w1 + PW - kw;
#else
            for (int kw = (w1 + PW) % SW; kw < KW; kw += SW) {
                const int w2 = (w1 + PW - kw) / SW;
#endif
                if (w2 < 0 || w2 >= W2) continue;

#if OPTIMIZE && C1_DIVIDABLE_BY_4
                col4[((((n * H2 + h2) * W2 + w2) * KH + kh) * KW + kw) * C1_4 + c1_4] = v4;
#else
                col[((((n * H2 + h2) * W2 + w2) * KH + kh) * KW + kw) * C1 + c1] = v;
#endif
            }
        }
    }


#undef SH_EQUAL_1
#undef SW_EQUAL_1
#undef C1_DIVIDABLE_BY_4
}


kernel void sgemm_1fe1c908701210c37c7cdbeb4c16e81d8255a13579dd647fe858078f(const device float *weight_buffer[[buffer(0)]],
                          device float *data_buffer[[buffer(1)]],
                          const device int * meta_buffer [[buffer(2)]],
                          ushort index[[thread_index_in_threadgroup]],
                          ushort2 group_position[[threadgroup_position_in_grid]])
{
#define TRANSPOSE_A 1
#define TRANSPOSE_B 1
#define M_DIVIDABLE_BY_64 1
#define N_DIVIDABLE_BY_64 0
#define K_DIVIDABLE_BY_8 0

#if TRANSPOSE_A
    #define A_STRIDE_K 1
    #define A_STRIDE_M K
#else
    #define A_STRIDE_K M
    #define A_STRIDE_M 1
#endif

#if TRANSPOSE_B
    #define B_STRIDE_K N
    #define B_STRIDE_N 1
#else
    #define B_STRIDE_K 1
    #define B_STRIDE_N K
#endif

#define WITH_BIAS 1
#define HAS_INLINE 1


#if K_DIVIDABLE_BY_8 && M_DIVIDABLE_BY_64  && N_DIVIDABLE_BY_64 && !TRANSPOSE_A && TRANSPOSE_B && OPTIMIZE
    const device float4 *load_target4 = (index & 32) 
        ? (const device float4 *)(weight_buffer + meta_buffer[1]) 
        : (const device float4 *)(data_buffer + meta_buffer[0]);
#else
    const device float *load_target = (index & 32) 
        ? (weight_buffer + meta_buffer[1]) 
        : (data_buffer + meta_buffer[0]);
#endif

    const int M = meta_buffer[4];
    const int N = meta_buffer[5];

    const int K = meta_buffer[6];

    threadgroup float4 shared4[32 * 8 * 2];

    const int stride_k = (index & 32) ? B_STRIDE_K : A_STRIDE_K;
    const int stride_mn = (index & 32) ? B_STRIDE_N : A_STRIDE_M;

    const int m_offset = index & 7;
    const int n_offset = index >> 3;

    const int mn_load_offset = ((index & 32) ? group_position.y : group_position.x) * 64 + (index & 15) * 4;
    const int k_load_offset = ((index & 16) ? 1 : 0);

    int track0 = mn_load_offset * stride_mn + (k_load_offset + 0) * stride_k;
    int track2 = track0 + 2 * stride_k;
    int track4 = track0 + 4 * stride_k;
    int track6 = track0 + 6 * stride_k;

#if !OPTIMIZE || !M_DIVIDABLE_BY_64 || !N_DIVIDABLE_BY_64
    const int max_MN = (index & 32) ? N : M;
#endif

    int shared_offset4 = ((index & 32) ? 16 : 0) + k_load_offset * 32 + (index & 15);
    int read_A_offset4 = m_offset * 2;
    int read_B_offset4 = n_offset * 2 + 16;

    float4 result[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int k = 0;

    while (k < K)
    {
        {
#if OPTIMIZE && K_DIVIDABLE_BY_8
    #if OPTIMIZE && M_DIVIDABLE_BY_64 && N_DIVIDABLE_BY_64
        #if OPTIMIZE && !TRANSPOSE_A && TRANSPOSE_B
            shared4[shared_offset4 + 32 * 0] = load_target4[track0 >> 2];
            shared4[shared_offset4 + 32 * 2] = load_target4[track2 >> 2];
            shared4[shared_offset4 + 32 * 4] = load_target4[track4 >> 2];
            shared4[shared_offset4 + 32 * 6] = load_target4[track6 >> 2];
        #else
            shared4[shared_offset4 + 32 * 0] = float4(
                load_target[track0 + stride_mn * 0],
                load_target[track0 + stride_mn * 1],
                load_target[track0 + stride_mn * 2],
                load_target[track0 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 2] = float4(
                load_target[track2 + stride_mn * 0],
                load_target[track2 + stride_mn * 1],
                load_target[track2 + stride_mn * 2],
                load_target[track2 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 4] = float4(
                load_target[track4 + stride_mn * 0],
                load_target[track4 + stride_mn * 1],
                load_target[track4 + stride_mn * 2],
                load_target[track4 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 6] = float4(
                load_target[track6 + stride_mn * 0],
                load_target[track6 + stride_mn * 1],
                load_target[track6 + stride_mn * 2],
                load_target[track6 + stride_mn * 3]
            ); 
        #endif
    #else
            shared4[shared_offset4 + 32 * 0] = float4(
                (mn_load_offset + 0 >= max_MN) ? 0 : load_target[track0 + stride_mn * 0],
                (mn_load_offset + 1 >= max_MN) ? 0 : load_target[track0 + stride_mn * 1],
                (mn_load_offset + 2 >= max_MN) ? 0 : load_target[track0 + stride_mn * 2],
                (mn_load_offset + 3 >= max_MN) ? 0 : load_target[track0 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 2] = float4(
                (mn_load_offset + 0 >= max_MN) ? 0 : load_target[track2 + stride_mn * 0],
                (mn_load_offset + 1 >= max_MN) ? 0 : load_target[track2 + stride_mn * 1],
                (mn_load_offset + 2 >= max_MN) ? 0 : load_target[track2 + stride_mn * 2],
                (mn_load_offset + 3 >= max_MN) ? 0 : load_target[track2 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 4] = float4(
                (mn_load_offset + 0 >= max_MN) ? 0 : load_target[track4 + stride_mn * 0],
                (mn_load_offset + 1 >= max_MN) ? 0 : load_target[track4 + stride_mn * 1],
                (mn_load_offset + 2 >= max_MN) ? 0 : load_target[track4 + stride_mn * 2],
                (mn_load_offset + 3 >= max_MN) ? 0 : load_target[track4 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 6] = float4(
                (mn_load_offset + 0 >= max_MN) ? 0 : load_target[track6 + stride_mn * 0],
                (mn_load_offset + 1 >= max_MN) ? 0 : load_target[track6 + stride_mn * 1],
                (mn_load_offset + 2 >= max_MN) ? 0 : load_target[track6 + stride_mn * 2],
                (mn_load_offset + 3 >= max_MN) ? 0 : load_target[track6 + stride_mn * 3]
            ); 
    #endif

            k += 8;
#else
    #if OPTIMIZE && M_DIVIDABLE_BY_64 && N_DIVIDABLE_BY_64
            shared4[shared_offset4 + 32 * 0] = float4(
                (k + k_load_offset >= K) ? 0 : load_target[track0 + stride_mn * 0],
                (k + k_load_offset >= K) ? 0 : load_target[track0 + stride_mn * 1],
                (k + k_load_offset >= K) ? 0 : load_target[track0 + stride_mn * 2],
                (k + k_load_offset >= K) ? 0 : load_target[track0 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 2] = float4(
                (k + k_load_offset >= K) ? 0 : load_target[track2 + stride_mn * 0],
                (k + k_load_offset >= K) ? 0 : load_target[track2 + stride_mn * 1],
                (k + k_load_offset >= K) ? 0 : load_target[track2 + stride_mn * 2],
                (k + k_load_offset >= K) ? 0 : load_target[track2 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 4] = float4(
                (k + k_load_offset >= K) ? 0 : load_target[track4 + stride_mn * 0],
                (k + k_load_offset >= K) ? 0 : load_target[track4 + stride_mn * 1],
                (k + k_load_offset >= K) ? 0 : load_target[track4 + stride_mn * 2],
                (k + k_load_offset >= K) ? 0 : load_target[track4 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 6] = float4(
                (k + k_load_offset >= K) ? 0 : load_target[track6 + stride_mn * 0],
                (k + k_load_offset >= K) ? 0 : load_target[track6 + stride_mn * 1],
                (k + k_load_offset >= K) ? 0 : load_target[track6 + stride_mn * 2],
                (k + k_load_offset >= K) ? 0 : load_target[track6 + stride_mn * 3]
            ); 
            k += 2;
    #else
            shared4[shared_offset4 + 32 * 0] = float4(
                (k + k_load_offset >= K || mn_load_offset + 0 >= max_MN) ? 0 : load_target[track0 + stride_mn * 0],
                (k + k_load_offset >= K || mn_load_offset + 1 >= max_MN) ? 0 : load_target[track0 + stride_mn * 1],
                (k + k_load_offset >= K || mn_load_offset + 2 >= max_MN) ? 0 : load_target[track0 + stride_mn * 2],
                (k + k_load_offset >= K || mn_load_offset + 3 >= max_MN) ? 0 : load_target[track0 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 2] = float4(
                (k + k_load_offset >= K || mn_load_offset + 0 >= max_MN) ? 0 : load_target[track2 + stride_mn * 0],
                (k + k_load_offset >= K || mn_load_offset + 1 >= max_MN) ? 0 : load_target[track2 + stride_mn * 1],
                (k + k_load_offset >= K || mn_load_offset + 2 >= max_MN) ? 0 : load_target[track2 + stride_mn * 2],
                (k + k_load_offset >= K || mn_load_offset + 3 >= max_MN) ? 0 : load_target[track2 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 4] = float4(
                (k + k_load_offset >= K || mn_load_offset + 0 >= max_MN) ? 0 : load_target[track4 + stride_mn * 0],
                (k + k_load_offset >= K || mn_load_offset + 1 >= max_MN) ? 0 : load_target[track4 + stride_mn * 1],
                (k + k_load_offset >= K || mn_load_offset + 2 >= max_MN) ? 0 : load_target[track4 + stride_mn * 2],
                (k + k_load_offset >= K || mn_load_offset + 3 >= max_MN) ? 0 : load_target[track4 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 6] = float4(
                (k + k_load_offset >= K || mn_load_offset + 0 >= max_MN) ? 0 : load_target[track6 + stride_mn * 0],
                (k + k_load_offset >= K || mn_load_offset + 1 >= max_MN) ? 0 : load_target[track6 + stride_mn * 1],
                (k + k_load_offset >= K || mn_load_offset + 2 >= max_MN) ? 0 : load_target[track6 + stride_mn * 2],
                (k + k_load_offset >= K || mn_load_offset + 3 >= max_MN) ? 0 : load_target[track6 + stride_mn * 3]
            ); 
            k += 2;
    #endif
#endif
        }

        threadgroup_barrier(mem_flags::mem_threadgroup);

        {
            for (int k_sub = 0; k_sub < 8; k_sub++)
            {
                float4 a00 = shared4[32 * k_sub + read_A_offset4 + 0];
                float4 a01 = shared4[32 * k_sub + read_A_offset4 + 1];
                float4 b00 = shared4[32 * k_sub + read_B_offset4 + 0];
                float4 b01 = shared4[32 * k_sub + read_B_offset4 + 1];

                result[4][0]  += b00[0] * a00[2];
                result[4][1]  += b00[1] * a00[2];
                result[0][1]  += b00[1] * a00[0];
                result[0][0]  += b00[0] * a00[0];
                result[6][0]  += b00[0] * a00[3];
                result[6][1]  += b00[1] * a00[3];
                result[2][1]  += b00[1] * a00[1];
                result[2][0]  += b00[0] * a00[1];
                result[12][0] += b00[0] * a01[2];
                result[12][1] += b00[1] * a01[2];
                result[8][1]  += b00[1] * a01[0];
                result[8][0]  += b00[0] * a01[0];
                result[14][0] += b00[0] * a01[3];
                result[14][1] += b00[1] * a01[3];
                result[10][1] += b00[1] * a01[1];
                result[10][0] += b00[0] * a01[1];

                result[14][2] += b00[2] * a01[3];
                result[14][3] += b00[3] * a01[3];
                result[10][3] += b00[3] * a01[1];
                result[10][2] += b00[2] * a01[1];
                result[12][2] += b00[2] * a01[2];
                result[12][3] += b00[3] * a01[2];
                result[8][3]  += b00[3] * a01[0];
                result[8][2]  += b00[2] * a01[0];
                result[6][2]  += b00[2] * a00[3];
                result[6][3]  += b00[3] * a00[3];
                result[2][3]  += b00[3] * a00[1];
                result[2][2]  += b00[2] * a00[1];
                result[4][2]  += b00[2] * a00[2];
                result[4][3]  += b00[3] * a00[2];
                result[0][3]  += b00[3] * a00[0];
                result[0][2]  += b00[2] * a00[0];

                result[5][0]  += b01[0] * a00[2];
                result[5][1]  += b01[1] * a00[2];
                result[1][1]  += b01[1] * a00[0];
                result[1][0]  += b01[0] * a00[0];
                result[7][0]  += b01[0] * a00[3];
                result[7][1]  += b01[1] * a00[3];
                result[3][1]  += b01[1] * a00[1];
                result[3][0]  += b01[0] * a00[1];
                result[13][0] += b01[0] * a01[2];
                result[13][1] += b01[1] * a01[2];
                result[9][1]  += b01[1] * a01[0];
                result[9][0]  += b01[0] * a01[0];
                result[15][0] += b01[0] * a01[3];
                result[15][1] += b01[1] * a01[3];
                result[11][1] += b01[1] * a01[1];
                result[11][0] += b01[0] * a01[1];

                result[15][2] += b01[2] * a01[3];
                result[15][3] += b01[3] * a01[3];
                result[11][3] += b01[3] * a01[1];
                result[11][2] += b01[2] * a01[1];
                result[13][2] += b01[2] * a01[2];
                result[13][3] += b01[3] * a01[2];
                result[9][3]  += b01[3] * a01[0];
                result[9][2]  += b01[2] * a01[0];
                result[7][2]  += b01[2] * a00[3];
                result[7][3]  += b01[3] * a00[3];
                result[3][3]  += b01[3] * a00[1];
                result[3][2]  += b01[2] * a00[1];
                result[5][2]  += b01[2] * a00[2];
                result[5][3]  += b01[3] * a00[2];
                result[1][3]  += b01[3] * a00[0];
                result[1][2]  += b01[2] * a00[0];
            }
        }

        shared_offset4 ^= 32 * 8;
        read_A_offset4 ^= 32 * 8;
        read_B_offset4 ^= 32 * 8;
        track0 += stride_k * 8;
        track2 += stride_k * 8;
        track4 += stride_k * 8;
        track6 += stride_k * 8;
    }

    {
    
#if OPTIMIZE && N_DIVIDABLE_BY_64
    #if WITH_BIAS
        float4 b[2];
        const device float4 *bias4 = (const device float4 *)(weight_buffer + meta_buffer[3]);
        b[0] = bias4[group_position.y * 16 + n_offset * 2 + 0];
        b[1] = bias4[group_position.y * 16 + n_offset * 2 + 1];
    #endif
    
        device float4 *C4 = (device float4 *)(data_buffer + meta_buffer[2]);
        const int N4 = N >> 2;
        int m = group_position.x * 64 + m_offset * 8;
        for (int m_sub = 0; m_sub < 8; m_sub++)
        {

    #if !M_DIVIDABLE_BY_64
            if (m >= M) continue;
    #endif

            const int n = group_position.y * 16 + n_offset * 2;
            float4 result0 = result[m_sub * 2 + 0];
            float4 result1 = result[m_sub * 2 + 1];

    #if WITH_BIAS
            result0 += b[0];
            result1 += b[1];
    #endif

    #if HAS_INLINE
            result0[0] = (result0[0]>0?result0[0]:(exp(result0[0])-1));
            result0[1] = (result0[1]>0?result0[1]:(exp(result0[1])-1));
            result0[2] = (result0[2]>0?result0[2]:(exp(result0[2])-1));
            result0[3] = (result0[3]>0?result0[3]:(exp(result0[3])-1));
            result1[0] = (result1[0]>0?result1[0]:(exp(result1[0])-1));
            result1[1] = (result1[1]>0?result1[1]:(exp(result1[1])-1));
            result1[2] = (result1[2]>0?result1[2]:(exp(result1[2])-1));
            result1[3] = (result1[3]>0?result1[3]:(exp(result1[3])-1));
    #endif

            C4[m * N4 + n + 0] = result0;
            C4[m * N4 + n + 1] = result1;
            
            m++;
        }
#else
    #if WITH_BIAS
        const device float *bias = weight_buffer + meta_buffer[3];
        float b[8];
        for (int n_sub = 0; n_sub < 8; n_sub++)
        {
            b[n_sub] = (group_position.y * 64 + n_offset * 8 + n_sub < N)
                ? bias[group_position.y * 64 + n_offset * 8 + n_sub]
                : 0;
        }
    #endif

        device float *C = data_buffer + meta_buffer[2];
        int m = group_position.x * 64 + m_offset * 8;
        for (int m_sub = 0; m_sub < 8; m_sub++)
        {
            int n = group_position.y * 64 + n_offset * 8;

            for (int n_sub1 = 0; n_sub1 < 2; n_sub1++)
            {
                for (int n_sub2 = 0; n_sub2 < 4; n_sub2++)
                {

    #if WITH_BIAS
        #if OPTIMIZE && M_DIVIDABLE_BY_64
                    (         n < N) ? (C[m * N + n] = (result[m_sub * 2 + n_sub1][n_sub2] + b[n_sub1*4+n_sub2]>0?result[m_sub * 2 + n_sub1][n_sub2] + b[n_sub1*4+n_sub2]:(exp(result[m_sub * 2 + n_sub1][n_sub2] + b[n_sub1*4+n_sub2])-1))) : 0;
        #else
                    (m < M && n < N) ? (C[m * N + n] = (result[m_sub * 2 + n_sub1][n_sub2] + b[n_sub1*4+n_sub2]>0?result[m_sub * 2 + n_sub1][n_sub2] + b[n_sub1*4+n_sub2]:(exp(result[m_sub * 2 + n_sub1][n_sub2] + b[n_sub1*4+n_sub2])-1))) : 0;
        #endif
    #else
        #if OPTIMIZE && M_DIVIDABLE_BY_64
                    (         n < N) ? (C[m * N + n] = (result[m_sub * 2 + n_sub1][n_sub2]>0?result[m_sub * 2 + n_sub1][n_sub2]:(exp(result[m_sub * 2 + n_sub1][n_sub2])-1))) : 0;
        #else
                    (m < M && n < N) ? (C[m * N + n] = (result[m_sub * 2 + n_sub1][n_sub2]>0?result[m_sub * 2 + n_sub1][n_sub2]:(exp(result[m_sub * 2 + n_sub1][n_sub2])-1))) : 0;
        #endif
    #endif
                    n++;
                }
            }
            
            m++;
        }
#endif

    }


#undef M_DIVIDABLE_BY_64
#undef N_DIVIDABLE_BY_64
#undef K_DIVIDABLE_BY_8
#undef TRANSPOSE_A
#undef TRANSPOSE_B
#undef A_STRIDE_K
#undef B_STRIDE_K
#undef A_STRIDE_M
#undef A_STRIDE_M
#undef WITH_BIAS
#undef HAS_INLINE
}


kernel void axiswisescale_eb0047e7bd89ea7055f72ebafeb790e891bb2505194856bd8fe5e84d(const device float *weight_buffer[[buffer(0)]],
                          device float *data_buffer[[buffer(1)]],
                          const device int * meta_buffer [[buffer(2)]],
                          uint index[[thread_position_in_grid]],
                          uint num_threads[[threads_per_grid]])
{
    const device float *X = data_buffer + meta_buffer[0];
    device float *Y = data_buffer + meta_buffer[1];
    const device float *S = weight_buffer + meta_buffer[2];
    const int N = meta_buffer[3];
    const int C = meta_buffer[4];
  
    for (int gid = index; gid < N; gid += num_threads) {
        int c = gid % C;

        float result = X[gid] * S[c];

        Y[gid] = result;
    }
}


kernel void axiswisebias_6be2df9cf1526331046f94fd7a124df8004489b3081b9db2c7e0fa03(const device float *weight_buffer[[buffer(0)]],
                          device float *data_buffer[[buffer(1)]],
                          const device int * meta_buffer [[buffer(2)]],
                          uint index[[thread_position_in_grid]],
                          uint num_threads[[threads_per_grid]])
{
    const device float *X = data_buffer + meta_buffer[0];
    device float *Y = data_buffer + meta_buffer[1];
    const device float *B = weight_buffer + meta_buffer[2];
    const int N = meta_buffer[3];
    const int C = meta_buffer[4];
  
    for (int gid = index; gid < N * C; gid += num_threads) {
        int c = gid % C;
        int n = gid / C;

        float result = X[gid] + B[c];
        Y[n * C + c] = result;
    }
}


kernel void im2col_094f568353bf9de7c01e0b777df297bf4afc40425d719b9b676a1d9f(const device float *param_buffer[[buffer(0)]],
                          device float *data_buffer[[buffer(1)]],
                          const device int * meta_buffer [[buffer(2)]],
                          ushort index_thread[[thread_position_in_threadgroup]],
                          ushort index_group[[threadgroup_position_in_grid]])
{
#define SH_EQUAL_1 0
#define SW_EQUAL_1 0
#define C1_DIVIDABLE_BY_4 1


#if OPTIMIZE && C1_DIVIDABLE_BY_4
    const device float4 *im4 = (const device float4 *)(data_buffer + meta_buffer[0]);
    device float4 *col4 = (device float4 *)(data_buffer + meta_buffer[1]);
    const int C1_4 = (meta_buffer[3]) >> 2;
#else
    const device float *im = data_buffer + meta_buffer[0];
    device float *col = data_buffer + meta_buffer[1];
    const int C1 = meta_buffer[3];
#endif

    // const int N = meta_buffer[2];
    const int H1 = meta_buffer[4];
    const int W1 = meta_buffer[5];
    const int H2 = meta_buffer[6];
    const int W2 = meta_buffer[7];
    const int KH = meta_buffer[8];
    const int KW = meta_buffer[9];
    const int PH = meta_buffer[12];
    const int PW = meta_buffer[13];

#if !OPTIMIZE || !SH_EQUAL_1
    const int SH = meta_buffer[10];
#endif

#if !OPTIMIZE || !SW_EQUAL_1
    const int SW = meta_buffer[11];
#endif

    const int H1P = H1 + 2 * PH;
    const int W1P = W1 + 2 * PW;

    const int w1 = (index_group % W1P) - PW;
    const int h1 = (index_group / W1P % H1P) - PH;
    const int  n = index_group / W1P / H1P;

#if OPTIMIZE && C1_DIVIDABLE_BY_4
    for (int c1_4 = index_thread; c1_4 < C1_4; c1_4 += 64) {
        const float4 v4 = (h1 < 0 || h1 >= H1 || w1 < 0 || w1 >= W1) ? 0 : im4[((n * H1 + h1) * W1 + w1) * C1_4 + c1_4];
#else
    for (int c1 = index_thread; c1 < C1; c1 += 64) {
        const float v = (h1 < 0 || h1 >= H1 || w1 < 0 || w1 >= W1) ? 0 : im[((n * H1 + h1) * W1 + w1) * C1 + c1];
#endif

#if OPTIMIZE && SH_EQUAL_1
        for (int kh = 0; kh < KH; kh++) {
            const int h2 = h1 + PH - kh;
#else
        for (int kh = (h1 + PH) % SH; kh < KH; kh += SH) {
            const int h2 = (h1 + PH - kh) / SH;
#endif
            if (h2 < 0 || h2 >= H2) continue;

#if OPTIMIZE && SH_EQUAL_1
            for (int kw = 0; kw < KW; kw++) {
                const int w2 = w1 + PW - kw;
#else
            for (int kw = (w1 + PW) % SW; kw < KW; kw += SW) {
                const int w2 = (w1 + PW - kw) / SW;
#endif
                if (w2 < 0 || w2 >= W2) continue;

#if OPTIMIZE && C1_DIVIDABLE_BY_4
                col4[((((n * H2 + h2) * W2 + w2) * KH + kh) * KW + kw) * C1_4 + c1_4] = v4;
#else
                col[((((n * H2 + h2) * W2 + w2) * KH + kh) * KW + kw) * C1 + c1] = v;
#endif
            }
        }
    }


#undef SH_EQUAL_1
#undef SW_EQUAL_1
#undef C1_DIVIDABLE_BY_4
}


kernel void sgemm_7bb39dea9313ec8d8e4c2c50047d68337ba2b05833316191856f2f17(const device float *weight_buffer[[buffer(0)]],
                          device float *data_buffer[[buffer(1)]],
                          const device int * meta_buffer [[buffer(2)]],
                          ushort index[[thread_index_in_threadgroup]],
                          ushort2 group_position[[threadgroup_position_in_grid]])
{
#define TRANSPOSE_A 1
#define TRANSPOSE_B 1
#define M_DIVIDABLE_BY_64 1
#define N_DIVIDABLE_BY_64 1
#define K_DIVIDABLE_BY_8 1

#if TRANSPOSE_A
    #define A_STRIDE_K 1
    #define A_STRIDE_M K
#else
    #define A_STRIDE_K M
    #define A_STRIDE_M 1
#endif

#if TRANSPOSE_B
    #define B_STRIDE_K N
    #define B_STRIDE_N 1
#else
    #define B_STRIDE_K 1
    #define B_STRIDE_N K
#endif

#define WITH_BIAS 1
#define HAS_INLINE 1


#if K_DIVIDABLE_BY_8 && M_DIVIDABLE_BY_64  && N_DIVIDABLE_BY_64 && !TRANSPOSE_A && TRANSPOSE_B && OPTIMIZE
    const device float4 *load_target4 = (index & 32) 
        ? (const device float4 *)(weight_buffer + meta_buffer[1]) 
        : (const device float4 *)(data_buffer + meta_buffer[0]);
#else
    const device float *load_target = (index & 32) 
        ? (weight_buffer + meta_buffer[1]) 
        : (data_buffer + meta_buffer[0]);
#endif

    const int M = meta_buffer[4];
    const int N = meta_buffer[5];

    const int K = meta_buffer[6];

    threadgroup float4 shared4[32 * 8 * 2];

    const int stride_k = (index & 32) ? B_STRIDE_K : A_STRIDE_K;
    const int stride_mn = (index & 32) ? B_STRIDE_N : A_STRIDE_M;

    const int m_offset = index & 7;
    const int n_offset = index >> 3;

    const int mn_load_offset = ((index & 32) ? group_position.y : group_position.x) * 64 + (index & 15) * 4;
    const int k_load_offset = ((index & 16) ? 1 : 0);

    int track0 = mn_load_offset * stride_mn + (k_load_offset + 0) * stride_k;
    int track2 = track0 + 2 * stride_k;
    int track4 = track0 + 4 * stride_k;
    int track6 = track0 + 6 * stride_k;

#if !OPTIMIZE || !M_DIVIDABLE_BY_64 || !N_DIVIDABLE_BY_64
    const int max_MN = (index & 32) ? N : M;
#endif

    int shared_offset4 = ((index & 32) ? 16 : 0) + k_load_offset * 32 + (index & 15);
    int read_A_offset4 = m_offset * 2;
    int read_B_offset4 = n_offset * 2 + 16;

    float4 result[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int k = 0;

    while (k < K)
    {
        {
#if OPTIMIZE && K_DIVIDABLE_BY_8
    #if OPTIMIZE && M_DIVIDABLE_BY_64 && N_DIVIDABLE_BY_64
        #if OPTIMIZE && !TRANSPOSE_A && TRANSPOSE_B
            shared4[shared_offset4 + 32 * 0] = load_target4[track0 >> 2];
            shared4[shared_offset4 + 32 * 2] = load_target4[track2 >> 2];
            shared4[shared_offset4 + 32 * 4] = load_target4[track4 >> 2];
            shared4[shared_offset4 + 32 * 6] = load_target4[track6 >> 2];
        #else
            shared4[shared_offset4 + 32 * 0] = float4(
                load_target[track0 + stride_mn * 0],
                load_target[track0 + stride_mn * 1],
                load_target[track0 + stride_mn * 2],
                load_target[track0 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 2] = float4(
                load_target[track2 + stride_mn * 0],
                load_target[track2 + stride_mn * 1],
                load_target[track2 + stride_mn * 2],
                load_target[track2 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 4] = float4(
                load_target[track4 + stride_mn * 0],
                load_target[track4 + stride_mn * 1],
                load_target[track4 + stride_mn * 2],
                load_target[track4 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 6] = float4(
                load_target[track6 + stride_mn * 0],
                load_target[track6 + stride_mn * 1],
                load_target[track6 + stride_mn * 2],
                load_target[track6 + stride_mn * 3]
            ); 
        #endif
    #else
            shared4[shared_offset4 + 32 * 0] = float4(
                (mn_load_offset + 0 >= max_MN) ? 0 : load_target[track0 + stride_mn * 0],
                (mn_load_offset + 1 >= max_MN) ? 0 : load_target[track0 + stride_mn * 1],
                (mn_load_offset + 2 >= max_MN) ? 0 : load_target[track0 + stride_mn * 2],
                (mn_load_offset + 3 >= max_MN) ? 0 : load_target[track0 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 2] = float4(
                (mn_load_offset + 0 >= max_MN) ? 0 : load_target[track2 + stride_mn * 0],
                (mn_load_offset + 1 >= max_MN) ? 0 : load_target[track2 + stride_mn * 1],
                (mn_load_offset + 2 >= max_MN) ? 0 : load_target[track2 + stride_mn * 2],
                (mn_load_offset + 3 >= max_MN) ? 0 : load_target[track2 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 4] = float4(
                (mn_load_offset + 0 >= max_MN) ? 0 : load_target[track4 + stride_mn * 0],
                (mn_load_offset + 1 >= max_MN) ? 0 : load_target[track4 + stride_mn * 1],
                (mn_load_offset + 2 >= max_MN) ? 0 : load_target[track4 + stride_mn * 2],
                (mn_load_offset + 3 >= max_MN) ? 0 : load_target[track4 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 6] = float4(
                (mn_load_offset + 0 >= max_MN) ? 0 : load_target[track6 + stride_mn * 0],
                (mn_load_offset + 1 >= max_MN) ? 0 : load_target[track6 + stride_mn * 1],
                (mn_load_offset + 2 >= max_MN) ? 0 : load_target[track6 + stride_mn * 2],
                (mn_load_offset + 3 >= max_MN) ? 0 : load_target[track6 + stride_mn * 3]
            ); 
    #endif

            k += 8;
#else
    #if OPTIMIZE && M_DIVIDABLE_BY_64 && N_DIVIDABLE_BY_64
            shared4[shared_offset4 + 32 * 0] = float4(
                (k + k_load_offset >= K) ? 0 : load_target[track0 + stride_mn * 0],
                (k + k_load_offset >= K) ? 0 : load_target[track0 + stride_mn * 1],
                (k + k_load_offset >= K) ? 0 : load_target[track0 + stride_mn * 2],
                (k + k_load_offset >= K) ? 0 : load_target[track0 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 2] = float4(
                (k + k_load_offset >= K) ? 0 : load_target[track2 + stride_mn * 0],
                (k + k_load_offset >= K) ? 0 : load_target[track2 + stride_mn * 1],
                (k + k_load_offset >= K) ? 0 : load_target[track2 + stride_mn * 2],
                (k + k_load_offset >= K) ? 0 : load_target[track2 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 4] = float4(
                (k + k_load_offset >= K) ? 0 : load_target[track4 + stride_mn * 0],
                (k + k_load_offset >= K) ? 0 : load_target[track4 + stride_mn * 1],
                (k + k_load_offset >= K) ? 0 : load_target[track4 + stride_mn * 2],
                (k + k_load_offset >= K) ? 0 : load_target[track4 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 6] = float4(
                (k + k_load_offset >= K) ? 0 : load_target[track6 + stride_mn * 0],
                (k + k_load_offset >= K) ? 0 : load_target[track6 + stride_mn * 1],
                (k + k_load_offset >= K) ? 0 : load_target[track6 + stride_mn * 2],
                (k + k_load_offset >= K) ? 0 : load_target[track6 + stride_mn * 3]
            ); 
            k += 2;
    #else
            shared4[shared_offset4 + 32 * 0] = float4(
                (k + k_load_offset >= K || mn_load_offset + 0 >= max_MN) ? 0 : load_target[track0 + stride_mn * 0],
                (k + k_load_offset >= K || mn_load_offset + 1 >= max_MN) ? 0 : load_target[track0 + stride_mn * 1],
                (k + k_load_offset >= K || mn_load_offset + 2 >= max_MN) ? 0 : load_target[track0 + stride_mn * 2],
                (k + k_load_offset >= K || mn_load_offset + 3 >= max_MN) ? 0 : load_target[track0 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 2] = float4(
                (k + k_load_offset >= K || mn_load_offset + 0 >= max_MN) ? 0 : load_target[track2 + stride_mn * 0],
                (k + k_load_offset >= K || mn_load_offset + 1 >= max_MN) ? 0 : load_target[track2 + stride_mn * 1],
                (k + k_load_offset >= K || mn_load_offset + 2 >= max_MN) ? 0 : load_target[track2 + stride_mn * 2],
                (k + k_load_offset >= K || mn_load_offset + 3 >= max_MN) ? 0 : load_target[track2 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 4] = float4(
                (k + k_load_offset >= K || mn_load_offset + 0 >= max_MN) ? 0 : load_target[track4 + stride_mn * 0],
                (k + k_load_offset >= K || mn_load_offset + 1 >= max_MN) ? 0 : load_target[track4 + stride_mn * 1],
                (k + k_load_offset >= K || mn_load_offset + 2 >= max_MN) ? 0 : load_target[track4 + stride_mn * 2],
                (k + k_load_offset >= K || mn_load_offset + 3 >= max_MN) ? 0 : load_target[track4 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 6] = float4(
                (k + k_load_offset >= K || mn_load_offset + 0 >= max_MN) ? 0 : load_target[track6 + stride_mn * 0],
                (k + k_load_offset >= K || mn_load_offset + 1 >= max_MN) ? 0 : load_target[track6 + stride_mn * 1],
                (k + k_load_offset >= K || mn_load_offset + 2 >= max_MN) ? 0 : load_target[track6 + stride_mn * 2],
                (k + k_load_offset >= K || mn_load_offset + 3 >= max_MN) ? 0 : load_target[track6 + stride_mn * 3]
            ); 
            k += 2;
    #endif
#endif
        }

        threadgroup_barrier(mem_flags::mem_threadgroup);

        {
            for (int k_sub = 0; k_sub < 8; k_sub++)
            {
                float4 a00 = shared4[32 * k_sub + read_A_offset4 + 0];
                float4 a01 = shared4[32 * k_sub + read_A_offset4 + 1];
                float4 b00 = shared4[32 * k_sub + read_B_offset4 + 0];
                float4 b01 = shared4[32 * k_sub + read_B_offset4 + 1];

                result[4][0]  += b00[0] * a00[2];
                result[4][1]  += b00[1] * a00[2];
                result[0][1]  += b00[1] * a00[0];
                result[0][0]  += b00[0] * a00[0];
                result[6][0]  += b00[0] * a00[3];
                result[6][1]  += b00[1] * a00[3];
                result[2][1]  += b00[1] * a00[1];
                result[2][0]  += b00[0] * a00[1];
                result[12][0] += b00[0] * a01[2];
                result[12][1] += b00[1] * a01[2];
                result[8][1]  += b00[1] * a01[0];
                result[8][0]  += b00[0] * a01[0];
                result[14][0] += b00[0] * a01[3];
                result[14][1] += b00[1] * a01[3];
                result[10][1] += b00[1] * a01[1];
                result[10][0] += b00[0] * a01[1];

                result[14][2] += b00[2] * a01[3];
                result[14][3] += b00[3] * a01[3];
                result[10][3] += b00[3] * a01[1];
                result[10][2] += b00[2] * a01[1];
                result[12][2] += b00[2] * a01[2];
                result[12][3] += b00[3] * a01[2];
                result[8][3]  += b00[3] * a01[0];
                result[8][2]  += b00[2] * a01[0];
                result[6][2]  += b00[2] * a00[3];
                result[6][3]  += b00[3] * a00[3];
                result[2][3]  += b00[3] * a00[1];
                result[2][2]  += b00[2] * a00[1];
                result[4][2]  += b00[2] * a00[2];
                result[4][3]  += b00[3] * a00[2];
                result[0][3]  += b00[3] * a00[0];
                result[0][2]  += b00[2] * a00[0];

                result[5][0]  += b01[0] * a00[2];
                result[5][1]  += b01[1] * a00[2];
                result[1][1]  += b01[1] * a00[0];
                result[1][0]  += b01[0] * a00[0];
                result[7][0]  += b01[0] * a00[3];
                result[7][1]  += b01[1] * a00[3];
                result[3][1]  += b01[1] * a00[1];
                result[3][0]  += b01[0] * a00[1];
                result[13][0] += b01[0] * a01[2];
                result[13][1] += b01[1] * a01[2];
                result[9][1]  += b01[1] * a01[0];
                result[9][0]  += b01[0] * a01[0];
                result[15][0] += b01[0] * a01[3];
                result[15][1] += b01[1] * a01[3];
                result[11][1] += b01[1] * a01[1];
                result[11][0] += b01[0] * a01[1];

                result[15][2] += b01[2] * a01[3];
                result[15][3] += b01[3] * a01[3];
                result[11][3] += b01[3] * a01[1];
                result[11][2] += b01[2] * a01[1];
                result[13][2] += b01[2] * a01[2];
                result[13][3] += b01[3] * a01[2];
                result[9][3]  += b01[3] * a01[0];
                result[9][2]  += b01[2] * a01[0];
                result[7][2]  += b01[2] * a00[3];
                result[7][3]  += b01[3] * a00[3];
                result[3][3]  += b01[3] * a00[1];
                result[3][2]  += b01[2] * a00[1];
                result[5][2]  += b01[2] * a00[2];
                result[5][3]  += b01[3] * a00[2];
                result[1][3]  += b01[3] * a00[0];
                result[1][2]  += b01[2] * a00[0];
            }
        }

        shared_offset4 ^= 32 * 8;
        read_A_offset4 ^= 32 * 8;
        read_B_offset4 ^= 32 * 8;
        track0 += stride_k * 8;
        track2 += stride_k * 8;
        track4 += stride_k * 8;
        track6 += stride_k * 8;
    }

    {
    
#if OPTIMIZE && N_DIVIDABLE_BY_64
    #if WITH_BIAS
        float4 b[2];
        const device float4 *bias4 = (const device float4 *)(weight_buffer + meta_buffer[3]);
        b[0] = bias4[group_position.y * 16 + n_offset * 2 + 0];
        b[1] = bias4[group_position.y * 16 + n_offset * 2 + 1];
    #endif
    
        device float4 *C4 = (device float4 *)(data_buffer + meta_buffer[2]);
        const int N4 = N >> 2;
        int m = group_position.x * 64 + m_offset * 8;
        for (int m_sub = 0; m_sub < 8; m_sub++)
        {

    #if !M_DIVIDABLE_BY_64
            if (m >= M) continue;
    #endif

            const int n = group_position.y * 16 + n_offset * 2;
            float4 result0 = result[m_sub * 2 + 0];
            float4 result1 = result[m_sub * 2 + 1];

    #if WITH_BIAS
            result0 += b[0];
            result1 += b[1];
    #endif

    #if HAS_INLINE
            result0[0] = (result0[0]>0?result0[0]:(exp(result0[0])-1));
            result0[1] = (result0[1]>0?result0[1]:(exp(result0[1])-1));
            result0[2] = (result0[2]>0?result0[2]:(exp(result0[2])-1));
            result0[3] = (result0[3]>0?result0[3]:(exp(result0[3])-1));
            result1[0] = (result1[0]>0?result1[0]:(exp(result1[0])-1));
            result1[1] = (result1[1]>0?result1[1]:(exp(result1[1])-1));
            result1[2] = (result1[2]>0?result1[2]:(exp(result1[2])-1));
            result1[3] = (result1[3]>0?result1[3]:(exp(result1[3])-1));
    #endif

            C4[m * N4 + n + 0] = result0;
            C4[m * N4 + n + 1] = result1;
            
            m++;
        }
#else
    #if WITH_BIAS
        const device float *bias = weight_buffer + meta_buffer[3];
        float b[8];
        for (int n_sub = 0; n_sub < 8; n_sub++)
        {
            b[n_sub] = (group_position.y * 64 + n_offset * 8 + n_sub < N)
                ? bias[group_position.y * 64 + n_offset * 8 + n_sub]
                : 0;
        }
    #endif

        device float *C = data_buffer + meta_buffer[2];
        int m = group_position.x * 64 + m_offset * 8;
        for (int m_sub = 0; m_sub < 8; m_sub++)
        {
            int n = group_position.y * 64 + n_offset * 8;

            for (int n_sub1 = 0; n_sub1 < 2; n_sub1++)
            {
                for (int n_sub2 = 0; n_sub2 < 4; n_sub2++)
                {

    #if WITH_BIAS
        #if OPTIMIZE && M_DIVIDABLE_BY_64
                    (         n < N) ? (C[m * N + n] = (result[m_sub * 2 + n_sub1][n_sub2] + b[n_sub1*4+n_sub2]>0?result[m_sub * 2 + n_sub1][n_sub2] + b[n_sub1*4+n_sub2]:(exp(result[m_sub * 2 + n_sub1][n_sub2] + b[n_sub1*4+n_sub2])-1))) : 0;
        #else
                    (m < M && n < N) ? (C[m * N + n] = (result[m_sub * 2 + n_sub1][n_sub2] + b[n_sub1*4+n_sub2]>0?result[m_sub * 2 + n_sub1][n_sub2] + b[n_sub1*4+n_sub2]:(exp(result[m_sub * 2 + n_sub1][n_sub2] + b[n_sub1*4+n_sub2])-1))) : 0;
        #endif
    #else
        #if OPTIMIZE && M_DIVIDABLE_BY_64
                    (         n < N) ? (C[m * N + n] = (result[m_sub * 2 + n_sub1][n_sub2]>0?result[m_sub * 2 + n_sub1][n_sub2]:(exp(result[m_sub * 2 + n_sub1][n_sub2])-1))) : 0;
        #else
                    (m < M && n < N) ? (C[m * N + n] = (result[m_sub * 2 + n_sub1][n_sub2]>0?result[m_sub * 2 + n_sub1][n_sub2]:(exp(result[m_sub * 2 + n_sub1][n_sub2])-1))) : 0;
        #endif
    #endif
                    n++;
                }
            }
            
            m++;
        }
#endif

    }


#undef M_DIVIDABLE_BY_64
#undef N_DIVIDABLE_BY_64
#undef K_DIVIDABLE_BY_8
#undef TRANSPOSE_A
#undef TRANSPOSE_B
#undef A_STRIDE_K
#undef B_STRIDE_K
#undef A_STRIDE_M
#undef A_STRIDE_M
#undef WITH_BIAS
#undef HAS_INLINE
}


kernel void im2col_884fd9f360ddda90bf2cf46891d08c8646216406a876e59098baf6e0(const device float *param_buffer[[buffer(0)]],
                          device float *data_buffer[[buffer(1)]],
                          const device int * meta_buffer [[buffer(2)]],
                          ushort index_thread[[thread_position_in_threadgroup]],
                          ushort index_group[[threadgroup_position_in_grid]])
{
#define SH_EQUAL_1 1
#define SW_EQUAL_1 1
#define C1_DIVIDABLE_BY_4 1


#if OPTIMIZE && C1_DIVIDABLE_BY_4
    const device float4 *im4 = (const device float4 *)(data_buffer + meta_buffer[0]);
    device float4 *col4 = (device float4 *)(data_buffer + meta_buffer[1]);
    const int C1_4 = (meta_buffer[3]) >> 2;
#else
    const device float *im = data_buffer + meta_buffer[0];
    device float *col = data_buffer + meta_buffer[1];
    const int C1 = meta_buffer[3];
#endif

    // const int N = meta_buffer[2];
    const int H1 = meta_buffer[4];
    const int W1 = meta_buffer[5];
    const int H2 = meta_buffer[6];
    const int W2 = meta_buffer[7];
    const int KH = meta_buffer[8];
    const int KW = meta_buffer[9];
    const int PH = meta_buffer[12];
    const int PW = meta_buffer[13];

#if !OPTIMIZE || !SH_EQUAL_1
    const int SH = meta_buffer[10];
#endif

#if !OPTIMIZE || !SW_EQUAL_1
    const int SW = meta_buffer[11];
#endif

    const int H1P = H1 + 2 * PH;
    const int W1P = W1 + 2 * PW;

    const int w1 = (index_group % W1P) - PW;
    const int h1 = (index_group / W1P % H1P) - PH;
    const int  n = index_group / W1P / H1P;

#if OPTIMIZE && C1_DIVIDABLE_BY_4
    for (int c1_4 = index_thread; c1_4 < C1_4; c1_4 += 64) {
        const float4 v4 = (h1 < 0 || h1 >= H1 || w1 < 0 || w1 >= W1) ? 0 : im4[((n * H1 + h1) * W1 + w1) * C1_4 + c1_4];
#else
    for (int c1 = index_thread; c1 < C1; c1 += 64) {
        const float v = (h1 < 0 || h1 >= H1 || w1 < 0 || w1 >= W1) ? 0 : im[((n * H1 + h1) * W1 + w1) * C1 + c1];
#endif

#if OPTIMIZE && SH_EQUAL_1
        for (int kh = 0; kh < KH; kh++) {
            const int h2 = h1 + PH - kh;
#else
        for (int kh = (h1 + PH) % SH; kh < KH; kh += SH) {
            const int h2 = (h1 + PH - kh) / SH;
#endif
            if (h2 < 0 || h2 >= H2) continue;

#if OPTIMIZE && SH_EQUAL_1
            for (int kw = 0; kw < KW; kw++) {
                const int w2 = w1 + PW - kw;
#else
            for (int kw = (w1 + PW) % SW; kw < KW; kw += SW) {
                const int w2 = (w1 + PW - kw) / SW;
#endif
                if (w2 < 0 || w2 >= W2) continue;

#if OPTIMIZE && C1_DIVIDABLE_BY_4
                col4[((((n * H2 + h2) * W2 + w2) * KH + kh) * KW + kw) * C1_4 + c1_4] = v4;
#else
                col[((((n * H2 + h2) * W2 + w2) * KH + kh) * KW + kw) * C1 + c1] = v;
#endif
            }
        }
    }


#undef SH_EQUAL_1
#undef SW_EQUAL_1
#undef C1_DIVIDABLE_BY_4
}


kernel void sgemm_d2120e4054f06f47109f63cd135af4cf5af06fd36eb3fc1661c88712(const device float *weight_buffer[[buffer(0)]],
                          device float *data_buffer[[buffer(1)]],
                          const device int * meta_buffer [[buffer(2)]],
                          ushort index[[thread_index_in_threadgroup]],
                          ushort2 group_position[[threadgroup_position_in_grid]])
{
#define TRANSPOSE_A 1
#define TRANSPOSE_B 1
#define M_DIVIDABLE_BY_64 1
#define N_DIVIDABLE_BY_64 1
#define K_DIVIDABLE_BY_8 1

#if TRANSPOSE_A
    #define A_STRIDE_K 1
    #define A_STRIDE_M K
#else
    #define A_STRIDE_K M
    #define A_STRIDE_M 1
#endif

#if TRANSPOSE_B
    #define B_STRIDE_K N
    #define B_STRIDE_N 1
#else
    #define B_STRIDE_K 1
    #define B_STRIDE_N K
#endif

#define WITH_BIAS 1
#define HAS_INLINE 1


#if K_DIVIDABLE_BY_8 && M_DIVIDABLE_BY_64  && N_DIVIDABLE_BY_64 && !TRANSPOSE_A && TRANSPOSE_B && OPTIMIZE
    const device float4 *load_target4 = (index & 32) 
        ? (const device float4 *)(weight_buffer + meta_buffer[1]) 
        : (const device float4 *)(data_buffer + meta_buffer[0]);
#else
    const device float *load_target = (index & 32) 
        ? (weight_buffer + meta_buffer[1]) 
        : (data_buffer + meta_buffer[0]);
#endif

    const int M = meta_buffer[4];
    const int N = meta_buffer[5];

    const int K = meta_buffer[6];

    threadgroup float4 shared4[32 * 8 * 2];

    const int stride_k = (index & 32) ? B_STRIDE_K : A_STRIDE_K;
    const int stride_mn = (index & 32) ? B_STRIDE_N : A_STRIDE_M;

    const int m_offset = index & 7;
    const int n_offset = index >> 3;

    const int mn_load_offset = ((index & 32) ? group_position.y : group_position.x) * 64 + (index & 15) * 4;
    const int k_load_offset = ((index & 16) ? 1 : 0);

    int track0 = mn_load_offset * stride_mn + (k_load_offset + 0) * stride_k;
    int track2 = track0 + 2 * stride_k;
    int track4 = track0 + 4 * stride_k;
    int track6 = track0 + 6 * stride_k;

#if !OPTIMIZE || !M_DIVIDABLE_BY_64 || !N_DIVIDABLE_BY_64
    const int max_MN = (index & 32) ? N : M;
#endif

    int shared_offset4 = ((index & 32) ? 16 : 0) + k_load_offset * 32 + (index & 15);
    int read_A_offset4 = m_offset * 2;
    int read_B_offset4 = n_offset * 2 + 16;

    float4 result[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int k = 0;

    while (k < K)
    {
        {
#if OPTIMIZE && K_DIVIDABLE_BY_8
    #if OPTIMIZE && M_DIVIDABLE_BY_64 && N_DIVIDABLE_BY_64
        #if OPTIMIZE && !TRANSPOSE_A && TRANSPOSE_B
            shared4[shared_offset4 + 32 * 0] = load_target4[track0 >> 2];
            shared4[shared_offset4 + 32 * 2] = load_target4[track2 >> 2];
            shared4[shared_offset4 + 32 * 4] = load_target4[track4 >> 2];
            shared4[shared_offset4 + 32 * 6] = load_target4[track6 >> 2];
        #else
            shared4[shared_offset4 + 32 * 0] = float4(
                load_target[track0 + stride_mn * 0],
                load_target[track0 + stride_mn * 1],
                load_target[track0 + stride_mn * 2],
                load_target[track0 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 2] = float4(
                load_target[track2 + stride_mn * 0],
                load_target[track2 + stride_mn * 1],
                load_target[track2 + stride_mn * 2],
                load_target[track2 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 4] = float4(
                load_target[track4 + stride_mn * 0],
                load_target[track4 + stride_mn * 1],
                load_target[track4 + stride_mn * 2],
                load_target[track4 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 6] = float4(
                load_target[track6 + stride_mn * 0],
                load_target[track6 + stride_mn * 1],
                load_target[track6 + stride_mn * 2],
                load_target[track6 + stride_mn * 3]
            ); 
        #endif
    #else
            shared4[shared_offset4 + 32 * 0] = float4(
                (mn_load_offset + 0 >= max_MN) ? 0 : load_target[track0 + stride_mn * 0],
                (mn_load_offset + 1 >= max_MN) ? 0 : load_target[track0 + stride_mn * 1],
                (mn_load_offset + 2 >= max_MN) ? 0 : load_target[track0 + stride_mn * 2],
                (mn_load_offset + 3 >= max_MN) ? 0 : load_target[track0 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 2] = float4(
                (mn_load_offset + 0 >= max_MN) ? 0 : load_target[track2 + stride_mn * 0],
                (mn_load_offset + 1 >= max_MN) ? 0 : load_target[track2 + stride_mn * 1],
                (mn_load_offset + 2 >= max_MN) ? 0 : load_target[track2 + stride_mn * 2],
                (mn_load_offset + 3 >= max_MN) ? 0 : load_target[track2 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 4] = float4(
                (mn_load_offset + 0 >= max_MN) ? 0 : load_target[track4 + stride_mn * 0],
                (mn_load_offset + 1 >= max_MN) ? 0 : load_target[track4 + stride_mn * 1],
                (mn_load_offset + 2 >= max_MN) ? 0 : load_target[track4 + stride_mn * 2],
                (mn_load_offset + 3 >= max_MN) ? 0 : load_target[track4 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 6] = float4(
                (mn_load_offset + 0 >= max_MN) ? 0 : load_target[track6 + stride_mn * 0],
                (mn_load_offset + 1 >= max_MN) ? 0 : load_target[track6 + stride_mn * 1],
                (mn_load_offset + 2 >= max_MN) ? 0 : load_target[track6 + stride_mn * 2],
                (mn_load_offset + 3 >= max_MN) ? 0 : load_target[track6 + stride_mn * 3]
            ); 
    #endif

            k += 8;
#else
    #if OPTIMIZE && M_DIVIDABLE_BY_64 && N_DIVIDABLE_BY_64
            shared4[shared_offset4 + 32 * 0] = float4(
                (k + k_load_offset >= K) ? 0 : load_target[track0 + stride_mn * 0],
                (k + k_load_offset >= K) ? 0 : load_target[track0 + stride_mn * 1],
                (k + k_load_offset >= K) ? 0 : load_target[track0 + stride_mn * 2],
                (k + k_load_offset >= K) ? 0 : load_target[track0 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 2] = float4(
                (k + k_load_offset >= K) ? 0 : load_target[track2 + stride_mn * 0],
                (k + k_load_offset >= K) ? 0 : load_target[track2 + stride_mn * 1],
                (k + k_load_offset >= K) ? 0 : load_target[track2 + stride_mn * 2],
                (k + k_load_offset >= K) ? 0 : load_target[track2 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 4] = float4(
                (k + k_load_offset >= K) ? 0 : load_target[track4 + stride_mn * 0],
                (k + k_load_offset >= K) ? 0 : load_target[track4 + stride_mn * 1],
                (k + k_load_offset >= K) ? 0 : load_target[track4 + stride_mn * 2],
                (k + k_load_offset >= K) ? 0 : load_target[track4 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 6] = float4(
                (k + k_load_offset >= K) ? 0 : load_target[track6 + stride_mn * 0],
                (k + k_load_offset >= K) ? 0 : load_target[track6 + stride_mn * 1],
                (k + k_load_offset >= K) ? 0 : load_target[track6 + stride_mn * 2],
                (k + k_load_offset >= K) ? 0 : load_target[track6 + stride_mn * 3]
            ); 
            k += 2;
    #else
            shared4[shared_offset4 + 32 * 0] = float4(
                (k + k_load_offset >= K || mn_load_offset + 0 >= max_MN) ? 0 : load_target[track0 + stride_mn * 0],
                (k + k_load_offset >= K || mn_load_offset + 1 >= max_MN) ? 0 : load_target[track0 + stride_mn * 1],
                (k + k_load_offset >= K || mn_load_offset + 2 >= max_MN) ? 0 : load_target[track0 + stride_mn * 2],
                (k + k_load_offset >= K || mn_load_offset + 3 >= max_MN) ? 0 : load_target[track0 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 2] = float4(
                (k + k_load_offset >= K || mn_load_offset + 0 >= max_MN) ? 0 : load_target[track2 + stride_mn * 0],
                (k + k_load_offset >= K || mn_load_offset + 1 >= max_MN) ? 0 : load_target[track2 + stride_mn * 1],
                (k + k_load_offset >= K || mn_load_offset + 2 >= max_MN) ? 0 : load_target[track2 + stride_mn * 2],
                (k + k_load_offset >= K || mn_load_offset + 3 >= max_MN) ? 0 : load_target[track2 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 4] = float4(
                (k + k_load_offset >= K || mn_load_offset + 0 >= max_MN) ? 0 : load_target[track4 + stride_mn * 0],
                (k + k_load_offset >= K || mn_load_offset + 1 >= max_MN) ? 0 : load_target[track4 + stride_mn * 1],
                (k + k_load_offset >= K || mn_load_offset + 2 >= max_MN) ? 0 : load_target[track4 + stride_mn * 2],
                (k + k_load_offset >= K || mn_load_offset + 3 >= max_MN) ? 0 : load_target[track4 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 6] = float4(
                (k + k_load_offset >= K || mn_load_offset + 0 >= max_MN) ? 0 : load_target[track6 + stride_mn * 0],
                (k + k_load_offset >= K || mn_load_offset + 1 >= max_MN) ? 0 : load_target[track6 + stride_mn * 1],
                (k + k_load_offset >= K || mn_load_offset + 2 >= max_MN) ? 0 : load_target[track6 + stride_mn * 2],
                (k + k_load_offset >= K || mn_load_offset + 3 >= max_MN) ? 0 : load_target[track6 + stride_mn * 3]
            ); 
            k += 2;
    #endif
#endif
        }

        threadgroup_barrier(mem_flags::mem_threadgroup);

        {
            for (int k_sub = 0; k_sub < 8; k_sub++)
            {
                float4 a00 = shared4[32 * k_sub + read_A_offset4 + 0];
                float4 a01 = shared4[32 * k_sub + read_A_offset4 + 1];
                float4 b00 = shared4[32 * k_sub + read_B_offset4 + 0];
                float4 b01 = shared4[32 * k_sub + read_B_offset4 + 1];

                result[4][0]  += b00[0] * a00[2];
                result[4][1]  += b00[1] * a00[2];
                result[0][1]  += b00[1] * a00[0];
                result[0][0]  += b00[0] * a00[0];
                result[6][0]  += b00[0] * a00[3];
                result[6][1]  += b00[1] * a00[3];
                result[2][1]  += b00[1] * a00[1];
                result[2][0]  += b00[0] * a00[1];
                result[12][0] += b00[0] * a01[2];
                result[12][1] += b00[1] * a01[2];
                result[8][1]  += b00[1] * a01[0];
                result[8][0]  += b00[0] * a01[0];
                result[14][0] += b00[0] * a01[3];
                result[14][1] += b00[1] * a01[3];
                result[10][1] += b00[1] * a01[1];
                result[10][0] += b00[0] * a01[1];

                result[14][2] += b00[2] * a01[3];
                result[14][3] += b00[3] * a01[3];
                result[10][3] += b00[3] * a01[1];
                result[10][2] += b00[2] * a01[1];
                result[12][2] += b00[2] * a01[2];
                result[12][3] += b00[3] * a01[2];
                result[8][3]  += b00[3] * a01[0];
                result[8][2]  += b00[2] * a01[0];
                result[6][2]  += b00[2] * a00[3];
                result[6][3]  += b00[3] * a00[3];
                result[2][3]  += b00[3] * a00[1];
                result[2][2]  += b00[2] * a00[1];
                result[4][2]  += b00[2] * a00[2];
                result[4][3]  += b00[3] * a00[2];
                result[0][3]  += b00[3] * a00[0];
                result[0][2]  += b00[2] * a00[0];

                result[5][0]  += b01[0] * a00[2];
                result[5][1]  += b01[1] * a00[2];
                result[1][1]  += b01[1] * a00[0];
                result[1][0]  += b01[0] * a00[0];
                result[7][0]  += b01[0] * a00[3];
                result[7][1]  += b01[1] * a00[3];
                result[3][1]  += b01[1] * a00[1];
                result[3][0]  += b01[0] * a00[1];
                result[13][0] += b01[0] * a01[2];
                result[13][1] += b01[1] * a01[2];
                result[9][1]  += b01[1] * a01[0];
                result[9][0]  += b01[0] * a01[0];
                result[15][0] += b01[0] * a01[3];
                result[15][1] += b01[1] * a01[3];
                result[11][1] += b01[1] * a01[1];
                result[11][0] += b01[0] * a01[1];

                result[15][2] += b01[2] * a01[3];
                result[15][3] += b01[3] * a01[3];
                result[11][3] += b01[3] * a01[1];
                result[11][2] += b01[2] * a01[1];
                result[13][2] += b01[2] * a01[2];
                result[13][3] += b01[3] * a01[2];
                result[9][3]  += b01[3] * a01[0];
                result[9][2]  += b01[2] * a01[0];
                result[7][2]  += b01[2] * a00[3];
                result[7][3]  += b01[3] * a00[3];
                result[3][3]  += b01[3] * a00[1];
                result[3][2]  += b01[2] * a00[1];
                result[5][2]  += b01[2] * a00[2];
                result[5][3]  += b01[3] * a00[2];
                result[1][3]  += b01[3] * a00[0];
                result[1][2]  += b01[2] * a00[0];
            }
        }

        shared_offset4 ^= 32 * 8;
        read_A_offset4 ^= 32 * 8;
        read_B_offset4 ^= 32 * 8;
        track0 += stride_k * 8;
        track2 += stride_k * 8;
        track4 += stride_k * 8;
        track6 += stride_k * 8;
    }

    {
    
#if OPTIMIZE && N_DIVIDABLE_BY_64
    #if WITH_BIAS
        float4 b[2];
        const device float4 *bias4 = (const device float4 *)(weight_buffer + meta_buffer[3]);
        b[0] = bias4[group_position.y * 16 + n_offset * 2 + 0];
        b[1] = bias4[group_position.y * 16 + n_offset * 2 + 1];
    #endif
    
        device float4 *C4 = (device float4 *)(data_buffer + meta_buffer[2]);
        const int N4 = N >> 2;
        int m = group_position.x * 64 + m_offset * 8;
        for (int m_sub = 0; m_sub < 8; m_sub++)
        {

    #if !M_DIVIDABLE_BY_64
            if (m >= M) continue;
    #endif

            const int n = group_position.y * 16 + n_offset * 2;
            float4 result0 = result[m_sub * 2 + 0];
            float4 result1 = result[m_sub * 2 + 1];

    #if WITH_BIAS
            result0 += b[0];
            result1 += b[1];
    #endif

    #if HAS_INLINE
            result0[0] = (result0[0]>0?result0[0]:0);
            result0[1] = (result0[1]>0?result0[1]:0);
            result0[2] = (result0[2]>0?result0[2]:0);
            result0[3] = (result0[3]>0?result0[3]:0);
            result1[0] = (result1[0]>0?result1[0]:0);
            result1[1] = (result1[1]>0?result1[1]:0);
            result1[2] = (result1[2]>0?result1[2]:0);
            result1[3] = (result1[3]>0?result1[3]:0);
    #endif

            C4[m * N4 + n + 0] = result0;
            C4[m * N4 + n + 1] = result1;
            
            m++;
        }
#else
    #if WITH_BIAS
        const device float *bias = weight_buffer + meta_buffer[3];
        float b[8];
        for (int n_sub = 0; n_sub < 8; n_sub++)
        {
            b[n_sub] = (group_position.y * 64 + n_offset * 8 + n_sub < N)
                ? bias[group_position.y * 64 + n_offset * 8 + n_sub]
                : 0;
        }
    #endif

        device float *C = data_buffer + meta_buffer[2];
        int m = group_position.x * 64 + m_offset * 8;
        for (int m_sub = 0; m_sub < 8; m_sub++)
        {
            int n = group_position.y * 64 + n_offset * 8;

            for (int n_sub1 = 0; n_sub1 < 2; n_sub1++)
            {
                for (int n_sub2 = 0; n_sub2 < 4; n_sub2++)
                {

    #if WITH_BIAS
        #if OPTIMIZE && M_DIVIDABLE_BY_64
                    (         n < N) ? (C[m * N + n] = (result[m_sub * 2 + n_sub1][n_sub2] + b[n_sub1*4+n_sub2]>0?result[m_sub * 2 + n_sub1][n_sub2] + b[n_sub1*4+n_sub2]:0)) : 0;
        #else
                    (m < M && n < N) ? (C[m * N + n] = (result[m_sub * 2 + n_sub1][n_sub2] + b[n_sub1*4+n_sub2]>0?result[m_sub * 2 + n_sub1][n_sub2] + b[n_sub1*4+n_sub2]:0)) : 0;
        #endif
    #else
        #if OPTIMIZE && M_DIVIDABLE_BY_64
                    (         n < N) ? (C[m * N + n] = (result[m_sub * 2 + n_sub1][n_sub2]>0?result[m_sub * 2 + n_sub1][n_sub2]:0)) : 0;
        #else
                    (m < M && n < N) ? (C[m * N + n] = (result[m_sub * 2 + n_sub1][n_sub2]>0?result[m_sub * 2 + n_sub1][n_sub2]:0)) : 0;
        #endif
    #endif
                    n++;
                }
            }
            
            m++;
        }
#endif

    }


#undef M_DIVIDABLE_BY_64
#undef N_DIVIDABLE_BY_64
#undef K_DIVIDABLE_BY_8
#undef TRANSPOSE_A
#undef TRANSPOSE_B
#undef A_STRIDE_K
#undef B_STRIDE_K
#undef A_STRIDE_M
#undef A_STRIDE_M
#undef WITH_BIAS
#undef HAS_INLINE
}


kernel void sgemm_a9be2ab4b5e5fb537f0aa8f376574568f7fb0a432167d3f83543bf78(const device float *weight_buffer[[buffer(0)]],
                          device float *data_buffer[[buffer(1)]],
                          const device int * meta_buffer [[buffer(2)]],
                          ushort index[[thread_index_in_threadgroup]],
                          ushort2 group_position[[threadgroup_position_in_grid]])
{
#define TRANSPOSE_A 1
#define TRANSPOSE_B 1
#define M_DIVIDABLE_BY_64 1
#define N_DIVIDABLE_BY_64 1
#define K_DIVIDABLE_BY_8 1

#if TRANSPOSE_A
    #define A_STRIDE_K 1
    #define A_STRIDE_M K
#else
    #define A_STRIDE_K M
    #define A_STRIDE_M 1
#endif

#if TRANSPOSE_B
    #define B_STRIDE_K N
    #define B_STRIDE_N 1
#else
    #define B_STRIDE_K 1
    #define B_STRIDE_N K
#endif

#define WITH_BIAS 1
#define HAS_INLINE 1


#if K_DIVIDABLE_BY_8 && M_DIVIDABLE_BY_64  && N_DIVIDABLE_BY_64 && !TRANSPOSE_A && TRANSPOSE_B && OPTIMIZE
    const device float4 *load_target4 = (index & 32) 
        ? (const device float4 *)(weight_buffer + meta_buffer[1]) 
        : (const device float4 *)(data_buffer + meta_buffer[0]);
#else
    const device float *load_target = (index & 32) 
        ? (weight_buffer + meta_buffer[1]) 
        : (data_buffer + meta_buffer[0]);
#endif

    const int M = meta_buffer[4];
    const int N = meta_buffer[5];

    const int K = meta_buffer[6];

    threadgroup float4 shared4[32 * 8 * 2];

    const int stride_k = (index & 32) ? B_STRIDE_K : A_STRIDE_K;
    const int stride_mn = (index & 32) ? B_STRIDE_N : A_STRIDE_M;

    const int m_offset = index & 7;
    const int n_offset = index >> 3;

    const int mn_load_offset = ((index & 32) ? group_position.y : group_position.x) * 64 + (index & 15) * 4;
    const int k_load_offset = ((index & 16) ? 1 : 0);

    int track0 = mn_load_offset * stride_mn + (k_load_offset + 0) * stride_k;
    int track2 = track0 + 2 * stride_k;
    int track4 = track0 + 4 * stride_k;
    int track6 = track0 + 6 * stride_k;

#if !OPTIMIZE || !M_DIVIDABLE_BY_64 || !N_DIVIDABLE_BY_64
    const int max_MN = (index & 32) ? N : M;
#endif

    int shared_offset4 = ((index & 32) ? 16 : 0) + k_load_offset * 32 + (index & 15);
    int read_A_offset4 = m_offset * 2;
    int read_B_offset4 = n_offset * 2 + 16;

    float4 result[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int k = 0;

    while (k < K)
    {
        {
#if OPTIMIZE && K_DIVIDABLE_BY_8
    #if OPTIMIZE && M_DIVIDABLE_BY_64 && N_DIVIDABLE_BY_64
        #if OPTIMIZE && !TRANSPOSE_A && TRANSPOSE_B
            shared4[shared_offset4 + 32 * 0] = load_target4[track0 >> 2];
            shared4[shared_offset4 + 32 * 2] = load_target4[track2 >> 2];
            shared4[shared_offset4 + 32 * 4] = load_target4[track4 >> 2];
            shared4[shared_offset4 + 32 * 6] = load_target4[track6 >> 2];
        #else
            shared4[shared_offset4 + 32 * 0] = float4(
                load_target[track0 + stride_mn * 0],
                load_target[track0 + stride_mn * 1],
                load_target[track0 + stride_mn * 2],
                load_target[track0 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 2] = float4(
                load_target[track2 + stride_mn * 0],
                load_target[track2 + stride_mn * 1],
                load_target[track2 + stride_mn * 2],
                load_target[track2 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 4] = float4(
                load_target[track4 + stride_mn * 0],
                load_target[track4 + stride_mn * 1],
                load_target[track4 + stride_mn * 2],
                load_target[track4 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 6] = float4(
                load_target[track6 + stride_mn * 0],
                load_target[track6 + stride_mn * 1],
                load_target[track6 + stride_mn * 2],
                load_target[track6 + stride_mn * 3]
            ); 
        #endif
    #else
            shared4[shared_offset4 + 32 * 0] = float4(
                (mn_load_offset + 0 >= max_MN) ? 0 : load_target[track0 + stride_mn * 0],
                (mn_load_offset + 1 >= max_MN) ? 0 : load_target[track0 + stride_mn * 1],
                (mn_load_offset + 2 >= max_MN) ? 0 : load_target[track0 + stride_mn * 2],
                (mn_load_offset + 3 >= max_MN) ? 0 : load_target[track0 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 2] = float4(
                (mn_load_offset + 0 >= max_MN) ? 0 : load_target[track2 + stride_mn * 0],
                (mn_load_offset + 1 >= max_MN) ? 0 : load_target[track2 + stride_mn * 1],
                (mn_load_offset + 2 >= max_MN) ? 0 : load_target[track2 + stride_mn * 2],
                (mn_load_offset + 3 >= max_MN) ? 0 : load_target[track2 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 4] = float4(
                (mn_load_offset + 0 >= max_MN) ? 0 : load_target[track4 + stride_mn * 0],
                (mn_load_offset + 1 >= max_MN) ? 0 : load_target[track4 + stride_mn * 1],
                (mn_load_offset + 2 >= max_MN) ? 0 : load_target[track4 + stride_mn * 2],
                (mn_load_offset + 3 >= max_MN) ? 0 : load_target[track4 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 6] = float4(
                (mn_load_offset + 0 >= max_MN) ? 0 : load_target[track6 + stride_mn * 0],
                (mn_load_offset + 1 >= max_MN) ? 0 : load_target[track6 + stride_mn * 1],
                (mn_load_offset + 2 >= max_MN) ? 0 : load_target[track6 + stride_mn * 2],
                (mn_load_offset + 3 >= max_MN) ? 0 : load_target[track6 + stride_mn * 3]
            ); 
    #endif

            k += 8;
#else
    #if OPTIMIZE && M_DIVIDABLE_BY_64 && N_DIVIDABLE_BY_64
            shared4[shared_offset4 + 32 * 0] = float4(
                (k + k_load_offset >= K) ? 0 : load_target[track0 + stride_mn * 0],
                (k + k_load_offset >= K) ? 0 : load_target[track0 + stride_mn * 1],
                (k + k_load_offset >= K) ? 0 : load_target[track0 + stride_mn * 2],
                (k + k_load_offset >= K) ? 0 : load_target[track0 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 2] = float4(
                (k + k_load_offset >= K) ? 0 : load_target[track2 + stride_mn * 0],
                (k + k_load_offset >= K) ? 0 : load_target[track2 + stride_mn * 1],
                (k + k_load_offset >= K) ? 0 : load_target[track2 + stride_mn * 2],
                (k + k_load_offset >= K) ? 0 : load_target[track2 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 4] = float4(
                (k + k_load_offset >= K) ? 0 : load_target[track4 + stride_mn * 0],
                (k + k_load_offset >= K) ? 0 : load_target[track4 + stride_mn * 1],
                (k + k_load_offset >= K) ? 0 : load_target[track4 + stride_mn * 2],
                (k + k_load_offset >= K) ? 0 : load_target[track4 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 6] = float4(
                (k + k_load_offset >= K) ? 0 : load_target[track6 + stride_mn * 0],
                (k + k_load_offset >= K) ? 0 : load_target[track6 + stride_mn * 1],
                (k + k_load_offset >= K) ? 0 : load_target[track6 + stride_mn * 2],
                (k + k_load_offset >= K) ? 0 : load_target[track6 + stride_mn * 3]
            ); 
            k += 2;
    #else
            shared4[shared_offset4 + 32 * 0] = float4(
                (k + k_load_offset >= K || mn_load_offset + 0 >= max_MN) ? 0 : load_target[track0 + stride_mn * 0],
                (k + k_load_offset >= K || mn_load_offset + 1 >= max_MN) ? 0 : load_target[track0 + stride_mn * 1],
                (k + k_load_offset >= K || mn_load_offset + 2 >= max_MN) ? 0 : load_target[track0 + stride_mn * 2],
                (k + k_load_offset >= K || mn_load_offset + 3 >= max_MN) ? 0 : load_target[track0 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 2] = float4(
                (k + k_load_offset >= K || mn_load_offset + 0 >= max_MN) ? 0 : load_target[track2 + stride_mn * 0],
                (k + k_load_offset >= K || mn_load_offset + 1 >= max_MN) ? 0 : load_target[track2 + stride_mn * 1],
                (k + k_load_offset >= K || mn_load_offset + 2 >= max_MN) ? 0 : load_target[track2 + stride_mn * 2],
                (k + k_load_offset >= K || mn_load_offset + 3 >= max_MN) ? 0 : load_target[track2 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 4] = float4(
                (k + k_load_offset >= K || mn_load_offset + 0 >= max_MN) ? 0 : load_target[track4 + stride_mn * 0],
                (k + k_load_offset >= K || mn_load_offset + 1 >= max_MN) ? 0 : load_target[track4 + stride_mn * 1],
                (k + k_load_offset >= K || mn_load_offset + 2 >= max_MN) ? 0 : load_target[track4 + stride_mn * 2],
                (k + k_load_offset >= K || mn_load_offset + 3 >= max_MN) ? 0 : load_target[track4 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 6] = float4(
                (k + k_load_offset >= K || mn_load_offset + 0 >= max_MN) ? 0 : load_target[track6 + stride_mn * 0],
                (k + k_load_offset >= K || mn_load_offset + 1 >= max_MN) ? 0 : load_target[track6 + stride_mn * 1],
                (k + k_load_offset >= K || mn_load_offset + 2 >= max_MN) ? 0 : load_target[track6 + stride_mn * 2],
                (k + k_load_offset >= K || mn_load_offset + 3 >= max_MN) ? 0 : load_target[track6 + stride_mn * 3]
            ); 
            k += 2;
    #endif
#endif
        }

        threadgroup_barrier(mem_flags::mem_threadgroup);

        {
            for (int k_sub = 0; k_sub < 8; k_sub++)
            {
                float4 a00 = shared4[32 * k_sub + read_A_offset4 + 0];
                float4 a01 = shared4[32 * k_sub + read_A_offset4 + 1];
                float4 b00 = shared4[32 * k_sub + read_B_offset4 + 0];
                float4 b01 = shared4[32 * k_sub + read_B_offset4 + 1];

                result[4][0]  += b00[0] * a00[2];
                result[4][1]  += b00[1] * a00[2];
                result[0][1]  += b00[1] * a00[0];
                result[0][0]  += b00[0] * a00[0];
                result[6][0]  += b00[0] * a00[3];
                result[6][1]  += b00[1] * a00[3];
                result[2][1]  += b00[1] * a00[1];
                result[2][0]  += b00[0] * a00[1];
                result[12][0] += b00[0] * a01[2];
                result[12][1] += b00[1] * a01[2];
                result[8][1]  += b00[1] * a01[0];
                result[8][0]  += b00[0] * a01[0];
                result[14][0] += b00[0] * a01[3];
                result[14][1] += b00[1] * a01[3];
                result[10][1] += b00[1] * a01[1];
                result[10][0] += b00[0] * a01[1];

                result[14][2] += b00[2] * a01[3];
                result[14][3] += b00[3] * a01[3];
                result[10][3] += b00[3] * a01[1];
                result[10][2] += b00[2] * a01[1];
                result[12][2] += b00[2] * a01[2];
                result[12][3] += b00[3] * a01[2];
                result[8][3]  += b00[3] * a01[0];
                result[8][2]  += b00[2] * a01[0];
                result[6][2]  += b00[2] * a00[3];
                result[6][3]  += b00[3] * a00[3];
                result[2][3]  += b00[3] * a00[1];
                result[2][2]  += b00[2] * a00[1];
                result[4][2]  += b00[2] * a00[2];
                result[4][3]  += b00[3] * a00[2];
                result[0][3]  += b00[3] * a00[0];
                result[0][2]  += b00[2] * a00[0];

                result[5][0]  += b01[0] * a00[2];
                result[5][1]  += b01[1] * a00[2];
                result[1][1]  += b01[1] * a00[0];
                result[1][0]  += b01[0] * a00[0];
                result[7][0]  += b01[0] * a00[3];
                result[7][1]  += b01[1] * a00[3];
                result[3][1]  += b01[1] * a00[1];
                result[3][0]  += b01[0] * a00[1];
                result[13][0] += b01[0] * a01[2];
                result[13][1] += b01[1] * a01[2];
                result[9][1]  += b01[1] * a01[0];
                result[9][0]  += b01[0] * a01[0];
                result[15][0] += b01[0] * a01[3];
                result[15][1] += b01[1] * a01[3];
                result[11][1] += b01[1] * a01[1];
                result[11][0] += b01[0] * a01[1];

                result[15][2] += b01[2] * a01[3];
                result[15][3] += b01[3] * a01[3];
                result[11][3] += b01[3] * a01[1];
                result[11][2] += b01[2] * a01[1];
                result[13][2] += b01[2] * a01[2];
                result[13][3] += b01[3] * a01[2];
                result[9][3]  += b01[3] * a01[0];
                result[9][2]  += b01[2] * a01[0];
                result[7][2]  += b01[2] * a00[3];
                result[7][3]  += b01[3] * a00[3];
                result[3][3]  += b01[3] * a00[1];
                result[3][2]  += b01[2] * a00[1];
                result[5][2]  += b01[2] * a00[2];
                result[5][3]  += b01[3] * a00[2];
                result[1][3]  += b01[3] * a00[0];
                result[1][2]  += b01[2] * a00[0];
            }
        }

        shared_offset4 ^= 32 * 8;
        read_A_offset4 ^= 32 * 8;
        read_B_offset4 ^= 32 * 8;
        track0 += stride_k * 8;
        track2 += stride_k * 8;
        track4 += stride_k * 8;
        track6 += stride_k * 8;
    }

    {
    
#if OPTIMIZE && N_DIVIDABLE_BY_64
    #if WITH_BIAS
        float4 b[2];
        const device float4 *bias4 = (const device float4 *)(weight_buffer + meta_buffer[3]);
        b[0] = bias4[group_position.y * 16 + n_offset * 2 + 0];
        b[1] = bias4[group_position.y * 16 + n_offset * 2 + 1];
    #endif
    
        device float4 *C4 = (device float4 *)(data_buffer + meta_buffer[2]);
        const int N4 = N >> 2;
        int m = group_position.x * 64 + m_offset * 8;
        for (int m_sub = 0; m_sub < 8; m_sub++)
        {

    #if !M_DIVIDABLE_BY_64
            if (m >= M) continue;
    #endif

            const int n = group_position.y * 16 + n_offset * 2;
            float4 result0 = result[m_sub * 2 + 0];
            float4 result1 = result[m_sub * 2 + 1];

    #if WITH_BIAS
            result0 += b[0];
            result1 += b[1];
    #endif

    #if HAS_INLINE
            result0[0] = result0[0];
            result0[1] = result0[1];
            result0[2] = result0[2];
            result0[3] = result0[3];
            result1[0] = result1[0];
            result1[1] = result1[1];
            result1[2] = result1[2];
            result1[3] = result1[3];
    #endif

            C4[m * N4 + n + 0] = result0;
            C4[m * N4 + n + 1] = result1;
            
            m++;
        }
#else
    #if WITH_BIAS
        const device float *bias = weight_buffer + meta_buffer[3];
        float b[8];
        for (int n_sub = 0; n_sub < 8; n_sub++)
        {
            b[n_sub] = (group_position.y * 64 + n_offset * 8 + n_sub < N)
                ? bias[group_position.y * 64 + n_offset * 8 + n_sub]
                : 0;
        }
    #endif

        device float *C = data_buffer + meta_buffer[2];
        int m = group_position.x * 64 + m_offset * 8;
        for (int m_sub = 0; m_sub < 8; m_sub++)
        {
            int n = group_position.y * 64 + n_offset * 8;

            for (int n_sub1 = 0; n_sub1 < 2; n_sub1++)
            {
                for (int n_sub2 = 0; n_sub2 < 4; n_sub2++)
                {

    #if WITH_BIAS
        #if OPTIMIZE && M_DIVIDABLE_BY_64
                    (         n < N) ? (C[m * N + n] = result[m_sub * 2 + n_sub1][n_sub2] + b[n_sub1*4+n_sub2]) : 0;
        #else
                    (m < M && n < N) ? (C[m * N + n] = result[m_sub * 2 + n_sub1][n_sub2] + b[n_sub1*4+n_sub2]) : 0;
        #endif
    #else
        #if OPTIMIZE && M_DIVIDABLE_BY_64
                    (         n < N) ? (C[m * N + n] = result[m_sub * 2 + n_sub1][n_sub2]) : 0;
        #else
                    (m < M && n < N) ? (C[m * N + n] = result[m_sub * 2 + n_sub1][n_sub2]) : 0;
        #endif
    #endif
                    n++;
                }
            }
            
            m++;
        }
#endif

    }


#undef M_DIVIDABLE_BY_64
#undef N_DIVIDABLE_BY_64
#undef K_DIVIDABLE_BY_8
#undef TRANSPOSE_A
#undef TRANSPOSE_B
#undef A_STRIDE_K
#undef B_STRIDE_K
#undef A_STRIDE_M
#undef A_STRIDE_M
#undef WITH_BIAS
#undef HAS_INLINE
}


kernel void elementwisesum_ddc83841f1d13720fab1c95ae5207cf535216a48ae9d41fa157da605(const device float *weight_buffer[[buffer(0)]],
                          device float *data_buffer[[buffer(1)]],
                          const device int * meta_buffer [[buffer(2)]],
                          uint index[[thread_position_in_grid]],
                          uint num_threads[[threads_per_grid]])
{
#define N_DIVIDABLE_BY_4 1
#define HAS_INLINE 1

#if N_DIVIDABLE_BY_4
    #define T_VALUE float4
#else
    #define T_VALUE float
#endif


#if OPTIMIZE && N_DIVIDABLE_BY_4
    const device float4 *X0 = (const device float4 *)(data_buffer + meta_buffer[0]);
    const device float4 *X1 = (const device float4 *)(data_buffer + meta_buffer[1]);
    device float4 *Y = (device float4 *)(data_buffer + meta_buffer[2]);
    const int N = (meta_buffer[3]) >> 2;
#else
    const device float *X0 = data_buffer + meta_buffer[0];
    const device float *X1 = data_buffer + meta_buffer[1];
    device float *Y = data_buffer + meta_buffer[2];
    const int N = meta_buffer[3];
#endif
  
    for (int gid = index; gid < N; gid += num_threads) {

        T_VALUE result = X0[gid] + X1[gid];

#if OPTIMIZE && HAS_INLINE
    #if OPTIMIZE && N_DIVIDABLE_BY_4
            result[0] = result[0];
            result[1] = result[1];
            result[2] = result[2];
            result[3] = result[3];
    #else
            result = result;
    #endif
#endif

        Y[gid] = result;
    }


#undef N_DIVIDABLE_BY_4 
#undef HAS_INLINE
}


kernel void sgemm_0d98ffd842960aa3b91ffe4db8338ca6d8dff5ae17c5b2fe737a1019(const device float *weight_buffer[[buffer(0)]],
                          device float *data_buffer[[buffer(1)]],
                          const device int * meta_buffer [[buffer(2)]],
                          ushort index[[thread_index_in_threadgroup]],
                          ushort2 group_position[[threadgroup_position_in_grid]])
{
#define TRANSPOSE_A 1
#define TRANSPOSE_B 1
#define M_DIVIDABLE_BY_64 1
#define N_DIVIDABLE_BY_64 1
#define K_DIVIDABLE_BY_8 1

#if TRANSPOSE_A
    #define A_STRIDE_K 1
    #define A_STRIDE_M K
#else
    #define A_STRIDE_K M
    #define A_STRIDE_M 1
#endif

#if TRANSPOSE_B
    #define B_STRIDE_K N
    #define B_STRIDE_N 1
#else
    #define B_STRIDE_K 1
    #define B_STRIDE_N K
#endif

#define WITH_BIAS 0
#define HAS_INLINE 1


#if K_DIVIDABLE_BY_8 && M_DIVIDABLE_BY_64  && N_DIVIDABLE_BY_64 && !TRANSPOSE_A && TRANSPOSE_B && OPTIMIZE
    const device float4 *load_target4 = (index & 32) 
        ? (const device float4 *)(weight_buffer + meta_buffer[1]) 
        : (const device float4 *)(data_buffer + meta_buffer[0]);
#else
    const device float *load_target = (index & 32) 
        ? (weight_buffer + meta_buffer[1]) 
        : (data_buffer + meta_buffer[0]);
#endif

    const int M = meta_buffer[4];
    const int N = meta_buffer[5];

    const int K = meta_buffer[6];

    threadgroup float4 shared4[32 * 8 * 2];

    const int stride_k = (index & 32) ? B_STRIDE_K : A_STRIDE_K;
    const int stride_mn = (index & 32) ? B_STRIDE_N : A_STRIDE_M;

    const int m_offset = index & 7;
    const int n_offset = index >> 3;

    const int mn_load_offset = ((index & 32) ? group_position.y : group_position.x) * 64 + (index & 15) * 4;
    const int k_load_offset = ((index & 16) ? 1 : 0);

    int track0 = mn_load_offset * stride_mn + (k_load_offset + 0) * stride_k;
    int track2 = track0 + 2 * stride_k;
    int track4 = track0 + 4 * stride_k;
    int track6 = track0 + 6 * stride_k;

#if !OPTIMIZE || !M_DIVIDABLE_BY_64 || !N_DIVIDABLE_BY_64
    const int max_MN = (index & 32) ? N : M;
#endif

    int shared_offset4 = ((index & 32) ? 16 : 0) + k_load_offset * 32 + (index & 15);
    int read_A_offset4 = m_offset * 2;
    int read_B_offset4 = n_offset * 2 + 16;

    float4 result[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int k = 0;

    while (k < K)
    {
        {
#if OPTIMIZE && K_DIVIDABLE_BY_8
    #if OPTIMIZE && M_DIVIDABLE_BY_64 && N_DIVIDABLE_BY_64
        #if OPTIMIZE && !TRANSPOSE_A && TRANSPOSE_B
            shared4[shared_offset4 + 32 * 0] = load_target4[track0 >> 2];
            shared4[shared_offset4 + 32 * 2] = load_target4[track2 >> 2];
            shared4[shared_offset4 + 32 * 4] = load_target4[track4 >> 2];
            shared4[shared_offset4 + 32 * 6] = load_target4[track6 >> 2];
        #else
            shared4[shared_offset4 + 32 * 0] = float4(
                load_target[track0 + stride_mn * 0],
                load_target[track0 + stride_mn * 1],
                load_target[track0 + stride_mn * 2],
                load_target[track0 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 2] = float4(
                load_target[track2 + stride_mn * 0],
                load_target[track2 + stride_mn * 1],
                load_target[track2 + stride_mn * 2],
                load_target[track2 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 4] = float4(
                load_target[track4 + stride_mn * 0],
                load_target[track4 + stride_mn * 1],
                load_target[track4 + stride_mn * 2],
                load_target[track4 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 6] = float4(
                load_target[track6 + stride_mn * 0],
                load_target[track6 + stride_mn * 1],
                load_target[track6 + stride_mn * 2],
                load_target[track6 + stride_mn * 3]
            ); 
        #endif
    #else
            shared4[shared_offset4 + 32 * 0] = float4(
                (mn_load_offset + 0 >= max_MN) ? 0 : load_target[track0 + stride_mn * 0],
                (mn_load_offset + 1 >= max_MN) ? 0 : load_target[track0 + stride_mn * 1],
                (mn_load_offset + 2 >= max_MN) ? 0 : load_target[track0 + stride_mn * 2],
                (mn_load_offset + 3 >= max_MN) ? 0 : load_target[track0 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 2] = float4(
                (mn_load_offset + 0 >= max_MN) ? 0 : load_target[track2 + stride_mn * 0],
                (mn_load_offset + 1 >= max_MN) ? 0 : load_target[track2 + stride_mn * 1],
                (mn_load_offset + 2 >= max_MN) ? 0 : load_target[track2 + stride_mn * 2],
                (mn_load_offset + 3 >= max_MN) ? 0 : load_target[track2 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 4] = float4(
                (mn_load_offset + 0 >= max_MN) ? 0 : load_target[track4 + stride_mn * 0],
                (mn_load_offset + 1 >= max_MN) ? 0 : load_target[track4 + stride_mn * 1],
                (mn_load_offset + 2 >= max_MN) ? 0 : load_target[track4 + stride_mn * 2],
                (mn_load_offset + 3 >= max_MN) ? 0 : load_target[track4 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 6] = float4(
                (mn_load_offset + 0 >= max_MN) ? 0 : load_target[track6 + stride_mn * 0],
                (mn_load_offset + 1 >= max_MN) ? 0 : load_target[track6 + stride_mn * 1],
                (mn_load_offset + 2 >= max_MN) ? 0 : load_target[track6 + stride_mn * 2],
                (mn_load_offset + 3 >= max_MN) ? 0 : load_target[track6 + stride_mn * 3]
            ); 
    #endif

            k += 8;
#else
    #if OPTIMIZE && M_DIVIDABLE_BY_64 && N_DIVIDABLE_BY_64
            shared4[shared_offset4 + 32 * 0] = float4(
                (k + k_load_offset >= K) ? 0 : load_target[track0 + stride_mn * 0],
                (k + k_load_offset >= K) ? 0 : load_target[track0 + stride_mn * 1],
                (k + k_load_offset >= K) ? 0 : load_target[track0 + stride_mn * 2],
                (k + k_load_offset >= K) ? 0 : load_target[track0 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 2] = float4(
                (k + k_load_offset >= K) ? 0 : load_target[track2 + stride_mn * 0],
                (k + k_load_offset >= K) ? 0 : load_target[track2 + stride_mn * 1],
                (k + k_load_offset >= K) ? 0 : load_target[track2 + stride_mn * 2],
                (k + k_load_offset >= K) ? 0 : load_target[track2 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 4] = float4(
                (k + k_load_offset >= K) ? 0 : load_target[track4 + stride_mn * 0],
                (k + k_load_offset >= K) ? 0 : load_target[track4 + stride_mn * 1],
                (k + k_load_offset >= K) ? 0 : load_target[track4 + stride_mn * 2],
                (k + k_load_offset >= K) ? 0 : load_target[track4 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 6] = float4(
                (k + k_load_offset >= K) ? 0 : load_target[track6 + stride_mn * 0],
                (k + k_load_offset >= K) ? 0 : load_target[track6 + stride_mn * 1],
                (k + k_load_offset >= K) ? 0 : load_target[track6 + stride_mn * 2],
                (k + k_load_offset >= K) ? 0 : load_target[track6 + stride_mn * 3]
            ); 
            k += 2;
    #else
            shared4[shared_offset4 + 32 * 0] = float4(
                (k + k_load_offset >= K || mn_load_offset + 0 >= max_MN) ? 0 : load_target[track0 + stride_mn * 0],
                (k + k_load_offset >= K || mn_load_offset + 1 >= max_MN) ? 0 : load_target[track0 + stride_mn * 1],
                (k + k_load_offset >= K || mn_load_offset + 2 >= max_MN) ? 0 : load_target[track0 + stride_mn * 2],
                (k + k_load_offset >= K || mn_load_offset + 3 >= max_MN) ? 0 : load_target[track0 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 2] = float4(
                (k + k_load_offset >= K || mn_load_offset + 0 >= max_MN) ? 0 : load_target[track2 + stride_mn * 0],
                (k + k_load_offset >= K || mn_load_offset + 1 >= max_MN) ? 0 : load_target[track2 + stride_mn * 1],
                (k + k_load_offset >= K || mn_load_offset + 2 >= max_MN) ? 0 : load_target[track2 + stride_mn * 2],
                (k + k_load_offset >= K || mn_load_offset + 3 >= max_MN) ? 0 : load_target[track2 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 4] = float4(
                (k + k_load_offset >= K || mn_load_offset + 0 >= max_MN) ? 0 : load_target[track4 + stride_mn * 0],
                (k + k_load_offset >= K || mn_load_offset + 1 >= max_MN) ? 0 : load_target[track4 + stride_mn * 1],
                (k + k_load_offset >= K || mn_load_offset + 2 >= max_MN) ? 0 : load_target[track4 + stride_mn * 2],
                (k + k_load_offset >= K || mn_load_offset + 3 >= max_MN) ? 0 : load_target[track4 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 6] = float4(
                (k + k_load_offset >= K || mn_load_offset + 0 >= max_MN) ? 0 : load_target[track6 + stride_mn * 0],
                (k + k_load_offset >= K || mn_load_offset + 1 >= max_MN) ? 0 : load_target[track6 + stride_mn * 1],
                (k + k_load_offset >= K || mn_load_offset + 2 >= max_MN) ? 0 : load_target[track6 + stride_mn * 2],
                (k + k_load_offset >= K || mn_load_offset + 3 >= max_MN) ? 0 : load_target[track6 + stride_mn * 3]
            ); 
            k += 2;
    #endif
#endif
        }

        threadgroup_barrier(mem_flags::mem_threadgroup);

        {
            for (int k_sub = 0; k_sub < 8; k_sub++)
            {
                float4 a00 = shared4[32 * k_sub + read_A_offset4 + 0];
                float4 a01 = shared4[32 * k_sub + read_A_offset4 + 1];
                float4 b00 = shared4[32 * k_sub + read_B_offset4 + 0];
                float4 b01 = shared4[32 * k_sub + read_B_offset4 + 1];

                result[4][0]  += b00[0] * a00[2];
                result[4][1]  += b00[1] * a00[2];
                result[0][1]  += b00[1] * a00[0];
                result[0][0]  += b00[0] * a00[0];
                result[6][0]  += b00[0] * a00[3];
                result[6][1]  += b00[1] * a00[3];
                result[2][1]  += b00[1] * a00[1];
                result[2][0]  += b00[0] * a00[1];
                result[12][0] += b00[0] * a01[2];
                result[12][1] += b00[1] * a01[2];
                result[8][1]  += b00[1] * a01[0];
                result[8][0]  += b00[0] * a01[0];
                result[14][0] += b00[0] * a01[3];
                result[14][1] += b00[1] * a01[3];
                result[10][1] += b00[1] * a01[1];
                result[10][0] += b00[0] * a01[1];

                result[14][2] += b00[2] * a01[3];
                result[14][3] += b00[3] * a01[3];
                result[10][3] += b00[3] * a01[1];
                result[10][2] += b00[2] * a01[1];
                result[12][2] += b00[2] * a01[2];
                result[12][3] += b00[3] * a01[2];
                result[8][3]  += b00[3] * a01[0];
                result[8][2]  += b00[2] * a01[0];
                result[6][2]  += b00[2] * a00[3];
                result[6][3]  += b00[3] * a00[3];
                result[2][3]  += b00[3] * a00[1];
                result[2][2]  += b00[2] * a00[1];
                result[4][2]  += b00[2] * a00[2];
                result[4][3]  += b00[3] * a00[2];
                result[0][3]  += b00[3] * a00[0];
                result[0][2]  += b00[2] * a00[0];

                result[5][0]  += b01[0] * a00[2];
                result[5][1]  += b01[1] * a00[2];
                result[1][1]  += b01[1] * a00[0];
                result[1][0]  += b01[0] * a00[0];
                result[7][0]  += b01[0] * a00[3];
                result[7][1]  += b01[1] * a00[3];
                result[3][1]  += b01[1] * a00[1];
                result[3][0]  += b01[0] * a00[1];
                result[13][0] += b01[0] * a01[2];
                result[13][1] += b01[1] * a01[2];
                result[9][1]  += b01[1] * a01[0];
                result[9][0]  += b01[0] * a01[0];
                result[15][0] += b01[0] * a01[3];
                result[15][1] += b01[1] * a01[3];
                result[11][1] += b01[1] * a01[1];
                result[11][0] += b01[0] * a01[1];

                result[15][2] += b01[2] * a01[3];
                result[15][3] += b01[3] * a01[3];
                result[11][3] += b01[3] * a01[1];
                result[11][2] += b01[2] * a01[1];
                result[13][2] += b01[2] * a01[2];
                result[13][3] += b01[3] * a01[2];
                result[9][3]  += b01[3] * a01[0];
                result[9][2]  += b01[2] * a01[0];
                result[7][2]  += b01[2] * a00[3];
                result[7][3]  += b01[3] * a00[3];
                result[3][3]  += b01[3] * a00[1];
                result[3][2]  += b01[2] * a00[1];
                result[5][2]  += b01[2] * a00[2];
                result[5][3]  += b01[3] * a00[2];
                result[1][3]  += b01[3] * a00[0];
                result[1][2]  += b01[2] * a00[0];
            }
        }

        shared_offset4 ^= 32 * 8;
        read_A_offset4 ^= 32 * 8;
        read_B_offset4 ^= 32 * 8;
        track0 += stride_k * 8;
        track2 += stride_k * 8;
        track4 += stride_k * 8;
        track6 += stride_k * 8;
    }

    {
    
#if OPTIMIZE && N_DIVIDABLE_BY_64
    #if WITH_BIAS
        float4 b[2];
        const device float4 *bias4 = (const device float4 *)(weight_buffer + meta_buffer[3]);
        b[0] = bias4[group_position.y * 16 + n_offset * 2 + 0];
        b[1] = bias4[group_position.y * 16 + n_offset * 2 + 1];
    #endif
    
        device float4 *C4 = (device float4 *)(data_buffer + meta_buffer[2]);
        const int N4 = N >> 2;
        int m = group_position.x * 64 + m_offset * 8;
        for (int m_sub = 0; m_sub < 8; m_sub++)
        {

    #if !M_DIVIDABLE_BY_64
            if (m >= M) continue;
    #endif

            const int n = group_position.y * 16 + n_offset * 2;
            float4 result0 = result[m_sub * 2 + 0];
            float4 result1 = result[m_sub * 2 + 1];

    #if WITH_BIAS
            result0 += b[0];
            result1 += b[1];
    #endif

    #if HAS_INLINE
            result0[0] = result0[0];
            result0[1] = result0[1];
            result0[2] = result0[2];
            result0[3] = result0[3];
            result1[0] = result1[0];
            result1[1] = result1[1];
            result1[2] = result1[2];
            result1[3] = result1[3];
    #endif

            C4[m * N4 + n + 0] = result0;
            C4[m * N4 + n + 1] = result1;
            
            m++;
        }
#else
    #if WITH_BIAS
        const device float *bias = weight_buffer + meta_buffer[3];
        float b[8];
        for (int n_sub = 0; n_sub < 8; n_sub++)
        {
            b[n_sub] = (group_position.y * 64 + n_offset * 8 + n_sub < N)
                ? bias[group_position.y * 64 + n_offset * 8 + n_sub]
                : 0;
        }
    #endif

        device float *C = data_buffer + meta_buffer[2];
        int m = group_position.x * 64 + m_offset * 8;
        for (int m_sub = 0; m_sub < 8; m_sub++)
        {
            int n = group_position.y * 64 + n_offset * 8;

            for (int n_sub1 = 0; n_sub1 < 2; n_sub1++)
            {
                for (int n_sub2 = 0; n_sub2 < 4; n_sub2++)
                {

    #if WITH_BIAS
        #if OPTIMIZE && M_DIVIDABLE_BY_64
                    (         n < N) ? (C[m * N + n] = result[m_sub * 2 + n_sub1][n_sub2] + b[n_sub1*4+n_sub2]) : 0;
        #else
                    (m < M && n < N) ? (C[m * N + n] = result[m_sub * 2 + n_sub1][n_sub2] + b[n_sub1*4+n_sub2]) : 0;
        #endif
    #else
        #if OPTIMIZE && M_DIVIDABLE_BY_64
                    (         n < N) ? (C[m * N + n] = result[m_sub * 2 + n_sub1][n_sub2]) : 0;
        #else
                    (m < M && n < N) ? (C[m * N + n] = result[m_sub * 2 + n_sub1][n_sub2]) : 0;
        #endif
    #endif
                    n++;
                }
            }
            
            m++;
        }
#endif

    }


#undef M_DIVIDABLE_BY_64
#undef N_DIVIDABLE_BY_64
#undef K_DIVIDABLE_BY_8
#undef TRANSPOSE_A
#undef TRANSPOSE_B
#undef A_STRIDE_K
#undef B_STRIDE_K
#undef A_STRIDE_M
#undef A_STRIDE_M
#undef WITH_BIAS
#undef HAS_INLINE
}


kernel void col2im_4acb2698e5a17440fd70c40ff9ad52a79bf489ffdf7bf17453d6b928(const device float *param_buffer[[buffer(0)]],
                          device float *data_buffer[[buffer(1)]],
                          const device int * meta_buffer [[buffer(2)]],
                          uint index[[thread_position_in_grid]],
                          uint num_threads[[threads_per_grid]])
{
    const device float *col = data_buffer + meta_buffer[1];
    device float *im = data_buffer + meta_buffer[0];

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

    for (int gid = index; gid < N*H1*W1*C1; gid += num_threads) {
        const int c1 = gid % C1;
        const int w1 = gid / C1 % W1;
        const int h1 = gid / C1 / W1 % H1;
        const int n  = gid / C1 / W1 / H1;

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


kernel void elu_c5570868e1ef853205898365aac544f7d5b3e4f367a8cccdbf2beda5(const device float *param_buffer[[buffer(0)]],
                          device float *data_buffer[[buffer(1)]],
                          const device int * meta_buffer [[buffer(2)]],
                          uint index[[thread_position_in_grid]],
                          uint num_threads[[threads_per_grid]])
{
    const device float *X = data_buffer + meta_buffer[0];
    device float *Y = data_buffer + meta_buffer[1];

    const int N = meta_buffer[2];
  
    for (int gid = index; gid < N; gid += num_threads) {
        float result = X[gid];
        result = result < 0.0 ? (exp(result)-1) : result;      

        Y[gid] = result;
    }
}


kernel void sgemm_56f1c9d6ee847286248ba3248b17ad2991a0c71fea8d0f7cfc7fdb15(const device float *weight_buffer[[buffer(0)]],
                          device float *data_buffer[[buffer(1)]],
                          const device int * meta_buffer [[buffer(2)]],
                          ushort index[[thread_index_in_threadgroup]],
                          ushort2 group_position[[threadgroup_position_in_grid]])
{
#define TRANSPOSE_A 1
#define TRANSPOSE_B 1
#define M_DIVIDABLE_BY_64 1
#define N_DIVIDABLE_BY_64 0
#define K_DIVIDABLE_BY_8 1

#if TRANSPOSE_A
    #define A_STRIDE_K 1
    #define A_STRIDE_M K
#else
    #define A_STRIDE_K M
    #define A_STRIDE_M 1
#endif

#if TRANSPOSE_B
    #define B_STRIDE_K N
    #define B_STRIDE_N 1
#else
    #define B_STRIDE_K 1
    #define B_STRIDE_N K
#endif

#define WITH_BIAS 0
#define HAS_INLINE 1


#if K_DIVIDABLE_BY_8 && M_DIVIDABLE_BY_64  && N_DIVIDABLE_BY_64 && !TRANSPOSE_A && TRANSPOSE_B && OPTIMIZE
    const device float4 *load_target4 = (index & 32) 
        ? (const device float4 *)(weight_buffer + meta_buffer[1]) 
        : (const device float4 *)(data_buffer + meta_buffer[0]);
#else
    const device float *load_target = (index & 32) 
        ? (weight_buffer + meta_buffer[1]) 
        : (data_buffer + meta_buffer[0]);
#endif

    const int M = meta_buffer[4];
    const int N = meta_buffer[5];

    const int K = meta_buffer[6];

    threadgroup float4 shared4[32 * 8 * 2];

    const int stride_k = (index & 32) ? B_STRIDE_K : A_STRIDE_K;
    const int stride_mn = (index & 32) ? B_STRIDE_N : A_STRIDE_M;

    const int m_offset = index & 7;
    const int n_offset = index >> 3;

    const int mn_load_offset = ((index & 32) ? group_position.y : group_position.x) * 64 + (index & 15) * 4;
    const int k_load_offset = ((index & 16) ? 1 : 0);

    int track0 = mn_load_offset * stride_mn + (k_load_offset + 0) * stride_k;
    int track2 = track0 + 2 * stride_k;
    int track4 = track0 + 4 * stride_k;
    int track6 = track0 + 6 * stride_k;

#if !OPTIMIZE || !M_DIVIDABLE_BY_64 || !N_DIVIDABLE_BY_64
    const int max_MN = (index & 32) ? N : M;
#endif

    int shared_offset4 = ((index & 32) ? 16 : 0) + k_load_offset * 32 + (index & 15);
    int read_A_offset4 = m_offset * 2;
    int read_B_offset4 = n_offset * 2 + 16;

    float4 result[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int k = 0;

    while (k < K)
    {
        {
#if OPTIMIZE && K_DIVIDABLE_BY_8
    #if OPTIMIZE && M_DIVIDABLE_BY_64 && N_DIVIDABLE_BY_64
        #if OPTIMIZE && !TRANSPOSE_A && TRANSPOSE_B
            shared4[shared_offset4 + 32 * 0] = load_target4[track0 >> 2];
            shared4[shared_offset4 + 32 * 2] = load_target4[track2 >> 2];
            shared4[shared_offset4 + 32 * 4] = load_target4[track4 >> 2];
            shared4[shared_offset4 + 32 * 6] = load_target4[track6 >> 2];
        #else
            shared4[shared_offset4 + 32 * 0] = float4(
                load_target[track0 + stride_mn * 0],
                load_target[track0 + stride_mn * 1],
                load_target[track0 + stride_mn * 2],
                load_target[track0 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 2] = float4(
                load_target[track2 + stride_mn * 0],
                load_target[track2 + stride_mn * 1],
                load_target[track2 + stride_mn * 2],
                load_target[track2 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 4] = float4(
                load_target[track4 + stride_mn * 0],
                load_target[track4 + stride_mn * 1],
                load_target[track4 + stride_mn * 2],
                load_target[track4 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 6] = float4(
                load_target[track6 + stride_mn * 0],
                load_target[track6 + stride_mn * 1],
                load_target[track6 + stride_mn * 2],
                load_target[track6 + stride_mn * 3]
            ); 
        #endif
    #else
            shared4[shared_offset4 + 32 * 0] = float4(
                (mn_load_offset + 0 >= max_MN) ? 0 : load_target[track0 + stride_mn * 0],
                (mn_load_offset + 1 >= max_MN) ? 0 : load_target[track0 + stride_mn * 1],
                (mn_load_offset + 2 >= max_MN) ? 0 : load_target[track0 + stride_mn * 2],
                (mn_load_offset + 3 >= max_MN) ? 0 : load_target[track0 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 2] = float4(
                (mn_load_offset + 0 >= max_MN) ? 0 : load_target[track2 + stride_mn * 0],
                (mn_load_offset + 1 >= max_MN) ? 0 : load_target[track2 + stride_mn * 1],
                (mn_load_offset + 2 >= max_MN) ? 0 : load_target[track2 + stride_mn * 2],
                (mn_load_offset + 3 >= max_MN) ? 0 : load_target[track2 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 4] = float4(
                (mn_load_offset + 0 >= max_MN) ? 0 : load_target[track4 + stride_mn * 0],
                (mn_load_offset + 1 >= max_MN) ? 0 : load_target[track4 + stride_mn * 1],
                (mn_load_offset + 2 >= max_MN) ? 0 : load_target[track4 + stride_mn * 2],
                (mn_load_offset + 3 >= max_MN) ? 0 : load_target[track4 + stride_mn * 3]
            ); 
            shared4[shared_offset4 + 32 * 6] = float4(
                (mn_load_offset + 0 >= max_MN) ? 0 : load_target[track6 + stride_mn * 0],
                (mn_load_offset + 1 >= max_MN) ? 0 : load_target[track6 + stride_mn * 1],
                (mn_load_offset + 2 >= max_MN) ? 0 : load_target[track6 + stride_mn * 2],
                (mn_load_offset + 3 >= max_MN) ? 0 : load_target[track6 + stride_mn * 3]
            ); 
    #endif

            k += 8;
#else
    #if OPTIMIZE && M_DIVIDABLE_BY_64 && N_DIVIDABLE_BY_64
            shared4[shared_offset4 + 32 * 0] = float4(
                (k + k_load_offset >= K) ? 0 : load_target[track0 + stride_mn * 0],
                (k + k_load_offset >= K) ? 0 : load_target[track0 + stride_mn * 1],
                (k + k_load_offset >= K) ? 0 : load_target[track0 + stride_mn * 2],
                (k + k_load_offset >= K) ? 0 : load_target[track0 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 2] = float4(
                (k + k_load_offset >= K) ? 0 : load_target[track2 + stride_mn * 0],
                (k + k_load_offset >= K) ? 0 : load_target[track2 + stride_mn * 1],
                (k + k_load_offset >= K) ? 0 : load_target[track2 + stride_mn * 2],
                (k + k_load_offset >= K) ? 0 : load_target[track2 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 4] = float4(
                (k + k_load_offset >= K) ? 0 : load_target[track4 + stride_mn * 0],
                (k + k_load_offset >= K) ? 0 : load_target[track4 + stride_mn * 1],
                (k + k_load_offset >= K) ? 0 : load_target[track4 + stride_mn * 2],
                (k + k_load_offset >= K) ? 0 : load_target[track4 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 6] = float4(
                (k + k_load_offset >= K) ? 0 : load_target[track6 + stride_mn * 0],
                (k + k_load_offset >= K) ? 0 : load_target[track6 + stride_mn * 1],
                (k + k_load_offset >= K) ? 0 : load_target[track6 + stride_mn * 2],
                (k + k_load_offset >= K) ? 0 : load_target[track6 + stride_mn * 3]
            ); 
            k += 2;
    #else
            shared4[shared_offset4 + 32 * 0] = float4(
                (k + k_load_offset >= K || mn_load_offset + 0 >= max_MN) ? 0 : load_target[track0 + stride_mn * 0],
                (k + k_load_offset >= K || mn_load_offset + 1 >= max_MN) ? 0 : load_target[track0 + stride_mn * 1],
                (k + k_load_offset >= K || mn_load_offset + 2 >= max_MN) ? 0 : load_target[track0 + stride_mn * 2],
                (k + k_load_offset >= K || mn_load_offset + 3 >= max_MN) ? 0 : load_target[track0 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 2] = float4(
                (k + k_load_offset >= K || mn_load_offset + 0 >= max_MN) ? 0 : load_target[track2 + stride_mn * 0],
                (k + k_load_offset >= K || mn_load_offset + 1 >= max_MN) ? 0 : load_target[track2 + stride_mn * 1],
                (k + k_load_offset >= K || mn_load_offset + 2 >= max_MN) ? 0 : load_target[track2 + stride_mn * 2],
                (k + k_load_offset >= K || mn_load_offset + 3 >= max_MN) ? 0 : load_target[track2 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 4] = float4(
                (k + k_load_offset >= K || mn_load_offset + 0 >= max_MN) ? 0 : load_target[track4 + stride_mn * 0],
                (k + k_load_offset >= K || mn_load_offset + 1 >= max_MN) ? 0 : load_target[track4 + stride_mn * 1],
                (k + k_load_offset >= K || mn_load_offset + 2 >= max_MN) ? 0 : load_target[track4 + stride_mn * 2],
                (k + k_load_offset >= K || mn_load_offset + 3 >= max_MN) ? 0 : load_target[track4 + stride_mn * 3]
            ); 
            k += 2;

            shared4[shared_offset4 + 32 * 6] = float4(
                (k + k_load_offset >= K || mn_load_offset + 0 >= max_MN) ? 0 : load_target[track6 + stride_mn * 0],
                (k + k_load_offset >= K || mn_load_offset + 1 >= max_MN) ? 0 : load_target[track6 + stride_mn * 1],
                (k + k_load_offset >= K || mn_load_offset + 2 >= max_MN) ? 0 : load_target[track6 + stride_mn * 2],
                (k + k_load_offset >= K || mn_load_offset + 3 >= max_MN) ? 0 : load_target[track6 + stride_mn * 3]
            ); 
            k += 2;
    #endif
#endif
        }

        threadgroup_barrier(mem_flags::mem_threadgroup);

        {
            for (int k_sub = 0; k_sub < 8; k_sub++)
            {
                float4 a00 = shared4[32 * k_sub + read_A_offset4 + 0];
                float4 a01 = shared4[32 * k_sub + read_A_offset4 + 1];
                float4 b00 = shared4[32 * k_sub + read_B_offset4 + 0];
                float4 b01 = shared4[32 * k_sub + read_B_offset4 + 1];

                result[4][0]  += b00[0] * a00[2];
                result[4][1]  += b00[1] * a00[2];
                result[0][1]  += b00[1] * a00[0];
                result[0][0]  += b00[0] * a00[0];
                result[6][0]  += b00[0] * a00[3];
                result[6][1]  += b00[1] * a00[3];
                result[2][1]  += b00[1] * a00[1];
                result[2][0]  += b00[0] * a00[1];
                result[12][0] += b00[0] * a01[2];
                result[12][1] += b00[1] * a01[2];
                result[8][1]  += b00[1] * a01[0];
                result[8][0]  += b00[0] * a01[0];
                result[14][0] += b00[0] * a01[3];
                result[14][1] += b00[1] * a01[3];
                result[10][1] += b00[1] * a01[1];
                result[10][0] += b00[0] * a01[1];

                result[14][2] += b00[2] * a01[3];
                result[14][3] += b00[3] * a01[3];
                result[10][3] += b00[3] * a01[1];
                result[10][2] += b00[2] * a01[1];
                result[12][2] += b00[2] * a01[2];
                result[12][3] += b00[3] * a01[2];
                result[8][3]  += b00[3] * a01[0];
                result[8][2]  += b00[2] * a01[0];
                result[6][2]  += b00[2] * a00[3];
                result[6][3]  += b00[3] * a00[3];
                result[2][3]  += b00[3] * a00[1];
                result[2][2]  += b00[2] * a00[1];
                result[4][2]  += b00[2] * a00[2];
                result[4][3]  += b00[3] * a00[2];
                result[0][3]  += b00[3] * a00[0];
                result[0][2]  += b00[2] * a00[0];

                result[5][0]  += b01[0] * a00[2];
                result[5][1]  += b01[1] * a00[2];
                result[1][1]  += b01[1] * a00[0];
                result[1][0]  += b01[0] * a00[0];
                result[7][0]  += b01[0] * a00[3];
                result[7][1]  += b01[1] * a00[3];
                result[3][1]  += b01[1] * a00[1];
                result[3][0]  += b01[0] * a00[1];
                result[13][0] += b01[0] * a01[2];
                result[13][1] += b01[1] * a01[2];
                result[9][1]  += b01[1] * a01[0];
                result[9][0]  += b01[0] * a01[0];
                result[15][0] += b01[0] * a01[3];
                result[15][1] += b01[1] * a01[3];
                result[11][1] += b01[1] * a01[1];
                result[11][0] += b01[0] * a01[1];

                result[15][2] += b01[2] * a01[3];
                result[15][3] += b01[3] * a01[3];
                result[11][3] += b01[3] * a01[1];
                result[11][2] += b01[2] * a01[1];
                result[13][2] += b01[2] * a01[2];
                result[13][3] += b01[3] * a01[2];
                result[9][3]  += b01[3] * a01[0];
                result[9][2]  += b01[2] * a01[0];
                result[7][2]  += b01[2] * a00[3];
                result[7][3]  += b01[3] * a00[3];
                result[3][3]  += b01[3] * a00[1];
                result[3][2]  += b01[2] * a00[1];
                result[5][2]  += b01[2] * a00[2];
                result[5][3]  += b01[3] * a00[2];
                result[1][3]  += b01[3] * a00[0];
                result[1][2]  += b01[2] * a00[0];
            }
        }

        shared_offset4 ^= 32 * 8;
        read_A_offset4 ^= 32 * 8;
        read_B_offset4 ^= 32 * 8;
        track0 += stride_k * 8;
        track2 += stride_k * 8;
        track4 += stride_k * 8;
        track6 += stride_k * 8;
    }

    {
    
#if OPTIMIZE && N_DIVIDABLE_BY_64
    #if WITH_BIAS
        float4 b[2];
        const device float4 *bias4 = (const device float4 *)(weight_buffer + meta_buffer[3]);
        b[0] = bias4[group_position.y * 16 + n_offset * 2 + 0];
        b[1] = bias4[group_position.y * 16 + n_offset * 2 + 1];
    #endif
    
        device float4 *C4 = (device float4 *)(data_buffer + meta_buffer[2]);
        const int N4 = N >> 2;
        int m = group_position.x * 64 + m_offset * 8;
        for (int m_sub = 0; m_sub < 8; m_sub++)
        {

    #if !M_DIVIDABLE_BY_64
            if (m >= M) continue;
    #endif

            const int n = group_position.y * 16 + n_offset * 2;
            float4 result0 = result[m_sub * 2 + 0];
            float4 result1 = result[m_sub * 2 + 1];

    #if WITH_BIAS
            result0 += b[0];
            result1 += b[1];
    #endif

    #if HAS_INLINE
            result0[0] = result0[0];
            result0[1] = result0[1];
            result0[2] = result0[2];
            result0[3] = result0[3];
            result1[0] = result1[0];
            result1[1] = result1[1];
            result1[2] = result1[2];
            result1[3] = result1[3];
    #endif

            C4[m * N4 + n + 0] = result0;
            C4[m * N4 + n + 1] = result1;
            
            m++;
        }
#else
    #if WITH_BIAS
        const device float *bias = weight_buffer + meta_buffer[3];
        float b[8];
        for (int n_sub = 0; n_sub < 8; n_sub++)
        {
            b[n_sub] = (group_position.y * 64 + n_offset * 8 + n_sub < N)
                ? bias[group_position.y * 64 + n_offset * 8 + n_sub]
                : 0;
        }
    #endif

        device float *C = data_buffer + meta_buffer[2];
        int m = group_position.x * 64 + m_offset * 8;
        for (int m_sub = 0; m_sub < 8; m_sub++)
        {
            int n = group_position.y * 64 + n_offset * 8;

            for (int n_sub1 = 0; n_sub1 < 2; n_sub1++)
            {
                for (int n_sub2 = 0; n_sub2 < 4; n_sub2++)
                {

    #if WITH_BIAS
        #if OPTIMIZE && M_DIVIDABLE_BY_64
                    (         n < N) ? (C[m * N + n] = result[m_sub * 2 + n_sub1][n_sub2] + b[n_sub1*4+n_sub2]) : 0;
        #else
                    (m < M && n < N) ? (C[m * N + n] = result[m_sub * 2 + n_sub1][n_sub2] + b[n_sub1*4+n_sub2]) : 0;
        #endif
    #else
        #if OPTIMIZE && M_DIVIDABLE_BY_64
                    (         n < N) ? (C[m * N + n] = result[m_sub * 2 + n_sub1][n_sub2]) : 0;
        #else
                    (m < M && n < N) ? (C[m * N + n] = result[m_sub * 2 + n_sub1][n_sub2]) : 0;
        #endif
    #endif
                    n++;
                }
            }
            
            m++;
        }
#endif

    }


#undef M_DIVIDABLE_BY_64
#undef N_DIVIDABLE_BY_64
#undef K_DIVIDABLE_BY_8
#undef TRANSPOSE_A
#undef TRANSPOSE_B
#undef A_STRIDE_K
#undef B_STRIDE_K
#undef A_STRIDE_M
#undef A_STRIDE_M
#undef WITH_BIAS
#undef HAS_INLINE
}


kernel void tanh_4d093f262c150685bb3ba8649b1aedd8a91b1d4e42056449061f72a2(const device float *param_buffer[[buffer(0)]],
                          device float *data_buffer[[buffer(1)]],
                          const device int * meta_buffer [[buffer(2)]],
                          uint index[[thread_position_in_grid]],
                          uint num_threads[[threads_per_grid]])
{
    const device float *X = data_buffer + meta_buffer[0];
    device float *Y = data_buffer + meta_buffer[1];

    const int N = meta_buffer[2];
  
    for (int gid = index; gid < N; gid += num_threads) {
        Y[gid] = tanh(X[gid]);
    }
}


kernel void scalaraffine_3daf34a046eb693a4a93f810672675819980cfa4239b993889317ec5(const device float *param_buffer[[buffer(0)]],
                          device float *data_buffer[[buffer(1)]],
                          const device int * meta_buffer [[buffer(2)]],
                          uint index[[thread_position_in_grid]],
                          uint num_threads[[threads_per_grid]])
{
    const device float *X = data_buffer + meta_buffer[0];
    device float *Y = data_buffer + meta_buffer[1];

    const float scale = *((const device float *)(& meta_buffer[3]));
    const float bias = *((const device float *)(& meta_buffer[4]));
    const int N = meta_buffer[2];

    for (int gid = index; gid < N; gid += num_threads) {
        float result = X[gid];
        result = result * scale + bias;
        Y[gid] = result;
    }
}
