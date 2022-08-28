#include "dgemm.h"
#include <cstdio>
#include <cstdlib>
#include <immintrin.h>
#include <omp.h>

#pragma GCC optimize ("-Ofast")

void dgemm(float alpha, const float *a, const float *b, float beta, float *c)
{
    for (int i = 0; i < MATRIX_SIZE; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            c[i * MATRIX_SIZE + j] *= beta;
            for (int k = 0; k < MATRIX_SIZE; k++)
            {
                c[i * MATRIX_SIZE + j] += alpha * a[i * MATRIX_SIZE + k] * b[j * MATRIX_SIZE + k];
            }
        }
    }
}

void dgemm_simd(float alpha, const float *a, const float *b, float beta, float *c, float *d)
{
    // #pragma omp parallel for
    // for (int n = 0; n < MATRIX_SIZE*MATRIX_SIZE; ++n)
    // {
    //     /* code */
    //     int i = n/(MATRIX_SIZE);
    //     int j = n%MATRIX_SIZE;

    //     d[n] = b[MATRIX_SIZE*j + i];
    // }
    
    for (int i = 0; i < MATRIX_SIZE; ++i)
    {
        /* code */
        for (int j = 0; j < MATRIX_SIZE; ++j)
        {
            /* code */
            int ub = MATRIX_SIZE - (MATRIX_SIZE%8);
            c[i*MATRIX_SIZE + j] *= beta;
            float mid_vector[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
            float mid_value = 0.0;
            for (int k = 0; k < ub; k+=8)
            {
                /* code */
                __m256 a_vector = _mm256_loadu_ps(a+i*MATRIX_SIZE + k);
                __m256 b_vector = _mm256_loadu_ps(b+j*MATRIX_SIZE + k); //need to transpose b matrix and then do this
                __m256 result = _mm256_mul_ps(a_vector, b_vector);

                // result = _mm256_mul_ps(alpha_vector, result);

                _mm256_store_ps(mid_vector, result);

                for (int l = 0; l < 8; ++l)
                {
                    /* code */
                    mid_value += mid_vector[l];
                }
            }
            
            c[i*MATRIX_SIZE + j] += alpha*mid_value;
            for (int k = ub; k < MATRIX_SIZE; ++k){
                c[i * MATRIX_SIZE + j] += alpha * a[i * MATRIX_SIZE + k] * b[j * MATRIX_SIZE + k];
            }
        }
        
    }
}

int main(int, char **)
{
    float alpha, beta;

    // mem allocations
    int mem_size = MATRIX_SIZE * MATRIX_SIZE * sizeof(float);
    auto a = (float *)malloc(mem_size);
    auto b = (float *)malloc(mem_size);
    auto c = (float *)malloc(mem_size);
    auto d = (float *)malloc(mem_size);

    // check if allocated
    if (nullptr == a || nullptr == b || nullptr == c)
    {
        printf("Memory allocation failed\n");
        if (nullptr != a)
            free(a);
        if (nullptr != b)
            free(b);
        if (nullptr != c)
            free(c);
        return 0;
    }

    generateProblemFromInput(alpha, a, b, beta, c);

    std::cerr << "Launching dgemm step." << std::endl;
    // matrix-multiplication
    dgemm_simd(alpha, a, b, beta, c, d);

    outputSolution(c);

    free(a);
    free(b);
    free(c);
    return 0;
}