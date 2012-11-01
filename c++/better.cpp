#include <complex>
#include <math.h>
#include <algorithm>
#include <alloca.h>
#include "tick.h"

#include <iostream>
using std::cout;
using std::cin;
using std::endl;

typedef std::complex<double> complex;

const complex I = complex(0,1);

void dft(complex *in, complex *out, int N)
{
	for (int i=0; i<N; ++i) {
		complex sum(0,0);
		for (int j=0; j<N; ++j) {
			double v = -j*M_PI*2*i/N;
			sum += in[j] * exp(-I*v);
		}
		out[i] = sum;
	}
}

void bit_reverse(const long n, complex * const x)
{
	int j=0;
	for (int i=0; i<n-1; ++i) { // count i up in normal order
		if (j > i) {
			std::swap(x[j], x[i]);
		}

		// count j up in bit-reversed order
		int k = n >> 1;
		while (k <= j) {
			j -= k;
			k >>= 1;
		}
		j += k;
	}
}

void myfft52(const long N, complex *x)
{
	// shuffle x[] to bit-reversed order
	bit_reverse(N, x);

	// do log2(n) iterations
	int fft_len = 2;
	int fft_num = N/2;
	int pos = 0;
	while (fft_len <= N) {
		for (int m=0; m<fft_len/2; ++m) {
			complex W_k = exp(I*complex(2.0*M_PI*double(m)/double(fft_len)));
			for (int i=m; i<N; i+=fft_len) {
				int j=i+fft_len/2;
				complex a =       x[i];
				complex b = W_k * x[j];
				x[i] = a + b;
				x[j] = a - b;
			}
		}

		// next outer iteration
		fft_len <<= 1;
		fft_num >>= 1;
	}
}

complex foo[65536];

void myfft_init(const long N)
{
	// find out M with 2^M==N
	int M = 0;
	while ((1<<M) < N) ++M;

	int pos = 0;

	// do log2(n) iterations
	int fft_len = 2;
	while (fft_len <= N) {
		foo[pos++] = exp(I*complex(1.0*M_PI/double(fft_len)));
		fft_len <<= 1;
	}

	cout << "foo len: " << pos << endl;
}

inline void fft_sub_2(int N, complex *x)
{
	for (int i=0; i<N; i+=2) {
		const complex a = x[i];
		const complex b = x[i+1];
		x[i+0] = a + b;
		x[i+1] = a - b;
	}
}

inline void fft_sub_4(int N, complex *x)
{
	for (int i=0; i<N; i+=4) {
		const complex a = x[i];
		const complex c = x[i+1];
		const complex b = x[i+2];
		const complex d(-imag(x[i+3]), real(x[i+3]));
		x[i+0] = a + b;
		x[i+1] = c + d;
		x[i+2] = a - b;
		x[i+3] = c - d;
	}
}

void fft_sub_8(int N, complex *x)
{
	const complex r = exp(I*complex(2.0*M_PI/8.0));

	const complex W_0(1.0);
	const complex W_1 = W_0 * r;
	const complex W_2 = W_1 * r;
	const complex W_3 = W_2 * r;

	for (int i=0; i<N; i+=8) {
		const complex a = x[i+0];
		const complex b = x[i+4];
		x[i+0] = a + b;
		x[i+4] = a - b;

		const complex c =       x[i+1];
		const complex d = W_1 * x[i+5];
		x[i+1] = c + d;
		x[i+5] = c - d;

		const complex e =       x[i+2];
		const complex f = W_2 * x[i+6];
		x[i+2] = e + f;
		x[i+6] = e - f;

		const complex g =       x[i+3];
		const complex h = W_3 * x[i+7];
		x[i+3] = g + h;
		x[i+7] = g - h;
	}
}

inline void fft_sub_2_4(int N, complex *x)
{
	for (int i=0; i<N; i+=4) {
		const complex a = x[i+0];
		const complex b = x[i+1];
		const complex c = x[i+2];
		const complex d = x[i+3];
		const complex h(imag(d) - imag(c), real(c) - real(d));
		x[i+0] = a + b + c + d;
		x[i+1] = a - b + h;
		x[i+2] = a + b - c - d;
		x[i+3] = a - b - h;
	}
}

void fft_sub_2_4_8(int N, complex *x)
{
	const complex r = exp(I*complex(2.0*M_PI/8.0));
	const complex W_1 = r;
	const complex W_3 = I * r;
	const complex W_2 = I;

	for (int i=0; i<N; i+=8) {
		const complex a = x[i+0];
		const complex b = x[i+1];
		const complex c = x[i+2];
		const complex d = x[i+3];
		const complex e = x[i+4];
		const complex f = x[i+5];
		const complex g = x[i+6];
		const complex h = x[i+7];
		const complex t0(imag(d) - imag(c), real(c) - real(d));
		const complex t1(imag(h) - imag(g), real(g) - real(h));
		const complex k = a + b + c + d;
		const complex l = e + f + g + h;
		const complex m = a - b + t0;
		const complex n = W_1 * (e - f + t1);
		const complex o = a + b - c - d;
		const complex p = W_2 * (e + f - g - h);
		const complex q = a - b - t0;
		const complex r = W_3 * (e - f - t1);
		x[i+0] = k + l;
		x[i+1] = m + n;
		x[i+2] = o + p;
		x[i+3] = q + r;
		x[i+4] = k - l;
		x[i+5] = m - n;
		x[i+6] = o - p;
		x[i+7] = q - r;
	}
}

void fft_sub(int fft_len, int fft_m, int N, complex *x)
{
	const complex r = foo[fft_m-1];

	complex W_k(1.0);
	for (int k=0; k<fft_len/2; ++k) {
		for (int i=k; i<N; i+=fft_len) {
			const int j=i+fft_len/2;
			const complex a =       x[i];
			const complex b = W_k * x[j];
			x[i] = a + b;
			x[j] = a - b;
		}
		W_k *= r;
		//cout << "k=" << k << ", " << W_k << endl;
	}
}

void myfft(const long N, complex *x)
{
	// shuffle x[] to bit-reversed order
	bit_reverse(N, x);

	// do log2(n) iterations
	int fft_m   = 1;
	int fft_len = 2;
	while (fft_len <= N) {
		switch (fft_m) {
			/*case 1: break;
			case 2: break;
			case 3: fft_sub_2_4_8(N, x); break;*/
			default: fft_sub(fft_len, fft_m, N, x); break;
		}

		// next outer iteration
		fft_len <<= 1;
		fft_m += 1;
	}
}

void myfft2(const long N, complex *x)
{
	// shuffle x[] to bit-reversed order
	bit_reverse(N, x);

	// do log2(n) iterations
	int fft_len = 2;
	int fft_num = N/2;
	while (fft_len <= N) {
		for (int m=0; m<fft_len/2; ++m) {
			complex W_m = exp(I*complex(M_PI*m/double(fft_len/2.0)));
			for (int i=m; i<N; i+=fft_len) {
				int j=i+fft_len/2;
				complex b = W_m * x[j];
				x[i] += b;
				x[j] = x[i] - 2.0*b;
			}
		}

		// next outer iteration
		fft_len <<= 1;
		fft_num >>= 1;
	}
}

int test2()
{
	const int T = 100;
	const int N = 65536;
	complex *in = new complex[N];
	
	cout << N << " points fft" << endl;
	for (int i=0; i<N; ++i) in[i] = complex(rand()/double(RAND_MAX),rand()/double(RAND_MAX));

	tic();
	for (int i=0; i<T; ++i) myfft(N, in);
	int dt = toc();
	cout << "Elapsed time per fft: " << double(dt) / T << endl;
}

int main(int argc, char *argv[])
{
	const int TIMES = 1;
	const long M = 3;
	const long N = 1<<M;

	cout << N << " points fft" << endl;

	complex *in = new complex [N];
	complex *out = new complex [N];

	/*srand(time(0));
	for (int i=0; i<N; ++i) in[i] = complex(rand()/double(RAND_MAX),rand()/double(RAND_MAX));*/
	for (int i=0; i<N; ++i) in[i] = complex(i,0);

	dft(in, out, N);
	myfft_init(65536);
	myfft2(N, in);

	/*tic();
	FFT(1, M, in);
	toc(true);*/
#if 1
	for (int i=0; i<N; ++i) {
		cout << out[i] << endl;
	}
	cout << endl;
	for (int i=0; i<N; ++i) {
		cout << in[i] << endl;
	}
	cout << endl;
#endif
	//test2();
}


