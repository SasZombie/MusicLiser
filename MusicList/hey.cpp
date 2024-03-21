#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <iomanip>
#include <valarray>

using cd = std::complex<double>;
const double PI = std::acos(-1);
std::vector<float> in;
using namespace std::complex_literals;
using Complex = std::complex<double> ;
using CArray = std::valarray<Complex> ;

void dft(const std::vector<double>in, std::vector<cd>&out)
{

    int n = in.size() ;
    
   for (size_t f = 0; f < n; ++f)
    {
        out[f] = 0;

        for (size_t i = 0; i < n; ++i)
        {
            double t = static_cast<double>(i)/n;
            out[f]+= in[i]*std::exp(2*PI*t*f*1i);

        }
    }
}

void fft(CArray& x)
{
    const size_t N = x.size();
    if (N <= 1) return;

    // divide
    CArray even = x[std::slice(0, N/2, 2)];
    CArray  odd = x[std::slice(1, N/2, 2)];

    // conquer
    fft(even);
    fft(odd);

    // combine
    for (size_t k = 0; k < N/2; ++k)
    {
        double l = static_cast<double>(k)/N;
        double theta = static_cast<double>(-2) * PI * l;
        Complex t = std::polar(1.0, theta) * odd[k];
        x[k    ] = even[k] + t;
        x[k+N/2] = even[k] - t;
    }
}


void fft2(double in[], size_t stride, std::complex<double>out[], size_t n)
{
    if( n == 1)
    {
        out[0] = in[0];
        return;
    }
    fft2(in, stride*2, out, n/2);

    fft2(in + stride, stride*2, out + n/2, n/2);

    for (size_t k = 0; k < n/2; ++k)
    {
        double t = (double) k/n;
        std::complex<double> v = std::exp((2) * PI * t * 1i) * out[k+n/2];
        std::complex<double> e = out[k];
        out[k] = e + v;
        out[k+n/2] = e-v;
    }
    

}
int main()
{
    size_t n = 8;
    std::vector<cd> out(8);
    std::vector<double> in;
    double in2[n];
    std::complex<double> out2[n];

   // Complex in2[n];

    for (size_t i = 0; i < n; ++i)
    {
        double t = static_cast<double>(i)/n;

        in.emplace_back(std::sin(2*PI*t) + std::sin(2*PI*t*2));

        in2[i] = std::sin(2*PI*t) + std::sin(2*PI*t*2) + std::cos(2*PI*t*4);
      
    }
    
   // CArray data(in2, n);


    // fft(data);
    //dft(in, out);
    fft2(in2, 1, out2, n);


    for (size_t f = 0; f < n; f++)
    {
        // std::cout <<std::setprecision(2)<< out[f].real() << " " << out[f].imag()<<'\n';
        
         std::cout <<std::fixed << std::setprecision(2)<< out2[f] << " " <<'\n';

    }
    
//********************************************
    // size_t n = 9;
    // double in[n];
    // std::complex<double> out[n];

    // for(size_t i = 0; i < n; ++i)
    // {
    //     float t = (float)i/n;
    //     in[i] = sinf(2*PI * t) + sinf(2 * PI * t * 2);
    // }

    // for (size_t f = 0; f < n; ++f)
    // {
    //     out[f] = 0;

    //     for(size_t i = 0; i < n; ++i)
    //     {
    //         float t = (float)i/n;
    //         out[f]+= in[i] * std::exp(2*PI*t*f*1i);
    //     }

    // }
    
    // for (size_t f = 0; f < n; f++)
    // {
    //     std::cout <<std::setprecision(2)<< out[f].real() << " " << out[f].imag()<<'\n';
    // }
    
//**********************************
    // for (size_t i = 0; i < n; ++i)
    // {
    //     //double t = static_cast<double>(i/n);
    //     numbersFloat.emplace_back(std::cos(2*PI*i/n) + std::sin(2*PI*i/n*2));

    //     std::cout << std::cos(2*PI*i/n) + std::sin(2*PI*i/n*2) << '\n';
        
    // }


    // std::puts("\n");
    // numbers.resize(9);
    // dft(numbersFloat, numbers);

    // for (auto &&i : numbers)
    // {
    //     std::cout << i.real() << " " << i.imag() <<'\n';
    // }
    
   
}

// (2.2e-16,0) 
// (6.7e-16,4) 
// (-8.7e-16,4) 
// (4,1.3e-15) 
// (1.1e-15,0) 
// (4,-2.2e-15) 
// (-1.4e-15,-4) 
// (-2.2e-16,-4) 


// (1.3e-06,0) 
// (0,4) 
// (-5.6e-08,4) 
// (4,-9.5e-07) 
// (-1.3e-06,0) 
// (4,1.1e-06) 
// (2.9e-07,-4) 
// (2.4e-07,-4) 