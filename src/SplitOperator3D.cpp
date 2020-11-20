// Alexander Schünemann
// 26-01-2019
// SplitOperator3D - Tool to solve Schrödinger equation on a cubic orthonormal grid

#define EIGEN_WARNINGS_DISABLED

#include <complex>
#include <vector>
#include <iostream>
#include <cstring>
#include <cmath>
#include <eigen3/Eigen/Dense>                       // using eigen
#include <eigen3/unsupported/Eigen/CXX11/Tensor>    // using eigen tensor
#include <fftw3.h>                                  // using fft

using namespace Eigen;

using complex = std::complex<double>;
using vector_real = std::vector<double>;
using vector_complex = std::vector<complex>;

struct Params 
{
    Params(double _xmax, unsigned int _res, double _dt, unsigned int _timesteps, bool im) 
    {
        xmax = _xmax;
        res = _res;
        dt = _dt;
        timesteps = _timesteps;
        dx = 2.0 * xmax / res;
        x.reserve(res);
        dk = M_PI / xmax;
        k.reserve(res);
        im_time = im;

        for (size_t i = 0; i < res; ++i) {
            x.emplace_back(xmax / res - xmax + i * (2.0 * xmax / res));
            if (i < res / 2) {
                k.push_back(i * M_PI / xmax);
            } else {
                k.push_back((static_cast<double>(i) - res) * M_PI / xmax);
            }
        }
    }

    double xmax;
    unsigned int res;
    double dt;
    unsigned int timesteps;
    double dx;
    vector_real x;
    double dk;
    vector_real k;
    bool im_time;
};

struct Operators 
{
    public:
    Operators(Params &par, double wfcXoffset, double wfcYoffset, double wfcZoffset) 
    {
        size = par.res;

        v   = Tensor<complex, 3>(size, size, size);
        pe  = Tensor<complex, 3>(size, size, size);
        ke  = Tensor<complex, 3>(size, size, size);
        wfc = Tensor<complex, 3>(size, size, size);

        for (size_t ix = 0; ix < size; ix++)
        {
            for (size_t iy = 0; iy < size; iy++)
            {
                for (size_t iz = 0; iz < size; iz++)
                {
                    // Setup potential as harmonic
                    v(ix,iy,iz) = 0.5 * ( pow( par.x[ix], 2) + pow( par.x[iy], 2) + pow( par.x[iz], 2) );
                    // Setup wfc as gaussian
                    wfc(ix,iy,iz) = 0.5 * exp(0.5 * ( pow( par.x[ix]-wfcXoffset, 2) + pow( par.x[iy]-wfcYoffset, 2) + pow( par.x[iz]-wfcZoffset, 2) ) );
                    
                    // Setup Hamiltonian
                    if (par.im_time)
                    {
                        ke(ix,iy,iz) = exp(-0.5 * par.dt * (pow(par.k[ix] , 2.0) + pow(par.k[iy] , 2.0) + pow(par.k[iz] , 2.0)) );
                        pe(ix,iy,iz) = exp(-0.5 * par.dt * v(ix,iy,iz));
                    }
                    else
                    {
                        ke(ix,iy,iz) = exp(-0.5 * par.dt * (pow(par.k[ix] , 2.0) + pow(par.k[iy] , 2.0) + pow(par.k[iz] , 2.0)) * complex(0.0,1.0) );
                        pe(ix,iy,iz) = exp(-0.5 * par.dt * v(ix,iy,iz) * complex(0.0,1.0) );
                    }
                }
            }
        }

    }

    size_t size;
    Tensor<complex, 3> v;
    Tensor<complex, 3> pe;
    Tensor<complex, 3> ke;
    Tensor<complex, 3> wfc;
};

void fft(Tensor<complex, 3> &x, bool inverse)
{
    const auto& d = x.dimensions();
    complex y[d[0]][d[1]][d[2]];
    memset(y, 0, sizeof(y));
    fftw_plan p;

    fftw_complex *in = reinterpret_cast<fftw_complex*>(x.data());
    fftw_complex *out = reinterpret_cast<fftw_complex*>(y);
    p = fftw_plan_dft_3d(d[0], d[1], d[2], in, out, (inverse ? FFTW_BACKWARD : FFTW_FORWARD), FFTW_ESTIMATE);

    fftw_execute(p);
    fftw_destroy_plan(p);

    for (size_t ix = 0; ix < d[0]; ix++) 
    {
        for (size_t iy = 0; iy < d[1]; iy++)
        {
            for (size_t iz = 0; iz < d[2]; iz++)
            {
                x(ix,iy,iz) = y[ix][iy][iz] / sqrt( static_cast<double>( d[0] * d[1] * d[2] ) );
            }
        }
    }
}

void split_op(Params &par, Operators &opr)
{
    Tensor<double, 3> density(opr.size,opr.size,opr.size);
    float progress = 0.0f;

    for (size_t t = 0; t < par.timesteps; t++)
    {
        // Half step in real space
        opr.wfc *= opr.pe;

        // Full step in momentum space
        fft(opr.wfc, false);
        opr.wfc *= opr.ke;

        // Final step in real space
        fft(opr.wfc, true);
        opr.wfc *= opr.pe;

        // get density
        for (size_t ix = 0; ix < opr.size; ix++)
        {
            for (size_t iy = 0; iy < opr.size; iy++)
            {
                for (size_t iz = 0; iz < opr.size; iz++)
                {
                    density(ix,iy,iz) = std::real(pow(abs(opr.wfc(ix,iy,iz)), 2.0));
                }
            }
        }

        // renorm
        if (par.im_time) 
        {
            double sum = 0;
            for (size_t ix = 0; ix < opr.size; ix++)
            {
                for (size_t iy = 0; iy < opr.size; iy++)
                {
                    for (size_t iz = 0; iz < opr.size; iz++)
                    {
                        sum += density(ix,iy,iz);
                    }
                }
            }
            sum *= par.dx * par.dx * par.dx;
            for (size_t ix = 0; ix < opr.size; ix++)
            {
                for (size_t iy = 0; iy < opr.size; iy++)
                {
                    for (size_t iz = 0; iz < opr.size; iz++)
                    {
                        opr.wfc(ix,iy,iz) = opr.wfc(ix,iy,iz) / sqrt(sum);
                    }
                }
            }
        }

        
        // progress bar
        int barWidth = 70;

        std::cout << "[";
        int pos = barWidth * progress;
        for (int i = 0; i < barWidth; ++i) 
        {
            if (i < pos) {std::cout << "=";}
            else if (i == pos) {std::cout << ">";}
            else {std::cout << " ";}
        }

        std::cout << "] " << int(progress * 100.0) << " %\r";
        std::cout.flush();

        progress += (float)(1.0 / par.timesteps);
    }
}

double calculate_energy(Params &par, Operators &opr) 
{
    Tensor<complex, 3> wfc_r(opr.wfc);
    Tensor<complex, 3> wfc_k(opr.wfc);
    Tensor<complex, 3> wfc_c(opr.size,opr.size,opr.size);
    fft(wfc_k, false);

    for (size_t ix = 0; ix < opr.size; ix++) 
    {
        for (size_t iy = 0; iy < opr.size; iy++)
        {
            for (size_t iz = 0; iz < opr.size; iz++)
            {
                wfc_c(ix,iy,iz) = conj(wfc_r(ix,iy,iz));
            }
        }
    }

    Tensor<complex, 3> energy_k(opr.size,opr.size,opr.size);
    Tensor<complex, 3> energy_r(opr.size,opr.size,opr.size);

    for (size_t ix = 0; ix < opr.size; ix++) 
    {
        for (size_t iy = 0; iy < opr.size; iy++)
        {
            for (size_t iz = 0; iz < opr.size; iz++)
            {
                energy_k(ix,iy,iz) = wfc_k(ix,iy,iz) * ( pow( complex(par.k[ix], 0.0), 2) + pow(complex(par.k[iy], 0.0), 2) + pow(complex(par.k[iz], 0.0), 2) );
            }
        }
    }

    fft(energy_k, true);

    energy_k *= 0.5 * wfc_c;
    energy_r = wfc_c * opr.v * wfc_r;

    double energy_final = 0.0;

    for (size_t ix = 0; ix < opr.size; ix++) 
    {
        for (size_t iy = 0; iy < opr.size; iy++)
        {
            for (size_t iz = 0; iz < opr.size; iz++)
            {
                energy_final += std::real(energy_k(ix,iy,iz) + energy_r(ix,iy,iz));
            }
        }
    }

    return energy_final * par.dx * par.dx * par.dx;
}

int main() 
{
    // 80 is upper limit for grid size
    Params par = Params(5.0, 80, 0.05, 500, true);
    Operators opr = Operators(par, -1.0, -1.0, -1.0);

    split_op(par, opr);

    std::cout << "\nThe ground state energy is " << calculate_energy(par, opr) << " ħω\n";

    return 0;
}
