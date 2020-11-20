// Alexander Schünemann
// 26-01-2019
// SplitOperator2D - Tool to solve Schrödinger equation on a square grid

#include <complex>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cmath>
#include <eigen3/Eigen/Dense>   // using eigen
#include <fftw3.h>              // using fftw

using namespace Eigen;

using complex = std::complex<double>;
using vector_real = std::vector<double>;
using vector_complex = std::vector<complex>;

// Container for simulation parameters
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
    Operators(Params &par, double vXoffset, double vYoffset, double wfcXoffset, double wfcYoffset) 
    {
        size = par.res;
        
        v   = ArrayXXcd::Zero(size, size);
        pe  = ArrayXXcd::Zero(size, size);
        ke  = ArrayXXcd::Zero(size, size);
        wfc = ArrayXXcd::Zero(size, size);

        for (size_t ix = 0; ix < size; ix++)
        {
            for (size_t iy = 0; iy < size; iy++)
            {
                v(ix, iy)   = 0.5 * (pow(par.x[ix] - vXoffset, 2.0) + pow(par.x[iy] - vYoffset, 2.0) );
                wfc(ix, iy) = exp(-0.5 * ( pow(par.x[ix] - wfcXoffset, 2.0) + pow(par.x[iy] - wfcYoffset, 2.0) ) );

                if (par.im_time)
                {
                    ke(ix,iy) = exp(-0.5 * par.dt * (pow(par.k[ix] , 2.0) + pow(par.k[iy] , 2.0)) );
                    pe(ix,iy) = exp(-0.5 * par.dt * v(ix,iy));
                }
                else
                {
                    ke(ix,iy) = exp(-0.5 * par.dt * (pow(par.k[ix] , 2.0) + pow(par.k[iy] , 2.0)) * complex(0.0,1.0) );
                    pe(ix,iy) = exp(-0.5 * par.dt * v(ix,iy) * complex(0.0,1.0) );
                }
                
            }
        }

    }

    size_t size;
    ArrayXXcd v;
    ArrayXXcd pe;
    ArrayXXcd ke;
    ArrayXXcd wfc;
};

void fft(ArrayXXcd &x, bool inverse)
{
    complex y[x.cols()][x.rows()];
    memset(y, 0, sizeof(y));

    fftw_plan p;

    fftw_complex *in = reinterpret_cast<fftw_complex*>(x.data());
    fftw_complex *out = reinterpret_cast<fftw_complex*>(y);
    p = fftw_plan_dft_2d(x.cols(), x.rows(), in, out, (inverse ? FFTW_BACKWARD : FFTW_FORWARD), FFTW_ESTIMATE);

    fftw_execute(p);
    fftw_destroy_plan(p);

    for (size_t ix = 0; ix < x.cols(); ix++) 
    {
        for (size_t iy = 0; iy < x.rows(); iy++)
        {
            x(ix,iy) = y[ix][iy] / sqrt( static_cast<double>(x.rows() * x.cols()) );
        }
    }
}

void split_op(Params &par, Operators &opr, bool output)
{
    ArrayXXd density = ArrayXXd::Zero(opr.size, opr.size);
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
        density = real(pow(abs(opr.wfc), 2.0));

        // renorm
        if (par.im_time) 
        {
            double sum = density.sum();
            sum *= par.dx * par.dx;
            opr.wfc /= sqrt(sum);
        }

        // output
        if (output)
        {
            std::string filename = "output" + std::to_string(t) + ".dat";
            std::ofstream off;
            off.open(filename);
            off << density << std::endl;
            off.close();
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

    std::cout << "[";
    for (int i = 0; i < 70; ++i) 
    {
        std::cout << "=";
    }
    std::cout << "] " << "100" << " %\r";
    std::cout.flush();
}

double calculate_energy(Params &par, Operators &opr) 
{
    ArrayXXcd wfc_r(opr.wfc);
    ArrayXXcd wfc_k(opr.wfc);
    ArrayXXcd wfc_c(opr.size, opr.size);
    fft(wfc_k, false);

    wfc_c = conj(wfc_r);

    ArrayXXcd energy_k(opr.size, opr.size);
    ArrayXXcd energy_r(opr.size, opr.size);

    for (size_t ix = 0; ix < opr.size; ix++) 
    {
        for (size_t iy = 0; iy < opr.size; iy++)
        {
            energy_k(ix, iy) = wfc_k(ix, iy) * ( pow( complex(par.k[ix], 0.0), 2) + pow(complex(par.k[iy], 0.0), 2) );
        }
    }

    fft(energy_k, true);

    energy_k *= 0.5 * wfc_c;
    energy_r = wfc_c * opr.v * wfc_r;

    double energy_final = real(energy_k.sum() + energy_r.sum());

    return energy_final * par.dx * par.dx;
}

int main() 
{
    Params par = Params(5.0, 512, 0.05, 500, true);
    Operators opr = Operators(par, 0.0, 0.0, -1.0, -1.0);

    split_op(par, opr, true);

    std::cout << "\nThe ground state energy is " << calculate_energy(par, opr) << " ħω\n";

    return 0;
}
