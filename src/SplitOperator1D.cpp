#include <complex>
#include <vector>
#include <iostream>
#include <cstring>
#include <fstream>

// Using fftw3 library.
#include <fftw3.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using complex = std::complex<double>;
using vector_real = std::vector<double>;
using vector_complex = std::vector<complex>;

struct Params {
    Params(double _xmax, unsigned int _res, double _dt, unsigned int _timesteps, bool im) {
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

struct Operators {
public:
    Operators(Params &par, double voffset,
              double wfcoffset) {
        size = par.res;
        v.reserve(size);
        pe.reserve(size);
        ke.reserve(size);
        wfc.reserve(size);

        for (size_t i = 0; i < size; ++i) {

            v.push_back(0.5 * pow(par.x[i] - voffset, 2));
            wfc.push_back(exp(-pow(par.x[i] - wfcoffset, 2) / 2.0));

            if (par.im_time) {
                ke.push_back(exp(-0.5 * par.dt * pow(par.k[i], 2)));
                pe.push_back(exp(-0.5 * par.dt * v[i]));
            } else {
                ke.push_back(exp(-0.5 * par.dt * pow(par.k[i], 2) * complex(0.0, 1.0)));
                pe.push_back(exp(-0.5 * par.dt * v[i] * complex(0.0, 1.0)));
            }
        }
    }

    size_t size;
    vector_complex v;
    vector_complex pe;
    vector_complex ke;
    vector_complex wfc;
};

void fft(vector_complex &x, int n, bool inverse) {
    complex y[n];
    memset(y, 0, sizeof(y));
    fftw_plan p;

    fftw_complex *in = reinterpret_cast<fftw_complex*>(x.data());
    fftw_complex *out = reinterpret_cast<fftw_complex*>(y);
    p = fftw_plan_dft_1d(n, in, out, (inverse ? FFTW_BACKWARD : FFTW_FORWARD), FFTW_ESTIMATE);

    fftw_execute(p);
    fftw_destroy_plan(p);

    for (size_t i = 0; i < n; ++i) {
        x[i] = y[i] / sqrt(static_cast<double>(n));
    }
}

void split_op(Params &par, Operators &opr) {
    double density[opr.size];

    float progress = 0.0f;

    for (size_t i = 0; i < par.timesteps; ++i) 
    {
        // Half step in real space
        for (size_t j = 0; j < opr.size; ++j) { opr.wfc[j] *= opr.pe[j]; }

        // Full step in k space
        fft(opr.wfc, opr.size, false);
        for (size_t j = 0; j < opr.size; ++j) { opr.wfc[j] *= opr.ke[j]; }

        // Half step in real space
        fft(opr.wfc, opr.size, true);
        for (size_t j = 0; j < opr.size; ++j) { opr.wfc[j] *= opr.pe[j]; }

        // Get |psi(x)|^2
        for (size_t j = 0; j < opr.size; ++j) { density[j] = pow(abs(opr.wfc[j]), 2); }

        // renorm
        if (par.im_time) 
        {
            double sum = 0;
            for (size_t j = 0; j < opr.size; ++j) { sum += density[j]; }
            sum *= par.dx;
            for (size_t j = 0; j < opr.size; ++j) { opr.wfc[j] /= sqrt(sum); }
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

        /*
        // Writing data into a file in the format of:
        // index, density, real potential.
        std::stringstream filename_stream;
        filename_stream << "output" << i << ".dat";

        std::ofstream fstream = std::ofstream(filename_stream.str());

        if (fstream) {
            for (int i = 0; i < opr.size; ++i) {
                std::stringstream data_stream;

                data_stream << i << "\t" << density[i] << "\t" << real(opr.v[i]) << "\n";

                fstream.write(data_stream.str().c_str(), data_stream.str().length());
            }
        }

        fstream.close();
        */
    }
}

double calculate_energy(Params &par, Operators &opr) {
    vector_complex wfc_r(opr.wfc);
    vector_complex wfc_k(opr.wfc);
    vector_complex wfc_c(opr.size);
    fft(wfc_k, opr.size, false);

    for (size_t i = 0; i < opr.size; ++i) {
        wfc_c[i] = conj(wfc_r[i]);
    }

    vector_complex energy_k(opr.size);
    vector_complex energy_r(opr.size);

    for (size_t i = 0; i < opr.size; ++i) {
        energy_k[i] = wfc_k[i] * pow(complex(par.k[i], 0.0), 2);
    }

    fft(energy_k, opr.size, true);

    for (size_t i = 0; i < opr.size; ++i) {
        energy_k[i] *= 0.5 * wfc_c[i];
        energy_r[i] = wfc_c[i] * opr.v[i] * wfc_r[i];
    }

    double energy_final = 0;

    for (size_t i = 0; i < opr.size; ++i) {
        energy_final += real(energy_k[i] + energy_r[i]);
    }

    return energy_final * par.dx;
}

int main() {
    Params par = Params(5.0, 512, 0.05, 1000, true);
    Operators opr = Operators(par, 0.0, -1.0);

    split_op(par, opr);

    std::cout << "\nThe energy is " << calculate_energy(par, opr) << "\n";

    return 0;
}
