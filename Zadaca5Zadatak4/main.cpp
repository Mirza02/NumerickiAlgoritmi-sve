#include <iostream>
#include <complex>
#include <cmath>
#include <limits>

const double PI = 4 * atan(1.);

void FFT(std::vector<double> &x, std::vector<std::complex<double>> &x_kon, int N, int s = 0, int d = 0, int t = 1){
    if(N == 1){
        x_kon[d] = x[s];
    }
    else{
        FFT(x, x_kon, N / 2, s, d, 2 * t);
        FFT(x, x_kon, N / 2, s + t, d + N / 2, 2 * t);
        std::complex<double> mi(1);
        std::complex<double> w(std::cos((2 * PI) / N), std::sin((2 * PI) / N));
        w = std::pow(w, -1);
        for(int k = d; k < d + N / 2; k++){
            std::complex<double> u = x_kon[k];
            std::complex<double> v = mi * x_kon[k + N / 2];
            x_kon[k] = u + v;
            x_kon[k + N / 2] = u - v;
            mi *= w;
        }
    }
}

std::vector<double> LossyCompress(std::vector<double> data, int new_size){
    if(data.size() - std::pow(2, std::log2(data.size())) > std::numeric_limits<double>::epsilon()){
        throw std::range_error("Data size must be a power of two");
    }
    if(new_size < 0 || new_size > data.size()){
        throw std::range_error("Bad new size");
    }
    std::vector<double> y(new_size);
    for(int i = 0; i < new_size / 2; i++){
        y[i] = data[2 * i];
    }
    for(int i = new_size / 2; i < new_size; i++){
        y[i] = data[2 * (new_size - i) - 1];
    }
    std::vector<std::complex<double>> y_kon(new_size);
    FFT(y, y_kon, new_size);
    std::vector<double> x_kon(new_size);
    std::complex<double> w(std::cos((2 * PI) / (2 * new_size)), std::sin((2 * PI) / (2 * new_size)));
    for(int k = 0; k < new_size; k++){
        x_kon[k] = std::real(std::pow(w, -k / 2.0) * y_kon[k]);
    }
    x_kon[new_size - 1] = data.size();
    return x_kon;
}

void InvFFT(std::vector<std::complex<double>> &x_kon, std::vector<std::complex<double>> &x, int N, int s = 0, int d = 0, int t = 1){
    if(N == 1){
        x[d] = x_kon[s];
    }
    else{
        InvFFT(x_kon, x, N / 2, s, d, 2 * t);
        InvFFT(x_kon, x, N / 2, s + t, d + N / 2, 2 * t);
        std::complex<double> mi = 1;
        std::complex<double> w(std::cos((2 * PI) / N), std::sin((2 * PI) / N));
        for(int k = d; k < d + N / 2; k++){
            std::complex<double> u = x[k];
            std::complex<double> v = mi * x[k + N / 2];
            x[k] = (u + v) / 2.;
            x[k + N / 2] = (u - v) / 2.;
            mi *= w;
        }
    }
}

std::vector<double> LossyDecompress(std::vector<double> decompressed){
    double N = decompressed[decompressed.size() - 1];
    if(N < decompressed.size() || N - std::pow(2, std::log2(N)) > std::numeric_limits<double>::epsilon()){
        throw std::logic_error("Bad compressed sequence");
    }
    decompressed[decompressed.size() - 1] = 0;
    decompressed.resize(decompressed.size() + (N - decompressed.size()));
    std::vector<std::complex<double>> y_kon(decompressed.size());
    std::complex<double> w(std::cos((2 * PI) / (2 * N)), std::sin((2 * PI) / (2 * N)));
    y_kon[0] = decompressed[0];
    for(int k = 1; k < N; k++){
        y_kon[k] = 2. * std::pow(w, (k / 2.)) * decompressed[k];
    }
    std::vector<std::complex<double>> y(N);
    InvFFT(y_kon, y, N);
    std::vector<double> x(N);
    for(int n = 0; n < N; n++){
        if(n % 2 == 0){
            x[n] = std::real(y[n / 2]);
        }
        else{
            x[n] = std::real(y[N - (n + 1) / 2.]);
        }
    }
    return x;
}

int main() {
    std::vector<double> test1(8);                   // testiranje kompresije kada se zadrzava pocetna velicina niza
    for(int i = 0; i < test1.size(); i++){
        test1[i] = i * 1.5;
    }
    auto test1_comp = LossyCompress(test1, 8);
    std::cout << "Kompresovani podaci: ";
    for(auto x : test1_comp){
        std::cout << x << " ";
    }
    std::cout << std::endl << "Dekompresovani podaci: ";
    auto test1_decomp = LossyDecompress(test1_comp);
    for(auto x : test1_decomp){
        std::cout << x << " ";
    }

    try{                                            // testiranje kompresije za preveliku vrijednost nove velicine
        auto test2_komp = LossyCompress(test1, 10);
    }
    catch(std::range_error e){
        std::cout << std::endl << std::endl << e.what() << std::endl << std::endl;
    }
}
