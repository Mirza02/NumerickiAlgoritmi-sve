#include <iostream>
#include <vector>
#include <cmath>


const double PI = 4 * std::atan(1);

template <typename FunType>
std::pair<double, bool> RombergIntegration(FunType f, double a, double b, double eps = 1e-8,
                                           int nmax = 1000000, int nmin = 50){
    if(eps < 0 || nmax < 0 || nmin < 0 || nmax < nmin){
        throw std::domain_error("Bad parameter");
    }
    std::pair<double, bool> result;
    double N = 2;
    double h = (b - a) / 2;
    double s = (f(a) + f(b)) / 2;
    double iOld = s;
    double i = 0;
    std::vector<double> I(nmax);
    while(N < nmax){
        for(int j = 1; j <= N/2; j++){
            s = s + f(a + (2 * j - 1) * h);
        }
        I[i] = h * s;
        double p = 4;
        for(int k = i - 1; k >= 0; k--){
            I[k] = (p * I[k + 1] - I[k]) / (p - 1);
            p = 4 * p;
        }
        if(std::fabs(I[0] - iOld) <= eps && N > nmin){
            result.first = I[0];
            result.second = true;
            return result;
        }
        iOld = I[0];
        h = h / 2;
        N = N * 2;
        i++;
        result.first = I[0];
    }
    result.second = false;
    return result;
}

template <typename FunType>
std::pair<double, bool> TanhSinhIntegration(FunType f, double a, double b, double eps = 1e-8,
                                            int nmax = 1000000, int nmin = 20, double range = 3.5){
    if(eps < 0 || nmax < 0 || nmin < 0 || nmax < nmin || range < 0){
        throw std::domain_error("Bad parameter");
    }
    std::pair<double, bool> result;
    double N = 2;
    double h = (2 * range) / N;
    double p = (b + a) / 2;
    double q = (b - a) / 2;
    double s = 0;
    double iOld = s;
    double I = 0;
    while(N < nmax){
        for(int i = 1; i <= N/2; i++){
            double t = -range + (2 * i - 1) * h;
            double u = PI * sinh(t) / 2;
            double v = f(p + q * tanh(u));
            if(std::isfinite(v)){
                s = s + q * PI * cosh(t) * v / (2 * (cosh(u) * cosh(u)));
            }
            else v = 0;
        }
        I = h * s;
        if(N >= nmin && std::fabs(I - iOld) <= eps){
            result.first = I;
            result.second = true;
            return result;
        }
        iOld = I;
        N = 2 * N;
        h = h / 2;
        result.first = I;
    }
    result.second = false;
    return result;
}

template <typename FunType>
std::pair<double, bool> AdaptiveAux(FunType f, double a, double b, double eps, double f1, double f2,
                                    double f3, int maxdepth){
    double c = (a + b) / 2;
    double I1 = (b - a) * (f1 + 4 * f3 + f2) / 6;
    double f4 =f((a + c) / 2);
    double f5 = f((c + b) / 2);
    double I2 = (b - a) * (f1 + 4 * f4 + 2 * f3 + 4 * f5 + f2) / 12;
    if(!std::isfinite(I2)) I2 = 0;
    if(std::fabs(I1 - I2) <= eps){
        return {I2, true};
    }
    if(maxdepth <= 0){
        return {I2, false};
    }
    auto x = AdaptiveAux(f, a, c, eps, f1, f3, f4, maxdepth - 1);
    auto y = AdaptiveAux(f, c, b, eps, f3, f2, f5, maxdepth - 1);
    if(x.second == false || y.second == false) return {x.first + y.first, false};
    return {x.first + y.first, true};
}

template <typename FunType>
std::pair<double, bool> AdaptiveIntegration(FunType f, double a, double b, double eps = 1e-10,
                                            int maxdepth = 30, int nmin = 1){
    if(eps < 0 || maxdepth < 0 || nmin < 0){
        throw std::domain_error("Bad parameter");
    }
    double s = 0;
    double h = (b - a) / nmin;
    bool flag = false;
    for(int i = 0; i < nmin; i++){
        auto x = AdaptiveAux(f, a, a + h, eps, f(a), f(a + h), f(a + h / 2), maxdepth);
        s = s + x.first;
        flag = x.second;
        a = a + h;
    }
    return {s, flag};
}

template <typename FunType>
void RombergIntegrationOutput(FunType f, double a, double b, double eps = 1e-8, int nmax = 1000000, int nmin = 50){
    try{
        auto x = RombergIntegration(f, a, b, eps, nmax, nmin);
        std::cout << "Rombergova inegracija: " << x.first << " , konvergira: " << x.second << std::endl;
    }
    catch(std::domain_error &e){
        std::cout << e.what() << std::endl;
    }
}

template <typename FunType>
void TanhSinhIntegrationOutput(FunType f, double a, double b, double eps = 1e-8, int nmax = 1000000, int nmin = 20, double range = 3.5){
    try{
        auto x = TanhSinhIntegration(f, a, b, eps, nmax, nmin, range);
        std::cout << "TanhSinh inegracija: " << x.first << " , konvergira: " << x.second << std::endl;
    }
    catch(std::domain_error &e){
        std::cout << e.what() << std::endl;
    }
}

template <typename FunType>
void AdaptiveIntegrationOutput(FunType f, double a, double b, double eps = 1e-10, int maxdepth = 30, int nmin = 1){
    try{
        auto x = TanhSinhIntegration(f, a, b, eps, maxdepth, nmin);
        std::cout << "Adaptivna inegracija: " << x.first << " , konvergira: " << x.second << std::endl;
    }
    catch(std::domain_error &e){
        std::cout << e.what() << std::endl;
    }
}


int main() {
    /*std::cout << "Integracija funkcije sin(x) na intervalu [0, PI] koristeci sve tri metode:" << std::endl;
    auto sinfun = [](double x){ return std::sin(x); };
    RombergIntegrationOutput(sinfun, 0, PI);
    TanhSinhIntegrationOutput(sinfun, 0, PI);
    AdaptiveIntegrationOutput(sinfun, 0, PI);

    std::cout << std::endl << "Primjeri pogresnih parametara za Rombergovu integraciju: " << std::endl;
    RombergIntegrationOutput(sinfun, 1, 2, -4);
    RombergIntegrationOutput(sinfun, 1, 2, 1e-8, -5);
    RombergIntegrationOutput(sinfun, 1, 2, 1e-8, 100, -21);
    RombergIntegrationOutput(sinfun, 1, 2, 1e-8, 100, 123);

    std::cout << std::endl << "Primjeri pogresnih parametara za TanhSinh integraciju: " << std::endl;
    TanhSinhIntegrationOutput(sinfun, 3, 6, -1e-8);
    TanhSinhIntegrationOutput(sinfun, 3, 6, 1e-8, -42, 0);
    TanhSinhIntegrationOutput(sinfun, 3, 6, 1e-8, 42, -2);
    TanhSinhIntegrationOutput(sinfun, 3, 6, 1e-8, 42, 2, -12);
    TanhSinhIntegrationOutput(sinfun, 3, 6, 1e-8, 42, 50);

    std::cout << std::endl << "Primjeri pogresnih parametara za Adaptivnu integraciju: " << std::endl;
    AdaptiveIntegrationOutput(sinfun, 5, 7, -1e-8);
    AdaptiveIntegrationOutput(sinfun, 5, 7, 1e-8, -3);
    AdaptiveIntegrationOutput(sinfun, 5, 7, 1e-8, 40, -5);

    std::cout << std::endl << "Integracija funkcije 1 / sqrt(x) u intervalu [0, 1]: " << std::endl;
    auto sqfun = [](double x){ return 1 / sqrt(x); };
    TanhSinhIntegrationOutput(sqfun, 0, 1);
    AdaptiveIntegrationOutput(sqfun, 0, 1);

    std::cout << std::endl << "Integracija funkcije cos(x) u intervalu [0, 2PI]: " << std::endl;
    auto cosfun = [](double x){ return std::cos(x); };
    RombergIntegrationOutput(cosfun, 0, 2 * PI);
    TanhSinhIntegrationOutput(cosfun, 0, 2 * PI);
    AdaptiveIntegrationOutput(cosfun, 0, 2 * PI);

    std::cout << std::endl << "Integracija funkcije 1 / sin(x) u intervalu [0, PI]: " << std::endl;
    auto sinfun2 = [](double x) { return 1 / std::sin(x); };
    TanhSinhIntegrationOutput(sinfun2, 0, PI);
    AdaptiveIntegrationOutput(sinfun2, 0, PI);

    std::cout << std::endl << "Integracija funkcije e^x u intervalu [0, 1]: " << std::endl;
    auto expfun = [](double x){return std::exp(x);};
    RombergIntegrationOutput(expfun, 0, 1);
    TanhSinhIntegrationOutput(expfun, 0, 1);
    AdaptiveIntegrationOutput(expfun, 0, 1);

    std::cout << std::endl << "Integracija funkcije x^2 + 2 * x u intervalu [-8, 5]: " << std::endl;
    auto xfun = [](double x){return x * x + 2 * x;};
    RombergIntegrationOutput(xfun, -8, 5);
    TanhSinhIntegrationOutput(xfun, -8, 5);
    AdaptiveIntegrationOutput(xfun, -8, 5);

    std::cout << std::endl << "Integracija funkcije lnx u intervalu [1, 3]: " << std::endl;
    auto lnfun = [](double x){return std::log(x);};
    RombergIntegrationOutput(lnfun, 1, 3);
    TanhSinhIntegrationOutput(lnfun, 1, 3);
    AdaptiveIntegrationOutput(lnfun, 1, 3);

    std::cout << std::endl << "Testiranje razlicitih vrijednosti nmin i nmax u Rombergovoj integraciji za razlicite funkcije: " << std::endl;
    RombergIntegrationOutput(sinfun, 0, PI, 1e-8, 5000, 5);
    RombergIntegrationOutput(sinfun, 0, PI, 1e-8, 500, 5);
    RombergIntegrationOutput(sinfun, 0, PI, 1e-8, 5, 1);

    RombergIntegrationOutput(xfun, -5, 5, 1e-8, 10000, 5);
    RombergIntegrationOutput(xfun, -5, 5, 1e-8, 5, 1);

    RombergIntegrationOutput(lnfun, 1, 3, 1e-8, 1000, 2);
    RombergIntegrationOutput(lnfun, 1, 3, 1e-8, 100, 2);
    RombergIntegrationOutput(lnfun, 1, 3, 1e-8, 3, 2);

    std::cout << std::endl << "Testiranje razlicitih nmax, nmin i range vrijednosti u  TanhSinh integraciji: " << std::endl;
    TanhSinhIntegrationOutput(sinfun, 0, PI, 1e-8, 5000, 5, 10);
    TanhSinhIntegrationOutput(sinfun, 0, PI, 1e-8, 5000, 5, 1);
    TanhSinhIntegrationOutput(sinfun, 0, PI, 1e-8, 500, 40, 1);

    TanhSinhIntegrationOutput(xfun, -5, 5, 1e-8, 10000, 55, 10);
    TanhSinhIntegrationOutput(xfun, -5, 5, 1e-8, 10000, 55, 2);
    TanhSinhIntegrationOutput(xfun, -5, 5, 1e-8, 1000000, 55, 2);
    TanhSinhIntegrationOutput(xfun, -5, 5, 1e-8, 1000, 55, 3.5);

    TanhSinhIntegrationOutput(lnfun, 1, 3, 1e-8, 100, 2, 1);
    TanhSinhIntegrationOutput(lnfun, 1, 3, 1e-8, 100000, 2, 1);
    TanhSinhIntegrationOutput(lnfun, 1, 3, 1e-8, 100000, 2, 3);

    std::cout << std::endl << "Testiranje razlicitih maxdepth i nmin vrijednosti za adaptivnu integraciju: " << std::endl;
    AdaptiveIntegrationOutput(sinfun, 0, PI, 1e-8, 10, 2);
    AdaptiveIntegrationOutput(sinfun, 0, PI, 1e-8, 10, 1);
    AdaptiveIntegrationOutput(sinfun, 0, PI, 1e-8, 15, 1);
    AdaptiveIntegrationOutput(sinfun, 0, PI, 1e-8, 35, 1);

    AdaptiveIntegrationOutput(lnfun, 0, 5, 1e-8, 5);
    AdaptiveIntegrationOutput(lnfun, 0, 5, 1e-8, 100);
    AdaptiveIntegrationOutput(lnfun, 0, 5, 1e-8, 100, 3);

    AdaptiveIntegrationOutput(cosfun, 0, PI, 1e-8, 2);
    AdaptiveIntegrationOutput(cosfun, 0, PI, 1e-8, 30);*/

    //AT15 - Adaptive integration - 6
    auto aig =  AdaptiveIntegration([](double x) { return 1 / std::sqrt(x); }, 0, 1, 1e-6, 40);
    std::cout << aig.first << " " << aig.second;





}