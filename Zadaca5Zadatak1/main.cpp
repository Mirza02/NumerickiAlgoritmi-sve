#include <iostream>
#include <cmath>
#include <functional>
#include <limits>
#include <complex>
#include <random>


enum RegulaFalsiMode {Unmodified, Illinois, Slavic, IllinoisSlavic};

template <typename FunType>
bool BracketRoot(FunType f, double x0, double &a, double &b, double hinit = 1e-5, double hmax = 1e10, double lambda = 1.4, int depth = 0){
    if(hinit < 0 || hmax < 0 || lambda < 0) throw std::domain_error("Invalid parameters");
    double h = hinit;
    a = x0;
    double f1 = f(a);
    while(std::fabs(h) < hmax){
        b = a + h;
        double f2 = f(b);
        while(!std::isfinite(f2)){
            h = h / (2 * (1 + lambda));
            if(std::fabs(h) <= (std::fabs(a) * std::numeric_limits<double>::epsilon())){
                if(depth == 0){
                    return BracketRoot(f, -x0, a, b, hinit, hmax, lambda, ++depth);
                }
                return false;
            }
            b = a + h;
            f2 = f(b);
        }
        if(f1 * f2 <= 0){
            return true;
        }
        h = lambda * h;
        a = b;
        f1 = f2;
    }
    if(depth == 0){
        return BracketRoot(f, -x0, a, b, hinit, hmax, lambda, ++depth);
    }
    return false;
}

template <typename FunType>
double RegulaFalsiSolve(FunType f, double a, double b, RegulaFalsiMode mode = Slavic, double eps = 1e-10, int maxiter = 100){
    if(eps < 0 || maxiter < 0) throw std::domain_error("Invalid parameters");
    if((f(a) * f(b)) > 0) throw std::range_error("Root must be bracketed");
    std::function<double(double)> fi = [](double x){ return x / (1 + std::fabs(x)); };
    double f1 = f(a);
    double f2 = f(b);
    if(mode == Slavic || mode == IllinoisSlavic){
        f1 = fi(f(a));
        f2 = fi(f(b));
    }
    double c = a;
    double cold = b;
    double briter = 0;
    while(std::fabs(c - cold) > std::numeric_limits<double>::epsilon()){
        if(briter == maxiter) throw std::domain_error("Given accuracy has not achieved");
        cold = c;
        c = (a * f2 - b * f1) / (f2 - f1);
        double f3 = f(c);
        if(mode == Slavic){
            f3 = fi(f(c));
        }
        if(std::fabs(f3 - 0) < std::numeric_limits<double>::epsilon()) return c;
        if(mode == Illinois || mode == IllinoisSlavic){
            if((f1 * f3) < 0){
                b = a;
                f2 = f1;
            }
            else f2 = f2 / 2;
        }
        else{
            if((f1 * f3) < 0){
                b = a;
                f2 = f1;
            }
        }
        a = c;
        f1 = f3;
        briter++;
    }
    return c;
}

template <typename FunType>
double RiddersSolve(FunType f, double a, double b, double eps = 1e-10, int maxiter = 100){
    if(eps < 0 || maxiter < 0) throw std::domain_error("Invalid parameters");
    if((f(a) * f(b)) > 0) throw std::range_error("Root must be bracketed");
    double f1 = f(a);
    double f2 = f(b);
    double briter = 0;
    while(std::fabs(b - a) > std::numeric_limits<double>::epsilon()){
        if(briter == maxiter) throw std::domain_error("Given accuracy has not achieved");
        double c = (a + b) / 2;
        double f3 = f(c);
        if(std::fabs(f3 - 0) < std::numeric_limits<double>::epsilon()) return c;
        double sgn;
        if((f1 - f2) < 0) sgn = -1;
        else sgn = 1;
        double d = c + (f3 * (c - a) * sgn) / (std::sqrt((f3 * f3) - (f1 * f2)));
        double f4 = f(d);
        if(std::fabs(f4 - 0) < std::numeric_limits<double>::epsilon()) return d;
        if ((f3 * f4) <= 0){
            a = c;
            b = d;
            f1 = f3;
            f2 = f4;
        }
        else if((f1 * f4) <= 0){
            b = d;
            f2 = f4;
        }
        else{
            a = d;
            f1 = f4;
        }
        briter++;
    }
    return (a + b) / 2;
}

template <typename FunType1, typename FunType2>
double NewtonRaphsonSolve(FunType1 f, FunType2 fprim, double x0, double eps = 1e-10, double damping = 0, int maxiter = 100){
    if(eps < 0 || maxiter < 0 || damping < 0 || damping >= 1) throw std::domain_error("Invalid parameters");
    int briter = 0;
    double delta = std::numeric_limits<double>::infinity();
    double d = fprim(x0);
    while(std::fabs(delta) > std::numeric_limits<double>::epsilon()){
        if(briter == maxiter) throw std::logic_error("Convergence has not achieved");
        double v = f(x0);
        if(std::fabs(v) <= std::numeric_limits<double>::epsilon()){
            return x0;
        }
        delta = v / d;
        if(std::fabs(damping) > std::numeric_limits<double>::epsilon()){
            double w = v;
            v = f(x0 - delta);
            d = fprim(x0 - delta);
            while(std::fabs(v) > std::fabs(w) || !std::isfinite(v) || d == 0){
                delta = damping * delta;
                v = f(x0 - delta);
                d = fprim(x0 - delta);
            }
        }
        x0 = x0 - delta;
        if(std::fabs(fprim(x0) < std::numeric_limits<double>::epsilon())) throw std::logic_error("Convergence has not achieved");
        if(!std::isfinite(f(x0))) throw std::logic_error("Convergence has not achieved");
        briter++;
    }
    return x0;
}


std::complex<double> RandomComplex(double re1, double re2, double im1, double im2){
    std::uniform_real_distribution<double> unifRe(re1, re2);
    std::uniform_real_distribution<double> unifIm(im1, im2);
    std::default_random_engine re;

    double randomRe = unifRe(re);
    double randomIm = unifIm(re);

    std::complex<double> result(randomRe, randomIm);
    return result;

}

std::pair<std::complex<double>, bool> Laguerre(std::vector<std::complex<double>> coefficients, int stepen, std::complex<double> x0, double maxiters = 100){
    std::complex<double> delta = std::numeric_limits<double>::infinity();
    int k = 1;
    while((std::abs(delta) > std::numeric_limits<double>::epsilon() && (k < maxiters))){
        std::complex<double> f = coefficients[stepen - 1];
        std::complex<double> d = 0;
        std::complex<double> s = 0;
        for(int i = stepen - 1; i >= 0; i--){
            s = s * x0 + 2. * d;
            d = d * x0 + f;
            f = f * x0 + coefficients[i];
        }
        if(std::fabs(f.real()) < std::numeric_limits<double>::epsilon() && std::fabs(f.imag()) < std::numeric_limits<double>::epsilon()){
            return {x0, true};
        }
        std::complex<double> r = std::sqrt((stepen - 1.) * ((stepen - 1.) * d * d  - static_cast<double>(stepen) * f * s));
        if(std::abs(d + r) > std::abs(d - r)){
            delta = static_cast<double>(stepen) * f / (d + r);
        }
        else{
            delta = static_cast<double>(stepen) * f / (d - r);
        }
        x0 = x0 - delta;
        k++;
    }
    if(std::abs(delta) < std::numeric_limits<double>::epsilon()){
        return {x0, true};
    }
    return {x0, false};
}

std::pair<double, bool> Laguerre(std::vector<double> coefficients, int stepen, double x0, double maxiters = 100){
    double delta = std::numeric_limits<double>::infinity();
    int k = 1;
    while((std::abs(delta) > std::numeric_limits<double>::epsilon() && (k < maxiters))){
        double f = coefficients[stepen - 1];
        double d = 0;
        double s = 0;
        for(int i = stepen - 1; i >= 0; i--){
            s = s * x0 + 2. * d;
            d = d * x0 + f;
            f = f * x0 + coefficients[i];
        }
        if(std::fabs(f) < std::numeric_limits<double>::epsilon()){
            return {x0, true};
        }
        double r = std::sqrt((stepen - 1.) * ((stepen - 1.) * d * d  - static_cast<double>(stepen) * f * s));
        if(std::abs(d + r) > std::abs(d - r)){
            delta = static_cast<double>(stepen) * f / (d + r);
        }
        else{
            delta = static_cast<double>(stepen) * f / (d - r);
        }
        x0 = x0 - delta;
        k++;
    }
    if(std::abs(delta) < std::numeric_limits<double>::epsilon()){
        return {x0, true};
    }
    return {x0, false};
}

std::vector<std::complex<double>> PolyRoots(std::vector<std::complex<double>> coefficients, double eps = 1e-10, int maxiters = 100, int maxtrials = 10){
    if(eps < 0 || maxiters < 0 || maxtrials < 0) throw std::domain_error("Invalid parameters");
    std::vector<std::complex<double>> nule(coefficients.size() - 1);
    for(int i = coefficients.size() - 2; i >= 0; i--){
        std::complex<double> x;
        int briters = 1;
        bool c = false;
        while(!c || (briters < maxiters)){
            x = RandomComplex(-10, 10, -10, 10);
            auto pair = Laguerre(coefficients, i, x);
            x = pair.first;
            c = pair.second;
            briters++;
        }
        if(!c){
            throw std::range_error("Convergence has not achieved");
        }
        if(std::fabs(x.imag()) <= std::numeric_limits<double>::epsilon()){
            x = x.real();
        }
        nule[i] = x;
        auto v = coefficients[i];
        for(int j = i - 1; j >= 0; j--){
            auto w = coefficients[j];
            coefficients[j] = v;
            v = w + (x * v);
        }
    }
    return nule;
}

double RandomDouble(double re1, double re2){
    std::uniform_real_distribution<double> unifRe(re1, re2);
    std::default_random_engine re;

    double randomRe = unifRe(re);

    return randomRe;

}

std::vector<double> PolyRoots(std::vector<double> coefficients, double eps = 1e-10, int maxiters = 100, int maxtrials = 10){
    if(eps < 0 || maxiters < 0 || maxtrials < 0) throw std::domain_error("Invalid parameters");
    std::vector<double> nule(coefficients.size() - 1);
    for(int i = coefficients.size() - 2; i >= 0; i--){
        double x;
        int briters = 1;
        bool c = false;
        while(!c || (briters < maxiters)){
            x = RandomDouble(-10, 10);
            auto pair = Laguerre(coefficients, i, x);
            x = pair.first;
            c = pair.second;
            briters++;
        }
        if(!c){
            throw std::range_error("Convergence has not achieved");
        }
        nule[i] = x;
        auto v = coefficients[i];
        for(int j = i - 1; j >= 0; j--){
            auto w = coefficients[j];
            coefficients[j] = v;
            v = w + (x * v);
        }
    }
    return nule;
}



int main() {
   std::function<double(double)> f = [](double x){ return x - 5;};

   try{             // testiranje BracketRoot za funkciju y = x - 5
       double a;
       double b;
       std::cout << "Ogradjena nula za funkciju y = x - 5: " << BracketRoot(f, 1, a, b) << std::endl << std::endl;
   }
   catch(...){

   }
   try{             // testiranje izuzetaka za BrackeRoot
       double a;
       double b;
       BracketRoot(f, 1, a, b, -2);
   }
   catch(std::domain_error e){
       std::cout << e.what() << std::endl << std::endl;
   }

   try{

   }
   catch(...){

   }
}
