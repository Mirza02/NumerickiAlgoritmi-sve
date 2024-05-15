#include <iostream>
#include <cmath>
#include <limits>

template <typename FunType>
double FindMinimum(FunType f, double x0, double eps = 1e-8, double hinit = 1e-5, double hmax = 1e10, double lambda = 1.4){
    if(eps < 0 || hinit < 0 || hmax < 0 || lambda < 0){
        throw std::domain_error("Invalid parameters");
    }
    double a, b, c, d;
    while(true){
        if(std::fabs(hinit) >= hmax) throw std::logic_error("Minimum has not found");
        if((f(x0 + std::fabs(hinit)) < f(x0))){
            hinit *= lambda;
            x0 += hinit;
        }
        else if(f(x0 - std::fabs(hinit)) < f(x0)){
            hinit *= -1;
            hinit *= lambda;
            x0 += hinit;
        }
        else{
            break;
        }
    }
    a = x0 - std::fabs(hinit);
    b = x0 + std::fabs(hinit);
    c = x0;
    double k = (1 + std::sqrt(5)) / 2;
    if(std::fabs(c - a) < std::fabs(b - c)){
        d = b - (b - c) / k;
    }
    else{
        d = c;
        c = a + (c - a) / k;
    }
    double u = f(c);
    double v = f(d);
    while(std::fabs(b - a) > eps){
        if(u < v){
            b = d;
            d = c;
            c = a + (c - a) / k;
            v = u;
            u = f(c);
        }
        else{
            a = c;
            c = d;
            d = b - (b - d) / k;
            u = v;
            v = f(d);
        }
    }
    return (a + b) / 2;
}

int main() {
    try{
        std::function<double(double)> f = [](double x){ return (5 * x * x + 2 * x); };
        std::cout << "Minimum funkcije 5x^2 + 2x : " << FindMinimum(f, -1) << std::endl << std::endl;
    }
    catch(...){

    }
    try{
        std::function<double(double)> f = [] (double x){ return 3 * x * x * x + 5 * x * x + 2 * x; };
        std::cout << "Minimum funkcije 3x^3 + 5x^2 + 2x : " << FindMinimum(f, -0.5);
    }
    catch(...){

    }
    try{
        std::function<double(double)> f = [] (double x) { return x; };
        FindMinimum(f, 1, -2);
    }
    catch(std::domain_error &e){
        std::cout << e.what() << std::endl << std::endl;
    }

}
