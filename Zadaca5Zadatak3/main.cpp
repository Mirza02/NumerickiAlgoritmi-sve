#include <iostream>
#include <limits>
#include <utility>
#include <vector>
#include <algorithm>
#include <cmath>

template <typename FunType>
double RK4Step(FunType f, double x, double y, double h){
    double k1 = f(x, y);
    double k2 = f(x + h / 2, y + h * k1 / 2);
    double k3 = f(x + h / 2, y + h * k2 / 2);
    double k4 = f(x + h, y + h * k3);
    return y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
}

template <typename FunType>
std::vector<std::pair<double, double>> RK4Integrator(FunType f, double x0, double y0,
                                                     double xmax, double h, double eps = 1e-8,
                                                     bool adaptive = false){
    std::vector<std::pair<double, double>> result;
    if((h > 0 && xmax <= x0) || (h < 0 && xmax >= x0)){
        result.push_back({x0, y0});
        return result;
    }
    if(!adaptive){
        double x = x0;
        double y = y0;
        double delta = h / 1000;
        if(h < 0){
            while(x >= (xmax + delta)){
                result.push_back({x, y});
                double k1 = f(x, y);
                double k2 = f(x + h / 2, y + h * k1 / 2);
                double k3 = f(x + h / 2, y + h * k2 / 2);
                double k4 = f(x + h, y + h * k3);
                y = y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
                x = x + h;
            }
        }
        else{
            while(x <= (xmax + delta)){
                result.push_back({x, y});
                double k1 = f(x, y);
                double k2 = f(x + h / 2, y + h * k1 / 2);
                double k3 = f(x + h / 2, y + h * k2 / 2);
                double k4 = f(x + h, y + h * k3);
                y = y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
                x = x + h;
            }
        }

    }
    else{
        double x = x0;
        double y = y0;
        result.push_back({x, y});
        if(h < 0){
            while(x >= xmax){
                double u = RK4Step(f, x, y, h / 2);
                double v = RK4Step(f, x + h / 2, u, h / 2);
                double w = RK4Step(f, x, y, h);
                double delta = std::fabs(w - v) / std::fabs(h);
                if(delta <= std::numeric_limits<double>::epsilon()){
                    x = x + h;
                    y = v;
                    result.push_back({x, y});
                }
                h = h * std::min(5., 0.9 * std::pow((std::numeric_limits<double>::epsilon() / delta), 1. / 4.));
            }
        }
        else{
            while(x <= xmax){
                double u = RK4Step(f, x, y, h / 2);
                double v = RK4Step(f, x + h / 2, u, h / 2);
                double w = RK4Step(f, x, y, h);
                double delta = std::fabs(w - v) / std::fabs(h);
                if(delta <= std::numeric_limits<double>::epsilon()){
                    x = x + h;
                    y = v;
                    result.push_back({x, y});
                }
                double a = 0.9 * std::pow((std::numeric_limits<double>::epsilon() / delta), 1. / 4.);
                h = h * std::min(5., 0.9 * std::pow((std::numeric_limits<double>::epsilon() / delta), 1. / 4.));
            }
        }

    }
    return result;
}

/*template <typename FunType>
std::vector<std::pair<double, std::vector<double>>> RK4SystemIntegrator(
        FunType f, double x0, std::vector<double> y0, double xmax, double h
        ){
    std::vector<double> test = FunType(x0, y0);
    if(y0.size() != test.size()){
        throw std::range_error("Incompatible formats");
    }
    std::vector<std::pair<double, std::vector<double>>> result;
    if((xmax < x0 && h > 0) || (xmax > x0 && h < 0)){
        result.push_back({x0, y0});
        return result;
    }
    double x = x0;
    std::vector<double> y(y0.size());
    for(int k = 0; k < y0.size(); k++){
        y[k] = y0[k];
    }
    if(h > 0){
        while(x >= xmax){
            result.push_back({x, y});
            std::vector<std::vector<double>> K1(y.size(), std::vector<double>(y.size()));
            std::vector<std::vector<double>> K2(y.size(), std::vector<double>(y.size()));
            std::vector<std::vector<double>> K3(y.size(), std::vector<double>(y.size()));
            std::vector<std::vector<double>> K4(y.size(), std::vector<double>(y.size()));
            std::vector<double> u(y.size());
            std::vector<double> v(y.size());
            for(int k = 0; k < y.size(); k++){
                K1[k] = f(x, y);
            }
        }
    }
}*/



int main() {
 //   std::function<double(double, double)> f = [] (double x, double y){ return 2 * x + y; };
 //   try{
 //       double
 //   }
}
