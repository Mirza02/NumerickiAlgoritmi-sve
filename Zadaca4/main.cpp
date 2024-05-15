#include <iostream>
#include <vector>
#include <cmath>

const double PI = 4 * atan(1.);

class ChebyshevApproximation{
    std::vector<double> koef;
    double m;
    double xmin, xmax;
    ChebyshevApproximation(std::vector<double> c, double xmin, double xmax, double m){
        this->koef.resize(m);
        this->xmin = xmin;
        this->xmax = xmax;
        this->m = m;
        double mi = 4 / (xmax - xmin);
        koef[m - 1] = mi * m * c[m];
        koef[m - 2] = mi * (m - 1) * c[m - 1];
        for(int k = m - 3; k >= 0; k--){
            koef[k] = koef[k + 2] + mi * (k + 1) * c[k + 1];
        }
    }
public:
    template <typename FunType>
    ChebyshevApproximation(FunType f, double xmin, double xmax, int n){
        if(xmin > xmax || n < 1) throw std::domain_error("Bad parameters");
        this->m = n;
        this->xmin = xmin;
        this->xmax = xmax;
        std::vector<double> v(n + 1);
        koef.resize(n + 1);
        for(int i = 0; i <= n; i++){
            v[i] = f((xmin + xmax + (xmax - xmin) * cos(M_PI * (2 * i + 1) / (2 * (n) + 2))) / 2);
        }
        for(int k = 0; k <= n; k++){
            double s = 0;
            for(int i = 0; i <= n; i++){
                s = s + v[i] * cos(k * M_PI * (2 * i + 1) / (2 * (n) + 2));
            }
            koef[k] = 2 * s / ((n) + 1);
        }

    }
    void set_m(int m){
        if(m <= 1 || m >= koef.size()) throw std::domain_error("Bad order");
        this->m = m;
    }
    void trunc(double eps){
        int tempM = m;
        while(tempM > 1 && std::fabs(koef[tempM - 1]) < eps){
            tempM--;
        }
        if(tempM < 1){
            throw std::domain_error("Bad tolerance");
        }
        m = tempM;
    }
    double operator()(double x) const{
        if(x < xmin || x > xmax) throw std::domain_error("Bad argument");
        double t = (2 * x - xmin - xmax) / (xmax - xmin);
        double p = 1;
        double q = t;
        double s = koef[0] / 2 + koef[1] * t;
        for(int k = 2; k <= m; k++){
            double r = 2 * t * q - p;
            s = s + koef[k] * r;
            p = q;
            q = r;
        }
        return s;
    }
    double derivative(double x) const{
        if(x < xmin || x > xmax) throw std::domain_error("Bad argument");
        double t = (2 * x - xmin - xmax) / (xmax - xmin);
        double p = 1;
        double q = 4 * t;
        double s = koef[1] + q * koef[2];
        for(int k = 3; k < m; k++){
            double r = k * (2 * t * q / (k - 1) - p / (k - 2));
            s = s + koef[k] * r;
            p = q;
            q = r;
        }
        return 2 * s / (xmax - xmin);
    }
    ChebyshevApproximation derivative() const{
        ChebyshevApproximation chebDeriv(koef, xmin, xmax, m);
        return chebDeriv;
    }
    ChebyshevApproximation antiderivative() const{

    }
    //---------------------------------------------------------------POPRAVITI!!!!
    double integrate(double a, double b) const{
        if(a < b){
            if(a < xmin || b > xmax) throw std::domain_error("Bad interval");
        }
        else{
            if(b < xmin || a > xmax) throw std::domain_error("Bad interval");
        }
        std::vector<double> koefzvijezda(m + 2);
        koefzvijezda[0] = 0;
        for(int k = 1; k < m + 2; k++){
            if(k == m || k == (m + 1)){
                koefzvijezda[k] = ((xmax - xmin) / (4 * k)) * (koef[k - 1] - 0);
            }
            else{
                koefzvijezda[k] = ((xmax - xmin) / (4 * k)) * (koef[k - 1] - 0);
            }
        }

        double t = (2 * a - xmin - xmax) / (xmax - xmin);
        double p = 1;
        double q = t;
        double fA = koefzvijezda[0] / 2 + koefzvijezda[1] * t;
        for(int k = 2; k < m + 2; k++){
            double r = 2 * t * q - p;
            fA = fA + koefzvijezda[k] * r;
            p = q;
            q = r;
        }

        t = (2 * b - xmin - xmax) / (xmax - xmin);
        p = 1;
        q = t;
        double fB = koefzvijezda[0] / 2 + koefzvijezda[1] * t;
        for(int k = 2; k < m + 2; k++){
            double r = 2 * t * q - p;
            fB = fB + koefzvijezda[k] * r;
            p = q;
            q = r;
        }

        return fB - fA;
    }
    double integrate() const{
        double loopresult = 0;
        for(int k = 1; k <= std::fabs((m - 1) / 2); k++){
            loopresult += ((2 * koef[2 * k]) / (1 - (4 * (k * k))));
        }
        double result = ((xmax - xmin) / 2) * koef[0] + ((xmax - xmin) / 2) * loopresult;
        return result;
    }
};

int main() {
    auto x = [](double x){return sin(x);};
    std::cout << "Razne operacije za funkciju sinx na intervalu [0, PI]: " << std::endl;
    auto cheb = ChebyshevApproximation(x, 0, PI, 10);
    std::cout << "U tacki 1: " << cheb(1) << std::endl;
    std::cout << "Izvod u tacki 1: " << cheb.derivative(1) << std::endl;
    auto izvodcheb = cheb.derivative().derivative();
    std::cout << "Drugi izvod u tacki 1: " << izvodcheb(1) << std::endl;
    cheb.set_m(4);
    std::cout << "Aproksimacija u tacki 1 za m = 4 : " << cheb(1) << std::endl;
    std::cout << "Integral funkcije u intervalu [0, PI]: " << cheb.integrate() << std::endl << std::endl;
    std::cout << "Testiranje izuzetaka metoda: " << std::endl;
    try{
        auto test = ChebyshevApproximation(x, 5, 2, 2);
    }
    catch(std::domain_error &e){
        std::cout << e.what() << std::endl;
    }
    try{
        auto test = ChebyshevApproximation(x, 1, 2, -3);
    }
    catch(std::domain_error &e){
        std::cout << e.what() << std::endl;
    }
    try{
        cheb(-4);
    }
    catch(std::domain_error &e){
        std::cout << e.what() << std::endl;
    }
    try{
        cheb(10);
    }
    catch(std::domain_error &e){
        std::cout << e.what() << std::endl;
    }
    try{
        cheb.set_m(15);
    }
    catch(std::domain_error &e){
        std::cout << e.what() << std::endl;
    }
    try{
        cheb.set_m(-4);
    }
    catch(std::domain_error &e){
        std::cout << e.what() << std::endl;
    }
    try{
        cheb.trunc(-4);
    }
    catch(std::domain_error &e){
        std::cout << e.what() << std::endl;
    }
    try{
        cheb.integrate(-5, 2);
    }
    catch(std::domain_error &e){
        std::cout << e.what() << std::endl;
    }
    try{
        cheb.integrate(1, 21);
    }
    catch(std::domain_error &e){
        std::cout << e.what() << std::endl;
    }
}
