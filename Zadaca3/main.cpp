#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

class AbstractInterpolator{
protected:
    std::vector<std::pair<double, double>> data;
    mutable int indeks;
    int Locate(double x) const{
        if(x < data[0].first){
            indeks++;
            return 0;
        }
        if(x > data[data.size() - 1].first){
            indeks = data.size();
            return indeks;
        }
        if(data[indeks - 1].first < x && (data[indeks].first > x)){
            return indeks;
        }
        auto it = std::lower_bound(data.begin(), data.end(), x,
                                   [](std::pair<double, double> podatak, double x){
                                       if(podatak.first < x) return true;
                                       return false;
                                   });
        indeks = static_cast<int>(it - data.begin());
        return indeks;
    }
public:
    AbstractInterpolator(const std::vector<std::pair<double, double>> &data) : data(data){
        std::sort(this->data.begin(), this->data.end(),
                  [](std::pair<double, double> x, std::pair<double, double> y){
            if(std::fabs(x.first - y.first) < std::numeric_limits<double>::epsilon()) throw std::domain_error("Invalid data set");
            if(x.first > y.first) return false;
            return true;
        });
        indeks = 0;
    }
    virtual double operator()(double x) const = 0;
};

class LinearInterpolator : public AbstractInterpolator{
public:
    LinearInterpolator(const std::vector<std::pair<double, double>> &data) : AbstractInterpolator(data){}
    double operator() (double x) const override{
        int lokacija = this->Locate(x);
        if(lokacija == 0) lokacija++;
        double y = (((data[lokacija].first - x) / (data[lokacija].first - data[lokacija - 1].first)) * data[lokacija - 1].second)
                + (((x - data[lokacija - 1].first) / (data[lokacija].first - data[lokacija - 1].first)) * data[lokacija].second);
        return y;
    }
};

class PolynomialInterpolator : public AbstractInterpolator{
    std::vector<double> newton_Koef;
    std::vector<double> pomocni;
public:
    PolynomialInterpolator(const std::vector<std::pair<double, double>> &data) : AbstractInterpolator(data){
        newton_Koef.resize(data.size());
        pomocni.resize(data.size());
        for(int i = 0; i < newton_Koef.size(); i++){
            newton_Koef[i] = this->data[i].second;
        }
        for(int j = 1; j < newton_Koef.size(); j++) {
            for(int i = this->data.size() - 1; i >= j; i--) {
                newton_Koef[i] = ((newton_Koef[i] - newton_Koef[i - 1]) / (this->data[i].first - this->data[i - j].first));
            }
        }
    }
    double operator() (double x) const override{
        double result = 0;
        for(int i = 0; i < newton_Koef.size(); i++){
            double jloop = 1;
            for(int j = 0; j < i; j++){
                jloop *= (x - data[j].first);
            }
            result += newton_Koef[i] * jloop;
        }
        return result;
    }
};

class PiecewisePolynomialInterpolator : public AbstractInterpolator{
    int order;
public:
    PiecewisePolynomialInterpolator(const std::vector<std::pair<double, double>> &data, int order): AbstractInterpolator(data){
        if(order < 1) throw std::domain_error("Invalid order");
        if(order >= data.size()) throw std::domain_error("Invalid order");
        this->order = order;
    }
    double operator() (double x) const override{
        int lokacija = this->Locate(x);
        double rangefloor;
        double rangeceiling;
        if(lokacija % 2 != 0){
            rangefloor = lokacija - ((order - 1) / 2.) - 1;
            rangeceiling = lokacija + ((order + 1) / 2.) - 1;
        }
        else{
            rangefloor = lokacija - (order / 2.) - 1;
            rangeceiling = lokacija + (order / 2.) - 1;
        }
        if(rangefloor < 0){
            rangefloor = 0;
            rangeceiling = order + 1;
        }
        if(rangeceiling >= data.size()){
            rangefloor = data.size() - order - 1;
            rangeceiling = data.size() - 1;
        }
        double result = 0;
        for(int i = rangefloor; i <= rangeceiling; i++){
            double jloop = 1;
            for(int j = rangefloor; j <= rangeceiling; j++){
                if(j == i) continue;
                jloop *= ((x - data[j].first) / (data[i].first - data[j].first));
            }
            result += (data[i].second * jloop);
        }
        return result;
    }
};

class SplineInterpolator : public AbstractInterpolator{
    std::vector<double> r;
    std::vector<double> q;
    std::vector<double> s;
public:
    SplineInterpolator(const std::vector<std::pair<double, double>> &data) : AbstractInterpolator(data){
        r.resize(data.size());
        q.resize(data.size() - 1);
        s.resize(data.size() - 1);
        std::vector<double> alpha(data.size());

        r[0] = 0;
        r[r.size() - 1] = 0;
        for(int i = 1; i < r.size() - 1; i++){
            alpha[i] = 2 * (data[i + 1].first - data[i - 1].first);
            r[i] = 3 * (((data[i + 1].second - data[i].second) / (data[i + 1].first - data[i].first)) - ((data[i].second - data[i - 1].second) / (data[i].first - data[i - 1].first)));
        }
        for(int i = 1; i < r.size() - 2; i++){
            double mi = (data[i + 1].first - data[i].first) / alpha[i];
            alpha[i + 1] = alpha[i + 1] - (mi * (data[i + 1].first - data[i].first));
            r[i + 1] = r[i + 1] - (mi * r[i]);
        }
        r[r.size() - 2] = r[r.size() - 2] / alpha[alpha.size() - 2];
        for(int i = r.size() - 3; i >= 1; i--){
            r[i] = (r[i] - ((data[i + 1].first - data[i].first) * r[i + 1])) / alpha[i];
        }

        for(int i = 0; i < s.size(); i++){
            double delta = data[i + 1].first - data[i].first;
            s[i] = (r[i + 1] - r[i]) / (3 * delta);
            q[i] = ((data[i + 1].second - data[i].second) / delta) - ((delta * (r[i + 1] + (2 * r[i]))) / 3);
        }
    }
    double operator() (double x) const override{
        int lokacija = this->Locate(x) - 1;
        double t = x - data[lokacija].first;
        double result = data[lokacija].second + (t * (q[lokacija] + (t * (r[lokacija] + (t * s[lokacija])))));
        return result;
    }
};

double min(double a, double b){
    if(a < b) return a;
    return b;
}

double max(double a, double b){
    if(a > b) return a;
    return b;
}

class BarycentricInterpolator : public AbstractInterpolator{
    std::vector<double> weights;
public:
    BarycentricInterpolator(const std::vector<std::pair<double, double>> &data, int order) : AbstractInterpolator(data){
        if(order < 0 || order > data.size()) throw std::domain_error("Invalid order");
        weights.resize(data.size());
        for(int i = 0; i < data.size(); i++){
            weights[i] = 0;
            for(int k = max(0, i - order); k <= min(i, weights.size() - order); k++){
                double p = 1;
                for(int j = k; j <= k + order; j++){
                    if(j != i){
                        p = p / (data[i].first - data[j].first);
                    }
                }
                if(k % 2 != 0){
                    p = -p;
                }
                weights[i] += p;
            }

        }
    }
    double operator () (double x) const override{
        double p = 0;
        double q = 0;
        for(int i = 0; i < data.size(); i++){
            if(x == data[i].first) return data[i].second;
            double u = weights[i] / (x - data[i].first);
            p = p + (u * data[i].second);
            q = q + u;
        }
        return p / q;
    }
    std::vector<double> GetWeights(){
        return weights;
    }
};

class TrigonometricInterpolator : public AbstractInterpolator{
public:
    TrigonometricInterpolator(const std::vector<std::pair<double, double>> &data) : AbstractInterpolator(data){
        if(this->data[0].second != this->data[data.size() - 1].second) throw std::domain_error("Function is not periodic");
    }
    double operator () (double x) const override{
        double result = 0;
        double omega = (2 * M_PI) / (data[data.size() - 1].first - data[0].first);
        if(data.size() % 2 == 0){
            for(int k = 0; k < data.size() - 1; k++){
                double jloop = 1;
                for(int j = 0; j < data.size() - 1; j++){
                    if(j != k){
                        jloop *= (std::sin(omega / 2) * (x - data[j].first)) / (std::sin(omega/2) * (data[k].first - data[j].first));
                    }
                }
                result += data[k].second * jloop;
            }
        }
        else{
            std::vector<double> alpha(data.size() - 1);
            for(int i = 0; i < alpha.size(); i++){
                alpha[i] = 0;
                for(int j = 0; j < alpha.size(); j++){
                    if(i != j){
                        alpha[i] += data[j].first;
                    }
                }
                alpha[i] *= -1;
            }

            for(int k = 0; k < data.size() - 1; k++){
                double jloop = 1;
                for(int j = 0; j < data.size() - 1; j++){
                    if(j != k){
                        jloop *= (std::sin(omega / 2) * (x - data[j].first)) / (std::sin(omega / 2) * (data[k].first - data[j].first));
                    }
                }
                result += data[k].second * jloop * ((std::sin(omega / 2) * (x - alpha[k])) / (std::sin(omega / 2) * (data[k].first - alpha[k])));
            }
        }
        return result;
    }
};

void LinearIntTest1(){          // testiranje LinearInterpolatora za sortirane i nesortirane vrijednosti
    std::cout << "LINEAR INTERPOLATOR TEST 1: " << std::endl;
    std::vector<std::pair<double, double>> data({{1, 1}, {3, 4}, {6, 6}, {8, 5}});
    LinearInterpolator lin1(data);
    std::vector<std::pair<double, double>> dataunsorted({{6, 6}, {8, 5}, {3, 4}, {1, 1}});
    LinearInterpolator lin2(data);
    for(int i = 1; i < 5; i++){
        std::cout << lin1(i * 1.5) << "  -  " << lin2(i * 1.5) << std::endl;
    }
    std::cout << std::endl << std::endl;
}

void LinearIntTest2(){          // testiranje izuzetaka za kada se unesu dvije jednake vrijednosti x[n]
    std::cout << "LINEAR INTERPOLATOR TEST 2: " << std::endl;
    std::vector<std::pair<double, double>> data({{2, 3}, {6, 4}, {8, 3}, {2, 4}});
    try{
        LinearInterpolator lin1(data);
    }
    catch(std::domain_error &e){
        std::cout << e.what();
    }
    std::cout << std::endl << std::endl;
}

void LinearIntTest3(){         // testiranje linearne interpolacija za vrijednosti izvan opsega podataka
    std::cout << "LINEAR INTERPOLATOR TEST 3: " << std::endl;
    std::vector<std::pair<double, double>> data({{6, 6}, {8, 5}, {3, 4}, {1, 1}});
    LinearInterpolator lin(data);
    std::cout << lin(0.5) << " " << lin(9) << std::endl << std::endl;
}

void PolyIntTest1(){          // testiranje polinomijalne interpolacije
    std::cout << "POLYNOMIAL INTERPOLATOR TEST 1: " << std::endl;
    std::vector<std::pair<double, double>> data({{1, 1}, {8, 5}, {6, 6}, {3, 4}});
    PolynomialInterpolator poli(data);
    std::cout << poli(2) << " " << poli(7.643) << " " << poli(3) << std::endl << std::endl;
}

void PolyIntTest2(){           // testiranje polinomijalne interpolacije na podacima sa jednakim x vrijednostima
    std::cout << "POLYNOMIAL INTERPOLATOR TEST 2: " << std::endl;
    std::vector<std::pair<double, double>> data({{1, 1}, {3, 4}, {6, 6}, {1, 6}});
    try{
        PolynomialInterpolator poli(data);
    }
    catch(std::domain_error &e){
        std::cout << e.what() << std::endl << std::endl;
    }
}

void PolyIntTest3(){            // testiranje polinomijalne interpolacije za vrijednosti polinoma x^2 + 3x + 6
    std::cout << "POLYNOMIAL INTERPOLATOR TEST 3: " << std::endl;
    std::vector<std::pair<double, double>> data({{0, 6}, {-1.5, 3.75}, {-4, 9.98}, {1, 9.99}});
    PolynomialInterpolator poli(data);
    std::cout << poli(0) << " " << poli(-3) << " " << poli(0.75) << " " << poli(1) << std::endl << std::endl;
}

void PolyIntTest4(){            // testiranje unosenja vrijednosti van opsega podataka za polinom x^2 + 3x + 6 i uporedjivanje sa vrijednostima na rubovima podataka
    std::cout << "POLYNOMIAL INTERPOLATOR TEST 4: " << std::endl;
    std::vector<std::pair<double, double>> data({{0, 6}, {-1.5, 3.75}, {-4, 9.98}, {1, 9.99}});
    PolynomialInterpolator poli(data);
    std::cout << poli(-5) << " " << poli(-4) << " " <<  poli(2) << " " <<  poli(1) << std::endl << std::endl;
}

void PieceIntTest1(){           // testiranje interpolacije po dijelovima razlicitih redova za polinom 2x^3 + x + 2
    std::cout << "PIECEWISE INTERPOLATION TEST 1: " << std::endl;
    std::vector<std::pair<double, double>> data({{0, 2}, {-0.83512, 0}, {1.11, 5.86}, {-1.05, -1.4}, {-0.51, 1.23}});
    PiecewisePolynomialInterpolator piece(data, 1);
    std::cout << piece(1) << " " << piece(-1) << " " << piece(0) << " " << piece(-0.3) << std::endl;
    PiecewisePolynomialInterpolator piece2(data, 2);
    std::cout << piece(1) << " " << piece(-1) << " " << piece(0) << " " << piece(-0.3) << std::endl;
    PiecewisePolynomialInterpolator piece3(data, 3);
    std::cout << piece(1) << " " << piece(-1) << " " << piece(0) << " " << piece(-0.3) << std::endl << std::endl;
}

void PieceIntTest2(){           // testiranje izuzetaka za pogresne vrijednosti reda
    std::cout << "PIECEWISE INTERPOLATION TEST 2: " << std::endl;
    std::vector<std::pair<double, double>> data({{0, 2}, {-0.83512, 0}, {1.11, 5.86}, {-1.05, -1.4}, {-0.51, 1.23}});
    try{
        PiecewisePolynomialInterpolator piece(data, -1);
    }
    catch(std::domain_error &e){
        std::cout << e.what() << std::endl;
    }
    try{
        PiecewisePolynomialInterpolator piece2(data, 10);
    }
    catch(std::domain_error &e){
        std::cout << e.what() << std::endl << std::endl;
    }
}

void PieceIntTest3(){           // testiranje vrijednosti van opsega podataka
    std::cout << "PIECEWISE INTERPOLATION TEST 3: " << std::endl;
    std::vector<std::pair<double, double>> data({{0, 2}, {-0.83512, 0}, {1.11, 5.86}, {-1.05, -1.4}, {-0.51, 1.23}});
    PiecewisePolynomialInterpolator piece(data, 2);
    std::cout << piece(-2) << " " << piece(3) << " " << piece(-5) << " " << piece(4) << std::endl << std::endl;
}

void SplineTest1(){             // testiranje spline interpolatora za polinom 4x^5 + 2x^3 + x^2 + 2
    std::cout << "SPLINE INTERPOLATION TEST 1: " << std::endl;
    std::vector<std::pair<double, double>> data({{0, 2}, {-0.26868, 2.0278}, {-0.82745, 0}, {0.87, 6.1}, {-0.7, 1.13}, {-0.94, -1.72}});
    SplineInterpolator spline(data);
    std::cout << spline(1) << " " << spline(0) << " " << spline(0.5) << " " << spline(-0.6) << " " << spline(-0.7) << std::endl << std::endl;
}

void SplineTest2(){             // testiranje bacanja izuzetka za nepravilne podatke
    std::cout << "SPLINE INTERPOLATOR TEST 2: " << std::endl;
    std::vector<std::pair<double, double>> data({{0, 2}, {-0.26868, 2.0278}, {-0.82745, 0}, {0.87, 6.1}, {-0.7, 1.13}, {-0.94, -1.72}, {0, 5}});
    try{
        SplineInterpolator spline(data);
    }
    catch(std::domain_error &e){
        std::cout << e.what() << std::endl << std::endl;
    }
}

void BarycentricTest1(){            // testiranje baricentricne interpolacije za polinom 6x^3 + 11x^2
    std::cout << "BARYCENTRIC INTERPOLATOR TEST 1: " << std::endl;
    std::vector<std::pair<double, double>> data({{0, 0}, {-8.722, 5.4774}, {-1.833, 0}, {0.72, 8.04}, {-1.9, -1.4}, {-1.73, 1.9}});
    BarycentricInterpolator bar(data, 1);
    std::cout << bar(-1.5) << " " << bar(0) << " " << bar(0.63) << " " << bar(-1.9) << " " << bar(-6) << " " << bar(-7.532) << std::endl;
    BarycentricInterpolator bar2(data, 3);
    std::cout << bar2(-1.5) << " " << bar2(0) << " " << bar2(0.63) << " " << bar2(-1.9) << " " << bar2(-6) << " " << bar2(-7.532) << std::endl;
    BarycentricInterpolator bar3(data, 4);
    std::cout << bar3(-1.5) << " " << bar3(0) << " " << bar3(0.63) << " " << bar3(-1.9) << " " << bar3(-6) << " " << bar3(-7.532) << std::endl << std::endl;
}

void BarycentricTest2(){        // testiranje izuzetaka za pogresne vrijednosti reda
    std::cout << "BARYCENTRIC INTERPOLATOR TEST 2: " << std::endl;
    std::vector<std::pair<double, double>> data({{0, 0}, {-8.722, 5.4774}, {-1.833, 0}, {0.72, 8.04}, {-1.9, -1.4}, {-1.73, 1.9}});
    try{
        BarycentricInterpolator bar(data, -5);
    }
    catch(std::domain_error &e){
        std::cout << e.what() << std::endl;
    }
    try{
        BarycentricInterpolator bar(data, 10);
    }
    catch(std::domain_error &e){
        std::cout << e.what() << std::endl << std::endl;
    }
}

void BarycentricTest3(){        // testiranje metode GetWeights
    std::cout << "BARYCENTRIC INTERPOLATOR TEST 3: " << std::endl;
    std::vector<std::pair<double, double>> data({{0, 0}, {-8.722, 5.4774}, {-1.833, 0}, {0.72, 8.04}, {-1.9, -1.4}, {-1.73, 1.9}});
    BarycentricInterpolator bar(data, 2);
    std::vector<double> weights = bar.GetWeights();
    for(auto x : weights) std::cout << x << " ";
    std::cout << std::endl << std::endl;
}

void TrigIntTest1(){            // testiranje trigonometrijske interpolacije za vrijednosti sinx
    std::cout << "TRIGONOMETRIC INTERPOLATOR TEST 1: " << std::endl;
    std::vector<std::pair<double, double>> data({{M_PI, 1}, {0, 0}, {0.92, 0.79}, {3.53, -0.38}, {2 * M_PI, 0}});
    TrigonometricInterpolator trig(data);
    std::cout << trig(1) << " " << trig(0) << " " << trig(0.42) << " " << trig(3.31) << " " << trig(0.97) << std::endl << std::endl;
}

void TrigIntTest2(){            // testiranja izuzetka kada se razlikuju y vrijednosti u prvom i posljednjem cvoru
    std::cout << "TRIGONOMETRIC INTERPOLATOR TEST 2: " << std::endl;
    std::vector<std::pair<double, double>> data({{M_PI, 1}, {0, 0}, {0.92, 0.79}, {3.53, -0.38}});
    try{
        TrigonometricInterpolator trig(data);
    }
    catch(std::domain_error &e){
        std::cout << e.what();
    }
}

int main() {
    LinearIntTest1();
    LinearIntTest2();
    LinearIntTest3();

    PolyIntTest1();
    PolyIntTest2();
    PolyIntTest3();
    PolyIntTest4();

    PieceIntTest1();
    PieceIntTest2();
    PieceIntTest3();

    SplineTest1();
    SplineTest2();

    BarycentricTest1();
    BarycentricTest2();
    BarycentricTest3();

    TrigIntTest1();
    TrigIntTest2();
}
