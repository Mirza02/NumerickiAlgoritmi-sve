#include <iostream>
#include <vector>
#include <stdexcept>
#include <initializer_list>
#include <cmath>
#include <limits>
#include <iomanip>

class Vector{
    std::vector<double> vektor;
public:
    explicit Vector(int n){
        if(n <= 0) throw std::range_error("Bad dimension");
        vektor.resize(n, 0);
    }
    Vector(std::initializer_list<double> l){
        if(l.size() == 0 ) throw std::range_error("Bad dimension");
        vektor.resize(l.size());
        int i = 0;
        for(auto x : l){
            vektor[i] = x;
            i++;
        }
    }
    int NElems() const{
        return vektor.size();
    }
    double &operator[](int i){
        return vektor[i];
    }
    double operator[](int i) const{
        return vektor[i];
    }
    double &operator()(int i){
        if(i < 1 || i > vektor.size()) throw std::range_error("Invalid index");
        return vektor[i - 1];
    }
    double operator()(int i) const{
        if(i < 1 || i > vektor.size()) throw std::range_error("Invalid index");
        return vektor[i - 1];
    }
    double Norm() const{
        double suma = 0;
        for(int i = 0; i < vektor.size(); i++){
            suma += (vektor[i] * vektor[i]);
        }
        return sqrt(suma);
    }
    friend double VectorNorm(const Vector &v);
    double GetEpsilon() const{
        return 10 * this->Norm() * std::numeric_limits<double>::epsilon();
    }
    void Print(char seperator = '\n', double eps = -1) const{
        if(eps < 0) eps = this->GetEpsilon();
        for(int i = 0; i < vektor.size() - 1; i++){
            if(std::abs(vektor[i]) < eps) std::cout << 0 << seperator;
            else std::cout << vektor[i] << seperator;
        }
        if(std::abs(vektor[vektor.size() - 1]) < eps) std::cout << 0;
        else std::cout << vektor[vektor.size() - 1];
        if(seperator == '\n') std::cout << seperator;
    }
    friend void PrintVector(const Vector &v, char seperator, double eps);
    friend Vector operator +(const Vector &v1, const Vector &v2);
    Vector &operator +=(const Vector &v){
        if(this->NElems() != v.NElems()) throw std::domain_error("Incompatible formats");
        *this = *this + v;
        return *this;
    }
    friend Vector operator -(const Vector &v1, const Vector &v2);
    Vector &operator -=(const Vector &v){
        if(this->NElems() != v.NElems()) throw std::domain_error("Incompatible formats");
        *this = *this - v;
        return *this;
    }
    friend Vector operator *(double s, const Vector &v);
    friend Vector operator *(const Vector &v, double s);
    Vector &operator *=(double s){
        *this = s * *this;
        return *this;
    }
    friend double operator *(const Vector &v1, const Vector &v2);
    friend Vector operator /(const Vector &v, double s);
    Vector &operator /=(double s){
        *this = *this / s;
        return *this;
    }
    void Chop(double eps = -1){
        if(eps < 0) eps = this->GetEpsilon();
        for(int i = 0; i < vektor.size(); i++){
            if(std::abs(vektor[i]) <= eps) vektor[i] = 0;
        }
    }
    bool EqualTo(const Vector &v, double eps = -1) const{
        if(eps < 0) eps = this->GetEpsilon();
        if(vektor.size() != v.NElems()) return false;
        for(int i = 0; i < vektor.size(); i++){
            if(std::abs(vektor[i] - v[i]) >= eps) return false;
        }
        return true;
    }
};

double VectorNorm(const Vector &v){
    double suma = 0;
    for(int i = 0; i < v.NElems(); i++){
        suma += (v[i] * v[i]);
    }
    return sqrt(suma);
}

void PrintVector(const Vector &v, char seperator = '\n', double eps = -1){
    v.Print(seperator, eps);
}

Vector operator+(const Vector &v1, const Vector &v2){
    if(v1.NElems() != v2.NElems()) throw std::domain_error("Incompatible formats");
    Vector result(v1.NElems());
    for(int i = 0; i < result.NElems(); i++){
        result[i] = v1[i] + v2[i];
    }
    return result;
}

Vector operator -(const Vector &v1, const Vector &v2){
    if(v1.NElems() != v2.NElems()) throw std::domain_error("Incompatible formats");
    Vector result(v1.NElems());
    for(int i = 0; i < result.NElems(); i++){
        result[i] = v1[i] - v2[i];
    }
    return result;
}

Vector operator *(double s, const Vector &v){
    Vector result(v.NElems());
    for(int i = 0; i < result.NElems(); i++){
        result[i] = s * v[i];
    }
    return result;
}

Vector operator *(const Vector &v, double s){
    return s * v;
}

double operator *(const Vector &v1, const Vector &v2){
    if(v1.NElems() != v2.NElems()) throw std::domain_error("Incompatible formats");
    double result = 0;
    for(int i = 0; i < v1.NElems(); i++){
        result += v1[i] * v2[i];
    }
    return result;
}

Vector operator /(const Vector &v, double s){
    if(s == 0) throw std::domain_error("Division by zero");
    Vector result(v.NElems());
    for(int i = 0; i < result.NElems(); i++){
        result[i] = v[i] / s;
    }
    return result;
}


class Matrix{
    std::vector<std::vector<double>> matrica;
public:
    Matrix(int m, int n){
        if(m <= 0 || n <= 0) throw std::range_error("Bad dimension");
        matrica.resize(m);
        for(int i = 0; i < m; i++){
            matrica[i].resize(n);
        }
    }
    Matrix(const Vector &v){
        matrica.resize(1);
        matrica[0].resize(v.NElems());
        for(int j = 0; j < v.NElems(); j++){
            matrica[0][j] = v[j];
        }
    }
    Matrix(std::initializer_list<std::vector<double>> l){
        if(l.size() <= 0 || l.begin()->size() <= 0) throw std::range_error("Bad dimension");
        int size = l.begin()->size();
        for(auto x : l){
            if(x.size() != size) throw std::logic_error("Bad matrix");
        }
        matrica.resize(l.size());
        int i = 0;
        for(auto x : l){
            matrica[i].resize(size);
            for(int j = 0; j < size; j++){
                matrica[i][j] = x[j];
            }
            i++;
        }
    }
    int NRows() const{
        return matrica.size();
    }
    int NCols() const{
        return matrica[0].size();
    }
    double *operator[](int i){
        if(i < 0 || i > matrica.size()) throw std::range_error("Invalid index");
        return &matrica[i][0];
    }
    const double *operator[](int i) const{
        if(i < 0 || i > matrica.size()) throw std::range_error("Invalid index");
        return &matrica[i][0];
    }
    double operator()(int i, int j) const{
        if(i < 1 || i >= matrica.size() || j < 1 || j >= matrica[0].size()) throw std::range_error("Invalid index");
        return matrica[i - 1][j - 1];
    }
    double &operator()(int i, int j){
        if(i < 1 || i >= matrica.size() || j < 1 || j >= matrica[0].size()) throw std::range_error("Invalid index");
        return matrica[i - 1][j - 1];
    }
    double Norm() const{
        double suma = 0;
        for(int i = 0; i < matrica.size(); i++){
            for(int j = 0; j < matrica[0].size(); j++){
                suma += matrica[i][j] * matrica[i][j];
            }
        }
        return std::sqrt(suma);
    }
    friend double MatrixNorm(const Matrix &m);
    double GetEpsilon() const{
        return (10 * this->Norm() * std::numeric_limits<double>::epsilon());
    }
    void Print(int width = 10, double eps = -1) const{
        if(eps < 0) eps = this->GetEpsilon();
        bool flag = false;
        for(int i = 0; i < matrica.size(); i++){
            for(int j = 0; j < matrica[0].size(); j++){
                if(matrica[i][j] < 0){
                    width = width + 1;
                    flag = true;
                    break;
                }
            }
            if(flag) break;
        }
        for(int i = 0; i < matrica.size(); i++){
            for(int j = 0; j < matrica[0].size(); j++){
                if(std::abs(matrica[i][j]) < eps){
                    std::cout << std::setw(width) << 0;
                }
                else std::cout << std::setw(width) << matrica[i][j];
            }
            std::cout << std::endl;
        }
    }
    friend void PrintMatrix(const Matrix &m, int width, double eps);
    friend Matrix operator +(const Matrix &m1, const Matrix &m2);
    Matrix &operator +=(const Matrix &m){
        *this = *this + m;
        return *this;
    }
    friend Matrix operator -(const Matrix &m1, const Matrix &m2);
    Matrix &operator -=(const Matrix &m){
        *this = *this - m;
        return *this;
    }
    friend Matrix operator *(double s, const Matrix &m);
    friend Matrix operator *(const Matrix &m, double s);
    Matrix &operator *=(double s){
        *this = *this * s;
        return *this;
    }
    friend Matrix operator *(const Matrix &m1, const Matrix &m2);
    Matrix &operator *=(const Matrix &m){
        *this = *this * m;
        return *this;
    }
    friend Vector operator *(const Matrix &m, const Vector &v);
    friend Matrix Transpose(const Matrix &m);
    void Transpose(){
        if(matrica.size() == matrica[0].size()){
            for(int i = 0; i < matrica.size(); i++){
                for(int j = i + 1; j < matrica.size(); j++){
                    double temp = matrica[i][j];
                    matrica[i][j] = matrica[j][i];
                    matrica[j][i] = temp;
                }
            }
        }
        else{
            std::vector<std::vector<double>> result(matrica[0].size(), std::vector<double>(matrica.size()));
            for(int i = 0; i < result.size(); i++){
                for(int j = 0; j < result[i].size(); j++){
                    result[i][j] = matrica[j][i];
                }
            }
            this->matrica = result;
        }
    }
    void Chop(double eps = -1){
        if(eps < 0) eps = this->GetEpsilon();
        for(int i = 0; i < matrica.size(); i++){
            for(int j = 0; j < matrica[0].size(); j++){
                if(std::abs(matrica[i][j]) <= eps) matrica[i][j] = 0;
            }
        }
    }
    bool EqualTo(const Matrix &m, double eps = -1) const{
        if((matrica.size() != m.NRows()) || (matrica[0].size() != m.NCols())) return false;
        if(eps <= 0) eps = this->GetEpsilon();
        for(int i = 0; i < matrica.size(); i++){
            for(int j = 0; j < matrica[0].size(); j++){
                if(std::abs(matrica[i][j] - m[i][j]) >= eps) return false;
            }
        }
        return true;
    }
    friend Matrix LeftDiv(Matrix m1, Matrix m2);
    friend Vector LeftDiv(Matrix m, Vector v);
    friend Matrix operator/(const Matrix &m, double s);
    Matrix &operator /=(double s){
        *this = *this / s;
        return *this;
    }
    friend Matrix operator /(Matrix m1, Matrix m2);
    Matrix &operator /=(Matrix m){
        Matrix x = *this / m;
        *this = x;
        return *this;
    }
    double Det() const{
        if(this->NRows() != this->NCols()) throw std::domain_error("Matrix is not square");
        double d = 1;
        Matrix test = *this;
        for(int k = 0; k < test.NRows(); k++){
            double p = k;
            for(int i = k + 1; i < test.NRows(); i++){
                if(std::fabs(test[i][k]) > test[p][k]) p = i;
            }
            if(std::fabs(test[p][k]) < this->GetEpsilon()) return 0;
            if(p != k) {
                for (int l = k; l < test.NCols(); l++) {
                    double temp;
                    temp = test[k][l];
                    test[k][l] = test[p][l];
                    test[p][l] = temp;
                }
                d *= -1;
            }
            d *= test[k][k];
            for(int i = k + 1; i < test.NRows(); i++){
                double mi = test[i][k] / test[k][k];
                for(int j = k + 1; j < test.NCols(); j++){
                    test[i][j] = test[i][j] - (mi * test[k][j]);
                }
            }
        }
        return d;
    }
    friend double Det(Matrix m);
    friend Matrix Inverse(Matrix m);
    void Invert(){
        *this = Inverse(*this);
    }
    friend Matrix RREF(Matrix m);
    void ReduceToRREF(){
        *this = RREF(*this);
    }
    friend int Rank(Matrix m);
    int Rank() const{
        Matrix r = RREF(*this);
        int rank = 0;
        bool flag = false;
        for(int i = 0; i < r.NCols(); i++){
            flag = false;
            for(int j = 0; j < r.NRows(); j++){
                if(flag){
                    if(r[j][i] == 0) continue;
                    flag = false;
                }
                else{
                    if(r[j][i] == 1) flag = true;
                }
            }
            if(flag) rank++;
        }
        return rank;
    }
};

double MatrixNorm(const Matrix &m){
    return m.Norm();
}

void PrintMatrix(const Matrix &m, int width = 10, double eps = -1){
    m.Print();
}

Matrix operator +(const Matrix &m1, const Matrix &m2){
    if(m1.NCols() != m2.NCols() || m1.NRows() != m2.NRows()) throw std::domain_error("Incompatible formats");
    Matrix result(m1.NRows(), m1.NCols());
    for(int i = 0; i < result.NRows(); i++){
        for(int j = 0; j < result.NCols(); j++){
            result[i][j] = m1[i][j] + m2[i][j];
        }
    }
    return result;
}

Matrix operator -(const Matrix &m1, const Matrix &m2){
    if(m1.NCols() != m2.NCols() || m1.NRows() != m2.NRows()) throw std::domain_error("Incompatible formats");
    Matrix result(m1.NRows(), m1.NCols());
    for(int i = 0; i < result.NRows(); i++){
        for(int j = 0; j < result.NCols(); j++){
            result[i][j] = m1[i][j] - m2[i][j];
        }
    }
    return result;
}

Matrix operator *(double s, const Matrix &m){
    Matrix result(m.NRows(), m.NCols());
    for(int i = 0; i < m.NRows(); i++){
        for(int j = 0; j < m.NCols(); j++){
            result[i][j] = m[i][j] * s;
        }
    }
    return result;
}

Matrix operator *(const Matrix &m, double s){
    Matrix result = s * m;
    return result;
}

Matrix operator *(const Matrix &m1, const Matrix &m2){
    if(m1.NCols() != m2.NRows()) throw std::domain_error("Incompatible formats");
    int m = m2.NCols();
    Matrix result(m1.NRows(), m2.NCols());
    for(int i = 0; i < result.NRows(); i++){
        for(int j = 0; j < result.NCols(); j++){
            double sum = 0;
            for(int k = 0; k < m1.NCols(); k++){
                sum += m1[i][k] * m2[k][j];
            }
            result[i][j] = sum;
        }
    }
    return result;
}

Vector operator*(const Matrix &m, const Vector &v){
    if(m.NCols() != v.NElems()) throw std::domain_error("Incompatible formats");
    Vector result(m.NRows());
    for(int i = 0; i < m.NRows(); i++){
        double suma = 0;
        for(int j = 0; j < m.NCols(); j++){
            suma += m[i][j] * v[j];
        }
        result[i] = suma;
    }
    return result;
}

Matrix Transpose(const Matrix &m){
    Matrix result(m.NCols(), m.NRows());
    for(int i = 0; i < result.NRows(); i++){
        for(int j = 0; j < result.NCols(); j++){
            result[i][j] = m[j][i];
        }
    }
    return result;
}

Matrix LeftDiv(Matrix m1, Matrix m2){
    if(m1.NRows() != m1.NCols()) throw std::domain_error("Divisor matrix is not square");
    if(m1.NRows() != m2.NRows()) throw std::domain_error("Incompatible formats");
    Matrix x(m2.NRows(), m2.NCols());
    for(int k = 0; k < m1.NRows(); k++){
        int p = k;
        for(int i = k + 1; i < m1.NRows(); i++){
            if(std::fabs(m1[i][k]) > std::fabs(m1[p][k])) p = i;
        }
        if(std::fabs(m1[p][k]) < m1.GetEpsilon()){
            throw std::domain_error("Divisor matrix is singular");
        }
        if(p != k){
            for(int l = 0; l < m1.NCols(); l++){
                double temp;
                temp = m1[k][l];
                m1[k][l] = m1[p][l];
                m1[p][l] = temp;
            }
            for(int l = 0; l < m2.NCols(); l++){
                double temp;
                temp = m2[k][l];
                m2[k][l] = m2[p][l];
                m2[p][l] = temp;
            }
        }
        for(int i = k + 1; i < m1.NRows(); i++){
            double mi = m1[i][k] / m1[k][k];
            for(int j = k + 1; j < m1.NRows(); j++){
                m1[i][j] = m1[i][j] - (mi * m1[k][j]);
            }
            for(int j = 0; j < m2.NCols(); j++){
                m2[i][j] = m2[i][j] - (mi * m2[k][j]);
            }
        }
    }

    for(int k = 0; k < m2.NCols(); k++){
        for(int i = m2.NRows() - 1; i >= 0; i--){
            double s = m2[i][k];
            for(int j = i + 1; j < m2.NRows(); j++){
                s = s - (m1[i][j] * x[j][k]);
            }
            x[i][k] = s / m1[i][i];
        }
    }
    return x;
}

Vector LeftDiv(Matrix m, Vector v){
    if(m.NRows() != m.NCols()) throw std::domain_error("Divisor matrix is not square");
    if(m.NRows() != v.NElems()) throw std::domain_error("Incompatible formats");
    Vector x(v.NElems());
    for(int k = 0; k < m.NRows(); k++){
        int p = k;
        for(int i = k + 1; i < m.NRows(); i++){
            if(std::fabs(m[i][k]) > std::fabs(m[p][k])) p = i;
        }
        if(std::fabs(m[p][k]) < m.GetEpsilon()){
            throw std::domain_error("Divisor matrix is singular");
        }
        if(p != k){
            for(int l = k; l < m.NCols(); l++){
                double temp;
                temp = m[k][l];
                m[k][l] = m[p][l];
                m[p][l] = temp;
            }
            double temp = v[k];
            v[k] = v[p];
            v[p] = temp;
        }
        for(int i = k + 1; i < m.NRows(); i++){
            double mi = m[i][k] / m[k][k];
            for(int j = k + 1; j < m.NRows(); j++){
                m[i][j] = m[i][j] - (mi * m[k][j]);
            }
            v[i] = v[i] - (mi * v[k]);
        }
    }

    for(int i = v.NElems() - 1; i >= 0; i--){
        double s = v[i];
        for(int j = i + 1; j < m.NRows(); j++){
            s = s - (m[i][j] * x[j]);
        }
        x[i] = s / m[i][i];
    }
    return x;
}

Matrix operator /(const Matrix &m, double s){
    if(s == 0) throw std::domain_error("Division by zero");
    Matrix result(m.NRows(), m.NCols());
    for(int i = 0; i < m.NRows(); i++){
        for(int j = 0; j < m.NCols(); j++){
            result[i][j] = m[i][j] / s;
        }
    }
    return result;
}

Matrix operator /(Matrix m1, Matrix m2){
    Matrix x(m1.NRows(), m1.NCols());
    Matrix temp = m1;
    m1 = m2;
    m2 = temp;
    if(m1.NCols() != m1.NRows()) throw std::domain_error("Divisor matrix is not square");
    if(m1.NCols() != m2.NCols()) throw std::domain_error("Incompatible formats");
    for(int k = 0; k < m1.NCols(); k++){
        int p = k;
        for(int i = k + 1; i < m1.NCols(); i++){
            if(std::fabs(m1[k][i]) > std::fabs(m1[k][p])) p = i;
        }
        if(std::fabs(m1[k][p]) < m1.GetEpsilon()){
            throw std::domain_error("Divisor matrix is singular");
        }
        if(p != k){
            for(int l = 0; l < m1.NRows(); l++){
                double temp;
                temp = m1[l][k];
                m1[l][k] = m1[l][p];
                m1[l][p] = temp;
            }
            for(int l = 0; l < m2.NRows(); l++){
                double temp;
                temp = m2[l][k];
                m2[l][k] = m2[l][p];
                m2[l][p] = temp;
            }
        }
        for(int i = k + 1; i < m1.NCols(); i++){
            double mi = m1[k][i] / m1[k][k];
            for(int j = k + 1; j < m1.NCols(); j++){
                m1[j][i] = m1[j][i] - (mi * m1[j][k]);
            }
            for(int j = 0; j < m2.NRows(); j++){
                m2[j][i] = m2[j][i] - (mi * m2[j][k]);
            }
        }
    }

    for(int k = 0; k < m2.NRows(); k++){
        for(int i = m2.NCols() - 1; i >= 0; i--){
            double s = m2[k][i];
            for(int j = i + 1; j < m2.NCols(); j++){
                s = s - (m1[j][i] * x[k][j]);
            }
            x[k][i] = s / m1[i][i];
        }
    }
    return x;
}

double Det(Matrix m){
    return m.Det();
}

Matrix Inverse(Matrix m){
    if (m.NRows() != m.NCols()) throw std::domain_error("Matrix is not square");

    Vector w(m.NRows());
    for (int k = 0; k < m.NRows(); k++) {
        double p = k;
        for (int i = k + 1; i < m.NRows(); i++) {
            if (std::fabs(m[i][k]) > std::fabs(m[p][k])) p = i;
        }
        if (std::fabs(m[p][k]) < m.GetEpsilon()) throw std::domain_error("Matrix is singular");
        if (p != k) {
            for (int l = 0; l < m.NCols(); l++) {
                std::swap(m[k][l], m[p][l]);
            }
            std::swap(w[k], w[p]);
        }
        w[k] = p;
        double mi = m[k][k];
        if (std::fabs(mi) < m.GetEpsilon()) throw std::domain_error("Matrix is singular");
        m[k][k] = 1;
        for (int j = 0; j < m.NRows(); j++) {
            m[k][j] = m[k][j] / mi;
        }
        for (int i = 0; i < m.NRows(); i++) {
            if (i != k) {
                mi = m[i][k];
                m[i][k] = 0;
                for (int j = 0; j < m.NRows(); j++) {
                    m[i][j] = m[i][j] - (mi * m[k][j]);
                }
            }
        }
    }

    for (int j = m.NRows() - 1; j >= 0; j--) {
        int p = static_cast<int>(w[j]);
        if (p != j) {
            for (int i = 0; i < m.NRows(); i++) {
                std::swap(m[i][p], m[i][j]);
            }
        }
    }

    return m;
}

Matrix RREF(Matrix m){
    Matrix result(m);
    int k = 0;
    int l = 0;
    double v = 0;
    int p = 0;
    double mi = 0;
    std::vector<bool> w(result.NCols());
    for(int j = 0; j < result.NCols(); j++){
        w[j] = false;
    }
    while(k < result.NRows() && l < result.NCols()){
        v = 0;
        while(v < result.GetEpsilon() && l < result.NCols()){
            p = k;
            for(int i = k; i < result.NRows(); i++){
                if(std::fabs(result[i][l]) > v){
                    v = std::fabs(result[i][l]);
                    p = i;
                }
            }
            if(v < result.GetEpsilon()) l++;
        }
        if(l < result.NCols()){
            w[l] = true;
            if(p != k) {
                for(int s = 0; s < result.NCols(); s++) {
                    double temp;
                    temp = result[k][s];
                    result[k][s] = result[p][s];
                    result[p][s] = temp;
                }
                bool temp = w[k];
                w[k] = w[p];
                w[p] = temp;
            }
            mi = result[k][l];

            for(int j = l; j < result.NCols(); j++){
                result[k][j] = result[k][j] / mi;
            }
            for(int i = 0; i < result.NRows(); i++){
                if(i != k){
                    mi = result[i][l];
                    for(int j = l; j < result.NCols(); j++){
                        result[i][j] = result[i][j] - (mi * result[k][j]);
                    }
                }
            }
        }
        l = l + 1;
        k = k + 1;
    }
    return result;
}

int Rank(Matrix m){
    Matrix r = RREF(m);
    int rank = 0;
    bool flag = false;
    for(int i = 0; i < r.NCols(); i++){
        flag = false;
        for(int j = 0; j < r.NRows(); j++){
            if(flag){
                if(r[j][i] == 0) continue;
                flag = false;
            }
            else{
                if(r[j][i] == 1) flag = true;
            }
        }
        if(flag) rank++;
    }
    return rank;
}


class LUDecomposer{
    Vector w;
    Matrix L;
    Matrix U;
public:
    LUDecomposer(Matrix m): L(m.NRows(), m.NCols()), U(m.NRows(), m.NCols()), w(m.NRows()){
        if(m.NRows() != m.NCols()) throw std::domain_error("Matrix is not square");

        double s = 0;
        int p = 0;
        double mi = 0;
        for(int j = 0; j < m.NRows(); j++){
            for(int i = 0; i <= j; i++){
                s = m[i][j];
                for(int k = 0; k < i ; k++){
                    s = s - (m[i][k] * m[k][j]);
                }
                m[i][j] = s;
            }
            p = j;
            for(int i = j + 1; i < m.NRows(); i++){
                s = m[i][j];
                for(int k = 0; k < j ; k++){
                    s = s - (m[i][k] * m[k][j]);
                }
                m[i][j] = s;
                if(std::fabs(s) > std::fabs(m[p][j])){
                    p = i;
                }
            }
            if(std::fabs(m[p][j]) < m.GetEpsilon()) throw std::domain_error("Matrix is singular");
            if(p != j){
                for(int d = 0; d < m.NCols(); d++) {
                    double temp;
                    temp = m[j][d];
                    m[j][d] = m[p][d];
                    m[p][d] = temp;
                }
            }
            this->w[j] = p;
            mi = m[j][j];
            for(int i = j + 1; i < m.NRows(); i++){
                m[i][j] = m[i][j] / mi;
            }
        }

        for(int i = 0; i < m.NCols(); i++){
            for(int j = 0; j < m.NRows(); j++){
                if(i > j){
                    L[i][j] = m[i][j];
                }
                if(i == j){
                    L[i][j] = 1;
                    U[i][j] = m[i][j];
                }
                else if(i < j){
                    U[i][j] = m[i][j];
                }
            }
        }
    }
    void Solve(const Vector &b, Vector &x) const{
        if(L.NRows() != b.NElems()) throw std::domain_error("Incompatible formats");
        if(b.NElems() != x.NElems()) throw std::domain_error("Incompatible formats");
        Vector y(b.NElems());
        for(int i = 0; i < L.NCols(); i++){
            double s = b[i];
            for(int j = 0; j <= i - 1; j++){
                s = s - (L[i][j] * y[j]);
            }
            y[i] = s;
        }

        for(int i = U.NCols() - 1; i >= 0; i--){
            double s = y[i];
            for(int j = i + 1; j < U.NCols(); j++){
                s = s - (U[i][j] * x[j]);
            }
            x[i] = s / U[i][i];
        }
    }
    Vector Solve(Vector b) const{
        Vector x(b.NElems());
        Solve(b, x);
        return x;
    }
    void Solve(const Matrix &b, Matrix &x) const{
        if(L.NRows() != b.NRows()) throw std::domain_error("Incompatible formats");
        if(x.NRows() != b.NRows() || x.NCols() != b.NCols()) throw std::domain_error("Incompatible formats");
        Matrix y(b.NRows(), b.NCols());

        for(int d = 0; d < b.NCols(); d++){
            for(int i = 0; i < L.NRows(); i++){
                double s = b[i][d];
                for(int j = 0; j < i - 1; j++){
                    s = s - (L[i][j] * y[j][d]);
                }
                y[i][d] = s;
            }
        }

        for(int d = 0; d < b.NCols(); d++){
            for(int i = U.NRows() - 1; i >= 0; i--){
                double s = y[i][d];
                for(int j = i + 1; i < U.NRows(); j++){
                    s = s - (U[i][j] * x[j][d]);
                }
                x[i][d] = s /  U[i][i];
            }
        }

    }
    Matrix Solve(Matrix b) const{
        Matrix x(b.NRows(), b.NCols());
        Solve(b, x);
        return x;
    }
    Matrix GetCompactLU() const{
        Matrix result(L.NRows(), L.NCols());
        for(int i = 0; i < result.NRows(); i++){
            for(int j = 0; j < result.NCols(); j++){
                if(i > j){
                    result[i][j] = L[i][j];
                }
                else{
                    result[i][j] = U[i][j];
                }
            }
        }
        return result;
    }
    Matrix GetL() const{
        return L;
    }
    Matrix GetU() const{
        return U;
    }
    Vector GetPermutation() const{
        return w;
    }
};


class QRDecomposer{
    Matrix R;
    Matrix V;
public:
    explicit QRDecomposer(Matrix m) : R(m), V(m.NRows(), m.NCols()){
        if(m.NRows() < m.NCols()) throw std::domain_error("Invalid matrix format");
        for(int k = 0; k < R.NCols(); k++){
            double s = 0;
            for(int i = k; i < R.NRows(); i++){
                s = s + (R[i][k] * R[i][k]);
            }
            s = std::sqrt(s);
            double mi = std::sqrt(s * (s + std::fabs(R[k][k])));
            if(R[k][k] < 0){
                s = -s;
            }
            if(std::fabs(mi) < R.GetEpsilon()) throw std::domain_error("Matrix is singular");
            V[k][k] = (R[k][k] + s) / mi;
            for(int i = k + 1; i < R.NRows(); i++){
                V[i][k] = R[i][k] / mi;
            }
            R[k][k] = -s;
            for(int j = k + 1; j < R.NCols(); j++){
                s = 0;
                for(int i = k; i < R.NRows(); i++){
                    s = s + (V[i][k] * R[i][j]);
                }
                for(int i = k; i < R.NRows(); i++){
                    R[i][j] = R[i][j] - (s * V[i][k]);
                }
            }
        }
    }
    void Solve(const Vector &b, Vector &x) const{
        Matrix q = this->GetQ();
        Matrix r = this->GetR();
        Matrix a = r * q;
        if(a.NRows() != a.NCols()) throw std::domain_error("Matrix is not square");
        if(q.NRows() != b.NElems()) throw std::domain_error("Incompatible formats");
        if(b.NElems() != x.NElems()) throw std::domain_error("Incompatible formats");
        Vector y(b.NElems());
        for(int i = 0; i < q.NCols(); i++){
            double s = b[i];
            for(int j = 0; j <= i - 1; j++){
                s = s - (q[i][j] * y[j]);
            }
            y[i] = s;
        }

        for(int i = r.NCols() - 1; i >= 0; i--){
            double s = y[i];
            for(int j = i + 1; j < r.NCols(); j++){
                s = s - (r[i][j] * x[j]);
            }
            x[i] = s / r[i][i];
        }

    }
    Vector MulQWith(Vector v) const{
        if(V.NCols() != v.NElems()) throw std::domain_error("Incompatible formats");
        Vector res(v);
        for(int k = V.NCols() - 1; k >= 0; k--){
            double s = 0;
            for(int i = k; i < V.NRows(); i++){
                s = s + (V[i][k] * res[i]);
            }
            for(int i = k; i < V.NRows(); i++){
                res[i] = res[i] - (s * V[i][k]);
            }
        }
        return res;
    }
    Matrix MulQWith(Matrix m) const{
        if(V.NCols() != m.NRows()) throw std::domain_error("Incompatible formats");
        Matrix res(m);
        for(int d = 0; d < res.NCols(); d++){
            for(int k = V.NCols() - 1; k >= 0; k--){
                double s = 0;
                for(int i = k; i < V.NRows(); i++){
                    s = s + (V[i][k] * res[i][d]);
                }
                for(int i = k; i < V.NRows(); i++){
                    res[i][d] = res[i][d] - (s * V[i][k]);
                }
            }
        }
        return res;
    }
    Vector MulQTWith(Vector v) const{
        if(V.NCols() != v.NElems()) throw std::domain_error("Incompatible formats");
        Vector res(v);
        for(int k = 0; k < V.NCols(); k++){
            double s = 0;
            for(int i = k; i < V.NRows(); i++){
                s = s + (V[i][k] * res[i]);
            }
            for(int i = k; i < V.NRows(); i++){
                res[i] = res[i] - (s * V[i][k]);
            }
        }
        return res;
    }
    Matrix MulQTWith(Matrix m) const{
        if(V.NRows() != m.NRows()) throw std::domain_error("Incompatible formats");
        Matrix res(m);
        for(int d = 0; d < res.NCols(); d++){
            for(int k = 0; k < V.NRows(); k++){
                double s = 0;
                for(int i = k; i < V.NCols(); i++){
                    s = s + (V[k][i] * res[i][d]);
                }
                for(int i = k; i < V.NCols(); i++){
                    res[i][d] = res[i][d] - (s * V[k][i]);
                }
            }
        }
        return res;
    }
    Matrix GetQ() const{
        Matrix Q(R.NRows(), R.NRows());
        for(int j = 0; j < Q.NRows(); j++){
            for(int i = 0; i < Q.NRows(); i++){
                Q[i][j] = 0;
            }
            Q[j][j] = 1;
            for(int k = V.NCols(); k >= 0; k--){
                double s = 0;
                for(int i = k; i < Q.NRows(); i++){
                    s = s + (V[i][k] * Q[i][j]);
                }
                for(int i = k; i < Q.NRows(); i++){
                    Q[i][j] = Q[i][j] - (s * V[i][k]);
                }
            }
        }
        return Q;
    }
    Matrix GetR() const{
        Matrix res(R);
        for(int i = 0; i < R.NRows(); i++){
            for(int j = 0; j < R.NCols(); j++){
                if(j < i) res[i][j] = 0;
            }
        }
        return res;
    }
};


int main() {
    //------------------------------------------------------------------------------------------------------------------
    std::cout << "Testiranje metoda klase Vector" << std::endl;
    std::cout << "Testiranje metode Chop:" << std::endl;
    try{
        Vector vtest({4, 3, 5, 2, 5, 0.00000000000000000004});
        vtest.Chop();       // Chop sa defaultnom vrijednoscu eps
        vtest.Print();

        std::cout << std::endl;

        vtest = {5, 4, 2, 4, 1};
        vtest.Chop(2);      // Chop sa proizvoljnom vrijednoscu eps
        vtest.Print();
    }
    catch(...){
    }

    std::cout << std::endl << "Testiranje metode EqualTo: " << std::endl;
    try{
        Vector vtest({4, 5, 3, 7, 3});
        Vector vtestcopy(vtest);

        std::cout << vtest.EqualTo(vtestcopy) << std::endl;         // EqualTo kada jesu jednaki

        Vector vtest2({5, 32, 6, 2});
        std::cout << vtest2.EqualTo(vtest);         // EqualTo kada nisu jednaki
    }
    catch(...){}

    //------------------------------------------------------------------------------------------------------------------
    std::cout << std::endl << "Testiranje Metoda klase Matrix" << std::endl;

    std::cout << "Testiranje metode Chop:" << std::endl;
    try{
        Matrix mtest1({{1, 3, 4}, {5, 6, 8}, {1, 2, 5}});
        mtest1.Chop();          // Chop bez parametra
        mtest1.Print();
        std::cout << std::endl;
        mtest1.Chop(2);         // Chop sa parametrom;
        mtest1.Print();
    }
    catch(...){}

    std::cout << std::endl << "Testiranje metode EqualTo: " << std::endl;
    try{
        Matrix mtest1({{5, 7, 3}, {6, 8, 3}});
        Matrix mtest2(mtest1);
        std::cout << mtest1.EqualTo(mtest2) << std::endl;           // EqualTo kada su jednake matrice
        Matrix mtest3({{4, 3}, {4, 3}});
        std::cout << mtest3.EqualTo(mtest1) << std::endl;           // EqualTo kada su matrice razlicitih dimenzija
        Matrix mtest4({{5, 7, 3}, {6, 8, 4}});
        std::cout << mtest4.EqualTo(mtest2) << std::endl;           // EqualTo kada se razlikuju samo za jednu vrijednost
    }
    catch(...){}

    std::cout << "Testiranje matricnog lijevog kolicnika: " << std::endl;
    try{
        try{
            Matrix mtest1({{5, 2, 5}, {4, 5, 2}});
            Matrix mtest2({{2, 4}, {7, 8}});
            Matrix res = LeftDiv(mtest1, mtest2);           // LeftDiv kada prva matrica nije kvadratna
            res.Print();
        }
        catch(std::domain_error izuzetak){
            std::cout << izuzetak.what() << std::endl;
        }
        try{
            Matrix mtest1({{6, 8}, {5, 2}});
            Matrix mtest2({{2, 6}, {9, 7}});
            Matrix res = LeftDiv(mtest1, mtest2);       // Pravilno korisenje LeftDiv-a
            res.Print();
        }
        catch(...){}
        try{
            Matrix mtest1({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
            Matrix mtest2({{4, 5, 3}, {6, 3, 2}, {3, 63, 12}});
            Matrix res = LeftDiv(mtest1, mtest2);       // LeftDiv kada je prva matrica singularna
        }
        catch(std::domain_error e){
            std::cout << e.what() << std::endl;
        }
        try{
            Matrix mtest1({{4, 5}, {7, 3}});
            Vector vtest1({5, 3});
            Vector res = LeftDiv(mtest1, vtest1);
            res.Print();
        }
        catch(...){}
    }
    catch(...){}
    std::cout << std::endl << "Testiranje desnog matricnog dijeljenja i dijeljenja skalarom: " << std::endl;
    try{
        try{
            Matrix mtest1({{2, 4}, {6, 8}});
            double s = 2;
            Matrix res = mtest1/s;          // Dijeljenje matrice skalarom
            res.Print();
            std::cout << std::endl;
        }
        catch(...){}
        try{
            Matrix mtest1({{2, 4}, {6, 8}});
            double s = 0;
            Matrix res = mtest1 / s;        // Dijeljenje matrice nulom
        }
        catch(std::domain_error e){
            std::cout << e.what() << std::endl;
        }
        try{
            Matrix mtest1({{2, 4}, {6, 8}});
            double s = 2;
            mtest1 /= s;            // Testiranje operatora /=
            mtest1.Print();
        }
        catch(...){}
        try{
            Matrix mtest1({{0, 3, 2}, {4, 6, 1}, {3, 1, 7}});
            Matrix mtest2({{4, 1, 5}, {1, 2, 1}});
            Matrix res = mtest1/mtest2;         // Testiranje desnog matricnog dijeljenja
            res.Print();
        }
        catch(...){}
        try{
            Matrix mtest1({{6, 8, 8}, {4, 2, 1}, {9, 2, 3}});
            Matrix mtest2{{2, 5, 2}, {8, 2, 1}, {3, 7, 1}};
            mtest1/=mtest2;         // Testiranje operatora /=
            mtest1.Print();
        }
        catch(...){}
    }
    catch(...){}
    std::cout << std::endl << "Testiranje racunanja determinante matrice: " << std::endl;
    try{
        try{
            Matrix mtest1({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
            double s = Det(mtest1);         // Testiranje funkcije Det na singularnoj matrici
            std::cout << s << std::endl;
        }
        catch(...){}
        try{
            Matrix mtest1({{5, 3, 2}, {4, 2, 1}, {6, 2, 5}});
            std::cout << mtest1.Det() << std::endl;         // Testiranje metode Det
        }
        catch(...){}
        try{
            Matrix mtest1({{5, 4, 2}, {5, 5, 2}});
            double s = Det(mtest1);             // Testiranje funkcije Det na matrici koja nije kvadratna
        }
        catch(std::domain_error e){
            std::cout << e.what() << std::endl;
        }
    }
    catch(...){}
    std::cout << std::endl << "Testiranje pronalazenja inverzne matrice: " << std::endl;
    try{
        try{
            Matrix mtest1({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
            Matrix inv = Inverse(mtest1);           // testiranje funckije Inverse kada joj je data singulana matrica
        }
        catch(std::domain_error e){
            std::cout << e.what() << std::endl;
        }
        try{
            Matrix mtest1({{1, 2}, {3, 4}});
            Matrix inv = Inverse(mtest1);           // Pravilno koristenje funkcije Inverse
            inv.Print();
        }
        catch(...){}
        try{
            Matrix mtest1({{1, 2}, {3, 4}});
            mtest1.Invert();                    // Pravilno koristenje metode Invert
            mtest1.Print();
        }
        catch(...){}
        try{
            Matrix mtest1({{2, 3, 4}, {7, 3, 2}});
            mtest1.Invert();            // Metoda Invert kada joj je poslana nekvadratna matrica
        }
        catch(std::domain_error e){
            std::cout << e.what() << std::endl;
        }
    }
    catch(...){}
    std::cout << std::endl << "Testiranje svodjenja na RREF formu: " << std::endl;
    try{
        try{
            Matrix mtest1({{5, 6, 7, 8, 1}, {7, 0, 8, 2, 9}, {9, 2, 6, 1, 5}});
            Matrix rref = RREF(mtest1);         // Testiranje funkcije RREF
            rref.Print();
        }
        catch(...){}
        try{
            Matrix mtest1({{5, 0, 6}, {3, 2, 1}, {5, 1, 6}});
            mtest1.ReduceToRREF();          // Testiranje kada je proslijedjena kvadratna matrica
            mtest1.Print();
        }
        catch(...){}
    }
    catch(...){}
    std::cout << std::endl << "Racunanje ranga matrice: " << std::endl;
    try{
        try{
            Matrix mtest1({{5, 6, 7, 8, 1}, {7, 0, 8, 2, 9}, {9, 2, 6, 1, 5}});
            double s = Rank(mtest1);        // testiranje funkcije Rank
            std::cout << s << std::endl;
        }
        catch(...){}
        try{
            Matrix mtest1({{5, 6, 3, 2}, {6, 3, 4, 2}, {4, 2, 6, 1}, {4, 1, 6, 1}});
            std::cout << mtest1.Rank() << std::endl;
        }
        catch(...){}
    }
    catch(...){}

    //------------------------------------------------------------------------------------------------------------------

    std::cout << std::endl << "Testiranje metoda klase LUDecomposer" << std::endl;
    std::cout << "Testiranje LU dekompozicije: " << std::endl;
    try{
        try{
            Matrix mtest1({{5, 1, 5}, {8, 3, 7}, {2, 5, 3}});
            LUDecomposer lu1(mtest1);           // Testiranje konstruktora LUDecomposer klase
            lu1.GetL().Print();             // L matrica
            lu1.GetU().Print();             // U matrica
            lu1.GetCompactLU().Print();     // CompactLU matrica
        }
        catch(...){}
        try{
            Matrix mtest1({{5, 1, 2}, {5, 1, 6}});
            LUDecomposer lu1(mtest1);           // Testiranje LUDecomposer konstruktora kada se proslijedi nekvadratna matrica
        }
        catch(std::domain_error e){
            std::cout << e.what() << std::endl;
        }
        try{
            Matrix mtest1({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
            LUDecomposer lu1(mtest1);           // testiranje LUDecomposer konstruktora kada se proslijedi singularna matrica
        }
        catch(std::domain_error e){
            std::cout << e.what() << std::endl;
        }
    }
    catch(...){}

    std::cout << "Testiranje Solve metoda za LU faktorizaciju: " << std::endl;
    try{
        try{
            Matrix mtest1({{5, 6, 7}, {2, 3, 4}, {8, 4, 2}});
            LUDecomposer lu1(mtest1);
            Vector b({4, 6, 3});
            lu1.Solve(b, b);         // Testiranje Solve kada se rjesenje upisuje preko drugog proslijedjenog vektora
            b.Print();
        }
        catch(...){}
        try{
            Matrix mtest1({{5, 6, 7}, {2, 3, 4}, {8, 4, 2}});
            LUDecomposer lu1(mtest1);
            Vector b({4, 6, 3});
            Vector res = lu1.Solve(b);
            res.Print();
        }
        catch(...){}
        try{
            Matrix mtest1({{5, 6, 7}, {2, 3, 4}, {8, 4, 2}});
            LUDecomposer lu1(mtest1);
            Vector b({2, 3, 4, 5});
            lu1.Solve(b, b);
            b.Print();
        }
        catch(std::domain_error izuzetak){
            std::cout << izuzetak.what() << std::endl;
        }
        try{
            Matrix mtest1({{5, 6, 7}, {2, 3, 4}, {8, 4, 2}});
            LUDecomposer lu1(mtest1);
            Vector b({3, 2, 4});
            Vector x({2, 3, 4, 1});
            lu1.Solve(b, x);
        }
        catch(std::domain_error e){
            std::cout << e.what() << std::endl;
        }

        try{
            Matrix mtest1({{5, 6, 7}, {2, 3, 4}, {8, 4, 2}});
            LUDecomposer lu1(mtest1);
            Matrix b({{2, 4}, {3, 1}});
            lu1.Solve(b, b);
            b.Print();
        }
        catch(std::domain_error e){
            std::cout << e.what() << std::endl;
        }
    }
    catch(...){}

    //------------------------------------------------------------------------------------------------------------------
    std::cout << "Testiranje QR faktorizacije: " << std::endl;
    try{
        try{
            Matrix mtest1({{2, 4, 5}, {6, 2, 7}, {5, 7, 2}});
            QRDecomposer qr1(mtest1);
            qr1.GetR().Print();
            qr1.GetQ().Print();
        }
        catch(...){}
        try{
            Matrix mtest1({{4, 2, 5, 1}, {4, 2, 3, 4}, {5, 2, 5, 1}});
            QRDecomposer qr1(mtest1);
        }
        catch(std::domain_error e){
            std::cout << e.what() << std::endl;
        }
        try{
            Matrix mtest1({{4, 5, 2}, {6, 3, 1}, {7, 3, 1}, {6, 8, 2}});
            QRDecomposer qr1(mtest1);
            qr1.GetR().Print();
            qr1.GetQ().Print();
        }
        catch(...){}
    }
    catch(...){}
    std::cout << std::endl << "Testiranje mnozenja Q matrice: " << std::endl;
    try{
        try{
            Matrix mtest1({{4, 5, 2}, {6, 3, 1}, {7, 3, 1}, {6, 8, 2}});
            QRDecomposer qr1(mtest1);
            Vector v1({4, 5, 3});
            Matrix mtest2 = qr1.MulQWith(v1);
            mtest2.Print();
        }
        catch(...){}
    }
    catch(...){}
    return 0;
}
