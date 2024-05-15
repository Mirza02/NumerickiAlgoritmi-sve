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

int main() {
    // Testiranje klase Vector
    //__________________________________________________________________________________________________________________
    std::cout << "Konstruktor kada se unese negativna vrijednost: " << std::endl;
    try{
        Vector vektor_negativna_vel(-1);    // testiranje konstruktora kada se unese negativna vrijednost
    }
    catch(std::range_error izuzetak){
        std::cout << izuzetak.what() << std::endl;
    }
    std::cout << "Sekvencijalni kontruktor kada se unese prazna lista: " << std::endl;
    try{
        Vector vektor_inic_lista_0({});     // testiranje sekvencijalnog konstruktora kada se posalje prazna inicijalizacijska lista
    }
    catch(std::range_error izuzetak){
        std::cout << izuzetak.what() << std::endl;
    }
    Vector vektor_test_1(5);                // pravilna inicijalizcija
    Vector vektor_test_2({5, 4, 3, 2, 1});      // pravilna upotreba sekvencijalnog konstruktora
    try{
        for(int i = 0; i < vektor_test_1.NElems(); i++){
            vektor_test_2(i) = i;               // testiranje operatora () za indeks koji pocinje od nule
        }
    }
    catch(std::range_error izuzetak){
        std::cout << izuzetak.what() << std::endl;
    }
    for(int i = 1; i <= vektor_test_1.NElems(); i++) {
        vektor_test_1(i) = i + 10;
    }

    std::cout << vektor_test_1.Norm() << " = " << VectorNorm(vektor_test_1) << std::endl;        // testiranje jednakosti obje varijante funkcije za racunanje euklidske norme vektora
    std::cout << "Vrijednost epsilona za oba vektora: " << std::endl;
    std::cout << vektor_test_1.GetEpsilon() << " , " << vektor_test_2.GetEpsilon() << std::endl;         // testiranje funkcija za racunanje epsilona

    vektor_test_1 += vektor_test_2;         // testiranje operatora +=
    vektor_test_1.Print();          // testiranje metode Print
    PrintVector(vektor_test_1);     // testiranje da li ce prijateljska funkcija PrintVector imati isti ispis

    double pomnozeni_vektori = (vektor_test_1 * vektor_test_2);         // testiranje operatora mnozenja dva vektora
    std::cout << "Rezultat mnozenja dva vektora: " << std::endl << pomnozeni_vektori << std::endl;

    std::cout << "Sabiranje dva vektora razlicitih duzina: " << std::endl;
    try{
        Vector duzi_vektor({1, 2, 3, 4, 5, 6, 7, 8, 9});
        Vector sabrani_vektori_test1 = (vektor_test_1 + duzi_vektor);       // testiranje sabiranja kada nisu vektori jednake duzine
    }
    catch(std::domain_error izuzetak){
        std::cout << izuzetak.what() << std::endl;
    }

    std::cout << "Sabiranje dva vektora jednakih duzina: " << std::endl;
    PrintVector((vektor_test_1 + vektor_test_2), '/');      // testiranje drugacijeg seperatora u funkciji PrintVector
    vektor_test_1 += vektor_test_2;         // testiranje operatora +=
    vektor_test_1.Print(',', 20);      // testiranje veceg epsilona od defaultnog

    std::cout << "Oduzimanje dva vektora razlicitih duzina: " << std::endl;
    try{
        Vector vektor_test_oduzimanje({1, 2, 3});       // testiranje operatora - kada vektori nisu iste duzine
        Vector rezultat = vektor_test_1 - vektor_test_oduzimanje;
    }
    catch(std::domain_error izuzetak){
        std::cout << izuzetak.what() << std::endl;
    }

    std::cout << "Oduzimanje dva vektora: " << std::endl;
    PrintVector(vektor_test_1 - vektor_test_2, ',');        // testiranje operatora - za vektore iste duzine
    vektor_test_1 -= vektor_test_2;         // testiranje operatora -=
    vektor_test_1.Print('/');

    std::cout << "Mnozenje vektora skalarom na oba nacina: " << std::endl;
    PrintVector(vektor_test_1 * 5, ',');        // mnozenje kada je skalar sa desne strane
    PrintVector(5 * vektor_test_1, ',');        // mnozenje kada je skalar sa lijeve strane

    std::cout << "Dijeljenje vektora nulom: " << std::endl;
    try{
        Vector rezultat_dijeljenje = (vektor_test_1 / 0);       // dijeljenje vektora nulom
    }
    catch(std::domain_error izuzetak){
        std::cout << izuzetak.what() << std::endl;
    }

    std::cout << "Dijeljenje vektora: " << std::endl;
    vektor_test_1 /= 5;
    vektor_test_1.Print(',');


    // Testiranje klase Matrix
    // _________________________________________________________________________________________________________________
    std::cout << "Testiranje matrice: " << std::endl << std::endl;
    std::cout << "Inicijalizacija negativnim vrijednostima: " << std::endl;
    try{
        Matrix matrica_test_negativni(-1, 5);       // testiranje konstruktora sa jednom negativnom dimenzijom
    }
    catch(std::range_error izuzetak){
        std::cout << izuzetak.what() << std::endl;
    }

    std::cout << "Inicijalizacija matrice kada je bar jedna od vrijednosti nula: " << std::endl;
    try{
        Matrix matrica_test_nula(0, 2);         // testiranje konstruktora kada je jedna dimenzija nula
    }
    catch(std::range_error izuzetak){
        std::cout << izuzetak.what() << std::endl;
    }

    std::cout << "Pravilna inicijalizacija matrice dimenzija 4 x 3" << std::endl;
    Matrix matrica_test_1(4, 3);            // pravilno koristenje konstruktora

    std::cout << "Inicijalizacija matrice praznom listom: " << std::endl;
    try{
        Matrix matrica_test_prazna_lista({});        // inicijalizacija matrice praznom listom
    }
    catch(std::range_error izuzetak){
        std::cout << izuzetak.what() << std::endl;
    }

    std::cout << "Inicijalizacija matrice grbavom matricom: " << std::endl;
    try{
        Matrix matrica_test_grbava({{1, 2, 3}, {2, 7}, {6, 35, 7, 2}});     // inicijalizacija grbavom matricom
    }
    catch(std::logic_error izuzetak){
        std::cout << izuzetak.what() << std::endl;
    }

    Matrix matrica_test_2({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {10, 11, 12}});         // inicijalizacija listom istih dimenzija kao matrica_test_1

    Matrix matrica_test_vektor({1, 2, 3, 4, 5, 6});     // inicijalizacija matrice vektorom 1
    Matrix matrica_test_vektor_2(vektor_test_1);        // inicijalizacija matrice vektorom 2


    int k = 1;
    for(int i = 0; i < matrica_test_1.NRows(); i++){        // testiranje metode NRows()
        for(int j = 0; j < matrica_test_1.NCols(); j++){        // testiranje metode NCOls()
            matrica_test_1[i][j] = k;               // postavljanje vrijednosti matrice koristenjem operatora []
            k *= 4;
        }
    }

    std::cout << "Testiranje metode .Print(): " << std::endl;
    matrica_test_1.Print(5);            // testiranje metode .Print() sa proizvoljnom sirinom ispisa
    std::cout << "Testiranje funkcije PrintMatrix(): " << std::endl;
    PrintMatrix(matrica_test_1, 5, 10);     // testiranje prijateljske funkcije PrintMatrix koristenjem proizvoljne duzine i proizvoljnog epsilona
    PrintMatrix(matrica_test_1);            // testiranje prijateljske funkcije sa podrazumijevanim vrijednostima

    std::cout << "Nepravilno koristenje operatora []: " << std::endl;
    try{
        std::cout << matrica_test_2[7][0];          // testiranje operatora [] kada je vrijednost veca od dimenzija matrice
    }
    catch(std::range_error izuzetak){
        std::cout << izuzetak.what() << std::endl;
    }
    try{
        std::cout << matrica_test_1[2][-1];         // testiranje operatora [] kada je vrijednost negativna
    }
    catch(std::range_error izuzetak){
        std::cout << izuzetak.what() << std::endl;
    }

    std::cout << "Testiranje operatora (): " << std::endl;
    try{
        std::cout << matrica_test_1(2, 0);       // testiranje kada je jedna od vrijednosti 0
    }
    catch(std::range_error izuzetak){
        std::cout << izuzetak.what() << std::endl;
    }
    try{
        std::cout << matrica_test_2(3, 4);       //testiranje operatora () kada su vrijednosti jednake dimenzijama matrice (ovo radi)
    }
    catch(std::range_error izuzetak){
        std::cout << izuzetak.what() << std::endl;
    }

    std::cout << "Testiranje jednakosti funkcije MatrixNorm() i metode .Norm() za obje matrice: " << std::endl;
    std::cout << MatrixNorm(matrica_test_1) << " " << matrica_test_1.Norm() << std::endl;
    std::cout << MatrixNorm(matrica_test_2) << " " << matrica_test_2.Norm() << std::endl;

    std::cout << "Testiranje sabiranja matrica: " << std::endl;
    try{
        Matrix matrica(matrica_test_1 + Matrix({{1, 2}, {2, 3}}));      // sabiranje matrica razlicitih dimenzija
    }
    catch(std::domain_error izuzetak){
        std::cout << izuzetak.what() << std::endl;
    }
    PrintMatrix(matrica_test_1 + matrica_test_2);       // pravilno sabiranje matrica
    matrica_test_1 += Matrix({{1, 2, 3}, {7, 32, 2}, {32, 63, 123}, {234, 532, 12}});       // pravilno sabiranje matrica koristenjem operatora +=
    matrica_test_1.Print(20, 10);

    std::cout << "Testiranje oduzimanja matrica" << std::endl;
    try{
        Matrix matrica(matrica_test_2 - Matrix({{1, 2}, {6, 32}}));     // oduzimanje matrica razlicitih dimenzija
    }
    catch(std::domain_error izuzetak){
        std::cout << izuzetak.what() << std::endl;
    }
    PrintMatrix(Matrix({{23, 432, 213, 532}, {5, 653, 23, 523}}) - Matrix({{2, 5, 1, 43}, {84, 2, 7, 3}}), 20, 30);     // pravilno oduzimanje matrica
    matrica_test_1 -= Matrix({{1, 2, 3}, {7, 32, 2}, {32, 63, 123}, {234, 532, 12}});       // testiranje operatora -=
    matrica_test_1.Print(20, 10);

    std::cout << "Testiranje mnozenja matrice skalarom: " << std::endl;
    PrintMatrix(matrica_test_2 * 10, 15);       // testiranje mnozenja kada je skalar desno
    PrintMatrix(10 * matrica_test_2, 15);       // testiranje mnozenja kada je skalar lijevo
    matrica_test_2 *= 10;       // testiranje operatora *= za mnozenje sa skalarom
    PrintMatrix(matrica_test_2);

    std::cout << "Testiranje mnozenja dvije matrice: " << std::endl;
    try{
        Matrix mnozenje_test(matrica_test_1 * Matrix({{2, 23}, {231, 21}}));        // testiranje nekompatibilnih matrica
    }
    catch(std::domain_error izuzetak){
        std::cout << izuzetak.what() << std::endl;
    }
    Matrix mnozenje_test(matrica_test_1 * Transpose(matrica_test_2));      // mnozenje dvije matrice
    mnozenje_test.Print(15);

    matrica_test_1 *= Transpose(matrica_test_2);       // testiranje operatora *=
    matrica_test_1.Print(15);

    PrintMatrix(matrica_test_2 * Matrix({{4321, 432, 32, 432}, {432, 32, 76, 23}, {543, 764, 32, 432}}), 15, 200);       // testiranje mnozenja dvije matrice razlicitih dimenzija

    std::cout << "Testiranje mnozenja vektorom: " << std::endl;
    try{
        Matrix mnozenje_vektorom(matrica_test_2 * Vector({5, 3, 2, 53, 64, 32}));           // mnozenje nekompatibilnim vektorom
    }
    catch(std::domain_error izuzetak){
        std::cout << izuzetak.what() << std::endl;
    }
    PrintVector(matrica_test_1 * Vector({2, 4, 5, 2}));     // pravilno mnozenje matrice vektorom

    std::cout << "Testiranje transpozicije matrice: " << std::endl;
    PrintMatrix(Transpose(matrica_test_2));     // testiranje funkcije Transpose;
    matrica_test_2.Transpose();             // testiranje metode Transpose;
    matrica_test_2.Print();
    Matrix matrica_kvadratna({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
    matrica_kvadratna.Transpose();          // metoda Transpose uradjena "u mjestu" (iskoristena na kvadratnoj matrici)
    return 0;
}
