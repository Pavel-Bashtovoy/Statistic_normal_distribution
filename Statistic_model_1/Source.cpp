#include <iostream>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;
const double PI = 3.141592653589793238463;
const int N = 500;
const double lambda = 1.0;
const double tetta = 3.07528535261;
const double epsilon = 0.3;
//const double tetta = 2.991;

const double delta = 0.1;//0,1 0,5 1
const double eps = 1.E-4;


//U  равномерно распределенные на отрезке
double randomRangeU(double min, double max) {
   return min + rand() / (double)RAND_MAX * (max - min);
}

//плотность распределения
double DistributionDensity(double x) {

    double f;

    f = (1 / sqrt(2 * PI)) * exp((- x * x) / 2);

    return f;
}

//производная плотности распределения 
double DistributionDensityDerivative(double x) {

    double f;

    f = -x * sqrt(2) * (1 /(2 * sqrt( PI))) * exp((-x * x) / 2);

    return f;
}

//Генерация случайной величины
double GeneratingRandomVariable(double u1, double u2) {

    double x;

    x = sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
    if (x > 3.7) {
        return 3.04;
    } else if (x < -3.7) {
        return -3.04;
    } else {

        return x;  
    }
}

//выборочная медиана;
double median(std::vector<double> &y) {

    auto m = y.begin() + y.size() / 2;
    std::nth_element(y.begin(), m, y.end());
    return *m;
}
//выборочная дисперсия
double dispersion(double  yArithmeticMean, std::vector<double>& y) {
    double D;
    double sum = 0.0;

    for (int i = 0; i < N;i++) {
        sum += pow((y[i] - yArithmeticMean), 2);
    }

    D = sum / N;
    return D;
 }
//коэффициентов асимметрии 
double betta_3(double  yArithmeticMean, double D, std::vector<double>& y){
    double sum = 0.0, betta;

    for (int i = 0; i < N;i++) {
        sum += pow((y[i] - yArithmeticMean),3);
    }

    betta = sum / (N * sqrt(pow(D, 3)));

    return betta;
}
 //эксцесса
double betta_4(double  yArithmeticMean, double D, std::vector<double>& y) {
    double sum = 0.0, betta;

    for (int i = 0; i < N;i++) {
        sum += pow((y[i] - yArithmeticMean), 4);
    }

    betta = sum / (N * pow(D, 2));

    return betta;
}

//оценка максимального правдоподобия;
double lossMaximumLikelihoodEstimation (double y) {
    double lossFunction = 0;

    lossFunction = - log(DistributionDensity((y - tetta) / lambda));
  
    return lossFunction;
}
//обобщенные радикальные оценки с разными значениями параметра(как минимум три обязательных значения 0.1, 0.5, 1).
//  функция потерь 
double lossGeneralizedRadicalAssessments( double y ) {
    double lossFunction = 0;

    lossFunction = (- 1 / pow( DistributionDensity(0) , delta)) * pow( DistributionDensity((y - tetta) / lambda), delta);

    return lossFunction;
}

//Оценочная функция 
double evaluationFunction(double y) {
    double evaluation = 0;
    double c = 1;

    evaluation = c * DistributionDensityDerivative( (y - tetta) / lambda ) * pow( DistributionDensity((y - tetta) / lambda), delta - 1 );

    return evaluation;
}
/**Использовать выборки, имеющие следующие виды распределений :
    чистое распределение;
    засоренное распределение с асимметричным засорением
    засоренное распределение с симметричным засорением(равные сдвиги у чистого и засоряющего распределений, масштаб у засоряющего больше в 2 - 3 раза, чем у чистого);
*/
// метод золотого сечения (для минимизации функций и нахождения ОМП)
void golden_ratio_for_omp(double a, double b, double eps, ofstream& out)
{
    double f1, f2, x1, x2;
    double length, PREVlength;
    double div;
    length = fabs(b - a);
    int iter = 0;//номер текущей итерации 
    x1 = a + 0.381966011 * (b - a);
    x2 = a + 0.618003399 * (b - a);
    f1 = lossMaximumLikelihoodEstimation(x1);
    f2 = lossMaximumLikelihoodEstimation(x2);

    out << "n|\t" << "a\t\t  " << "|\t" << "b\t\t   " << "|" << "b(i)-a(i) \t|" << "x1\t\t  |\t" << "x2\t\t   |\t" << "f1\t\t|\t\t" << "f2" << endl;
    while (length > eps)
    {
        iter++;
        if (f1 <= f2)
        {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + 0.381966011 * (b - a);
            f1 = lossMaximumLikelihoodEstimation(x1);
        }
        else
        {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + 0.618003399 * (b - a);
            f2 = lossMaximumLikelihoodEstimation(x2);
        }
        PREVlength = length;
        length = abs(b - a);
        div = PREVlength / length;
        out.setf(ios::scientific);
        out << iter << "|" << a << "|" << b << "|" << length << "|"  << x1 << "|" << x2 << "|" << f1 << "|" << f2 << endl;
    }
}

// метод золотого сечения (для минимизации функций и нахождения ОРО)
void golden_ratio_for_ORO(double a, double b, double eps, ofstream& out)
{
    double f1, f2, x1, x2;
    double length, PREVlength;
    double div;
    length = fabs(b - a);
    int iter = 0;//номер текущей итерации 
    x1 = a + 0.381966011 * (b - a);
    x2 = a + 0.618003399 * (b - a);
    f1 = lossGeneralizedRadicalAssessments(x1);
    f2 = lossGeneralizedRadicalAssessments(x2);

    out << "n|\t" << "a\t\t  " << "|\t" << "b\t\t   " << "|" << "b(i)-a(i) \t|" << "x1\t\t  |\t" << "x2\t\t   |\t" << "f1\t\t|\t\t" << "f2" << endl;
    while (length > eps)
    {
        iter++;
        if (f1 < f2)
        {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + 0.381966011 * (b - a);
            f1 = lossGeneralizedRadicalAssessments(x1);
        }
        else
        {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + 0.618003399 * (b - a);
            f2 = lossGeneralizedRadicalAssessments(x2);
        }
        PREVlength = length;
        length = abs(b - a);
        div = PREVlength / length;
        out.setf(ios::scientific);
        out << iter << "|" << a << "|" << b << "|" << length << "|" << x1 << "|" << x2 << "|" << f1 << "|" << f2 << endl;
    }
}

vector<double> xRandom (std::vector<double>& u1, std::vector<double>& u2 , double min, double max) {

    vector<double> x;
    x.resize(N);

    //Получение U1  и U2
    for (int i = 0; i < N; i++)
    {
        u1[i] = randomRangeU(min, max);
        u2[i] = randomRangeU(min, max);
        std::cout << i << '|' << u1[i] << '|' << u2[i] << endl;

    }

    // Найдем случайную величину х
    for (int i = 0; i < N; i++)
    {
        x[i] = GeneratingRandomVariable(u1[i], u2[i]);
    }

    return  x;
}

int main()
{
    double min = 0.0;
    double max = 1.0;
    double sum_y = 0;
    double yArithmeticMean, M, D, b_3 , b_4;
    vector<double> u1;
    vector<double> u2;
    u1.resize(N);
    u2.resize(N);

    vector<double> x;
    x.resize(N);
    x = xRandom(u1, u2, min, max);
    vector<double> y;

    ofstream fout;




    y.resize(N);

    setlocale(LC_ALL, "Russian");


    //Найдем случайную величину y  и y`- среднее y , найдем значения функции плотности
    fout.open("DistributionDensity.txt");

    double a = -7, b = 3.7;
    double step = (b - a) / N;
    double point = 0;
    //std::cout << "Найдем y:" << endl;
   // std::cout << 'i' << '|' << "Y[i]    " << '|' << "Func" << endl;
    fout << "i    " << '|' << "y(i)" << endl;

    for (int i = 0; i < N; i++)
    {
        y[i] = x[i] * lambda + tetta;

        sum_y += y[i];

      // std::cout << i << '|' << y[i] <<'|'<< DistributionDensity(y[i]) << endl;

        //fout << y[i] << ','<< DistributionDensity(y[i]) << endl;
      // fout <<i<<'|'<< y[i]  << endl;
       point =a + i * step;
     //  point = point * lambda + tetta;
       
      fout << point << ',' << DistributionDensity(point + tetta) << endl;


    }
    fout.close();

    /// Среднее арифметическое 
    yArithmeticMean = sum_y / N;
    std::cout << " Среднее арифметическое  y:" << yArithmeticMean<< endl;
    //медианы
     M = median(y);
     std::cout << " Медиана  y:" <<M<< endl;
    //дисперсии 
     D = dispersion(yArithmeticMean, y);
     std::cout << "Дисперсия y:" << D << endl;
    //коэффициентов асимметрии 
     b_3 = betta_3(yArithmeticMean, D, y);
     std::cout << "Коэффициент асимметрии:" << b_3 << endl;
     //эксцесса;
     b_4 = betta_4(yArithmeticMean, D, y);
     std::cout << "Коэффициент Эксцесса:" << b_4 << endl;

     //минимизация для омп
     fout.open("OMP.txt");
     golden_ratio_for_omp(min, max, eps, fout);
     fout.close();

     //минимизация для ORO
     fout.open("ORO.txt");
     golden_ratio_for_ORO(min, max, eps, fout);
     fout.close();
     
     return 0;
}