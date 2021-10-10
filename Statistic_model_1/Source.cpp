#include <iostream>
#include <stdio.h>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;
const double PI = 3.141592653589793238463;
const int N = 200;
const float lambda = 1.0;
const float tetta = 0.0;


//U  ���������� �������������� �� �������
float randomRangeU(float min, float max) {
   return min + rand() / (float)RAND_MAX * (max - min);
}

//��������� �������������
float DistributionDensity(float x) {

    float f;

    f = (1 / sqrt(2 * PI)) * exp((- x * x) / 2);

    return f;
}

//����������� ��������� ������������� 
float DistributionDensityDerivative(float x) {

    float f;

    f = -x * sqrt(2) * (1 /(2 * sqrt( PI))) * exp((-x * x) / 2);

    return f;
}

//��������� ��������� ��������
float GeneratingRandomVariable(float u1, float u2) {

    float x;

    x = sqrt(-2 * log(u1)) * cos(2 * PI * u2);

    return x;
}

//���������� �������;
float median(std::vector<float> &y) {

    auto m = y.begin() + y.size() / 2;
    std::nth_element(y.begin(), m, y.end());
    return y[y.size() / 2];
}
//���������� ���������
float dispersion(float  yArithmeticMean, std::vector<float>& y) {
    float D;
    float sum = 0.0;

    for (int i = 0; i < N;i++) {
        sum = pow((y[i] - yArithmeticMean), 2);
    }

    D = sum / N;
    return D;
 }
//������������� ���������� 
float betta_3(float  yArithmeticMean, float D, std::vector<float>& y){
    float sum = 0.0, betta;

    for (int i = 0; i < N;i++) {
        sum = pow((y[i] - yArithmeticMean),3);
    }

    betta = sum / (N * sqrt(pow(D, 3)));

    return betta;
}
 //��������
float betta_4(float  yArithmeticMean, float D, std::vector<float>& y) {
    float sum = 0.0, betta;

    for (int i = 0; i < N;i++) {
        sum = pow((y[i] - yArithmeticMean), 4);
    }

    betta = sum / (N * pow(D, 2));

    return betta;
}

//������ ������������� �������������;
float lossMaximumLikelihoodEstimation (float y, float tetta) {
    float lossFunction = 0;
    
    lossFunction = - log(DistributionDensity((y - tetta) / lambda));
  
    return lossFunction;
}
//���������� ����������� ������ � ������� ���������� ���������(��� ������� ��� ������������ �������� 0.1, 0.5, 1).
//  ������� ������ 
float lossGeneralizedRadicalAssessments( float y, float tetta ) {
    float lossFunction = 0;
    float  delta = 1;

    lossFunction = (- 1 / pow( DistributionDensity(0) , delta)) * pow( DistributionDensity((y - tetta) / lambda), delta);

    return lossFunction;
}

//��������� ������� 
float evaluationFunction(float y, float tetta) {
    float evaluation = 0;
    float  delta = 1;
    float c = 1;

    evaluation = c * DistributionDensityDerivative( (y - tetta) / lambda ) * pow( DistributionDensity((y - tetta) / lambda), delta - 1 );

    return evaluation;
}
/**������������ �������, ������� ��������� ���� ������������� :
    ������ �������������;
    ���������� ������������� � ������������ ����������(������ ������ � ������� � ����������� �������������, ������� � ����������� ������ � 2 - 3 ����, ��� � �������);
    ���������� ������������� � ������������� ����������
*/
int main()
{
    float min = 0.0;
    float max = 1.0;
    float sum_y = 0;
    float yArithmeticMean, M, D, b_3 , b_4;
    vector<float> u1;
    vector<float> u2;
    vector<float> x;
    vector<float> y;


    u1.resize(N);
    u2.resize(N);
    x.resize(N);
    y.resize(N);

    setlocale(LC_ALL, "Russian");

    //��������� U1  � U2
    for (int i = 0; i < N; i++)
    {
        u1[i] = randomRangeU(min, max);
        u2[i] = randomRangeU(min, max);
        std::cout << i << '|' << u1[i] << '|' << u2[i] << endl;
    }

    // ������ ��������� �������� �
    for (int i = 0; i < N; i++)
    {
        x[i] = GeneratingRandomVariable(u1[i], u2[i]);
    }

    std::cout << "������ y:" << endl;
    //������ ��������� �������� y  � y`- ������� y 
    for (int i = 0; i < N; i++)
    {
        y[i] = x[i] * lambda + tetta;
        std::cout << i << '|' << y[i] << endl;
        sum_y += y[i];
    }

    /// ������� �������������� 
    yArithmeticMean = sum_y / N;
    std::cout << " ������� ��������������  y:" << yArithmeticMean<< endl;
    //�������
     M = median(y);
     std::cout << " �������  y:" <<M<< endl;
    //��������� 
     D = dispersion(yArithmeticMean, y);
     std::cout << "��������� y:" << D << endl;
    //������������� ���������� 
     b_3 = betta_3(yArithmeticMean, D, y);
     std::cout << "����������� ����������:" << b_3 << endl;
     //��������;
     b_4 = betta_4(yArithmeticMean, D, y);
     std::cout << "����������� ��������:" << b_4 << endl;

    return 0;
}