/* Kurkov Ivan, 22.B05-MM, 08.05.2024 */
#include <iostream>
#include <conio.h>
#include <vector>
#include <cmath>
#include <algorithm>

#include "gmqr.h"
#include "nleq.h"
#include "gqr.h"
#include "fort.hpp"
#include "int.h"

const size_t MAX_N = 9, STEPS_N = 100, ITER = 10000000;
const double EPS = 1e-12;

std::string Variant7Str = "(x + 0.8) / sqrt(x^2 + 1.2)";
double Variant7( double x)
{
  return (x + 0.8) / sqrt(x * x + 1.2);
}

int main( void )
{
  int key;
  size_t N[3];
  double roots[MAX_N], I[3], a = 1.6, b = 2.7, q, p;
  bool run = true;
  GQR Gauss[MAX_N];
  AcmFVals acm;
  std::vector<std::pair<double, double>> sgms;
  auto test_poly = [&N]( double x ) { return pow(x, 2 * N[0] - 1) + 30; };
  auto test_2 = [&N](double x) { return cos(x); };
  auto test_weight = [&N]( double x ) { return cos(x) / sqrt(1 - x * x); };

  do
  {
    std::cout << "Choose quadrature rule (g - Gaussian, m - Melers)\n";
    key = _getch();
  } while (key != 'g' && key != 'm');

  if (key == 'g')
  {
    for (size_t i = 0; i < MAX_N; i++)
    {
      auto Poly = [&i](double x) { return LegendrePoly(x, i); };
      fort::char_table Gauss_params;

      Gauss_params << fort::header << "N = " + std::to_string(i) << "x_i" << "A_i" << fort::endr;

      sgms = RootSeparation(Poly, -1, 1, STEPS_N);
      for (size_t i = 0; i < sgms.size(); i++)
        roots[i] = Secant(Poly, sgms[i].first, sgms[i].second, EPS);
      Gauss[i] = GQR(roots, i);
      
      for (size_t j = 0; j < i; j++)
        Gauss_params << j + 1 << Gauss[i].Nodes[j] << Gauss[i].Coeffs[j] << fort::endr;
      std::cout << Gauss_params.to_string();
    }

    do
    {
      std::cout << "\nInput number of nodes < " << MAX_N << " for which to check Gaussian QR: ";
      std::cin >> N[0];
    } while (N[0] >= MAX_N);
    
    acm = Accumulate(test_poly, ITER, -1, 1);
    I[0] = QRMiddleRectangles(acm);
    I[1] = Gauss[N[0] - 1].Integrate(test_poly, -1, 1);
    std::cout << "I_[ * ](x^" << 2 * N[0] - 1 << " + 30, -1, 1) = " << std::setprecision(15) << I[0] << '\n';
    std::cout << "I_Gauss(x^" << 2 * N[0] - 1 << " + 30, -1, 1) = " << std::setprecision(15) << I[1] << '\n';
    std::cout << "|I_[ * ] - I_Gauss| = " << std::setprecision(15) << fabs(I[0] - I[1]) << '\n';

    while (run)
    {
      fort::char_table integrals;

      std::cout << "\nGaussian quadrature rule menu [variant #7]:\n"
        "0 - exit\n"
        "1 - calculate integral\n";
      switch (_getch())
      {
      case '0':
        run = false;
        break;
      case '1':
        do
        {
          std::cout << "Input [a, b]: ";
          std::cin >> a >> b;
        } while (a >= b);
        q = (b - a) / 2;
        p = (b + a) / 2;
        do
        {
          std::cout << "Input 3 number of nodes < " << MAX_N << " for which to check Gaussian QR: ";
          std::cin >> N[0] >> N[1] >> N[2];
        } while (N[0] >= MAX_N || N[1] >= MAX_N || N[2] >= MAX_N);
        std::sort(N, N + 3);
        for (size_t i = 0; i < 3; i++)
        {
          fort::char_table similar_params;

          similar_params << fort::header << "N = " + std::to_string(N[i]) << "x_i" << "A_i" << fort::endr;
          I[i] = Gauss[N[i]].Integrate(Variant7, a, b);
          for (size_t j = 0; j < N[i]; j++)
            similar_params << j + 1 << p + q * Gauss[N[i]].Nodes[j] << q * Gauss[N[i]].Coeffs[j] << fort::endr;
          std::cout << similar_params.to_string();
        }
        integrals << fort::header << '#' << N[0] << N[1] << N[2] << fort::endr;
        integrals << "I(" + Variant7Str + ", " + std::to_string(a) + ", " + std::to_string(b) + ')'
          << std::setprecision(15) << I[0] << I[1] << I[2] << fort::endr;
        integrals.column(0).set_cell_text_align(fort::text_align::right);
        std::cout << integrals.to_string();
        std::cout << "Delta I: " << fabs(I[1] - I[0]) << " -> " << fabs(I[2] - I[1]) << "\n";
        break;
      default:
        std::cout << "[Error]: Incorrect choice!";
        break;
      }
    }
  }
  else
  {
    while (run)
    {
      std::cout << "\nMelers quadrature rule menu [variant #7]:\n"
        "0 - exit\n"
        "1 - calculate integral\n";
      switch (_getch())
      {
      case '0':
        run = false;
        break;
      case '1':
        std::cout << "\nInput number of nodes for which to check Melers QR: ";
        std::cin >> N[0];

        acm = Accumulate(test_weight, ITER, -1, 1);
        I[0] = QRMiddleRectangles(acm);
        I[1] = QRMeller(test_2, N[0]);
        std::cout << "I_[ * ](cos(x) / sqrt(1 - x^2), -1, 1) = " << std::setprecision(15) << I[0] << '\n';
        std::cout << "I_Meler(cos(x) / sqrt(1 - x^2), -1, 1) = " << std::setprecision(15) << I[1] << '\n';
        std::cout << "|I_[ * ] - I_Meler| = " << std::setprecision(15) << fabs(I[0] - I[1]) << '\n';
        break;
      default:
        std::cout << "[Error]: Incorrect choice!";
        break;
      }
    }
  }
  return 0;
}