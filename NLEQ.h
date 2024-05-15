/* Kurkov Ivan, 22.B05-MM, 03.03.2024 */
#ifndef NLEQ_H
#define NLEQ_H

#include <vector>
#include <functional>

double Secant( std::function<double(double)> f, double prev, double cur, double eps );
std::vector<std::pair<double, double>>
  RootSeparation(std::function<double(double)> f, double left, double right, size_t N );

#endif // !NLEQ_H