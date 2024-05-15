/* Kurkov Ivan, 22.B05-MM, 03.03.2024 */
#include <vector>
#include <functional>

#include "nleq.h"

double Secant( std::function<double(double)> f, double prev, double cur, double eps )
{
  double next;

  while (fabs(cur - prev) >= eps)
  {
    next = cur - f(cur) / (f(cur) - f(prev)) * (cur - prev);
    prev = cur;
    cur = next;
  }
  return cur;
}

std::vector<std::pair<double, double>>
  RootSeparation( std::function<double(double)> f, double left, double right, size_t N )
{
  std::vector<std::pair<double, double>> segments;
  double h = (right - left) / N, a = left, b = a + h;

  for (size_t i = 0; i < N; i++)
  {
    if (f(a) * f(b) < 0)
      segments.push_back({a, b});
    if (f(a) == 0)
      segments.push_back({a - h, a + h});
    a += h;
    b += h;
  }
  if (f(a) == 0)
    segments.push_back({a - h, a + h});
  return segments;
}