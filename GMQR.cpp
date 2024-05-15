/* Kurkov Ivan, 22.B05-MM, 08.05.2024 */
#define _USE_MATH_DEFINES

#include <vector>
#include <map>
#include <functional>
#include <cmath>

std::map<std::pair<size_t, double>, double> Cache;

double LegendrePoly( double x, size_t n )
{
  if (n == 0)
    return 1;
  if (n == 1)
    return x;

  if (Cache.contains({n, x}))
    return Cache[{n, x}];
  else
    return Cache[{n , x}] = (double)(2 * n - 1) / n * LegendrePoly(x, n - 1) * x
      - (double)(n - 1) / n * LegendrePoly(x, n - 2);
}

double QRMeller( std::function<double(double)> f, size_t N )
{
  double sum = 0;

  for (size_t i = 1; i <= N; i++)
    sum += f(cos(M_PI_2 * (2 * i - 1) / N));
  return sum * M_PI / N;
}