/* Kurkov Ivan, 22.B05-MM, 14.04.2024 */
#include <cstring>
#include <utility>
#include <functional>

#include "gqr.h"
#include "gmqr.h"

GQR::GQR( const double *x, size_t n ) : Nodes(new double[n]), Coeffs(new double[n]), Size(n)
{
  auto P = [n]( double x ){ return LegendrePoly(x, n - 1); };

  memcpy(Nodes, x, n * sizeof(double));
  for (size_t i = 0; i < n; i++)
  {
    double tmp = n * P(Nodes[i]);

    Coeffs[i] = 2 * (1 - Nodes[i] * Nodes[i]) / (tmp * tmp);
  }
}

void GQR::swap( GQR &R )
{
  std::swap(Nodes, R.Nodes);
  std::swap(Coeffs, R.Coeffs);
  std::swap(Size, R.Size);
}

GQR::GQR( const GQR &R ) : Nodes(new double[R.Size]), Coeffs(new double[R.Size]), Size(R.Size)
{
  memcpy(Nodes, R.Nodes, Size * sizeof(double));
  memcpy(Coeffs, R.Coeffs, Size * sizeof(double));
}

GQR & GQR::operator=( const GQR &R )
{
  if (this != &R)
    GQR(R).swap(*this);
  return *this;
}

GQR::GQR( GQR &&R ) { swap(R); }

GQR & GQR::operator=( GQR &&R )
{
  swap(R);
  return *this;
}

GQR::~GQR( void )
{
  delete[] Nodes;
  delete[] Coeffs;
  Size = 0;
}

double GQR::Integrate( std::function<double(double)> f, double left, double right )
{
  double sum = 0, q = (right - left) / 2, p = (right + left) / 2;

  for (size_t i = 0; i < Size; i++)
    sum += Coeffs[i] * f(p + q * Nodes[i]);
  return q * sum;
}