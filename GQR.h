/* Kurkov Ivan, 22.B05-MM, 14.04.2024 */
#ifndef GQR_H
#define GQR_H

#include <functional>

struct GQR
{
  double *Nodes = nullptr, *Coeffs = nullptr;
  size_t Size = 0;

  GQR( void ) = default;
  GQR( const double *x, size_t n);
  GQR( const GQR &Rule );
  GQR & operator=( const GQR &R );
  GQR( GQR &&R );
  GQR & operator=( GQR &&R );
  ~GQR( void );

  void swap( GQR &R );

  double Integrate( std::function<double(double)> f, double left, double right );
};

#endif // !GQR_H

