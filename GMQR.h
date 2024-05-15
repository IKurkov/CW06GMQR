/* Kurkov Ivan, 22.B05-MM, 08.05.2024 */
#ifndef GMQR_H
#define GMQR_H

#include <functional>

double LegendrePoly( double x, size_t n );
double QRMeller( std::function<double(double)> f, size_t N );
#endif // !GMQR_H

