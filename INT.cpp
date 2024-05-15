/* Kurkov Ivan, 22.B05-MM, 14.04.2024 */
#include <functional>

#include "int.h"

const QRUnit QRList[] = 
{
  {"[*  ]", QRLeftRectangles},
  {"[  *]", QRRightRectangles},
  {"[ * ]", QRMiddleRectangles},
  {"/Trapeze\\", QRTrapeze},
  {"Simpson", QRSimpson}
};
const size_t NumOfQR = Length(QRList);

AcmFVals Accumulate( std::function<double(double)> f, size_t n, double left, double right )
{
  double h = (right - left) / n, odd_offset = left + h / 2, even_offset = left + h;
  AcmFVals acm = {f(left), f(right), 0, f(odd_offset), h};

  odd_offset += h;
  for (size_t i = 1; i < n; i++)
  {
    acm.odd += f(odd_offset);
    acm.even += f(even_offset);
    odd_offset += h;
    even_offset += h;
  }
  return acm;
}

double QRLeftRectangles( AcmFVals acm )
{
  return acm.step * (acm.begin + acm.even);
}

double QRRightRectangles( AcmFVals acm )
{
  return acm.step * (acm.end + acm.even);
}

double QRMiddleRectangles( AcmFVals acm )
{
  return acm.step * acm.odd;
}

double QRTrapeze( AcmFVals acm )
{
  return acm.step * (acm.even + (acm.begin + acm.end) / 2);
}

double QRSimpson( AcmFVals acm )
{
  return acm.step * (acm.begin + acm.end + 2 * acm.even + 4 * acm.odd) / 6;
}

double QRThreeEighths( AcmFVals acm )
{
  return 3 * acm.step * (acm.begin + acm.end + 3 * acm.even) / 8;
}