#include <cmath>
#include <rose499/spline.hpp>

using namespace Eigen;

Spline::Spline(MatrixXd points)
  : mSplineCount(points.cols() - 1),
    mPoly(2 * mSplineCount, PolyOrder + 1),
    mDPoly(2 * mSplineCount, PolyOrder + 1),
    mDDPoly(2 * mSplineCount, PolyOrder + 1) { }

Matrix<Spline::ValueType, 2, 2> Spline::frame(ValueType parameter) const
{
    Matrix<ValueType, 2, 2> basis;
    Matrix<ValueType, 2, 1> tangent = ((*this)(parameter, 1)).matrix();

    // Form a right-handed frame from the tangent that spans R^2
    basis.col(0) = tangent / tangent.norm();
    basis.col(1) = tangent.reverse();
    basis(0, 1) = -basis(0, 1);

    return basis;
}

Array<Spline::ValueType, 2, 1> Spline::operator() (ValueType parameter, uint32_t derivative) const
{
    using Point = Array<Spline::ValueType, 2, 1>;

    // Parameter out of bounds
    if( parameter < 0 || parameter > 1 )
        throw new InvalidParameterException();

    // Identically zero derivative
    if( derivative > PolyOrder )
        return Point::Zero();

    // Calculate x^7, x^6, ..., x^0
    ArrayXd value(PolyOrder+1);
    ArrayXd powers(PolyOrder+1);
    for(int pi = PolyOrder; pi >= 0; --pi)
    {
        powers[PolyOrder - pi] = pi;
    }

    value.fill(parameter);
    value = pow(value,powers);

    // Choose spline based on parameter value.
    int indSpline = std::min(mSplineCount, (int)std::floor(mSplineCount * parameter));
    auto coeff = mPoly.block<2, PolyOrder + 1>(2 * indSpline, 0);

    // And evaluate!
    switch(derivative)
    {
    case 0:
        return (coeff * value.matrix()).array();

    default:
        throw new InvalidParameterException();
    };
}
