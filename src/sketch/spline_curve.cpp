#include "gector/sketch/spline_curve.h"
#include <stdexcept>

namespace gector {

SplineCurve::SplineCurve(NURBSCurve curve)
    : m_curve(std::move(curve)) {}

SplineCurve SplineCurve::interpolate(const std::vector<Point3D>& points, int degree) {
    return SplineCurve(NURBSCurve::makeInterpolated(points, degree));
}

Point3D SplineCurve::startPoint() const {
    return m_curve.evaluate(m_curve.paramStart());
}

Point3D SplineCurve::endPoint() const {
    return m_curve.evaluate(m_curve.paramEnd());
}

double SplineCurve::length() const {
    // Numerical approximation
    const int N = 64;
    double t0 = m_curve.paramStart();
    double t1 = m_curve.paramEnd();
    double len = 0.0;
    for (int i = 0; i < N; ++i) {
        double ta = t0 + (t1 - t0) *  i      / N;
        double tb = t0 + (t1 - t0) * (i + 1) / N;
        double tm = (ta + tb) * 0.5;
        double fa = m_curve.derivative(ta).length();
        double fm = m_curve.derivative(tm).length();
        double fb = m_curve.derivative(tb).length();
        len += (tb - ta) / 6.0 * (fa + 4.0*fm + fb);
    }
    return len;
}

std::shared_ptr<NURBSCurve> SplineCurve::toCurve() const {
    return std::make_shared<NURBSCurve>(m_curve);
}

} // namespace gector
