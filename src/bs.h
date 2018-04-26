double BS_call_cpp(
    const double V, const double D, const double T, const double r,
    const double std);

double BS_call_cpp_inv(
    const double S, const double D, const double T, const double r,
    const double std, const double tol,
    double V_min, double V_max, double V_mid);
