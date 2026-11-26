// threshold.cpp -- certified enclosures of m(3), m(23), m(24), m(25) (Section 5 of the paper)
//
//   m(y) = y(y-1) int_1^oo (-theta') (theta + 2t theta') [ (2 ttheta)^{y-2} + t (sqrt(t) theta - 1)^{y-2} ] dt
//
// Interval Riemann sums on [1, a] with step h = 2^-p, plus the tail bound T_a(y) on [a, oo).
// Ball arithmetic: Arb (FLINT >= 3.0).
//
//   g++ -O2 threshold.cpp -o threshold -lflint && ./threshold

#include <flint/arb.h>
#include <cstdio>
#include <cstdlib>

const slong PREC = 128;   // working precision (bits)
const long K = 6;         // terms of the theta series
const long A = 20;        // splitting point a > max{1, y/(2 pi)}

// A real ball with the usual arithmetic operators.
struct R {
    arb_t v;
    R() { arb_init(v); }
    R(long n) { arb_init(v); arb_set_si(v, n); }
    R(const R& x) { arb_init(v); arb_set(v, x.v); }
    R& operator=(const R& x) { arb_set(v, x.v); return *this; }
    ~R() { arb_clear(v); }
};
R operator+(const R& a, const R& b) { R r; arb_add(r.v, a.v, b.v, PREC); return r; }
R operator-(const R& a, const R& b) { R r; arb_sub(r.v, a.v, b.v, PREC); return r; }
R operator*(const R& a, const R& b) { R r; arb_mul(r.v, a.v, b.v, PREC); return r; }
R operator/(const R& a, const R& b) { R r; arb_div(r.v, a.v, b.v, PREC); return r; }
R operator-(const R& a) { R r; arb_neg(r.v, a.v); return r; }
R pow(const R& a, long n) { R r; arb_pow_ui(r.v, a.v, n, PREC); return r; }
R exp(const R& a) { R r; arb_exp(r.v, a.v, PREC); return r; }
R sqrt(const R& a) { R r; arb_sqrt(r.v, a.v, PREC); return r; }
R pi() { R r; arb_const_pi(r.v, PREC); return r; }

// Integrand F_y(t), for every t in the ball t.
// Series truncated at K; both remainders are <= M q^M / (1-q)^2 with M = (K+1)^2.
R F(const R& t, long y)
{
    R q = exp(-pi() * t), tt, s;                   // tt = ttheta, s = sum k^2 q^{k^2}
    for (long k = 1; k <= K; k++) {
        R qk = pow(q, k * k);
        tt = tt + qk;
        s = s + k * k * qk;
    }
    long M = (K + 1) * (K + 1);
    R rem = M * pow(q, M) / ((1 - q) * (1 - q));
    arb_add_error(tt.v, rem.v);
    arb_add_error(s.v, rem.v);

    R mdt = 2 * pi() * s;                          // -theta'
    R th = 1 + 2 * tt;                             // theta
    return y * (y - 1) * mdt * (th - 2 * t * mdt)
           * (pow(2 * tt, y - 2) + t * pow(sqrt(t) * th - 1, y - 2));
}

// Tail bound T_a(y), with c = 1/(1 - e^{-pi}).
R tail(long y, long a)
{
    R c = 1 / (1 - exp(-pi())), e = exp(-pi() * a);
    R lambda = pi() - R(y) / R(2 * a);
    if (!arb_is_positive(lambda.v)) { printf("need a > y/(2 pi)\n"); exit(1); }
    return 2 * pi() * c * c * y * (y - 1)
           * (pow(2 * c, y - 2) * pow(e, y - 1) / (pi() * (y - 1))
              + pow(1 + 2 * c * e, y - 2) * pow(sqrt(R(a)), y) * e / lambda);
}

// m(y) in h * sum_j F_y(I_j) + [-T, T],  I_j = [1 + jh, 1 + (j+1)h].
R m(long y, long a, slong p)
{
    R h(1), unit, sum;
    arb_mul_2exp_si(h.v, h.v, -p);                 // h = 2^-p
    arb_set_str(unit.v, "[0.5 +/- 0.5]", PREC);    // [0, 1]
    for (long j = 0; j < ((a - 1) << p); j++)
        sum = sum + F(1 + (j + unit) * h, y);
    R result = h * sum, T = tail(y, a);
    arb_add_error(result.v, T.v);
    return result;
}

int main(int argc, char** argv)
{
    slong p = (argc > 1) ? atol(argv[1]) : 10;
    printf("FLINT %s, %ld bits, K = %ld, a = %ld, h = 2^-%ld\n\n", flint_version, (long)PREC, K, A, (long)p);
    const long ys[] = {3, 23, 24, 25};
    for (long y : ys) {
        R r = m(y, A, p), one(1);
        printf("m(%2ld) in ", y);
        arb_printd(r.v, 10);
        printf("   %s\n", arb_lt(r.v, one.v) ? "< 1" : arb_gt(r.v, one.v) ? "> 1" : "undecided: increase p");
    }
    return 0;
}
