/*
  rk4.h / rk4.c (single-file example)

  Classical 4th-order Runge–Kutta (RK4) for first-order ODE systems:
      y' = f(t, y; params), y in R^n

  Usage:
    - Define an RHS function matching RHSFunc.
    - Call rk4_integrate(...) or rk4_step(...).

  Notes:
    - Fixed step size h (non-adaptive).
    - Workspace arrays are required (k1,k2,k3,k4,ytmp), each length n.
*/

#include <stddef.h>   // size_t
#include <math.h>     // fabs (optional)
#include <stdio.h>    // printf, fflush

/* RHS function signature: write dydt = f(t, y, params) */
typedef void (*RHSFunc)(double t, const double *y, double *dydt, void *params);

static void rk4_print_final(double t, const double *y, size_t n)
{
    printf("t_final = %.17g\n", t);
    printf("y_final = [");
    for (size_t i = 0; i < n; ++i) {
        printf("%.17g%s", y[i], (i + 1 < n) ? ", " : "");
    }
    printf("]\n");
    fflush(stdout);
}

/* One RK4 step: advances y(t) -> yout(t+h) */
void rk4_step(RHSFunc f,
              double t,
              const double *y,
              double h,
              size_t n,
              void *params,
              double *yout,
              double *k1, double *k2, double *k3, double *k4,
              double *ytmp)
{
    size_t i;

    /* k1 = f(t, y) */
    f(t, y, k1, params);

    /* k2 = f(t + h/2, y + (h/2) k1) */
    for (i = 0; i < n; ++i) ytmp[i] = y[i] + 0.5 * h * k1[i];
    f(t + 0.5 * h, ytmp, k2, params);

    /* k3 = f(t + h/2, y + (h/2) k2) */
    for (i = 0; i < n; ++i) ytmp[i] = y[i] + 0.5 * h * k2[i];
    f(t + 0.5 * h, ytmp, k3, params);

    /* k4 = f(t + h, y + h k3) */
    for (i = 0; i < n; ++i) ytmp[i] = y[i] + h * k3[i];
    f(t + h, ytmp, k4, params);

    /* yout = y + (h/6)(k1 + 2k2 + 2k3 + k4) */
    for (i = 0; i < n; ++i) {
        yout[i] = y[i] + (h / 6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    }
}

/*
  Integrate from t0 to t1 with fixed step h.
  - y is updated in-place to the final value at t1 (approximately).
  - If (t1 - t0) is not an integer multiple of h, the last step is shortened.
  Returns: number of steps taken (including the final shortened step if any).
*/
size_t rk4_integrate(RHSFunc f,
                     double t0,
                     double t1,
                     double h,
                     double *y,
                     size_t n,
                     void *params,
                     double *k1, double *k2, double *k3, double *k4,
                     double *ytmp)
{
    size_t steps = 0;
    double t = t0;

    if (h == 0.0) return 0;

    /* Ensure we step in the correct direction */
    if ((t1 > t0 && h < 0.0) || (t1 < t0 && h > 0.0)) {
        h = -h;
    }

    while ((h > 0.0 && t < t1) || (h < 0.0 && t > t1)) {
        double h_step = h;

        /* Shorten final step to land exactly on t1 (within floating limits) */
        if (h > 0.0 && t + h_step > t1) h_step = t1 - t;
        if (h < 0.0 && t + h_step < t1) h_step = t1 - t;

        /* Take one step: y -> y_next stored in ytmp */
        rk4_step(f, t, y, h_step, n, params, ytmp, k1, k2, k3, k4, k4);
        /* Copy ytmp (y_next) back into y */
        for (size_t i = 0; i < n; ++i) y[i] = ytmp[i];

        t += h_step;
        steps++;
        if (h_step == 0.0) break; /* safety */
    }

    rk4_print_final(t, y, n);
    return steps;
}

/* ------------------ Example RHS (optional) ------------------
   Simple harmonic oscillator as a first-order system:
     x' = v
     v' = -x
   y = [x, v]
*/
typedef struct {
    double omega;
} SHOParams;

void sho_rhs(double t, const double *y, double *dydt, void *params)
{
    (void)t; /* unused */
    SHOParams *p = (SHOParams*)params;
    double x = y[0];
    double v = y[1];
    dydt[0] = v;
    dydt[1] = -(p->omega * p->omega) * x;
}

/* If you want a standalone demo, compile with:
   cc -O2 rk4.c -lm -DDEMO && ./a.out
*/
#ifdef DEMO
#include <stdio.h>

int main(void)
{
    const size_t n = 2;
    double y[2] = {1.0, 0.0};   /* x(0)=1, v(0)=0 */
    SHOParams p = {.omega = 1.0};

    /* workspace */
    double k1[2], k2[2], k3[2], k4[2], ytmp[2];

    double t0 = 0.0, t1 = 10.0, h = 0.01;

    rk4_integrate(sho_rhs, t0, t1, h, y, n, &p, k1, k2, k3, k4, ytmp);

    printf("y(%.2f) = [x, v] = [%.15f, %.15f]\n", t1, y[0], y[1]);
    return 0;
}
#endif
