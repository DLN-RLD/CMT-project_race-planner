#include "model.h"
#include <stdio.h>
#include <math.h>

/* Minetti (2002) energy cost of running [J / (kg*m)] */
double minetti_run(double slope) {
    return 155.4 * pow(slope, 5) - 30.4 * pow(slope, 4) - 43.3 * pow(slope, 3) 
         + 46.3 * pow(slope, 2) + 19.5 * slope + 3.6;
}

/* Minetti (2002) energy cost of walking [J / (kg*m)] */
double minetti_walk(double slope) {
    return 280.5 * pow(slope, 5) - 58.7 * pow(slope, 4) - 76.8 * pow(slope, 3) 
         + 51.9 * pow(slope, 2) + 19.6 * slope + 2.5;
}

/* Computes net maximal aerobic power P_max [W = J/s].
 * Prioritizes VMA if provided; otherwise estimates VMA from VO2max (Leger & Mercier, 1984).
 */
double compute_Pmax(Runner r) {
    double pmax; 
    double vma_kmh = 0.0;

    if (r.vma > 0.0) {
        vma_kmh = r.vma;
    } else if (r.vo2max > 0.0) {
        vma_kmh = r.vo2max / 3.5;
    } else {
        fprintf(stderr, "Error: provide a valid VMA (> 0) or VO2max (> 0).\n");
        return -1.0;
    }

    // Convert to m/s and compute P_max
    double vma_ms = vma_kmh / 3.6;
    pmax = minetti_run(0.0) * r.mass * vma_ms;
    return pmax;
}

/* Computes individual run-to-walk transition speed [m/s] via the Froude number */
double compute_v_trans(double height) {
    double leg_length = LEG_LENGTH_RATIO * height;
    return sqrt(FROUDE_UPHILL * G_GRAVITY * leg_length);
}

/* Computes exogenous carbohydrate oxidation efficiency (eta)
 * from intake rate Rin_h [g/h] (calibrated on Podlogar et al., 2022; Hearris et al., 2022).
 */
double compute_carb_efficiency(double Rin_h) {
    if (Rin_h <= 60.0) {
        // Regime 1: moderate intake (<= 60 g/h), constant 95% efficiency
        return 0.95;
    } else if (Rin_h <= 120.0) {
        // Regime 2: linear decrease from 95% (at 60 g/h) to 75% (at 120 g/h)
        // y = a * x + b passing through (60, 0.95) and (120, 0.75)
        double a = (0.75 - 0.95) / (120.0 - 60.0);
        double b = 0.95 - a * 60.0;
        return a * Rin_h + b;
    } else {
        // Regime 3: intestinal saturation above 120 g/h, effective uptake capped at 90 g/h
        return 90.0 / Rin_h;
    }
}

/* Computes total race time [s] for a given target power P_target [W] */
double compute_total_time(double P_target, double mass, double v_trans,
                          const Segment course[], int n_segments) {
    double t_tot = 0.0;

    for (int k = 0; k < n_segments; k++) {
        double i_k = course[k].slope;

        // Minetti polynomial validity domain [-45%, +45%]
        if (i_k > MAX_SLOPE) {
            i_k = MAX_SLOPE;
        } else if (i_k < -MAX_SLOPE) {
            i_k = -MAX_SLOPE;
        }

        // 1. Theoretical running speed [m/s]
        double v_k = P_target / (minetti_run(i_k) * mass);

        // 2. Dynamic transition to walking on uphills below transition speed
        if (i_k > MIN_WALK_SLOPE && v_k <= v_trans) {
            v_k = P_target / (minetti_walk(i_k) * mass);
        }

        // 3. Accumulate segment time [s]
        t_tot += course[k].distance / v_k;
    }

    return t_tot;
}

/* Computes remaining glycogen stock at the finish line [g] */
double compute_carb_balance(double I, double P_target, double t_tot,
                            double G0, double Rin_h) {
    double Gfinal = G0;
    // 1. Exogenous carbohydrate inflow [g]
    double eta = compute_carb_efficiency(Rin_h);
    double Rin_s = Rin_h / 3600.0; // Convert g/h to g/s
    double total_intake = eta * Rin_s * t_tot;

    // 2. Carbohydrate oxidation outflow using Brooks crossover f(I) = I^2 [g]
    double P_carb = P_target * I * I;              // [W = J/s]
    double Rout_s = P_carb / GLYCOGEN_ENERGY;      // [g/s]
    double total_expenditure = Rout_s * t_tot;

    // 3. Final glycogen balance [g]
    Gfinal = G0 + total_intake - total_expenditure;

    return Gfinal;
}

/* Simulates the full race at a given relative intensity I.
 * Returns final glycogen [g] and optionally writes total time [s] to *t_tot_out.
 */
double simulate_race(double I, Runner r, const Segment course[],
                     int n_segments, double *t_tot_out) {
    double P_max = compute_Pmax(r);
    if (P_max <= 0.0) {
        return -1.0;
    }

    double P_target = P_max * I;
    double v_trans = compute_v_trans(r.height);

    double t_tot = compute_total_time(P_target, r.mass, v_trans, course, n_segments);
    if (t_tot_out != NULL) {
        *t_tot_out = t_tot;
    }

    double Gfinal = compute_carb_balance(I, P_target, t_tot, r.G0, r.Rin_h);
    return Gfinal;
}

/* Finds the optimal intensity I via binary search so that G_final = 10% of G0 */
Result find_I_opt(Runner r, const Segment course[], int n_segments) {
    Result res = {0.0, 0.0, 0.0};
    double G_target = WALL_THRESHOLD * r.G0;

    double I_low = 0.05;     // Lower bound (5% of P_max, avoids division by zero)
    double I_high = 1.00;    // Upper bound (100% of P_max)
    double tol = 1e-5;       // Convergence tolerance on I
    int max_iter = 100;      // Safety limit

    double t_test = 0.0;

    // Edge case 1: race is short enough that full intensity stays above the wall
    double G_max_effort = simulate_race(I_high, r, course, n_segments, &t_test);
    if (G_max_effort >= G_target) {
        res.I_opt = I_high;
        res.t_tot = t_test;
        res.G_final = G_max_effort;
        return res;
    }

    // Edge case 2: even at minimum intensity, the runner hits the wall
    double G_min_effort = simulate_race(I_low, r, course, n_segments, &t_test);
    if (G_min_effort <= G_target) {
        res.I_opt = I_low;
        res.t_tot = t_test;
        res.G_final = G_min_effort;
        return res;
    }

    // Main binary search loop
    double I_mid = 0.5;
    double G_mid = 0.0;
    int iter = 0;

    while ((I_high - I_low) > tol && iter < max_iter) {
        I_mid = 0.5 * (I_low + I_high);
        G_mid = simulate_race(I_mid, r, course, n_segments, &t_test);

        if (G_mid > G_target) {
            I_low = I_mid;   // Surplus glycogen: runner can push harder
        } else {
            I_high = I_mid;  // Below wall threshold: runner must slow down
        }
        iter++;
    }

    res.I_opt = I_mid;
    res.t_tot = t_test;
    res.G_final = G_mid;
    return res;
}