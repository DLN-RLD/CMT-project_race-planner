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

double downhill_pace_factor(double i) {
    /* Empirical downhill braking curve (min 0.88 at -9%, 1.0 at -18%, braking beyond) */
    double d = i + 0.09;
    return 0.88 + (0.12 / (0.09 * 0.09)) * d * d;
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

double peronnet_I_max(double t_s, double endurance_index) {
    if (t_s <= PERONNET_T_MAP) {
        return 1.0;
    }
    double val = 1.0 + (endurance_index / 100.0) * log(t_s / PERONNET_T_MAP);
    if (val > 1.0) {
        return 1.0;
    }
    if (val < 0.05) {
        return 0.05;
    }
    return val;
}


/* Simulates the full race at a given relative intensity I.
 * Returns final glycogen [g] and optionally writes total time [s] to *t_tot_out.
 */
void simulate_course(Runner r, const Segment course[], int n_segments,
                     double I_target, double *t_tot, double *G_final,
                     double *G_min, ProfileArrays *prof) {
    double P_max = compute_Pmax(r);
    double P_target = P_max * I_target;
    double v_trans = compute_v_trans(r.height);
    double eta = compute_carb_efficiency(r.Rin_h);
    double Rin_s = r.Rin_h / 3600.0; // Convert g
    double P_bmr = BMR_PER_KG * r.mass; // Resting metabolic rate [W]

    double cum_dist = 0.0;
    double cum_time = 0.0;
    double G_current = r.G0;
    double g_min_val = r.G0;
    
    for (int k = 0; k < n_segments; k++) {
        double d_k = course[k].distance;
        double i_k = course[k].slope;
        if (i_k > MAX_SLOPE) {
            i_k = MAX_SLOPE;
        } else if (i_k < -MAX_SLOPE) {
            i_k = -MAX_SLOPE;
        }

        /* Progressive fatigue drift evaluated at segment midpoint, capped at +15% */
        double dist_mid_km = (cum_dist + 0.5 * d_k) / 1000.0;
        double fatigue_add = r.fatigue_rate * dist_mid_km;
        if (fatigue_add > FATIGUE_CAP) {
            fatigue_add = FATIGUE_CAP;
        }
        double cost_mult = r.terrain_factor * (1.0 + fatigue_add);

        double c_run = minetti_run(i_k) * cost_mult;
        double c_walk = minetti_walk(i_k) * cost_mult;
        double c_flat = minetti_run(0.0) * cost_mult;

        /* 1. Metabolic speed & uphill walk-run transition*/
        double v_k = P_target / (c_run * r.mass);
        double is_walking = 0.0;
        double c_eff = c_run;

        if (i_k > MIN_WALK_SLOPE && v_k <= v_trans) {
            v_k = P_target / (c_walk * r.mass);
            is_walking = 1.0;
            c_eff = c_walk;
        }

        /* 2. Downhill biomechanical braking limit */
        if (i_k < 0.0) {
            double v_flat_I = P_target / (c_flat * r.mass);
            double v_downhill_max = v_flat_I * downhill_pace_factor(i_k);
            if (v_k > v_downhill_max) {
                v_k = v_downhill_max;
            }
        }
        
        /* 3. Local metabolic power and glycogen balance (Brooks crossover I^2) */
        double P_net = c_eff * r.mass * v_k;
        double I_k = P_net / P_max;
        if (I_k > 1.0) {
            I_k = 1.0;
        } else if (I_k < 0.0) {
            I_k = 0.0;
        }

        double P_tot = P_net + P_bmr;
        double f_cho = pow(I_k, CHO_FRAC_EXP); // I^1,5 instead of I^2

        double dt = d_k / v_k;
        double Rout_s = (P_tot * f_cho) / GLYCOGEN_ENERGY; // repalce I^2 by f_cho

        cum_dist += d_k;
        cum_time += dt;
        G_current += (eta * Rin_s - Rout_s) * dt;

        if (G_current > r.G0) {
            G_current = r.G0;
        }
        if (G_current < g_min_val) {
            g_min_val = G_current;
        }

        if (prof != NULL) {
            prof->cum_dist_m[k] = cum_dist;
            prof->cum_time_s[k] = cum_time;
            prof->speed_ms[k] = v_k;
            prof->is_walking[k] = is_walking;
            prof->local_intensity[k] = I_k;
            prof->glycogen_g[k] = G_current;
        }   
    }

    *t_tot = cum_time;
    *G_final = G_current;
    *G_min = g_min_val;
}


/* Finds the optimal intensity I via binary search so that G_final = 10% of G0 */
Result find_I_opt(Runner r, const Segment course[], int n_segments) {
    const double tol     = 1e-5;
    const int max_iter   = 100;
    const double G_floor = WALL_THRESHOLD * r.G0;

    double t_tmp = 0.0;
    double g_fin_tmp = 0.0;
    double g_min_tmp = 0.0;

    /* --- 1. Bisection for glycogen constraint (G_min >= 10% of G0) --- */
    double I_glyc = 1.0;
    simulate_course(r, course, n_segments, 1.0, &t_tmp, &g_fin_tmp, &g_min_tmp, NULL);
    if (g_min_tmp < G_floor) {
        double low = 0.05;
        double high = 1.0;
        for (int iter = 0; iter < max_iter && (high - low) > tol; iter++) {
            double mid = 0.5 * (low + high);
            simulate_course(r, course, n_segments, mid, &t_tmp, &g_fin_tmp, &g_min_tmp, NULL);
            if (g_min_tmp > G_floor) {
                low = mid;
            } else {
                high = mid;
            }
        }
        I_glyc = 0.5 * (low + high);
    }

    /* --- 2. Bisection for Peronnet-Thibault aerobic endurance limit --- */
    double I_aero = 1.0;
    simulate_course(r, course, n_segments, 1.0, &t_tmp, &g_fin_tmp, &g_min_tmp, NULL);
    if (1.0 > peronnet_I_max(t_tmp, r.endurance_index)) {
        /* Lower bound at |E|/100 guarantees strict monotonicity of I - I_max(t(I)) */
        double low = fabs(r.endurance_index) / 100.0;
        if (low < 0.05) {
            low = 0.05;
        }
        if (low > 0.99) {
            low = 0.99;
        }
        double high = 1.0;

        simulate_course(r, course, n_segments, low, &t_tmp, &g_fin_tmp, &g_min_tmp, NULL);
        if (low >= peronnet_I_max(t_tmp, r.endurance_index)) {
            I_aero = low;
        } else {
            for (int iter = 0; iter < max_iter && (high - low) > tol; iter++) {
                double mid = 0.5 * (low + high);
                simulate_course(r, course, n_segments, mid, &t_tmp, &g_fin_tmp, &g_min_tmp, NULL);
                if (mid < peronnet_I_max(t_tmp, r.endurance_index)) {
                    low = mid;
                } else {
                    high = mid;
                }
            }
            I_aero = 0.5 * (low + high);
        }
    }

    /* --- 3. Select dominant constraint & run final evaluation --- */
    Result res;
    res.I_glyc = I_glyc;
    res.I_aero = I_aero;

    if (I_glyc < I_aero) {
        res.I_opt   = I_glyc;
        res.limiter = LIMIT_GLYCOGEN;
    } else {
        res.I_opt   = I_aero;
        res.limiter = LIMIT_AEROBIC;
    }

    simulate_course(r, course, n_segments, res.I_opt,
                    &res.t_tot, &res.G_final, &res.G_min, NULL);

    return res;
}