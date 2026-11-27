#include "model.h"
#include <math.h>
#include <stddef.h>

/* ========================================================================= */
/* Internal numeric tolerances & helper                                      */
/* ========================================================================= */
#define BISECT_ITERS   30
#define MIN_SPEED_MS   1e-6

static double clamp_slope(double slope) {
    if (slope >  MAX_SLOPE) return  MAX_SLOPE;
    if (slope < -MAX_SLOPE) return -MAX_SLOPE;
    return slope;
}

/* ========================================================================= */
/* Layer 1: Elemental published equations                                    */
/* ========================================================================= */

/* Gross maximal aerobic metabolic power [W/kg] */
double compute_M_max(Runner r) {
    if (r.vo2max > 0.0) {
        return (r.vo2max * O2_ENERGY_EQ) / 60.0;
    }
    if (r.vma > 0.0) {
        double vma_ms = r.vma / 3.6;
        return re3_run_power(vma_ms, 0.0);
    }
    return 0.0;
}

/* RE3 running gross metabolic power [W/kg] (Looney, Hoogkamer & Kram, 2025/2026) */
double re3_run_power(double speed_ms, double slope) {
    slope = clamp_slope(slope);

    double G_pct = 100.0 * slope;
    double grade_coeff = 30.43 * slope * (1.0 - pow(1.133, 1.0 - pow(1.056, G_pct + 43.0)));
    return RE3_INTERCEPT + (1.51 + grade_coeff) * speed_ms + 0.37 * speed_ms * speed_ms;
}

/* LCDA walking gross metabolic power [W/kg] (Looney et al., 2019) */
double lcda_walk_power(double speed_ms, double slope) {
    slope = clamp_slope(slope);

    double G_pct = 100.0 * slope;
    double v2 = speed_ms * speed_ms;
    double v4 = v2 * v2;
    double grade_term = 0.34 * speed_ms * G_pct * (1.0 - pow(1.05, 1.0 - pow(1.1, G_pct + 32.0)));
    return LCDA_STANDING_M + 1.94 * pow(speed_ms, 0.43) + 0.24 * v4 + grade_term;
}

/* Carbohydrate oxidation fraction f_c(I) [-] (Rapoport, 2010 / Romijn et al., 1993) */
double rapoport_cho_fraction(double I_gross) {
    /* Romijn et al. (1993) reference points at 25%, 65%, 85%, and 100% VO2max */
    double I1 = 0.25, F1 = 0.11;
    double I2 = 0.65, F2 = 0.52;
    double I3 = 0.85, F3 = 0.77;
    double I4 = 1.00, F4 = 1.00;

    if (I_gross <= 0.0) return 0.0;
    if (I_gross <= I1)  return F1 * (I_gross / I1);
    if (I_gross <= I2)  return F1 + (F2 - F1) * (I_gross - I1) / (I2 - I1);
    if (I_gross <= I3)  return F2 + (F3 - F2) * (I_gross - I2) / (I3 - I2);
    if (I_gross <  I4)  return F3 + (F4 - F3) * (I_gross - I3) / (I4 - I3);
    return 1.0;
}

/* Exogenous CHO oxidation rate [g/s] (Jeukendrup, 2014) */
double jeukendrup_exo_rate(double Rin_h, double Rmax_h) {
    if (Rin_h == 0.0) return 0.0;
    if (Rmax_h > 0.0 && Rin_h > Rmax_h) {
        return Rmax_h / 3600.0;
    }
    return Rin_h / 3600.0;
}

/* Aerobic intensity ceiling I_max(T) [-] (Peronnet & Thibault, 1989) */
double peronnet_I_max(double t_s, double endurance_index) {
    if (t_s <= PERONNET_T_MAP) return 1.0;
    double I_max = 1.0 + (endurance_index / 100.0) * log(t_s / PERONNET_T_MAP);
    if (I_max < 0.0) return 0.0;
    if (I_max > 1.0) return 1.0;
    return I_max;
}

/* ========================================================================= */
/* Layer 2: Local segment physics                                            */
/* ========================================================================= */

/* Numerical inversion of RE3 for running speed v_run [m/s] */
double re3_run_speed(double M_target, double slope) {
    double lo = 0.0;
    double hi = 10.0; /* 36 km/h */
    for (int iter = 0; iter < BISECT_ITERS; iter++) {
        double mid = 0.5 * (lo + hi);
        if (re3_run_power(mid, slope) < M_target) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return 0.5 * (lo + hi);
}

/* Numerical inversion of LCDA for walking speed v_walk [m/s] */
double lcda_walk_speed(double M_target, double slope, double v_low) {
    double lo = v_low;
    double hi = 4.0; /* 14.4 km/h */
    for (int iter = 0; iter < BISECT_ITERS; iter++) {
        double mid = 0.5 * (lo + hi);
        if (lcda_walk_power(mid, slope) < M_target) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return 0.5 * (lo + hi);
}

/* Solve locomotion mode, speed, and actual power on a single segment */
SegmentState solve_segment(Runner r, double slope, double M_target) {
    SegmentState st;

    /* 1. Running speed from RE3 */
    double v_run = re3_run_speed(M_target, slope);
    st.speed_ms   = v_run;
    st.M_actual   = M_target;
    st.is_walking = 0.0;

    /* 2. Continuous walk/run transition (EOTS: v_k = max(v_run, v_walk)) */
    if (lcda_walk_power(v_run, slope) < M_target) {
        st.speed_ms   = lcda_walk_speed(M_target, slope, v_run);
        st.is_walking = 1.0;
    }

    /* 3. Optional biomechanical downhill speed cap */
    if (slope < 0.0 && r.v_downhill_max > 0.0 && st.speed_ms > r.v_downhill_max) {
        st.speed_ms = r.v_downhill_max;
        if (st.is_walking > 0.5) {
            st.M_actual = lcda_walk_power(st.speed_ms, slope);
        } else {
            st.M_actual = re3_run_power(st.speed_ms, slope);
        }
    }

    return st;
}

/* ========================================================================= */
/* Layer 3: Course simulation & global optimization                          */
/* ========================================================================= */

void simulate_course(Runner r, const Segment course[], int n_segments,
                     double I_target, double *t_tot, double *G_final,
                     double *G_min, ProfileArrays *prof) {
    double M_max    = compute_M_max(r);
    double M_target = M_max * I_target;
    double R_exo    = jeukendrup_exo_rate(r.Rin_h, r.Rmax_h);

    double t_acc     = 0.0;
    double d_acc     = 0.0;
    double G_current = r.G0;
    double G_min_val = r.G0;

    for (int k = 0; k < n_segments; k++) {
        SegmentState st = solve_segment(r, course[k].slope, M_target);
        if (st.speed_ms < MIN_SPEED_MS) {
            st.speed_ms = MIN_SPEED_MS;
        }

        double dt    = course[k].distance / st.speed_ms;
        double I_k   = st.M_actual / M_max;
        double fc    = rapoport_cho_fraction(I_k);

        double R_out = (r.mass * st.M_actual * fc) / GLYCOGEN_ENERGY;

        t_acc     += dt;
        d_acc     += course[k].distance;
        G_current += (R_exo - R_out) * dt;

        if (G_current > r.G0)      G_current = r.G0;
        if (G_current < G_min_val) G_min_val = G_current;

        if (prof != NULL) {
            prof->cum_dist_m[k]      = d_acc;
            prof->cum_time_s[k]      = t_acc;
            prof->speed_ms[k]        = st.speed_ms;
            prof->is_walking[k]      = st.is_walking;
            prof->local_intensity[k] = I_k;
            prof->glycogen_g[k]      = G_current;
        }
    }

    if (t_tot   != NULL) *t_tot   = t_acc;
    if (G_final != NULL) *G_final = G_current;
    if (G_min   != NULL) *G_min   = G_min_val;
}

double solve_I_glyc(Runner r, const Segment course[], int n_segments, double M_max) {
    double I_min = (LCDA_STANDING_M / M_max) + 0.01;
    double G_floor   = WALL_THRESHOLD * r.G0;
    double g_min_tmp = 0.0;

    double high  = 1.0;
    double low   = I_min;

    for (int iter = 0; iter < BISECT_ITERS; iter++) {
        double mid = 0.5 * (low + high);
        simulate_course(r, course, n_segments, mid, NULL, NULL, &g_min_tmp, NULL);
        if (g_min_tmp >= G_floor) {
            low = mid;
        } else {
            high = mid;
        }
    }
    return low;
}

double solve_I_aero(Runner r, const Segment course[], int n_segments, double M_max) {
    double E     = r.endurance_index;
    double low = (LCDA_STANDING_M / M_max) + 0.01;
    if (-E / 100.0 > low) low = -E / 100.0;
    double t_tmp = 0.0;
    double high  = 1.0;
   
    for (int iter = 0; iter < BISECT_ITERS; iter++) {
        double mid = 0.5 * (low + high);
        simulate_course(r, course, n_segments, mid, &t_tmp, NULL, NULL, NULL);
        if (mid < peronnet_I_max(t_tmp, E)) {
            low = mid;
        } else {
            high = mid;
        }
    }
    return low;
}

Result find_I_opt(Runner r, const Segment course[], int n_segments) {
    Result res;
    res.I_opt   = 0.0;
    res.I_glyc  = 0.0;
    res.I_aero  = 0.0;
    res.t_tot   = 0.0;
    res.G_final = r.G0;
    res.G_min   = r.G0;
    res.limiter = LIMIT_AEROBIC;

    double M_max = compute_M_max(r);

    res.I_glyc = solve_I_glyc(r, course, n_segments, M_max);
    res.I_aero = solve_I_aero(r, course, n_segments, M_max);

    if (res.I_glyc < res.I_aero) {
        res.I_opt   = res.I_glyc;
        res.limiter = LIMIT_GLYCOGEN;
    } else {
        res.I_opt   = res.I_aero;
        res.limiter = LIMIT_AEROBIC;
    }

    simulate_course(r, course, n_segments, res.I_opt,
                    &res.t_tot, &res.G_final, &res.G_min, NULL);
    return res;
}