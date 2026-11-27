#ifndef MODEL_H
#define MODEL_H

/* --- Locomotion constants (Looney et al., 2019 & 2026) --- */
#define MAX_SLOPE               0.45    /* RE3 & LCDA validated slope domain [-45%, +45%] [-] */
#define RE3_INTERCEPT           4.43    /* RE3 running metabolic rate intercept [W/kg] */
#define LCDA_STANDING_M         1.44    /* LCDA standing metabolic rate at v = 0 [W/kg] */

/* --- Bioenergetic & glycogen constants (Rapoport, 2010; Jeukendrup, 2014) --- */
#define O2_ENERGY_EQ            20.92   /* Caloric equivalent of O2 at RER ~ 1.0 (5.0 kcal/L) [J/mL O2] */
#define GLYCOGEN_ENERGY         16736.0 /* Energy released per gram of glycogen (4.0 kcal/g) [J/g] */
#define WALL_THRESHOLD          0.0     /* Depletion of accessible glycogen pool (Rapoport, 2010) [-] */
#define DEFAULT_RMAX_GLUCOSE    60.0    /* Exogenous oxidation cap: glucose only, SGLT1 [g/h] */
#define DEFAULT_RMAX_MIX        90.0    /* Exogenous oxidation cap: glucose + fructose [g/h] */

/* --- Aerobic endurance constants (Peronnet & Thibault, 1989) --- */
#define PERONNET_T_MAP          420.0   /* Reference sustain time at VMA/MAP (7 min) [s] */
#define DEFAULT_ENDURANCE_INDEX -7.5    /* Default Peronnet & Thibault endurance index E [-] */

/* --- Data types --- */
typedef enum {
    LIMIT_GLYCOGEN = 0,
    LIMIT_AEROBIC  = 1
} LimitingFactor;

typedef struct {
    double distance;         /* Segment length delta d_k [m] */
    double slope;            /* Dimensionless grade i_k (e.g., 0.12 for +12%) [-] */
} Segment;

typedef struct {
    double speed_ms;         /* Effective segment speed v_k [m/s] */
    double M_actual;         /* Actual gross specific metabolic power [W/kg] */
    double is_walking;       /* Locomotion mode indicator (1 = walk, 0 = run) [-] */
} SegmentState;

typedef struct {
    double mass;             /* Body mass [kg] */
    double vma;              /* Maximal Aerobic Speed (set to 0 if unknown) [km/h] */
    double vo2max;           /* VO2max (set to 0 if unknown) [mL/(kg*min)] */
    double endurance_index;  /* Peronnet & Thibault endurance index E (< 0) [-] */
    double G0;               /* Initial accessible glycogen: liver + leg muscles [g] */
    double Rin_h;            /* Exogenous carbohydrate intake rate [g/h] */
    double Rmax_h;           /* Intestinal oxidation cap (60 or 90, Jeukendrup) [g/h] */
    double v_downhill_max;   /* Biomechanical downhill speed cap (0 = disabled) [m/s] */
} Runner;

typedef struct {
    double I_opt;            /* Optimal sustainable gross intensity fraction [-] */
    double I_glyc;           /* Maximum intensity allowed by glycogen constraint [-] */
    double I_aero;           /* Maximum intensity allowed by aerobic endurance [-] */
    double t_tot;            /* Total predicted race time [s] */
    double G_final;          /* Remaining accessible glycogen at finish [g] */
    double G_min;            /* Minimum accessible glycogen reached during race [g] */
    LimitingFactor limiter;  /* Dominant physiological constraint */
} Result;

typedef struct {
    double *cum_dist_m;      /* Cumulative distance per segment [m] */
    double *cum_time_s;      /* Cumulative elapsed time per segment [s] */
    double *speed_ms;        /* Effective segment speed v_k [m/s] */
    double *is_walking;      /* Locomotion mode indicator (1 = walk, 0 = run) [-] */
    double *local_intensity; /* Local gross intensity I_k = M_k / M_max [-] */
    double *glycogen_g;      /* Remaining accessible glycogen at end of segment [g] */
} ProfileArrays;

/* --- Layer 1: Elemental published equations --- */
double compute_M_max(Runner r);
double re3_run_power(double speed_ms, double slope);
double lcda_walk_power(double speed_ms, double slope);
double rapoport_cho_fraction(double I_gross);
double jeukendrup_exo_rate(double Rin_h, double Rmax_h);
double peronnet_I_max(double t_s, double endurance_index);

/* --- Layer 2: Local segment physics --- */
double re3_run_speed(double M_target, double slope);
double lcda_walk_speed(double M_target, double slope, double v_low);
SegmentState solve_segment(Runner r, double slope, double M_target);

/* --- Layer 3: Course simulation & global optimization --- */
void simulate_course(Runner r, const Segment course[], int n_segments,
                     double I_target, double *t_tot, double *G_final,
                     double *G_min, ProfileArrays *prof);

double solve_I_glyc(Runner r, const Segment course[], int n_segments, double M_max);
double solve_I_aero(Runner r, const Segment course[], int n_segments, double M_max);
Result find_I_opt(Runner r, const Segment course[], int n_segments);

#endif /* MODEL_H */