#ifndef MODEL_H
#define MODEL_H

/* --- Physical & biomechanical constants --- */
#define MAX_SLOPE               0.45    /* Minetti (2002) polynomial validity domain [-45%, +45%] */
#define MIN_WALK_SLOPE          0.05    /* Minimum uphill slope to allow walk transition (5%) */
#define G_GRAVITY               9.81    /* Gravitational acceleration [m/s^2] */
#define FROUDE_UPHILL           0.28    /* Critical Froude number for uphill run-to-walk transition [-] */
#define LEG_LENGTH_RATIO        0.53    /* Leg length approximated as 53% of standing height [-] */

/* --- Metabolic & glycogen constants --- */
#define GLYCOGEN_ENERGY         16736.0 /* Energy released per gram of glycogen [J/g] */
#define WALL_THRESHOLD          0.10    /* Safety floor: 10% of initial glycogen stock G0 [-] */
#define BMR_PER_KG              1.2     /* Resting metabolic rate (~1 MET) [W/kg] */
#define CHO_FRAC_EXP            1.5     /* CHO crossover exponent fitted to Romijn et al. (1993) [-] */

/* --- Endurance, terrain & fatigue parameters --- */
#define PERONNET_T_MAP          420.0   /* Reference sustain time at VMA/MAP (7 min) [s] */
#define DEFAULT_ENDURANCE_INDEX -7.5    /* Default Peronnet & Thibault (1989) endurance index [-] */
#define DEFAULT_TERRAIN_FACTOR  1.10    /* Default Soule & Goldman (1972) trail cost multiplier [-] */
#define DEFAULT_FATIGUE_RATE    0.0015  /* Default Brueckner et al. (1991) cost drift (+0.15%/km) [km^-1] */
#define FATIGUE_CAP             0.15    /* Vernillo et al. (2016) maximum fatigue drift (+15%) [-] */

/* --- Data types --- */
typedef enum {
    LIMIT_GLYCOGEN = 0,
    LIMIT_AEROBIC  = 1
} LimitingFactor;

typedef struct {
    double distance;          /* Segment length delta d_k [m] */
    double slope;             /* Dimensionless grade i_k (e.g., 0.12 for +12%) [-] */
} Segment;

typedef struct {
    double mass;              /* Body mass [kg] */
    double vma;               /* Maximal Aerobic Speed (VMA) (set to 0 if unknown) [km/h] */
    double vo2max;            /* VO2max (set to 0 if unknown) [mL/(kg*min)] */
    double height;            /* Standing height (used for Froude transition) [m] */
    double G0;                /* Initial usable glycogen stock [g] */
    double Rin_h;             /* Exogenous carbohydrate intake rate [g/h] */
    double endurance_index;   /* Peronnet-Thibault index E (< 0, default -7.5) [-] */
    double terrain_factor;    /* Surface cost multiplier k_t (> 0, default 1.10) [-] */
    double fatigue_rate;      /* Progressive cost drift per km (>= 0, default 0.0015) [km^-1] */
} Runner;

typedef struct {
    double I_opt;             /* Optimal sustainable intensity fraction [-] */
    double I_glyc;            /* Maximum intensity allowed by glycogen constraint [-] */
    double I_aero;            /* Maximum intensity allowed by aerobic endurance [-] */
    double t_tot;             /* Total predicted race time [s] */
    double G_final;           /* Remaining glycogen at finish [g] */
    double G_min;             /* Minimum glycogen reached during race [g] */
    LimitingFactor limiter;   /* Dominant physiological constraint */
} Result;

typedef struct {
    double *cum_dist_m;       /* Cumulative distance per segment [m] */
    double *cum_time_s;       /* Cumulative elapsed time per segment [s] */
    double *speed_ms;         /* Segment speed [m/s] */
    double *is_walking;       /* Locomotion mode indicator (1 = walk, 0 = run) [-] */
    double *local_intensity;  /* Local relative intensity I_k = P_net / P_max [-] */
    double *glycogen_g;       /* Remaining glycogen at end of segment [g] */
} ProfileArrays;

/* --- Core bioenergetic functions --- */
double minetti_run(double slope);
double minetti_walk(double slope);
double downhill_pace_factor(double slope);
double compute_Pmax(Runner r);
double compute_v_trans(double height);
double compute_carb_efficiency(double Rin_h);
double peronnet_I_max(double t_s, double endurance_index);

/* --- Simulation & optimization functions --- */
void simulate_course(Runner r, const Segment course[], int n_segments,
                     double I_target, double *t_tot, double *G_final,
                     double *G_min, ProfileArrays *prof);

Result find_I_opt(Runner r, const Segment course[], int n_segments);

#endif /* MODEL_H */