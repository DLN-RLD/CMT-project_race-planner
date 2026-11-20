#ifndef MODEL_H
#define MODEL_H


#define MAX_SLOPE         0.45    // Minetti (2002) polynomial validity domain [-45%, +45%]
#define MIN_WALK_SLOPE    0.05    // Minimum uphill slope to allow walk transition (5%)
#define G_GRAVITY         9.81    // Gravitational acceleration [m/s^2]
#define FROUDE_UPHILL     0.28    // Critical Froude number for uphill run-to-walk transition
#define LEG_LENGTH_RATIO  0.53    // Leg length approximated as 53% of standing height
#define GLYCOGEN_ENERGY   16736.0 // Energy released per gram of glycogen [J/g]
#define WALL_THRESHOLD    0.10    // Stopping condition: 10% of initial glycogen stock G0
#define BMR_PER_KG              1.2       /* Resting metabolic rate (W/kg, ~1 MET) */
#define CHO_FRAC_EXP            1.5       /* Crossover exponent fitted to Romijn et al. (1993) 
                                             Instead of Brooks law (more accurate)  */

#define PERONNET_T_MAP          420.0     /* 7 min reference time at VMA/MAP (s) */
#define DEFAULT_ENDURANCE_INDEX -7.5      /* Peronnet & Thibault (1989) */
#define DEFAULT_TERRAIN_FACTOR  1.10      /* Soule & Goldman (1972) trail factor */
#define DEFAULT_FATIGUE_RATE    0.0015    /* Brueckner et al. (1991): +0.15%/km */
#define FATIGUE_CAP             0.15      /* Vernillo et al. (2016): +15% max drift */

typedef enum {
    LIMIT_GLYCOGEN = 0,
    LIMIT_AEROBIC  = 1
} LimitingFactor;

typedef struct {
    double distance; // Segment length delta d_k [m]
    double slope;    // Dimensionless grade i_k (e.g., 0.12 for +12%)
} Segment;

typedef struct {
    double mass;     // Body mass [kg]
    double vma;      // Maximal Aerobic Speed (VMA) [km/h] (set to 0 if unknown)
    double vo2max;   // VO2max [mL/(kg*min)] (set to 0 if unknown)
    double height;   // Standing height [m] (used for Froude transition)
    double G0;       // Initial glycogen stock [g]
    double Rin_h;    // Exogenous carbohydrate intake rate [g/h]
    double endurance_index;  /* Peronnet-Thibault E (< 0, default -7.5) */
    double terrain_factor;   /* Surface cost multiplier k_t (> 0, default 1.10) */
    double fatigue_rate;     /* Cost drift per km (>= 0, default 0.0015) */
} Runner;

typedef struct {
    double I_opt;            /* Optimal sustainable intensity fraction */
    double I_glyc;           /* Maximum intensity allowed by glycogen */
    double I_aero;           /* Maximum intensity allowed by aerobic endurance */
    double t_tot;            /* Total race time (s) */
    double G_final;          /* Remaining glycogen at finish (g) */
    double G_min;            /* Minimum glycogen reached during race (g) */
    LimitingFactor limiter;  /* Dominant physiological constraint */
} Result;

typedef struct {
    double *cum_dist_m;
    double *cum_time_s;
    double *speed_ms;
    double *is_walking;
    double *local_intensity;
    double *glycogen_g;
} ProfileArrays;

/* Core bioenergetic functions */
double minetti_run(double slope);
double minetti_walk(double slope);
double downhill_pace_factor(double slope);
double compute_Pmax(Runner r);
double compute_v_trans(double height);
double compute_carb_efficiency(double Rin_h);
double peronnet_I_max(double t_s, double endurance_index);

/* Simulation & optimization functions */
void simulate_course(Runner r, const Segment course[], int n_segments,
                   double I_target, double *t_tot, double *G_final, 
                   double *G_min, ProfileArrays *prf);

Result find_I_opt(Runner r, const Segment course[], int n_segments);

#endif