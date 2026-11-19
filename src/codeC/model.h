#ifndef MODEL_H
#define MODEL_H


#define MAX_SLOPE         0.45    // Minetti (2002) polynomial validity domain [-45%, +45%]
#define MIN_WALK_SLOPE    0.05    // Minimum uphill slope to allow walk transition (5%)
#define G_GRAVITY         9.81    // Gravitational acceleration [m/s^2]
#define FROUDE_UPHILL     0.28    // Critical Froude number for uphill run-to-walk transition
#define LEG_LENGTH_RATIO  0.53    // Leg length approximated as 53% of standing height
#define GLYCOGEN_ENERGY   16736.0 // Energy released per gram of glycogen [J/g]
#define WALL_THRESHOLD    0.10    // Stopping condition: 10% of initial glycogen stock G0


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
} Runner;

typedef struct {
    double I_opt;    // Optimal relative intensity [0, 1]
    double t_tot;    // Total predicted race time [s]
    double G_final;  // Remaining glycogen at finish [g]
} Result;


double minetti_run(double slope);
double minetti_walk(double slope);

double compute_Pmax(Runner r);
double compute_v_trans(double height);
double compute_carb_efficiency(double Rin_h);

double compute_total_time(double P_target, double mass, double v_trans,
                          const Segment course[], int n_segments);

double compute_carb_balance(double I, double P_target, double t_tot,
                            double G0, double Rin_h);

double simulate_race(double I, Runner r, const Segment course[],
                     int n_segments, double *t_tot_out);

Result find_I_opt(Runner r, const Segment course[], int n_segments);

#endif