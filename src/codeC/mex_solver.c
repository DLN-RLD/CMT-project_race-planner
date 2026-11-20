#include "mex.h"
#include "model.h"

/* Reads a scalar double field from a MATLAB struct */
static double get_struct_scalar(const mxArray *s, const char *field_name,
                                double default_val, int required);

/* Validates input and output arguments (same structure as arrayProduct.c) */
static void checkArguments(int nlhs, mxArray *plhs[],
                           int nrhs, const mxArray *prhs[]);

/* Fills segment-by-segment output arrays at the optimal intensity */
static void fill_profile(Runner r, Result res, const Segment course[], mwSize n,
                         double *out_dist, double *out_time, double *out_spd,
                         double *out_walk, double *out_glyc);


/*
 * Gateway function
 * Syntax: [I_opt, t_tot, G_final, profile] = mex_solver(distances, slopes, runner)
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    /* --- 1. Checks --- */
    checkArguments(nlhs, plhs, nrhs, prhs);

    /* --- 2. Extract MATLAB inputs into C types --- */
    mwSize n        = mxGetNumberOfElements(prhs[0]);
    double *dist_in = mxGetDoubles(prhs[0]);
    double *slp_in  = mxGetDoubles(prhs[1]);

    Runner r;
    r.mass   = get_struct_scalar(prhs[2], "mass",   0.0, 1);
    r.height = get_struct_scalar(prhs[2], "height", 0.0, 1);
    r.G0     = get_struct_scalar(prhs[2], "G0",     0.0, 1);
    r.Rin_h  = get_struct_scalar(prhs[2], "Rin_h",  0.0, 1);
    r.vma    = get_struct_scalar(prhs[2], "vma",    0.0, 0);
    r.vo2max = get_struct_scalar(prhs[2], "vo2max", 0.0, 0);

    if (r.mass <= 0.0 || r.height <= 0.0 || r.G0 <= 0.0 || r.Rin_h < 0.0) {
        mexErrMsgIdAndTxt("TrailSolver:invalidBiometrics",
                          "mass, height, and G0 must be > 0, and Rin_h must be >= 0.");
    }
    if (compute_Pmax(r) <= 0.0) {
        mexErrMsgIdAndTxt("TrailSolver:invalidAerobicCapacity",
                          "Provide a valid vma (> 0) or vo2max (> 0) in runner struct.");
    }

    Segment *course = mxMalloc(n * sizeof(Segment));
    for (mwSize k = 0; k < n; k++) {
        if (dist_in[k] <= 0.0) {
            mexErrMsgIdAndTxt("TrailSolver:invalidSegment",
                              "All segment distances must be strictly positive.");
        }
        course[k].distance = dist_in[k];
        course[k].slope    = slp_in[k];
    }

    /* --- 3. Call computational routine (model.c) --- */
    Result res = find_I_opt(r, course, (int)n);

    /* --- 4. Populate MATLAB outputs --- */
    plhs[0] = mxCreateDoubleScalar(res.I_opt);
    if (nlhs >= 2) {
        plhs[1] = mxCreateDoubleScalar(res.t_tot);
    }
    if (nlhs >= 3) {
        plhs[2] = mxCreateDoubleScalar(res.G_final);
    }

    if (nlhs >= 4) {
        const char *field_names[] = {
            "cum_dist_m", "cum_time_s", "speed_ms", "is_walking", "glycogen_g"
        };
        plhs[3] = mxCreateStructMatrix(1, 1, 5, field_names);

        mxArray *mx_dist = mxCreateDoubleMatrix(n, 1, mxREAL);
        mxArray *mx_time = mxCreateDoubleMatrix(n, 1, mxREAL);
        mxArray *mx_spd  = mxCreateDoubleMatrix(n, 1, mxREAL);
        mxArray *mx_walk = mxCreateDoubleMatrix(n, 1, mxREAL);
        mxArray *mx_glyc = mxCreateDoubleMatrix(n, 1, mxREAL);

        fill_profile(r, res, course, n,
                     mxGetDoubles(mx_dist), mxGetDoubles(mx_time),
                     mxGetDoubles(mx_spd),  mxGetDoubles(mx_walk),
                     mxGetDoubles(mx_glyc));

        mxSetField(plhs[3], 0, "cum_dist_m", mx_dist);
        mxSetField(plhs[3], 0, "cum_time_s", mx_time);
        mxSetField(plhs[3], 0, "speed_ms",   mx_spd);
        mxSetField(plhs[3], 0, "is_walking", mx_walk);
        mxSetField(plhs[3], 0, "glycogen_g", mx_glyc);
    }

    mxFree(course);
}





/* Functions definitions */

/* Fills segment-by-segment output arrays at the optimal intensity */
static void fill_profile(Runner r, Result res, const Segment course[], mwSize n,
                         double *out_dist, double *out_time, double *out_spd,
                         double *out_walk, double *out_glyc) {
    double P_target = compute_Pmax(r) * res.I_opt;
    double v_trans  = compute_v_trans(r.height);
    double eta      = compute_carb_efficiency(r.Rin_h);
    double Rin_s    = r.Rin_h / 3600.0;
    double Rout_s   = (P_target * res.I_opt * res.I_opt) / GLYCOGEN_ENERGY;

    double cum_dist  = 0.0;
    double cum_time  = 0.0;
    double G_current = r.G0;

    for (mwSize k = 0; k < n; k++) {
        double i_k = course[k].slope;
        if (i_k > MAX_SLOPE) {
            i_k = MAX_SLOPE;
        } else if (i_k < -MAX_SLOPE) {
            i_k = -MAX_SLOPE;
        }

        double v_k = P_target / (minetti_run(i_k) * r.mass);
        double is_walking = 0.0;

        if (i_k > MIN_WALK_SLOPE && v_k <= v_trans) {
            v_k = P_target / (minetti_walk(i_k) * r.mass);
            is_walking = 1.0;
        }

        double dt = course[k].distance / v_k;
        cum_dist  += course[k].distance;
        cum_time  += dt;
        G_current += (eta * Rin_s - Rout_s) * dt;
        if (G_current > r.G0) {
            G_current = r.G0;
        }

        out_dist[k] = cum_dist;
        out_time[k] = cum_time;
        out_spd[k]  = v_k;
        out_walk[k] = is_walking;
        out_glyc[k] = G_current;
    }
}


/* Validates input and output arguments (same structure as arrayProduct.c) */
static void checkArguments(int nlhs, mxArray *plhs[],
                           int nrhs, const mxArray *prhs[]) {
    (void)plhs; // Unused, kept to match standard MEX signature

    if (nrhs != 3) {
        mexErrMsgIdAndTxt("TrailSolver:nrhs",
                          "Expected 3 inputs: mex_solver(distances, slopes, runner).");
    }
    if (nlhs > 4) {
        mexErrMsgIdAndTxt("TrailSolver:nlhs",
                          "Too many output arguments (max 4).");
    }

    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
        !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("TrailSolver:notDouble",
                          "distances and slopes must be real double arrays.");
    }

    mwSize n = mxGetNumberOfElements(prhs[0]);
    if (n == 0 || mxGetNumberOfElements(prhs[1]) != n) {
        mexErrMsgIdAndTxt("TrailSolver:sizeMismatch",
                          "distances and slopes must be non-empty and have the same length.");
    }

    if (!mxIsStruct(prhs[2]) || mxGetNumberOfElements(prhs[2]) != 1) {
        mexErrMsgIdAndTxt("TrailSolver:notStruct",
                          "Third input 'runner' must be a scalar struct.");
    }
}


/* Reads a scalar double field from a MATLAB struct */
static double get_struct_scalar(const mxArray *s, const char *field_name,
                                double default_val, int required) {
    mxArray *field = mxGetField(s, 0, field_name);
    if (field == NULL) {
        if (required) {
            mexErrMsgIdAndTxt("TrailSolver:MissingField",
                              "Missing required runner field: '%s'.", field_name);
        }
        return default_val;
    }

    if (!mxIsDouble(field) || mxIsComplex(field) || mxGetNumberOfElements(field) != 1) {
        mexErrMsgIdAndTxt("TrailSolver:InvalidField",
                          "Field '%s' must be a real scalar double.", field_name);
    }

    return mxGetScalar(field);
}