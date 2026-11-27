#include "mex.h"
#include "model.h"

/* Reads a scalar double field from a MATLAB struct */
static double get_struct_scalar(const mxArray *s, const char *field_name,
                                double default_val, int required);

/* Validates input and output arguments (same structure as arrayProduct.c) */
static void checkArguments(int nlhs, mxArray *plhs[],
                           int nrhs, const mxArray *prhs[]);



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
    r.mass            = get_struct_scalar(prhs[2], "mass",            0.0, 1);
    r.G0              = get_struct_scalar(prhs[2], "G0",              0.0, 1);
    r.Rin_h           = get_struct_scalar(prhs[2], "Rin_h",           0.0, 1);
    r.vma             = get_struct_scalar(prhs[2], "vma",             0.0, 0);
    r.vo2max          = get_struct_scalar(prhs[2], "vo2max",          0.0, 0);
    r.endurance_index = get_struct_scalar(prhs[2], "endurance_index", DEFAULT_ENDURANCE_INDEX, 0);
    r.Rmax_h          = get_struct_scalar(prhs[2], "Rmax_h",          DEFAULT_RMAX_MIX,        0);
    r.v_downhill_max  = get_struct_scalar(prhs[2], "v_downhill_max",  0.0,                     0);

    if (r.mass <= 0.0 || r.G0 <= 0.0 || r.Rin_h < 0.0) {
        mexErrMsgIdAndTxt("TrailSolver:invalidBiometrics",
                          "mass and G0 must be > 0, and Rin_h must be >= 0.");
    }
    if (r.endurance_index >= 0.0) {
        mexErrMsgIdAndTxt("TrailSolver:invalidEndurance",
                          "endurance_index must be strictly negative (e.g. -7.5).");
    }
    if (r.Rmax_h < 0.0) {
        mexErrMsgIdAndTxt("TrailSolver:invalidRmax",
                          "Rmax_h must be >= 0 (e.g. 60 or 90 g/h, 0 = uncapped).");
    }
    if (r.v_downhill_max < 0.0) {
        mexErrMsgIdAndTxt("TrailSolver:invalidDownhillCap",
                          "v_downhill_max must be >= 0 (0 = disabled).");
    }
    if (compute_M_max(r) <= 0.0) {
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
            "cum_dist_m", "cum_time_s", "speed_ms", "is_walking",
            "local_intensity", "glycogen_g",
            "G_min_g", "I_glyc", "I_aero", "limiting_factor"
        };
        // create a struct with 10 fields for the runner profile output
        plhs[3] = mxCreateStructMatrix(1, 1, 10, field_names);

        mxArray *mx_dist = mxCreateDoubleMatrix(n, 1, mxREAL);
        mxArray *mx_time = mxCreateDoubleMatrix(n, 1, mxREAL);
        mxArray *mx_spd  = mxCreateDoubleMatrix(n, 1, mxREAL);
        mxArray *mx_walk = mxCreateDoubleMatrix(n, 1, mxREAL);
        mxArray *mx_glyc = mxCreateDoubleMatrix(n, 1, mxREAL);
        mxArray *mx_int  = mxCreateDoubleMatrix(n, 1, mxREAL);
       
        ProfileArrays prof;
        prof.cum_dist_m      = mxGetDoubles(mx_dist);
        prof.cum_time_s      = mxGetDoubles(mx_time);
        prof.speed_ms        = mxGetDoubles(mx_spd);
        prof.is_walking      = mxGetDoubles(mx_walk);
        prof.local_intensity = mxGetDoubles(mx_int);
        prof.glycogen_g      = mxGetDoubles(mx_glyc);

        double t_tmp = 0.0, g_fin_tmp = 0.0, g_min_tmp = 0.0;
        // call simulate_course to fill the runner profile arrays
        simulate_course(r, course, (int)n, res.I_opt,
                        &t_tmp, &g_fin_tmp, &g_min_tmp, &prof);

        mxSetField(plhs[3], 0, "cum_dist_m", mx_dist);
        mxSetField(plhs[3], 0, "cum_time_s", mx_time);
        mxSetField(plhs[3], 0, "speed_ms", mx_spd);
        mxSetField(plhs[3], 0, "is_walking", mx_walk);
        mxSetField(plhs[3], 0, "local_intensity", mx_int);
        mxSetField(plhs[3], 0, "glycogen_g", mx_glyc);
        mxSetField(plhs[3], 0, "G_min_g", mxCreateDoubleScalar(g_min_tmp));
        mxSetField(plhs[3], 0, "I_glyc", mxCreateDoubleScalar(res.I_glyc));
        mxSetField(plhs[3], 0, "I_aero", mxCreateDoubleScalar(res.I_aero));
        mxSetField(plhs[3], 0, "limiting_factor", 
            mxCreateString(res.limiter == LIMIT_GLYCOGEN ? "glycogen" : "aerobic"));
    }

    mxFree(course);
}





/* Functions definitions */

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