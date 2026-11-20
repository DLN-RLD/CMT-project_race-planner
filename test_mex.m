addpath("bin");

if ~isfile(fullfile("bin", "mex_solver." + mexext))
    mex -R2018a src/codeC/mex_solver.c src/codeC/model.c -Isrc/codeC -outdir bin
end

% Test parcours factice de 40 km
distances = [15000; 10000; 15000];
slopes    = [0.10;   0.00; -0.10];

runner.mass   = 70.0;
runner.vma    = 15.0;
runner.height = 1.75;
runner.G0     = 500.0;
runner.Rin_h  = 0.0;

% call C code 
[I_opt, t_tot, G_final, profile] = mex_solver(distances, slopes, runner);



fprintf("Intensité optimale : %.2f %%\n", I_opt * 100);
fprintf("Temps prédit       : %s (%.1f s)\n", char(duration(0, 0, t_tot)), t_tot);
fprintf("Glycogène final    : %.1f g / %.1f g\n", G_final, runner.G0);
disp(profile);