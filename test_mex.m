clear;

addpath("bin");
force_rebuild = false; % Passer a false une fois le code C stabilise
if force_rebuild || ~isfile(fullfile("bin", "mex_solver." + mexext))
    mex -R2018a src/codeC/mex_solver.c src/codeC/model.c -Isrc/codeC -outdir bin
end

% Parcours test (40 km, 1500 m D+)
distances = [15000; 10000; 15000];
slopes    = [0.10;   0.00; -0.10];

% Profil coureur
runner.mass            = 70.0;
runner.vma             = 17.0;
runner.G0              = 300.0;   % ou Rapoport : 100 + runner.mass * 0.214 * 20
runner.endurance_index = -7.5;    % -5/-6 : tres endurant | -7/-8 : entraine | -9/-10 : peu endurant
runner.Rmax_h          = 90.0;    % 60 : glucose seul | 90 : mix glucose + fructose [g/h]
runner.v_downhill_max  = 14.0/3.6; % 0 = desactive | ex. 18/3.6 pour brider a 18 km/h [m/s]
                                  % 10 km/h pour amateur, 12 km/h pour
                                  % entrainé, 14 km/h pour élite (déduit de
                                  % test_downhill

% Comparaison des strategies glucidiques
Rin_list  = [0; 30; 60; 90];
n_cases   = length(Rin_list);

temps_str = strings(n_cases, 1);
I_opt_pct = zeros(n_cases, 1);
G_fin_g   = zeros(n_cases, 1);
limite    = strings(n_cases, 1);

for k = 1:n_cases
    runner.Rin_h = Rin_list(k);
    [I_k, t_k, G_fin_k, prof_k] = mex_solver(distances, slopes, runner);

    temps_str(k) = string(duration(0, 0, round(t_k)));
    I_opt_pct(k) = round(100 * I_k, 1);
    G_fin_g(k)   = round(G_fin_k, 1);
    limite(k)    = string(prof_k.limiting_factor);
end

T = table(Rin_list, temps_str, I_opt_pct, G_fin_g, limite, ...
    'VariableNames', {'Rin_g_h', 'Temps', 'I_opt_pct', 'G_final_g', 'Limite'});
disp(T);