clear;

addpath("bin");
if ~isfile(fullfile("bin", "mex_solver." + mexext))
    mex -R2018a src/codeC/mex_solver.c src/codeC/model.c -Isrc/codeC -outdir bin
end

% Parcours test
distances = [15000; 10000; 15000];
slopes    = [0.10;    0.00; -0.10];

% Profil coureur
runner.mass   = 70.0;
runner.vma    = 20.0;
runner.height = 1.75;
runner.G0     = 400.0;
runner.endurance_index = -7.5;    % -5 à -6 : très endurant | -7 à -8 : entraîné | -9 à -10 : peu endurant
runner.terrain_factor  = 1.20;    % 1.00 : bitume | 1.10 : sentier roulant | 1.20 : technique
runner.fatigue_rate    = 0.010;  % mettre 0 pour désactiver la dérive de fatigue, max = 0.15


% Comparaison des stratégies glucidiques
Rin_list  = [0; 30; 60; 90; 120];
n_cases   = length(Rin_list);

temps_str = strings(n_cases, 1);
G_fin_g   = zeros(n_cases, 1);
limite    = strings(n_cases, 1);

for k = 1:n_cases
    runner.Rin_h = Rin_list(k);
    [~, t_k, G_fin_k, prof_k] = mex_solver(distances, slopes, runner);

    temps_str(k) = string(duration(0, 0, round(t_k)));
    G_fin_g(k)   = round(G_fin_k, 1);
    limite(k)    = string(prof_k.limiting_factor);
end

T = table(Rin_list, temps_str, G_fin_g, limite, ...
    'VariableNames', {'Rin_g_h', 'Temps', 'G_final_g', 'Limite'});
disp(T);