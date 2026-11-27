clear;
addpath("bin");

%% 1. Parcours fixe : Marathon trail (42.2 km, +2000 m D+ / -2000 m D-)
% 16 km de montee, 7.2 km de plat, 19 km de descente
distances = [3200; 6000; 5000;  6000; 4000; 5000;  8000;  5000];
slopes    = [0.00; 0.10; 0.18; -0.15; 0.00; 0.10; -0.10; -0.06];

%% 2. Trois profils de coureurs, caps de descente [km/h] et fourchettes reelles [s]
r_amateur = struct("mass", 75.0, "vma", 15.0, "vo2max", 0.0, ...
    "G0", 400.0, "endurance_index", -9.0, "Rin_h", 40.0, "Rmax_h", 60.0);

r_train = struct("mass", 70.0, "vma", 18.0, "vo2max", 0.0, ...
    "G0", 400.0, "endurance_index", -7.5, "Rin_h", 60.0, "Rmax_h", 90.0);

r_elite = struct("mass", 62.0, "vma", 0.0, "vo2max", 78.0, ...
    "G0", 460.0, "endurance_index", -5.5, "Rin_h", 80.0, "Rmax_h", 90.0);

% Helper pour convertir (heures, minutes) -> secondes
hms = @(h, m) h * 3600 + m * 60;

% Colonnes : {Nom, struct_runner, caps_kmh, [t_reel_min_s, t_reel_max_s]}
groups = {
    "Amateur (VMA 15)",  r_amateur, [0, 14, 12, 10],     [hms(5, 15), hms(5, 45)];
    "Entraine (VMA 18)", r_train,   [0, 18, 15, 13, 12], [hms(4, 00), hms(4, 30)];
    "Elite (VO2 78)",    r_elite,   [0, 20, 17, 15, 14], [hms(3, 10), hms(3, 25)];
    };

%% 3. Execution, decomposition et comparaison aux chronos reels
n_rows     = sum(cellfun(@numel, groups(:, 3)));
Profil     = strings(n_rows, 1);
Cap_kmh    = zeros(n_rows, 1);
Temps_mod  = strings(n_rows, 1);
Reel_plage = strings(n_rows, 1);
Ecart_pct  = zeros(n_rows, 1);
Verdict    = strings(n_rows, 1);
Montee_m   = zeros(n_rows, 1);
Plat_m     = zeros(n_rows, 1);
Desc_m     = zeros(n_rows, 1);
G_fin_g    = zeros(n_rows, 1);

row = 0;
for g = 1:size(groups, 1)
    nom_g   = groups{g, 1};
    r_base  = groups{g, 2};
    caps    = groups{g, 3};
    t_range = groups{g, 4};

    t_min_s = t_range(1);
    t_max_s = t_range(2);
    t_mid_s = 0.5 * (t_min_s + t_max_s);

    plage_str = sprintf("%s - %s", ...
        string(duration(0, 0, t_min_s, 'Format', 'hh:mm')), ...
        string(duration(0, 0, t_max_s, 'Format', 'hh:mm')));

    for c = 1:numel(caps)
        row = row + 1;
        r_curr = r_base;
        r_curr.v_downhill_max = caps(c) / 3.6;

        [~, t_k, G_k, prof] = mex_solver(distances, slopes, r_curr);
        dt_s = distances ./ prof.speed_ms;

        Profil(row)     = nom_g;
        Cap_kmh(row)    = caps(c);
        Temps_mod(row)  = string(duration(0, 0, round(t_k)));
        Reel_plage(row) = plage_str;
        Ecart_pct(row)  = round(100 * (t_k - t_mid_s) / t_mid_s, 1);

        if t_k < t_min_s
            Verdict(row) = "Trop rapide";
        elseif t_k > t_max_s
            Verdict(row) = "Trop lent";
        else
            Verdict(row) = "OK";
        end

        Montee_m(row) = round(sum(dt_s(slopes > 0))  / 60.0, 1);
        Plat_m(row)   = round(sum(dt_s(slopes == 0)) / 60.0, 1);
        Desc_m(row)   = round(sum(dt_s(slopes < 0))  / 60.0, 1);
        G_fin_g(row)  = round(G_k, 1);
    end
end

T = table(Profil, Cap_kmh, Temps_mod, Reel_plage, Ecart_pct, Verdict, ...
    Montee_m, Plat_m, Desc_m, G_fin_g);
disp(T);