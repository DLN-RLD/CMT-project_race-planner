#include <stdio.h>
#include <math.h>



typedef struct {
    double distance; // delta d_k [m]
    double pente;    // i_k [sans unité, ex: 0.12 pour 12 %]
} Segment;

typedef struct {
    double mass;     // m [kg]
    double vma;      // VMA [km/h] (mettre 0 si inconnue)
    double vo2max;   // VO2max [mL/(kg·min)] (mettre 0 si inconnue)
    double height;   // taille [m] (pour le nombre de Froude)
    double G0;       // Stock initial de glycogène [g]
    double Rin_h;    // Apport nutritionnel [g/h]
} Runner;

typedef struct {
    double I_opt;    // Intensité optimale trouvée [0, 1]
    double t_tot;    // Temps total de course correspondant [s]
    double G_final;  // Stock de glycogène à l'arrivée [g]
} Result;

// Calcule la puissance aérobie maximale nette P_max [W = J/s]
// Priorité à la VMA si elle est renseignée, sinon utilisation de la VO2max.
double calcul_Pmax(Runner c) {
    double vma_kmh = 0.0;

    if (c.vma > 0.0) {
        // Cas 1 (prioritaire) : VMA connue directement [km/h]
        vma_kmh = c.vma;
    } 
    else if (c.vo2max > 0.0) {
        // Cas 2 (repli) : estimation de la VMA [km/h] via Léger & Mercier (1984)
        vma_kmh = c.vo2max / 3.5;
    } 
    else {
        fprintf(stderr, "Erreur : fournir une VMA ou une VO2max valide.\n");
        return -1.0;
    }

    // Conversion en m/s et calcul de P_max cohérent avec Minetti sur le plat
    double vma_ms = vma_kmh / 3.6;
    return minetti_run(0.0) * c.mass * vma_ms;
}


/* simule_course :
   paramètre : data coureur, data segments, I
   utilise : calcul_ttot et bilan_glucidique
   return Gfinal et son temps total

   calcul_ttot :
   paramètre : I, data coureur, data segments
   utilise : minetti course / marche
   calcul : vk
   return : ttot

   bilan_glucidique :
   paramètre : data coureur, ttot
   calcul : Gfinal = G0 + (Rin - Rout) ttot
   utilise : calcul_ttot
   return Gfinal

   minetti course / marche :
   paramètre : pente i d'un segment (data coureur)
   return : coût du terrain
*/

double minetti_run(double i) {
    return 155.4 * pow(i, 5) - 30.4 * pow(i, 4) - 43.3 * pow(i, 3) 
         + 46.3 * pow(i, 2) + 19.5 * i + 3.6;
}

double minetti_walk(double i) {
    return 280.5 * pow(i, 5) - 58.7 * pow(i, 4) - 76.8 * pow(i, 3) 
         + 51.9 * pow(i, 2) + 19.6 * i + 2.5;
}


#define PENTE_LIMITE 0.45
#define PENTE_MIN_MARCHE 0.05

// Calcule le temps total de course (en secondes) pour une puissance cible donnée
double calcul_ttot(double P_cible, double masse, double v_trans,
                   const Segment parcours[], int nb_segments) {
    double t_tot = 0.0;

    for (int k = 0; k < nb_segments; k++) {
        double i_k = parcours[k].pente;

        // Sécurité : borner la pente dans le domaine de validité de Minetti [-45 %, +45 %]
        if (i_k > PENTE_LIMITE) {
            i_k = PENTE_LIMITE;
        } else if (i_k < -PENTE_LIMITE) {
            i_k = -PENTE_LIMITE;
        }

        // 1. Calcul de la vitesse théorique en courant (m/s)
        double v_k = P_cible / (minetti_run(i_k) * masse);

        // 2. Transition dynamique vers la marche si montée et vitesse sous le seuil
        if (i_k > PENTE_MIN_MARCHE && v_k <= v_trans) {
            v_k = P_cible / (minetti_walk(i_k) * masse);
        }

        // 3. Accumulation du temps passé sur le segment (dt_k = delta_d_k / v_k)
        t_tot += parcours[k].distance / v_k;
    }

    return t_tot;
}



#define G_GRAVITY 9.81
#define FROUDE 0.28   // Nombre de Froude critique en côte
#define LEG_LENGTH 0.53   // Longueur de jambe ≈ 53 % de la taille

// À appeler une seule fois dans simule_course ou à la lecture du coureur
double calcul_v_trans(double taille) {
    double L_jambe = LEG_LENGTH * taille;
    return sqrt(FROUDE * G_GRAVITY * L_jambe);
}


// Rendement d'oxydation exogène (eta) selon l'apport horaire [g/h]
// (Jeukendrup, 2014 ; Pfeiffer et al., 2011 ; Hearris et al., 2022)
static double calcul_rendement(double Rin_h) {
    if (Rin_h <= 60.0) {
        return 0.95;
    } else if (Rin_h <= 120.0) {
        return 0.95 - 0.20 * ((Rin_h - 60.0) / 60.0);
    } else {
        return 0.75;
    }
}
// Calcule le rendement d'absorption intestinale eta
// en fonction de l'apport en glucides ingéré (Rin_h, en g/h).
// Paramétrisation basée sur les données d'oxydation exogène (Podlogar et al., 2022 ; Hearris et al., 2022).
double calcul_rendement(double Rin_h) {
    if (Rin_h <= 60.0) {
        return 0.95; // Rendement maximal constant de 95 %
    } 
    else if (Rin_h <= 120.0) {
        // Régime 2 : Décroissance linéaire de 95 % (à 60 g/h) jusqu'à 75 % (à 120 g/h)
        // Droite y = b + a*x passant par (60 ; 0.95) et (120 ; 0.75)
        
        double a = 0.75 - 0.95 / (120.0 - 60.0); 
        double b = 0.95 - a * 60.0;
        return a * Rin_h + b;
    } 
    else {
       return 0.75; // Rendement minimal constant de 75 % au-delà de 120 g/h
    }
}

#define ENERGIE_GLYCOGENE 16736.0 // Énergie libérée par 1 g de glycogène [J/g]

// Calcule le stock de glycogène restant à l'arrivée
double bilan_glucidique(double I, double P_cible, double t_tot, double G0, double Rin_h) {
    double Gfinal = G0;

    // 1. Flux entrant : 
    double eta = calcul_rendement(Rin_h); // rendement 
    double Rin_s = Rin_h / 3600.0; // conversion de g/h en g/s 
    double input_total = eta * Rin_s * t_tot;

    // 2. Flux sortant : part glucidique selon Brooks f(I) = I^2
    double P_glucides = P_cible * I * I;            // [W] = [J/s]
    double Rout_s = P_glucides / ENERGIE_GLYCOGENE; // [g/s]
    double output_total = Rout_s * t_tot;

    // 3. Bilan final du réservoir [g]
    Gfinal = G0 + input_total - output_total;

    return Gfinal;
}

/* findIopt :
   paramètre : data coureur, data segments
   utilise : simule_course
   return : Iopt
*/

