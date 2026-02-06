#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <inttypes.h>
#include <string.h>
#include "polya_adj.h"

#define NARGC 5

#define N_max pow(10, 6)

int main(int argc, char *argv[])
{
    char fn1[100], fn2[100];
    FILE *dati1, *dati2;

    char *rinforzo_rho, *rinforzo_nu;

    int tau_max, N0, i, j, *freq, *times, estr, insize, D, T, step_print, step_fprint, seed, t, N_reinforcements, index;

    double ran, *prob, delta_w, w, w_actual, total, nu, rho, N_objects, N_total, N_objects_temp, N_total_temp, num_estr, *first_reinforcements, *prob_reinforcements, total_reinforcements, adj_possible_reinforced, adj_possible_normal;

    if (argc != NARGC)
    {
        printf("Inserire <N0> <T>\n");
        exit(EXIT_FAILURE);
    }

    N0 = atoi(argv[1]);
    tau_max = atoi(argv[2]);
    rho = atof(argv[3]);
    nu = atof(argv[4]);
    // prendi in entrata anche: a, b, w0

    seed = 1;

    freq = calloc(N_max, sizeof(*freq));
    times = calloc(N_max, sizeof(*times));
    prob = calloc(N_max, sizeof(*prob));
    
    first_reinforcements = calloc(N0, sizeof(*first_reinforcements)+1);
    prob_reinforcements = calloc(N0, sizeof(*prob_reinforcements)+1);

    N_total = 0;
    
    // READ FROM FILE INITIAL CONDITIONS
    
    FILE *fp = fopen("dati/frequencies_rho=1.0_nu=1.0_Dmax=15000.dat", "r");
    if (!fp) {
        perror("Errore apertura file");
        exit(EXIT_FAILURE);
    }    
    
    // substitute this part. instead of loading the file for initial conditions, it creates it now with this code. put some print at scree: creating initial conditions,...finisched... ecc.
    
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <inttypes.h>
#include <string.h>
#include "polya_adj.h"

#define NARGC 6

#define N_max pow(10, 7)

int main(int argc, char *argv[])
{
    char fn1[100];
    FILE *dati1;

    int T_max, N0, random,  i, j, *seq, insize, D, T, step_print, step_fprint, seed, t, estr, Dmax;

    double ran, *prob, total, nu, rho, N_objects, N_total, N_objects_temp, N_total_temp, num_estr, p_new, p_old, total_weight, total_weight_old, total_weight_new;

    if (argc != NARGC)
    {
        printf("Inserire <N0> <T>\n");
        exit(EXIT_FAILURE);
    }

    N0 = atoi(argv[1]);
    T_max = atoi(argv[2]);
    rho = atof(argv[3]);
    nu = atof(argv[4]);
    Dmax = atoi(argv[5]);

    seed = 1;

    seq = calloc(N_max, sizeof(*seq));
    
    total_weight_old = 0;
    total_weight_new = N0;
    total_weight = N0;

    p_old = total_weight_old / total_weight;
    p_new = total_weight_new / total_weight;


    sprintf(fn1, "dati/initial_conditions_rho=%.1lf_nu=%.1lf_Dmax=%d.dat", rho, nu, Dmax);
    dati1 = fopen(fn1, "w");

    srand48(seed);

    D = 0;
    

    step_print = 1000;
    step_fprint = 1;


    for (int t = 0; t < T_max; t++)
    {
        
            p_old = total_weight_old / total_weight;
            p_new = total_weight_new / total_weight;
            
            ran = (lrand48() / (RAND_MAX + 1.0));

            if (ran > p_old)
            {
                estr = D;
                D++;
                seq[t] = estr;
		total_weight_new += nu;
		total_weight_old += rho;
		total_weight += rho;
		total_weight += nu;

            }
            else
            {
		random = rand() % t;
		estr = seq[random];
            	seq[t] = estr;
		total_weight_old += rho;
		total_weight += rho;

            }
            if (t % step_print == 0)
            {
                printf("t = %d\tD = %d\testr = %d\n", t, D, estr);
            }
            if (t % step_fprint == 0)
            {
            
            	fprintf(dati1, "%d\t%d\t%d\n", t, D, estr);
            }
            
	    // --- qui il controllo su D ---
	    if (D == Dmax)
	    {
		break;
	    }

    }


    fclose(dati1);
    free(seq);
}
    
    
    
    // this finished the code of create initial conditions, merge with it, without changing it, minimally. you should only rewrite the block producing the initial conditions, all the aprameters you have (D0= N0 that is in this script), but do not creates the file of initial conditions, instead, it takes directly the distirbutions and use it. 
    
    // from this on, do not change anything (only dynamic arameters a,b,w0)
    
for (int i = 0; i < N0; i++) {
    double temp;
    if (fscanf(fp, "%lf", &temp) != 1) {
	fprintf(stderr, "Errore lettura valore alla riga %d\n", i+1);
	exit(EXIT_FAILURE);
    }
    first_reinforcements[i] = rho * temp;  // moltiplica per rho
    N_total += first_reinforcements[i];
}

    fclose(fp);

    // Stampa di controllo
    // for (int i = 0; i < N0; i++) {
    //     printf("first_reinforcements[%d] = %lf\n", i, first_reinforcements[i]);
    // }
    // printf("N_total = %.0f\n", N_total);



    total_reinforcements = N_total;
    N_reinforcements = N0; 
    adj_possible_reinforced = N_total;
    // no elementi già triggerati
    // adj_possible_normal = 0;
    // inserisci elementi triggerati
    adj_possible_normal = nu*N0;
    N_objects = 1;

    t = 0;
    w = 0;
    delta_w = 0.9;

    prob[0] = 1.;

    sprintf(fn1, "dati/UMT_papers_IN_sample_from_file_rho=%.1lf_nu=%.1lf_N0=%d.dat", rho, nu, N0);
    sprintf(fn2, "dati/n_UMT_papers_IN_sample_from_file_rho=%.1lf_nu=%.1lf_N0=%d.dat", rho, nu, N0);
    dati1 = fopen(fn1, "w");
    dati2 = fopen(fn2, "w");

    // srand48(time(0));
    srand48(seed);

    D = 0;

    step_print = 1000;
    step_fprint = 1;

    double prob_sat = 0;

    double total_old = 0, total_new = N0;

    for (int tau = 0; tau < tau_max; tau++)
    {
        num_estr = 0.177*(289000 + w);

        for (int j = 0; j < num_estr; j++){

            t += 1;

            ran = (lrand48() / (RAND_MAX + 1.0));
            estr = sample(prob, ran);

            if (t % step_print == 0)
            {
                printf("tau = %d\t t = %d\tw=%d\tD = %d\testr=%d\n", tau, t, (int)w, D, estr);
            }
            if (t % step_fprint == 0)
            {
            
            	fprintf(dati1, "%d\t%d\t%d\t%d\t%d\t%d\n", tau, t, D, estr, (int)num_estr, (int)w);
	
            }

            if (estr == N_objects - 1)
            {
                times[estr] = t;
                D++;
                w += delta_w;

                prob[estr + 1] = prob[estr];

                N_objects++;
                N_total_temp = N_total;

                freq[estr] = 1;


            if (N_reinforcements > 0)
            {

                ran = (lrand48() / (RAND_MAX + 1.0));

                if(ran <= (adj_possible_reinforced / (adj_possible_reinforced + adj_possible_normal)))

                {

                for(int i = 0; i < N_reinforcements; i++)
                {
                    prob_reinforcements[i] = first_reinforcements[i] / total_reinforcements;
                }

                ran = (lrand48() / (RAND_MAX + 1.0));
                index = sample(prob_reinforcements, ran);

                N_total += (first_reinforcements[index]+rho);
                N_total += nu;

                prob[estr] = (first_reinforcements[index]+rho) / N_total_temp;
                prob[estr + 1] += (nu) / N_total_temp;	

                total_reinforcements -= first_reinforcements[index];
                adj_possible_reinforced -= first_reinforcements[index];
                adj_possible_normal += (first_reinforcements[index]+nu);

                for(int i = index; i < N_reinforcements; i++)
                {
                    first_reinforcements[i] = first_reinforcements[i + 1];
                }

                N_reinforcements--;
                }
                else
                {
                    N_total += rho;
                    N_total += nu;

                    prob[estr] = rho / N_total_temp;
                    prob[estr + 1] += (nu) / N_total_temp;

                    adj_possible_normal += (nu+1);
                }

            }

		else
		{	
                N_total += rho;
                N_total += nu;

                prob[estr] = rho / N_total_temp;
                prob[estr + 1] += (nu) / N_total_temp;
                
                }

            }
            else
            {
                N_total_temp = N_total;
                N_total += rho;

                freq[estr] += 1;
                prob[estr] += rho / N_total_temp;

            }

            for (i = 0; i < N_objects; i++)
            {
                prob[i] *= (N_total_temp / N_total);
            }
        }
    }

    for (i = 0; i < N_objects; i++)
    {
        fprintf(dati2, "%d\t%d\t%lf\t%lf\n", times[i], freq[i], prob[i] * N_total, prob[i]);
    }

    fclose(dati1);
    fclose(dati2);
    free(prob);
    free(freq);
    free(times);
    free(first_reinforcements);
    free(prob_reinforcements);

}
