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

#define NARGC 8
#define N_MAX_INIT ((int)1e7)
#define N_MAX_SIM ((int)1e6)

static int cmp_desc_double(const void *a, const void *b)
{
    double da = *(const double *)a;
    double db = *(const double *)b;
    if (da < db)
        return 1;
    if (da > db)
        return -1;
    return 0;
}

int main(int argc, char *argv[])
{
    char fn1[256], fn2[256];
    FILE *dati1, *dati2;

    int tau_max, N0;
    double nu, rho, a, b, w0;

    if (argc != NARGC)
    {
        printf("Usage: %s <T> <nu> <rho> <N0> <a> <b> <w0>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    tau_max = atoi(argv[1]);
    nu = atof(argv[2]);
    rho = atof(argv[3]);
    N0 = atoi(argv[4]);
    a = atof(argv[5]);
    b = atof(argv[6]);
    w0 = atof(argv[7]);

    if (tau_max <= 0 || N0 <= 0)
    {
        fprintf(stderr, "Error: T and N0 must be positive.\n");
        return EXIT_FAILURE;
    }

    int seed = 1;
    srand48(seed);

    // ------------------------------
    // Stage 1: create initial conditions (UMT without explorers)
    // ------------------------------
    printf("Creating initial conditions (no explorers)...\n");

    int *seq = calloc(N_MAX_INIT, sizeof(*seq));
    double *freq_init = calloc(N0, sizeof(*freq_init));

    if (!seq || !freq_init)
    {
        perror("calloc initial conditions");
        free(seq);
        free(freq_init);
        return EXIT_FAILURE;
    }

    int D = 0;
    int t = 0;
    int step_print_init = 1000;

    // Match create_initial_conditions_UMT: initial urn weight N0_init, D grows until Dmax (= N0)
    int N0_init = 1;
    double total_weight_old = 0.0;
    double total_weight_new = (double)N0_init;
    double total_weight = (double)N0_init;

    while (D < N0 && t < N_MAX_INIT)
    {
        double p_old = total_weight_old / total_weight;
        double ran = (lrand48() / (RAND_MAX + 1.0));
        int estr;

        if (ran > p_old)
        {
            estr = D;
            D++;
            total_weight_new += nu;
            total_weight_old += rho;
            total_weight += rho;
            total_weight += nu;
        }
        else
        {
            if (t == 0)
            {
                continue;
            }
            int random = rand() % t;
            estr = seq[random];
            total_weight_old += rho;
            total_weight += rho;
        }

        seq[t] = estr;

        if (estr >= 0 && estr < N0)
        {
            freq_init[estr] += 1.0;
        }

        if (t % step_print_init == 0)
        {
            printf("init t = %d\tD = %d\testr = %d\n", t, D, estr);
        }

        t++;
    }

    if (D < N0)
    {
        fprintf(stderr, "Reached N_MAX_INIT before D reached N0. Increase N_MAX_INIT.\n");
        free(seq);
        free(freq_init);
        return EXIT_FAILURE;
    }

    printf("Initial conditions finished. Steps = %d, D = %d\n", t, D);

    free(seq);

    // ------------------------------
    // Stage 2: UMT-TAP dynamics with explorers
    // ------------------------------
    int *freq = calloc(N_MAX_SIM, sizeof(*freq));
    int *times = calloc(N_MAX_SIM, sizeof(*times));
    double *prob = calloc(N_MAX_SIM, sizeof(*prob));
    double *first_reinforcements = calloc(N0, sizeof(*first_reinforcements));
    double *prob_reinforcements = calloc(N0, sizeof(*prob_reinforcements));

    if (!freq || !times || !prob || !first_reinforcements || !prob_reinforcements)
    {
        perror("calloc arrays");
        free(freq);
        free(times);
        free(prob);
        free(first_reinforcements);
        free(prob_reinforcements);
        free(freq_init);
        return EXIT_FAILURE;
    }

    // Sort initial frequencies by rank (descending), as in frequencies_*.dat
    double *freq_sorted = malloc((size_t)N0 * sizeof(*freq_sorted));
    if (!freq_sorted)
    {
        perror("malloc freq_sorted");
        free(freq_init);
        return EXIT_FAILURE;
    }
    for (int i = 0; i < N0; i++)
    {
        freq_sorted[i] = freq_init[i];
    }

    qsort(freq_sorted, (size_t)N0, sizeof(*freq_sorted), cmp_desc_double);

    double N_total = 0.0;
    for (int i = 0; i < N0; i++)
    {
        first_reinforcements[i] = rho * freq_sorted[i];
        N_total += first_reinforcements[i];
    }

    free(freq_sorted);
    free(freq_init);

    double total_reinforcements = N_total;
    int N_reinforcements = N0;
    double adj_possible_reinforced = N_total;
    double adj_possible_normal = nu * N0;
    int N_objects = 1;

    t = 0;
    double w = 0.0;
    double delta_w = b;

    prob[0] = 1.0;

    // reset RNG to match a separate run of the simulation stage
    srand48(seed);

    make_dir("data_simulations");

    snprintf(fn1, sizeof(fn1),
             "data_simulations/model_T=%d_nu=%.6g_rho=%.6g_N0=%d_a=%.6g_b=%.6g_w0=%.6g.dat",
             tau_max, nu, rho, N0, a, b, w0);
    snprintf(fn2, sizeof(fn2),
             "data_simulations/n_model_T=%d_nu=%.6g_rho=%.6g_N0=%d_a=%.6g_b=%.6g_w0=%.6g.dat",
             tau_max, nu, rho, N0, a, b, w0);

    dati1 = fopen(fn1, "w");
    if (!dati1)
    {
        perror("fopen dati1");
        return EXIT_FAILURE;
    }

    dati2 = fopen(fn2, "w");
    if (!dati2)
    {
        perror("fopen dati2");
        fclose(dati1);
        return EXIT_FAILURE;
    }

    D = 0;

    int step_print = 1000;
    int step_fprint = 1;

    for (int tau = 0; tau < tau_max; tau++)
    {
        double num_estr = a * (w0 + w);

        for (int j = 0; j < num_estr; j++)
        {
            t += 1;

            double ran = (lrand48() / (RAND_MAX + 1.0));
            int estr = sample(prob, ran);

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
                double N_total_temp = N_total;

                freq[estr] = 1;

                if (N_reinforcements > 0)
                {
                    ran = (lrand48() / (RAND_MAX + 1.0));

                    if (ran <= (adj_possible_reinforced / (adj_possible_reinforced + adj_possible_normal)))
                    {
                        for (int i = 0; i < N_reinforcements; i++)
                        {
                            prob_reinforcements[i] = first_reinforcements[i] / total_reinforcements;
                        }

                        ran = (lrand48() / (RAND_MAX + 1.0));
                        int index = sample(prob_reinforcements, ran);

                        N_total += (first_reinforcements[index] + rho);
                        N_total += nu;

                        prob[estr] = (first_reinforcements[index] + rho) / N_total_temp;
                        prob[estr + 1] += (nu) / N_total_temp;

                        total_reinforcements -= first_reinforcements[index];
                        adj_possible_reinforced -= first_reinforcements[index];
                        adj_possible_normal += (first_reinforcements[index] + nu);

                        for (int i = index; i < N_reinforcements; i++)
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

                        adj_possible_normal += (nu + 1);
                    }
                }
                else
                {
                    N_total += rho;
                    N_total += nu;

                    prob[estr] = rho / N_total_temp;
                    prob[estr + 1] += (nu) / N_total_temp;
                }

                for (int i = 0; i < N_objects; i++)
                {
                    prob[i] *= (N_total_temp / N_total);
                }
            }
            else
            {
                double N_total_temp = N_total;
                N_total += rho;

                freq[estr] += 1;
                prob[estr] += rho / N_total_temp;

                for (int i = 0; i < N_objects; i++)
                {
                    prob[i] *= (N_total_temp / N_total);
                }
            }
        }
    }

    for (int i = 0; i < N_objects; i++)
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

    return 0;
}
