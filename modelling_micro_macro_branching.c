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
#define N_max ((int)1e6)

int main(int argc, char *argv[])
{
    char fn1[256], fn2[256];
    FILE *traj_file, *freq_file;

    int mode, tau_max, N0, i, j, *freq, *times, estr, D, step_print, step_fprint, seed, t;
    double ran, *prob, nu, rho, N_objects, N_total, N_total_temp, w, w_actual, w_prev, delta_w, tau_decimal;

    if (argc != NARGC)
    {
        printf("Usage: %s <mode> <T> <rho> <nu>\n", argv[0]);
        printf("  mode: 0=log, 1=simple, 2=exp(singularity)\n");
        exit(EXIT_FAILURE);
    }

    mode = atoi(argv[1]);
    tau_max = atoi(argv[2]);
    rho = atof(argv[3]);
    nu = atof(argv[4]);

    if (mode < 0 || mode > 2)
    {
        fprintf(stderr, "Invalid mode. Use 0 (log), 1 (simple), 2 (exp).\n");
        return EXIT_FAILURE;
    }

    N0 = 1;

    seed = 1;
    srand48(seed);

    freq = calloc(N_max, sizeof(*freq));
    times = calloc(N_max, sizeof(*times));
    prob = calloc(N_max, sizeof(*prob));

    if (!freq || !times || !prob)
    {
        perror("calloc");
        return EXIT_FAILURE;
    }

    N_total = N0;
    N_objects = 1;

    t = 0;
    w = 1;
    w_actual = 0;

    prob[0] = 1.;

    make_dir("data_simulations");

    if (mode == 0)
    {
        snprintf(fn1, sizeof(fn1), "data_simulations/model_log_mode=0_rho=%.1lf_nu=%.1lf.dat", rho, nu);
        snprintf(fn2, sizeof(fn2), "data_simulations/n_model_log_mode=0_rho=%.1lf_nu=%.1lf.dat", rho, nu);
    }
    else if (mode == 1)
    {
        snprintf(fn1, sizeof(fn1), "data_simulations/model_simple_branching_mode=1_rho=%.1lf_nu=%.1lf.dat", rho, nu);
        snprintf(fn2, sizeof(fn2), "data_simulations/n_model_simple_branching_mode=1_rho=%.1lf_nu=%.1lf.dat", rho, nu);
    }
    else
    {
        snprintf(fn1, sizeof(fn1), "data_simulations/model_singularity_mode=2_rho=%.1lf_nu=%.1lf.dat", rho, nu);
        snprintf(fn2, sizeof(fn2), "data_simulations/n_model_singularity_mode=2_rho=%.1lf_nu=%.1lf.dat", rho, nu);
    }

    traj_file = fopen(fn1, "w");
    freq_file = fopen(fn2, "w");
    if (!traj_file || !freq_file)
    {
        perror("fopen");
        return EXIT_FAILURE;
    }

    D = 0;
    step_print = 1000;
    step_fprint = 1;

    for (int tau = 0; tau < tau_max; tau++)
    {
        w_prev = w_actual;
        w_actual = w;

        if (mode == 0)
        {
            delta_w = 100 * 1. / (D + 1);
        }
        else if (mode == 1)
        {
            delta_w = 0.1;
        }
        else
        {
            delta_w = w / 1000 + 1;
        }

        for (int j = 0; j < w_actual; j++)
        {
            tau_decimal = tau + (double)j / w_actual;
            t += 1;

            ran = (lrand48() / (RAND_MAX + 1.0));
            estr = sample(prob, ran);

            if (t % step_print == 0)
            {
                printf("tau = %lf\t t = %d\tw=%lf\tD = %d\testr=%d\n", tau_decimal, t, w, D, estr);
            }
            if (t % step_fprint == 0)
            {
                fprintf(traj_file, "%lf\t%d\t%d\t%.0lf\t%d\n", tau_decimal, t, D, w, estr);
            }

            if (estr == N_objects - 1)
            {
                times[estr] = t;
                D++;

                if (mode == 0)
                {
                    float r = (float)rand() / RAND_MAX;
                    if (r < delta_w)
                    {
                        w += 1;
                    }
                }
                else
                {
                    w += delta_w;
                }

                prob[estr + 1] = prob[estr];

                N_objects++;
                N_total_temp = N_total;

                freq[estr] = 1;

                N_total += rho;
                N_total += nu;

                prob[estr] = rho / N_total_temp;
                prob[estr + 1] += (nu) / N_total_temp;
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
        fprintf(freq_file, "%d\t%d\t%lf\t%lf\n", times[i], freq[i], prob[i] * N_total, prob[i]);
    }

    fclose(traj_file);
    fclose(freq_file);
    free(prob);
    free(freq);
    free(times);

    return 0;
}
