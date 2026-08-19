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
#define N_MAX_SIM ((int)1e7)

static int cmp_desc_int(const void *a, const void *b)
{
    const int ia = *(const int *)a;
    const int ib = *(const int *)b;
    return (ib - ia);
}

static int is_close_to_integer(double x)
{
    double nearest = round(x);
    return fabs(x - nearest) < 1e-9;
}

static void fenwick_add(long long *tree, int n, int idx, long long delta)
{
    while (idx <= n)
    {
        tree[idx] += delta;
        idx += idx & -idx;
    }
}

static long long fenwick_sum(const long long *tree, int idx)
{
    long long out = 0;
    while (idx > 0)
    {
        out += tree[idx];
        idx -= idx & -idx;
    }
    return out;
}

static int fenwick_find_prefix(const long long *tree, int n, long long target)
{
    int idx = 0;
    int bit = 1;
    while ((bit << 1) <= n)
    {
        bit <<= 1;
    }
    while (bit > 0)
    {
        int next = idx + bit;
        if (next <= n && tree[next] < target)
        {
            idx = next;
            target -= tree[next];
        }
        bit >>= 1;
    }
    return idx + 1;
}

int main(int argc, char *argv[])
{
    char fn1[256], fn2[256];
    FILE *traj_file, *freq_file;

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
    if (!is_close_to_integer(nu) || !is_close_to_integer(rho))
    {
        fprintf(stderr, "Error: this implementation currently requires integer nu and rho.\n");
        return EXIT_FAILURE;
    }
    const int nu_i = (int)llround(nu);
    const int rho_i = (int)llround(rho);

    int seed = 1;
    srand48(seed);

    /* ------------------------------
       Stage 1: exact initial conditions
       ------------------------------ */
    int *seq_init = calloc(N_MAX_INIT, sizeof(*seq_init));
    int *freq_init = calloc((size_t)N0, sizeof(*freq_init));
    if (!seq_init || !freq_init)
    {
        perror("calloc initial conditions");
        free(seq_init);
        free(freq_init);
        return EXIT_FAILURE;
    }

    int D = 0;
    int t_init = 0;
    int N0_init = 1;
    long long total_weight_old = 0;
    long long total_weight_new = N0_init;
    long long total_weight = N0_init;

    while (D < N0 && t_init < N_MAX_INIT)
    {
        double p_old = (double)total_weight_old / (double)total_weight;
        double ran = (lrand48() / (RAND_MAX + 1.0));
        int estr = -1;

        if (ran > p_old)
        {
            estr = D;
            D += 1;
            total_weight_new += nu_i;
            total_weight_old += rho_i;
            total_weight += rho_i + nu_i;
        }
        else
        {
            if (t_init == 0)
            {
                continue;
            }
            int random = rand() % t_init;
            estr = seq_init[random];
            total_weight_old += rho_i;
            total_weight += rho_i;
        }

        seq_init[t_init] = estr;
        if (estr >= 0 && estr < N0)
        {
            freq_init[estr] += rho_i;
        }
        t_init += 1;
    }

    if (D < N0)
    {
        fprintf(stderr, "Reached N_MAX_INIT before D reached N0. Increase N_MAX_INIT.\n");
        free(seq_init);
        free(freq_init);
        return EXIT_FAILURE;
    }

    qsort(freq_init, (size_t)N0, sizeof(*freq_init), cmp_desc_int);

    /* ------------------------------
       Stage 2: exact dynamics with efficient structures
       ------------------------------ */
    int *times = calloc(N_MAX_SIM, sizeof(*times));
    int *freq_seen = calloc(N_MAX_SIM, sizeof(*freq_seen));
    int *old_weight = calloc(N_MAX_SIM, sizeof(*old_weight));
    int *reservoir_weight = calloc((size_t)N0, sizeof(*reservoir_weight));
    long long *fenwick = calloc((size_t)N0 + 1, sizeof(*fenwick));
    long long *old_fenwick = calloc((size_t)N_MAX_SIM + 2, sizeof(*old_fenwick));

    if (!times || !freq_seen || !old_weight || !reservoir_weight || !fenwick || !old_fenwick)
    {
        perror("calloc stage 2");
        free(seq_init);
        free(freq_init);
        free(times);
        free(freq_seen);
        free(old_weight);
        free(reservoir_weight);
        free(fenwick);
        free(old_fenwick);
        return EXIT_FAILURE;
    }

    long long total_reinforcements = 0;
    for (int i = 0; i < N0; i++)
    {
        reservoir_weight[i] = rho_i * freq_init[i];
        total_reinforcements += reservoir_weight[i];
        fenwick_add(fenwick, N0, i + 1, reservoir_weight[i]);
    }
    free(seq_init);
    free(freq_init);

    /* Match the full code: reset the stage-2 RNG stream after building
       the initial conditions. */
    srand48(seed);

    long long adj_possible_reinforced = total_reinforcements;
    long long adj_possible_normal = (long long)nu_i * (long long)N0;
    int N_objects = 1;

    long long t = 0;
    D = 0;
    int max_object_id = 0;
    double w = 0.0;
    long long frontier_mass = total_reinforcements;

    make_dir("data_simulations");

    snprintf(fn1, sizeof(fn1),
             "data_simulations/model_T=%d_nu=%.6g_rho=%.6g_N0=%d_a=%.6g_b=%.6g_w0=%.6g.dat",
             tau_max, nu, rho, N0, a, b, w0);
    snprintf(fn2, sizeof(fn2),
             "data_simulations/n_model_T=%d_nu=%.6g_rho=%.6g_N0=%d_a=%.6g_b=%.6g_w0=%.6g.dat",
             tau_max, nu, rho, N0, a, b, w0);

    traj_file = fopen(fn1, "w");
    if (!traj_file)
    {
        perror("fopen traj_file");
        free(times);
        free(freq_seen);
        free(old_weight);
        free(reservoir_weight);
        free(fenwick);
        free(old_fenwick);
        return EXIT_FAILURE;
    }

    freq_file = fopen(fn2, "w");
    if (!freq_file)
    {
        perror("fopen freq_file");
        fclose(traj_file);
        free(times);
        free(freq_seen);
        free(old_weight);
        free(reservoir_weight);
        free(fenwick);
        free(old_fenwick);
        return EXIT_FAILURE;
    }

    int step_print = 10000;
    int step_fprint = 1;

    for (int tau = 0; tau < tau_max; tau++)
    {
        double num_estr = a * (w0 + w);
        if (num_estr < 1.0)
        {
            num_estr = 1.0;
        }

        for (int j = 0; j < num_estr; j++)
        {
            t += 1;
            if (t >= N_MAX_SIM)
            {
                fprintf(stderr, "Reached N_MAX_SIM. Increase the buffer size.\n");
                fclose(traj_file);
                fclose(freq_file);
                free(times);
                free(freq_seen);
                free(old_weight);
                free(reservoir_weight);
                free(fenwick);
                free(old_fenwick);
                return EXIT_FAILURE;
            }

            int estr = -1;
            long long total_old_mass = fenwick_sum(old_fenwick, D);
            long long total_mass = total_old_mass + frontier_mass;
            long long draw = (long long)(lrand48() / (RAND_MAX + 1.0) * total_mass);
            if (draw < 0)
            {
                draw = 0;
            }
            if (draw >= total_mass)
            {
                draw = total_mass - 1;
            }

            if (draw < total_old_mass)
            {
                estr = fenwick_find_prefix(old_fenwick, D, draw + 1) - 1;
            }
            else
            {
                estr = N_objects - 1;
            }

            if (t % step_print == 0)
            {
                printf("tau = %d\t t = %lld\tw=%d\tD = %d\testr=%d\n", tau, t, (int)w, D, estr);
            }
            if (t % step_fprint == 0)
            {
                fprintf(traj_file, "%d\t%lld\t%d\t%d\t%d\t%d\n", tau, t, D, estr, (int)num_estr, (int)w);
            }

            if (estr == N_objects - 1)
            {
                int obj_id = estr;
                times[obj_id] = (int)t;
                freq_seen[obj_id] = 1;

                D += 1;
                w += b;
                frontier_mass += nu_i;
                N_objects += 1;
                if (N_objects > max_object_id)
                {
                    max_object_id = N_objects;
                }

                int added_old_mass = rho_i;
                if (adj_possible_reinforced > 0)
                {
                    double ran = (lrand48() / (RAND_MAX + 1.0));
                    double p_reinf = (double)adj_possible_reinforced /
                                     (double)(adj_possible_reinforced + adj_possible_normal);

                    if (ran <= p_reinf)
                    {
                        long long target = (long long)(lrand48() / (RAND_MAX + 1.0) * total_reinforcements) + 1;
                        if (target < 1)
                        {
                            target = 1;
                        }
                        if (target > total_reinforcements)
                        {
                            target = total_reinforcements;
                        }
                        int idx = fenwick_find_prefix(fenwick, N0, target) - 1;
                        int reinforcement = reservoir_weight[idx];

                        added_old_mass = reinforcement + rho_i;
                        total_reinforcements -= reinforcement;
                        adj_possible_reinforced -= reinforcement;
                        adj_possible_normal += reinforcement + nu_i;

                        fenwick_add(fenwick, N0, idx + 1, -reinforcement);
                        reservoir_weight[idx] = 0;
                    }
                    else
                    {
                        adj_possible_normal += nu_i + rho_i;
                    }
                }

                old_weight[obj_id] += added_old_mass;
                fenwick_add(old_fenwick, N_MAX_SIM + 1, obj_id + 1, added_old_mass);
            }
            else
            {
                freq_seen[estr] += 1;
                old_weight[estr] += rho_i;
                fenwick_add(old_fenwick, N_MAX_SIM + 1, estr + 1, rho_i);
            }
        }
    }

    long long total_mass = fenwick_sum(old_fenwick, D) + frontier_mass;
    if (total_mass <= 0)
    {
        total_mass = 1;
    }

    for (int i = 0; i < D; i++)
    {
        double prob = (double)old_weight[i] / (double)total_mass;
        fprintf(freq_file, "%d\t%d\t%d\t%.12g\n", times[i], freq_seen[i], old_weight[i], prob);
    }

    fclose(traj_file);
    fclose(freq_file);

    free(times);
    free(freq_seen);
    free(old_weight);
    free(reservoir_weight);
    free(fenwick);
    free(old_fenwick);

    return EXIT_SUCCESS;
}
