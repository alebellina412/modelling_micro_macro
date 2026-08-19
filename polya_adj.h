/* Small helpers shared by the three model programs. */

#ifndef POLYA_ADJ_H
#define POLYA_ADJ_H

#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>

/* Create the output directory, ignoring the error if it already exists. */
void make_dir(char *dirname)
{
    mkdir(dirname, 0777);
    errno = 0;
}

/* Draw an index from a cumulative distribution.

   prob[] holds the per-element probabilities and r is a uniform draw in
   [0, 1). Walks the array accumulating mass until it passes r, and returns
   the index where that happens. The caller guarantees that prob[] sums to at
   least r, so the walk always terminates inside the array. */
int sample(double *prob, double r)
{
    int tmp = 0;
    double cum = prob[tmp];
    while (r > cum)
        cum += prob[++tmp];
    return tmp;
}

#endif /* POLYA_ADJ_H */
