//
// Created by Nan on 7/23/2017.
//

#include <Python.h>
//#include <stdio.h>
//#include <stdlib.h>
#include <math.h>
//#include <string.h>
#include "ProbFunc.h"

long int statistical_bound_of_waiting_time1(double p, long int k, double alpha)
{
    double a1 = dpow(p, k);
    double a2 = (1 - p) * a1;
    long int x = 0;
    double Sum = 0.00;

    double *last_k_prob =
            (double *) MALLOC((k + 1) * sizeof(double));
    //ASSERT(last_k_prob, waiting_time_distrib_1);


    while (Sum < (1 - alpha)) {

        if (x < k) {
            last_k_prob[x % (k + 1)] = 0.00;
        } else if (x == k) {
            last_k_prob[x % (k + 1)] = a1;
        } else if (x <= 2 * k) {
            last_k_prob[x % (k + 1)] = a2;
        } else {
            last_k_prob[x % (k + 1)] =
                    last_k_prob[(x + k) % (k + 1)] -
                    a2 * last_k_prob[x % (k + 1)];
        }
        Sum += last_k_prob[x % (k + 1)];
        x++;
    }
    FREE(last_k_prob, (k + 1) * sizeof(double));

    return x;
}

double *randomwalk_probability_of_pos3(double pI, long int L)
{

    long int i, j;
    long int P = L;
    double a = pI * 0.50;
    double b = 1 - pI;
    double *u, *f, *t, *s;

    u = (double *) MALLOC((2 * L + 1) * sizeof(double));
    ASSERT(u, randomwalk_probability_of_pos3);
    f = (double *) MALLOC((2 * L + 1) * sizeof(double));
    ASSERT(f, randomwalk_probability_of_pos3);
    t = (double *) MALLOC((2 * L + 1) * sizeof(double));
    ASSERT(t, randomwalk_probability_of_pos3);

    /* (1) tables inits */

    for (i = 0; i < 2 * L + 1; i++)
        u[i] = 0;
    for (i = 0; i < 2 * L + 1; i++)
        f[i] = 0;

    f[0] = 1.00;
    u[0] = a;
    u[1] = b;
    u[2] = a;

    /* (2) qpow */
    while (P > 0) {
        if (P & 1) {
            /* f = f*u */
            for (i = 0; i < 2 * L + 1; i++) {
                t[i] = 0;
                for (j = 0; j <= i; j++)
                    t[i] += u[i - j] * f[j];
            }
            s = t;
            t = f;
            f = s;
        }
        /* u = u * u; */
        for (i = 0; i < 2 * L + 1; i++) {
            t[i] = 0;
            for (j = 0; j <= i; j++)
                t[i] += u[i - j] * u[j];
        }
        s = t;
        t = u;
        u = s;
        P >>= 1;
    }
    FREE(u, (2 * L + 1) * sizeof(double));
    FREE(t, (2 * L + 1) * sizeof(double));
    return f;

}

long int statistical_bound_of_randomwalk2(double pI, long int L, double alpha)
{
    double *RDW_Bound = randomwalk_probability_of_pos3(pI, L);
    double Sum = RDW_Bound[L];
    long int bound = 0;
    while (Sum < (1 - alpha) && bound < L) {
        bound++;
        Sum += RDW_Bound[L - bound] + RDW_Bound[L + bound];
    }

    FREE(RDW_Bound, (2 * L + 1) * sizeof(double));
    return bound;
}


static char module_docstring[] =
        "This module provides an interface for calculating chi-squared using C.";
static char stat_randomwalk_docstring[] =
        "Calculate the chi-squared of some data given a model.";

static PyObject *bound_randomwalk(PyObject *self, PyObject *args) {
    double p;
    long int l;
    double alpha;
    /* Do your stuff here. */
    if (!PyArg_ParseTuple(args, "dld", &p, &l, &alpha)) {
        return NULL;
    }

    long int value = statistical_bound_of_randomwalk2(p, l, alpha);

    /* Build the output tuple */
    PyObject *ret = Py_BuildValue("l", value);
    return ret;
}


static PyMethodDef module_methods[] = {
        { "statbound_randomwalk", bound_randomwalk, METH_NOARGS, stat_randomwalk_docstring },
        { NULL, NULL, 0, NULL }
};

PyMODINIT_FUNC init_prob() {
    Py_InitModule3("prob", module_methods, module_docstring);
}
