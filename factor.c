/* factor -- return prime factors of n.
   Copyright (C) 2012, Eric Willisson <ericwillisson@gmail.com>

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* Based on factor.c from "coreutils", with original information:
   
   Written by Paul Rubin <phr@ocf.berkeley.edu>.
   Adapted for GNU, fixed to factor UINT_MAX by Jim Meyering.
   Arbitrary-precision code adapted by James Youngman from Torbjörn
   Granlund's factorize.c, from GNU MP version 4.2.2.
*/

#include <Python.h>
#include <gmp.h>

#define DOCSTRING "This module defines a function that calculates the prime factors of an integer."

#define MIN(a,b) (((a) < (b)) ? (a) : (b))

/* The maximum number of factors, including -1, for negative numbers.  */
#define MAX_N_FACTORS (sizeof (unsigned long long) * CHAR_BIT)

#define ARRAY_CARDINALITY(Array) (sizeof (Array) / sizeof *(Array))

/* True if the arithmetic type T is signed.  */
# define TYPE_SIGNED(t) (! ((t) 0 < (t) -1))

/* Return zero if T can be determined to be an unsigned type.
   Otherwise, return 1.
   When compiling with GCC, INT_STRLEN_BOUND uses this macro to obtain a
   tighter bound.  Otherwise, it overestimates the true bound by one byte
   when applied to unsigned types of size 2, 4, 16, ... bytes.
   The symbol signed_type_or_expr__ is private to this header file.  */
# if __GNUC__ >= 2
#  define signed_type_or_expr__(t) TYPE_SIGNED (__typeof__ (t))
# else
#  define signed_type_or_expr__(t) 1
# endif

/* Bound on length of the string representing an integer type or expression T.
   Subtract 1 for the sign bit if T is signed; log10 (2.0) < 146/485;
   add 1 for integer division truncation; add 1 more for a minus sign
   if needed.  */
#define INT_STRLEN_BOUND(t) \
  ((sizeof (t) * CHAR_BIT - signed_type_or_expr__ (t)) * 146 / 485 \
   + signed_type_or_expr__ (t) + 1)

#include "wheel-size.h"  /* For the definition of WHEEL_SIZE.  */
static const unsigned char wheel_tab[] =
{
#include "wheel.h"
};
#define WHEEL_START (wheel_tab + WHEEL_SIZE)
#define WHEEL_END (wheel_tab + ARRAY_CARDINALITY (wheel_tab))

struct factor {
	struct factor *next;
	unsigned long long n;
};

struct factor *factors = NULL;
int num_factors;

void *
xcalloc (size_t nmemb, size_t size)
{
	void *mem;

	mem = calloc (nmemb, size);
	if (mem == NULL && nmemb != 0 && size != 0) {
		fprintf (stderr, "ERROR: calloc() returned NULL!\n");
		exit (-1);
	}

	return mem;
}

void *
xmalloc (size_t size)
{
	void *mem;

	mem = malloc (size);
	if (mem == NULL && size != 0) {
		fprintf (stderr, "ERROR: malloc() returned NULL!\n");
		exit (-1);
	}

	return mem;
}

void *
xrealloc (void *ptr, size_t size)
{
	void *mem;

	mem = realloc (ptr, size);
	if (mem == NULL && size != 0) {
		fprintf (stderr, "ERROR: realloc() returned NULL!\n");
		exit (-1);
	}

	return mem;
}

void
new_factor (unsigned long long n)
{
	struct factor *factor;

	factor = xcalloc (1, sizeof *factor);

	factor->next = factors;
	factors = factor;

	factor->n = n;

	num_factors++;
}		

int *
factors_array_ints (void)
{
	int *factor_array, *fp;
	struct factor *factor;

	factor_array = xcalloc (num_factors, sizeof *factor_array);

	for (factor = factors, fp = factor_array; factor != NULL;
	     factor = factor->next, fp++) {
		*fp = factor->n;
	}

	return factor_array;
}

unsigned long long *
factors_array (void)
{
	unsigned long long *factor_array, *fp;
	struct factor *factor;

	factor_array = xcalloc (num_factors, sizeof *factor_array);

	for (factor = factors, fp = factor_array; factor != NULL;
	     factor = factor->next, fp++) {
		*fp = factor->n;
	}

	return factor_array;
}

static size_t
factor_wheel (unsigned long long n0, size_t max_n_factors)
{
	unsigned long long n = n0, d, q;
	size_t n_factors = 0;
	unsigned char const *w = wheel_tab;

	if (n <= 1)
		return n_factors;

	/* The exit condition in the following loop is correct because
	   any time it is tested one of these 3 conditions holds:
	   (1) d divides n
	   (2) n is prime
	   (3) n is composite but has no factors less than d.
	   If (1) or (2) obviously the right thing happens.
	   If (3), then since n is composite it is >= d^2. */

	d = 2;
	do
	{
		q = n / d;
		while (n == q * d)
		{
			assert (n_factors < max_n_factors);
			new_factor (d);
			n = q;
			q = n / d;
		}
		d += *(w++);
		if (w == WHEEL_END)
			w = WHEEL_START;
	}
	while (d <= q);

	if (n != 1 || n0 == 1)
	{
		assert (n_factors < max_n_factors);
		new_factor (n);
	}

	return n_factors;
}

/* Single-precision factoring */
static unsigned long long *
get_factors_single (unsigned long long n)
{
	factor_wheel (n, MAX_N_FACTORS);

	return factors_array ();
}

static void
factor_using_division (mpz_t t, unsigned int limit)
{
	mpz_t q, r;
	unsigned long int f;
	int ai;
	static unsigned int const add[] = {4, 2, 4, 2, 4, 6, 2, 6};
	unsigned int const *addv = add;
	unsigned int failures;

	mpz_init (q);
	mpz_init (r);

	f = mpz_scan1 (t, 0);
	mpz_div_2exp (t, t, f);
	while (f)
	{
		new_factor (2);
		--f;
	}

	for (;;)
	{
		mpz_tdiv_qr_ui (q, r, t, 3);
		if (mpz_cmp_ui (r, 0) != 0)
			break;
		mpz_set (t, q);
		new_factor (3);
	}

	for (;;)
	{
		mpz_tdiv_qr_ui (q, r, t, 5);
		if (mpz_cmp_ui (r, 0) != 0)
			break;
		mpz_set (t, q);
		new_factor (5);
	}

	failures = 0;
	f = 7;
	ai = 0;
	while (mpz_cmp_ui (t, 1) != 0)
	{
		mpz_tdiv_qr_ui (q, r, t, f);
		if (mpz_cmp_ui (r, 0) != 0)
		{
			f += addv[ai];
			if (mpz_cmp_ui (q, f) < 0)
				break;
			ai = (ai + 1) & 7;
			failures++;
			if (failures > limit)
				break;
		}
		else
		{
			mpz_swap (t, q);
			new_factor (f);
			failures = 0;
		}
	}

	mpz_clear (q);
	mpz_clear (r);
}

static void
factor_using_pollard_rho (mpz_t n, int a_int)
{
	mpz_t x, x1, y, P;
	mpz_t a;
	mpz_t g;
	mpz_t t1, t2;
	int k, l, c, i;

	mpz_init (g);
	mpz_init (t1);
	mpz_init (t2);

	mpz_init_set_si (a, a_int);
	mpz_init_set_si (y, 2);
	mpz_init_set_si (x, 2);
	mpz_init_set_si (x1, 2);
	k = 1;
	l = 1;
	mpz_init_set_ui (P, 1);
	c = 0;

	while (mpz_cmp_ui (n, 1) != 0)
	{
	S2:
		mpz_mul (x, x, x); mpz_add (x, x, a); mpz_mod (x, x, n);

		mpz_sub (t1, x1, x); mpz_mul (t2, P, t1); mpz_mod (P, t2, n);
		c++;
		if (c == 20)
		{
			c = 0;
			mpz_gcd (g, P, n);
			if (mpz_cmp_ui (g, 1) != 0)
				goto S4;
			mpz_set (y, x);
		}

		k--;
		if (k > 0)
			goto S2;

		mpz_gcd (g, P, n);
		if (mpz_cmp_ui (g, 1) != 0)
			goto S4;

		mpz_set (x1, x);
		k = l;
		l = 2 * l;
		for (i = 0; i < k; i++)
		{
			mpz_mul (x, x, x); mpz_add (x, x, a); mpz_mod (x, x, n);
		}
		mpz_set (y, x);
		c = 0;
		goto S2;
	S4:
		do
		{
			mpz_mul (y, y, y); mpz_add (y, y, a); mpz_mod (y, y, n);
			mpz_sub (t1, x1, y); mpz_gcd (g, t1, n);
		}
		while (mpz_cmp_ui (g, 1) == 0);

		mpz_div (n, n, g);	/* divide by g, before g is overwritten */

		if (!mpz_probab_prime_p (g, 3))
		{
			do
			{
				mp_limb_t a_limb;
				mpn_random (&a_limb, (mp_size_t) 1);
				a_int = (int) a_limb;
			}
			while (a_int == -2 || a_int == 0);

			factor_using_pollard_rho (g, a_int);
		}
		else
		{
			new_factor (mpz_get_ui (g));
		}
		mpz_mod (x, x, n);
		mpz_mod (x1, x1, n);
		mpz_mod (y, y, n);
		if (mpz_probab_prime_p (n, 3))
		{
			new_factor (mpz_get_ui (n));
			break;
		}
	}

	mpz_clear (g);
	mpz_clear (P);
	mpz_clear (t2);
	mpz_clear (t1);
	mpz_clear (a);
	mpz_clear (x1);
	mpz_clear (x);
	mpz_clear (y);
}

static int
mpcompare (const void *av, const void *bv)
{
	unsigned long long *const *a = av;
	unsigned long long *const *b = bv;
	return a > b;
}

static unsigned long long *
sort_and_return_factors (void)
{
	unsigned long long *factor_array = factors_array ();

	qsort (factor_array, num_factors, sizeof *factor_array, mpcompare);

	return factor_array;
}

/* Arbitrary-precision factoring */
static unsigned long long *
get_factors_multi (mpz_t t)
{
	if (mpz_sgn (t) != 0)
	{
		/* Set the trial division limit according to the size of t.  */
		size_t n_bits = mpz_sizeinbase (t, 2);
		unsigned int division_limit = MIN (n_bits, 1000);
		division_limit *= division_limit;

		factor_using_division (t, division_limit);

		if (mpz_cmp_ui (t, 1) != 0)
		{
			if (mpz_probab_prime_p (t, 3))
				new_factor (mpz_get_ui (t));
			else
				factor_using_pollard_rho (t, 1);
		}
	}

	mpz_clear (t);
	return sort_and_return_factors ();
}

unsigned long long *
get_factors (unsigned long long n)
{
	enum { GMP_TURNOVER_POINT = 100000 };

	if (GMP_TURNOVER_POINT <= n) {
		mpz_t t;
		mpz_init (t);
		mpz_set_ui (t, n);
		return get_factors_multi (t);
	}

	return get_factors_single (n);
}

void
free_factors (void)
{
	struct factor *factor, *next_factor;

	for (factor = factors; factor != NULL; factor = next_factor) {
		next_factor = factor->next;
		free (factor);
	}
	factors = NULL;
	num_factors = 0;
}

static PyObject *
factor_factor (PyObject *self, PyObject *args)
{
	PyObject *factor_list;
	unsigned long long *factor_array;
        int n, i;

        if (!PyArg_ParseTuple (args, "i", &n)) {
                return NULL;
        }

	if (n < 0) {
		PyErr_SetString (PyExc_ValueError, "Can only factor "
				 "nonnegative numbers.");
		return NULL;
	}

	factor_array = get_factors (n);

	factor_list = PyList_New (num_factors);
	for (i = 0; i < num_factors; i++) {
		if (factor_array[i] < LONG_MAX) {
			PyList_SetItem (factor_list, num_factors - (i + 1),
					Py_BuildValue ("k", factor_array[i]));
		} else {
			PyList_SetItem (factor_list, num_factors - (i + 1),
					Py_BuildValue ("K", factor_array[i]));
		}
	}

	free_factors ();

        return factor_list;
}

static PyMethodDef FactorMethods[] = {
        {"factor", factor_factor, METH_VARARGS,
	 "factor(x): Return the prime factors of x. Currently fails if x >= 2^31"},
        {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initfactor (void)
{
        (void) Py_InitModule3 ("factor", FactorMethods, DOCSTRING);
}
