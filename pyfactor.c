#include <Python.h>
#include <gmp.h>

#define _(x) (x)
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

/* The maximum number of factors, including -1, for negative numbers.  */
#define MAX_N_FACTORS (sizeof (uintmax_t) * CHAR_BIT)

static int verbose = 0;

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

/* Bound on buffer size needed to represent an integer type or expression T,
   including the terminating null.  */
#define INT_BUFSIZE_BOUND(t) (INT_STRLEN_BOUND (t) + 1)

#include "wheel-size.h"  /* For the definition of WHEEL_SIZE.  */
static const unsigned char wheel_tab[] =
{
#include "wheel.h"
};
#define WHEEL_START (wheel_tab + WHEEL_SIZE)
#define WHEEL_END (wheel_tab + ARRAY_CARDINALITY (wheel_tab))

enum strtol_error
  {
    LONGINT_OK = 0,

    /* These two values can be ORed together, to indicate that both
       errors occurred.  */
    LONGINT_OVERFLOW = 1,
    LONGINT_INVALID_SUFFIX_CHAR = 2,

    LONGINT_INVALID_SUFFIX_CHAR_WITH_OVERFLOW = (LONGINT_INVALID_SUFFIX_CHAR
                                                 | LONGINT_OVERFLOW),
    LONGINT_INVALID = 4
  };
typedef enum strtol_error strtol_error;

#define _DECLARE_XSTRTOL(name, type) \
  strtol_error name (const char *, char **, int, type *, const char *);
_DECLARE_XSTRTOL (xstrtol, long int)
_DECLARE_XSTRTOL (xstrtoul, unsigned long int)
_DECLARE_XSTRTOL (xstrtoimax, intmax_t)
_DECLARE_XSTRTOL (xstrtoumax, uintmax_t)

static mpz_t *factor = NULL;
static size_t nfactors_found = 0;
static size_t nfactors_allocated = 0;

static void
debug (char const *fmt, ...)
{
	if (verbose)
	{
		va_list ap;
		va_start (ap, fmt);
		vfprintf (stderr, fmt, ap);
		va_end (ap);
	}
}

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
xalloc_die (void)
{
	fprintf (stderr, "ERROR: Memory exhausted\n");
	exit (-1);
  /* The `noreturn' cannot be given to error, since it may return if
     its first argument is 0.  To help compilers understand the
     xalloc_die does not return, call abort.  Also, the abort is a
     safety feature if exit_failure is 0 (which shouldn't happen).  */
	abort ();
}

/* Convert I to a printable string in BUF, which must be at least
   INT_BUFSIZE_BOUND (INTTYPE) bytes long.  Return the address of the
   printable string, which need not start at BUF.  */

char *
umaxtostr (uintmax_t i, char *buf)
{
	char *p = buf + INT_STRLEN_BOUND (uintmax_t);
	*p = 0;

#if inttype_is_signed
	if (i < 0)
	{
		do
			*--p = '0' - i % 10;
		while ((i /= 10) != 0);

		*--p = '-';
	}
	else
#endif
	{
		do
			*--p = '0' + i % 10;
		while ((i /= 10) != 0);
	}

	return p;
}

static size_t
factor_wheel (uintmax_t n0, size_t max_n_factors, uintmax_t *factors)
{
	uintmax_t n = n0, d, q;
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
			factors[n_factors++] = d;
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
		factors[n_factors++] = n;
	}

	return n_factors;
}

/* Single-precision factoring */
static void
print_factors_single (uintmax_t n)
{
	uintmax_t factors[MAX_N_FACTORS];
	size_t n_factors = factor_wheel (n, MAX_N_FACTORS, factors);
	size_t i;
	char buf[INT_BUFSIZE_BOUND (uintmax_t)];

	printf ("%s:", umaxtostr (n, buf));
	for (i = 0; i < n_factors; i++)
		printf (" %s", umaxtostr (factors[i], buf));
	putchar ('\n');
}

static void
emit_factor (mpz_t n)
{
	if (nfactors_found == nfactors_allocated)
		factor = xrealloc(factor, nfactors_allocated * sizeof(*factor));
//		factor = X2NREALLOC (factor, &nfactors_allocated);
	mpz_init (factor[nfactors_found]);
	mpz_set (factor[nfactors_found], n);
	++nfactors_found;
}

static void
emit_ul_factor (unsigned long int i)
{
	mpz_t t;
	mpz_init (t);
	mpz_set_ui (t, i);
	emit_factor (t);
	mpz_clear (t);
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

	debug ("[trial division (%u)] ", limit);

	mpz_init (q);
	mpz_init (r);

	f = mpz_scan1 (t, 0);
	mpz_div_2exp (t, t, f);
	while (f)
	{
		emit_ul_factor (2);
		--f;
	}

	for (;;)
	{
		mpz_tdiv_qr_ui (q, r, t, 3);
		if (mpz_cmp_ui (r, 0) != 0)
			break;
		mpz_set (t, q);
		emit_ul_factor (3);
	}

	for (;;)
	{
		mpz_tdiv_qr_ui (q, r, t, 5);
		if (mpz_cmp_ui (r, 0) != 0)
			break;
		mpz_set (t, q);
		emit_ul_factor (5);
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
			emit_ul_factor (f);
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

	debug ("[pollard-rho (%d)] ", a_int);

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

			debug ("[composite factor--restarting pollard-rho] ");
			factor_using_pollard_rho (g, a_int);
		}
		else
		{
			emit_factor (g);
		}
		mpz_mod (x, x, n);
		mpz_mod (x1, x1, n);
		mpz_mod (y, y, n);
		if (mpz_probab_prime_p (n, 3))
		{
			emit_factor (n);
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
	mpz_t *const *a = av;
	mpz_t *const *b = bv;
	return mpz_cmp (**a, **b);
}

static void
sort_and_print_factors (void)
{
	mpz_t **faclist;
	size_t i;

	faclist = xcalloc (nfactors_found, sizeof *faclist);
	for (i = 0; i < nfactors_found; ++i)
	{
		faclist[i] = &factor[i];
	}
	qsort (faclist, nfactors_found, sizeof *faclist, mpcompare);

	for (i = 0; i < nfactors_found; ++i)
	{
		fputc (' ', stdout);
		mpz_out_str (stdout, 10, *faclist[i]);
	}
	putchar ('\n');
	free (faclist);
}

static void
free_factors (void)
{
	size_t i;

	for (i = 0; i < nfactors_found; ++i)
	{
		mpz_clear (factor[i]);
	}
	/* Don't actually free factor[] because in the case where we are
	   reading numbers from stdin, we may be about to use it again.  */
	nfactors_found = 0;
}

/* Arbitrary-precision factoring */
static void
print_factors_multi (mpz_t t)
{
	mpz_out_str (stdout, 10, t);
	putchar (':');

	if (mpz_sgn (t) != 0)
	{
		/* Set the trial division limit according to the size of t.  */
		size_t n_bits = mpz_sizeinbase (t, 2);
		unsigned int division_limit = MIN (n_bits, 1000);
		division_limit *= division_limit;

		factor_using_division (t, division_limit);

		if (mpz_cmp_ui (t, 1) != 0)
		{
			debug ("[is number prime?] ");
			if (mpz_probab_prime_p (t, 3))
				emit_factor (t);
			else
				factor_using_pollard_rho (t, 1);
		}
	}

	mpz_clear (t);
	sort_and_print_factors ();
	free_factors ();
}

/* This is derived from the function used in the factor program provided by
   coreutils.

   Emit the factors of the indicated number.  If we have the option of using
   either algorithm, we select on the basis of the length of the number.
   For longer numbers, we prefer the MP algorithm even if the native algorithm
   has enough digits, because the algorithm is better.  The turnover point
   depends on the value.  */
static int
print_factors (char const *s)
{
	uintmax_t n;
	strtol_error err = xstrtoumax (s, NULL, 10, &n, "");

	enum { GMP_TURNOVER_POINT = 100000 };

	if (err == LONGINT_OVERFLOW
	    || (err == LONGINT_OK && GMP_TURNOVER_POINT <= n)) {
		mpz_t t;
		mpz_init (t);
		if (gmp_sscanf (s, "%Zd", t) == 1) {
			debug ("[%s]", _("using arbitrary-precision arithmetic"));
			print_factors_multi (t);
			return 1;
		}
		err = LONGINT_INVALID;
	}

	switch (err)
	{
	case LONGINT_OK:
		debug ("[%s]", _("using single-precision arithmetic"));
		print_factors_single (n);
		return 1;

	case LONGINT_OVERFLOW:
		fprintf (stderr, "ERROR: \"%s\" is too large\n", s);
		exit (0);
		return 0;

	default:
		fprintf (stderr, "ERROR: \"%s\" is not a valid positive "
			 "integer\n", s);
		exit (0);
		return 0;
	}
}


static PyObject *
pyfactor_factor (PyObject *self, PyObject *args)
{
	PyObject *factors;
        int n;

        if (!PyArg_ParseTuple (args, "i", &n)) {
                return NULL;
        }

	factors = PyList_New (1);

	PyList_SetItem (factors, 0, Py_BuildValue ("i", n));

	PyList_Append (factors, Py_BuildValue ("i", 1));

        return factors;
}

static PyMethodDef PyfactorMethods[] = {
        {"factor", pyfactor_factor, METH_VARARGS,
	 "Return the prime factors of x."},
        {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initpyfactor (void)
{
        (void) Py_InitModule ("pyfactor", PyfactorMethods);
}
