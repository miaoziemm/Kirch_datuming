

#ifndef se_sort_type
#error se_sort_type should be defined
#endif
#ifndef se_sort_name
#error se_sort_type should be defined
#endif


#define QS_SMALL_SIZE  7
#define QS_MEDIUM_SIZE 40


/** Swaps x[a] with x[b]. */
static void se_sort_name(_se_sort_swap)(se_sort_type* x, size_t a, size_t b)
{
    se_sort_type t = x[a];
    x[a] = x[b];
    x[b] = t;
}

/** Swaps x[a .. (a+n-1)] with x[b .. (b+n-1)]. */
static void se_sort_name(_se_sort_vecswap)(se_sort_type* x, size_t a, size_t b, size_t n)
{
    size_t i;
    for(i = 0; i < n; ++i, ++a, ++b) {
        se_sort_type t = x[a];
        x[a] = x[b];
        x[b] = t;
    }
}

/** Returns the index of the median of the three indexed elements. */
static size_t se_sort_name(_se_sort_med3)(se_sort_type* x, size_t a, size_t b, size_t c)
{
    se_sort_type xa = x[a], xb = x[b], xc = x[c];
    if(se_sort_lt(xa, xb)) {
        if(se_sort_lt(xb, xc)) return b;
        if(se_sort_lt(xa, xc)) return c;
        return a;
    } else {
        if(se_sort_lt(xc, xb)) return b;
        if(se_sort_lt(xc, xa)) return c;
        return a;
    }
}


void se_sort_name(se_insert_sort)(se_sort_type* x, size_t off, size_t len)
{
    size_t i, j;
    size_t last = off + len;
    for(i = off+1; i < last; i++) {
        for(j = i; j > off && se_sort_gt(x[j-1], x[j]); --j) {
            se_sort_name(_se_sort_swap)(x, j, j-1);
        }
    }
}


void se_sort_name(se_quick_sort)(se_sort_type* x, size_t off, size_t len)
{
    ssize_t m, l, n;
    ssize_t a, b, c, d;
    ssize_t s;
    se_sort_type v;

    /* ASSERT( (off+len) <= ((unsigned long long int)1<<(sizeof(size_t)*8 - 1)) - 1); */

    /* Insertion sort on smallest arrays */
    if (len < QS_SMALL_SIZE) {
        se_sort_name(se_insert_sort)(x, off, len);
        return;
    }

    /* Choose a partition element, v */
    l = off;
    n = off + len - 1;
    m = off + (len >> 1); /* middle element */
    if (len > QS_MEDIUM_SIZE) { /* Big arrays, pseudomedian of 9 */
        s = len/8;
        l = se_sort_name(_se_sort_med3)(x, l,     l+s, l+2*s);
        m = se_sort_name(_se_sort_med3)(x, m-s,   m,   m+s);
        n = se_sort_name(_se_sort_med3)(x, n-2*s, n-s, n);
    }
    m = se_sort_name(_se_sort_med3)(x, l, m, n); /* Mid-size, med of 3 */
  
    v = x[m];


    /*  Establish Invariant: (==v)* (<v)* (>v)* (==v)*  */
    a = off;
    b = a;
    c = off + len - 1;
    d = c;

    for(;;) {
        while (b <= c && (!(se_sort_lt(v, x[b])))) { /* x[b] <= v <=> !(v < x[b]) */
            if (se_sort_eq(x[b], v))
                se_sort_name(_se_sort_swap)(x, a++, b);
            ++b;
        }

        while (c >= b && (!(se_sort_lt(x[c], v)))) { /* x[c] >= v <=> !(x[c] < v)*/
            if (se_sort_eq(x[c], v))
                se_sort_name(_se_sort_swap)(x, c, d--);
            --c;
        }

        if (b > c)
            break;
    
        se_sort_name(_se_sort_swap)(x, b++, c--);
    }


    /* Swap partition elements back to middle */
    n = off + len;
    s = minssz(a-off, b-a  );  se_sort_name(_se_sort_vecswap)(x, off, b-s, s);
    s = minssz(d-c,   n-d-1);  se_sort_name(_se_sort_vecswap)(x, b,   n-s, s);

  
    /* Recursively sort non-partition-elements */
    if ( (s = b-a) > 1 )
        se_sort_name(se_quick_sort)(x, off, s);

    if ( (s = d-c) > 1 )
        se_sort_name(se_quick_sort)(x, n-s, s);
}


void se_sort_name(se_bubble_sort)(se_sort_type* x, size_t off, size_t len)
{
    int sorted = 0;
    size_t last = off + len;
  
    for(; !sorted; --last) {
        size_t i;
        sorted = 1;
        for(i = off+1; i < last; ++i) {
            if( se_sort_lt(x[i], x[i-1]) ) {
                se_sort_name(_se_sort_swap)(x, i, i-1);
                sorted = 0;
            }
        }
    }
}

/**
    a - the array to be searched
    from_index - the index of the first element (inclusive) to be searched
    to_index - the index of the last element (exclusive) to be searched
    key - the value to be searched for
*/
ssize_t se_sort_name(se_binary_search)(se_sort_type key,
                                         const se_sort_type* a,
                                         size_t from, size_t to)
{
    ssize_t low = from;
    ssize_t high = to - 1;
    ASSERT(from <= to);
    /* ASSERT(from <= ((unsigned long long int)1<<(sizeof(size_t)*8 - 1)) - 1); */
    /* ASSERT(to <= ((unsigned long long int)1<<(sizeof(size_t)*8 - 1)) - 1); */

    while (low <= high) {
        ssize_t mid = (low + high) / 2;
        se_sort_type mid_val = a[mid];
    
        if( se_sort_lt(mid_val, key) ) {
            low = mid + 1;
        } else if( se_sort_gt(mid_val, key) ) {
            high = mid - 1;
        } else {
            return mid; /* key found */
        }
    }
    return -(low + 1);  /* key not found. */
}

