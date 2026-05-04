#include "../include/se_hash.h"


/**
 * The default initial capacity - MUST be a power of two.
 */
#define HT_INITIAL_CAPACITY 16      //Hash_table_Initial_capacity
#define HT_LOAD_FACTOR 0.75f         //LOAD FACTOR
#define HT_MAXIMUM_CAPACITY (1 << 30)

typedef struct ht_entry_s {
    const char* key;
    size_t keylen;
    void* value;
    int hash;
    struct ht_entry_s* next;
} ht_entry;               //链表

typedef struct se_hash_s {
    /**
     * The table, resized as necessary. Length MUST Always be a
     * power of two.
     */
    ht_entry** table;
    int table_length;

    /**
     * The number of key-value mappings contained in this identity hash map. 关键字
     */
    int size;

    /**
     * The table is rehashed when its size exceeds this threshold.  (The
     * value of this field is (int)(capacity * loadFactor).)
     */
    int threshold;

    /**
     * The load factor for the hashtable.  转载因子
     */
    float load_factor;

} _hash;

/**
   Basic character hash algoztrhm, taken from Python's string hash:
   h = h * 1000003 ^ character, the constant being a prime number.
*/
#define CHAR_HASH(h, c) ((((unsigned int)(h)) * 0xF4243) ^ (unsigned char)(c))

/**
 * Value representing null keys inside tables. NULL keys
 */
static const char* NULL_KEY = "nullkey";
#define MASK_NULL(key)   ( (key) == NULL ? NULL_KEY : (key) )
#define UNMASK_NULL(key) ( (key) == NULL_KEY ? NULL : (key) )

/*********************************************************************/
/******* Internal functions ******************************************/
static int hash( const char* x );
static int eq( const char* x, const char* y, size_t xlen, size_t ylen );
static int index_for( int h, int length );     //调整大小
static void resize( _hash* ht, int new_capacity );
static void add_entry( _hash* ht,
                       int hash,
                       const char* key,
                       void* value,
                       size_t keylen,
                       int index ); //添加


/**
 * Constructs an empty hash table with the default initial capacity
 * (16) and the default load factor (0.75).
 */
se_hash create_hash()
{
    return create_hash_args( HT_INITIAL_CAPACITY, HT_LOAD_FACTOR );
}

/**
 * Constructs an empty hash table with the specified initial //空的hash表
 * capacity and load factor.
 */
se_hash create_hash_args(int initial_capacity, float load_factor)
{
    _hash* hash;
    int capacity;

    if( initial_capacity < 0 )
        ERROR(( "Negative hash initial capacity: %d", initial_capacity ));
    if( initial_capacity > HT_MAXIMUM_CAPACITY )
        initial_capacity = HT_MAXIMUM_CAPACITY;

    if( load_factor <= 0 || !std::isfinite((double)load_factor) )
        ERROR(( "Illegal load factor: %f", load_factor ));


    hash = (_hash*)malloc(sizeof(*hash));
    memset(hash, 0, sizeof(*hash)); //清零

    /* Find a power of 2 >= initial_capacity */ //非2的指数的话 补充 即移位
    capacity = 1;
    while( capacity < initial_capacity )
        capacity <<= 1;

    hash->load_factor = load_factor;
    hash->threshold = (int)( capacity * load_factor );
    hash->table = (ht_entry**) malloc( capacity * sizeof(ht_entry*) );
    memset(hash->table, 0, capacity * sizeof(ht_entry*)); //清零
    hash->table_length = capacity;

    return (se_hash)hash;
}

/**
 * Destroys the hash map;
 */
void destroy_hash( se_hash _ht )
{
    destroy_hash_and_entries(_ht, 0);
}

void destroy_hash_and_entries( se_hash _ht, int free_entries )
{
    destroy_hash_and_keyval(_ht, free_entries, free_entries);
}

void destroy_hash_and_keyval( se_hash _ht, int free_keys, int free_values )
{
    _hash* ht = (_hash*)_ht;
    int i;

    for( i = 0; i < ht->table_length; i++) {
        ht_entry* e = ht->table[i];
        while(e != NULL) {
            ht_entry* next = e->next;

            if( free_values && e->value != NULL ) free(e->value);
            if( free_keys   && e->key   != NULL ) free((void*)e->key);

            free(e);
            e = next;
        }
    }

    free(ht->table);
    free(ht);
}

/**
 * Returns the value to which the specified key is mapped in the h identity
 * hash map, or NULL if the map contains no mapping for this key.
 * A return value of NULL does not necessarily indicate that the map
 * contains no mapping for the key; it is also possible that the map
 * explicitly maps the key to NULL. The ht_contains_key function may be
 * used to distinguish these two cases.
 */
void* ht_get(se_hash _ht, const char* key)
{
    const char* k = MASK_NULL(key);
    int h = hash(k);
    return ht_get_hash(_ht, k, strlen(k), h);
}

void* ht_get_hash(se_hash _ht, const char* key, size_t keylen, int h)
{
    _hash* ht = (_hash*)_ht;
    int i = index_for( h, ht->table_length );

    ht_entry* e = ht->table[i];

    for(;;) {
        if( e == NULL )
            return e;
        if( e->hash == h && eq(e->key, key, e->keylen, keylen) )
            return e->value;
        e = e->next;
    }
}


/**
 * Returns true (1) if this map contains a mapping for the
 * specified key.
 */
int ht_contains_key(se_hash _ht, const char* key)
{
    const char* k = MASK_NULL(key);
    int h = hash(k);
    return ht_contains_key_hash(_ht, k, strlen(k), h);
}

int ht_contains_key_hash(se_hash _ht, const char* key, size_t keylen, int h)
{
    _hash* ht = (_hash*)_ht;
    int i = index_for( h, ht->table_length );

    ht_entry* e = ht->table[i];

    while( e != NULL ) {
        if( e->hash == h && eq(e->key, key, e->keylen, keylen) )
            return 1;
        e = e->next;
    }

    return 0;
}

/**
 * Associates the specified value with the specified key in this map.
 * If the map previously contained a mapping for this key, the old
 * value is replaced.
 *
 * Returns previous value associated with specified key, or NULL
 * if there was no mapping for key. A NULL return can also indicate
 * that the hash map previously associated NULL with the specified key.
 */
char* ht_put( se_hash _ht, const char* key, void* value )
{
    const char* k = MASK_NULL(key);
    int h = hash(k);
    return ht_put_hash( _ht, k, strlen(k), h, value );
}

char* ht_put_hash( se_hash _ht, const char* key, size_t keylen, int h, void* value )
{
    _hash* ht = (_hash*)_ht;
    int i = index_for( h, ht->table_length );
    ht_entry* e;

    for( e = ht->table[i]; e != NULL; e = e->next ) {
        if( e->hash == h && eq(e->key, key, e->keylen, keylen) ) {
            void* old_value = e->value;
            e->value = value;
            return (char*)old_value;
        }
    }

    add_entry(ht, h, key, value, keylen, i);

    return NULL;
}

/**
 * Removes the mapping for this key from this map if present.
 * Returns previous value associated with specified key, or NULL
 * if there was no mapping for key. A NULL return can also indicate
 * that the map previously associated NULL with the specified key.
 */
void* ht_remove( se_hash _ht, const char* key )
{
    const char* k = MASK_NULL(key);
    int h = hash(k);
    return ht_remove_hash(_ht, k, strlen(k), h);
}

void* ht_remove_hash( se_hash _ht, const char* key, size_t keylen, int h)
{
    _hash* ht = (_hash*)_ht;
    int i = index_for( h, ht->table_length );
    ht_entry* prev = ht->table[i];
    ht_entry* e = prev;

    while( e != NULL ) {
        ht_entry* next = e->next;
        if( e->hash == h && eq(e->key, key, e->keylen, keylen) ) {
            ht->size--;
            if( prev == e ) ht->table[i] = next;
            else            prev->next = next;
            break;
        }

        prev = e;
        e = next;
    }

    if( e != NULL ) {
        void* val = e->value;
        free(e);
        return val;
    }
    else {
        return NULL;
    }
}

double ht_get_dispersion( se_hash _ht )
{
    _hash* ht = (_hash*)_ht;
    int i;
    int nnull = 0;

    for( i = 0; i < ht->table_length; i++) {
        if( ht->table[i] != NULL ) nnull++;
    }

    return (double)nnull / (double)ht->size ;
}

/**
 * Returns the number of entries from the hash table
 **/
int ht_size( se_hash _ht )
{
    _hash* ht = (_hash*)_ht;
    return ht->size;
}

/**
 * Returns a hash value for the specified string.
 */
static int hash(const char* x)
{
    unsigned h = 0;

    while( *x )
        h = CHAR_HASH(h, *x++);

//#define CHAR_HASH(h, c) ((((unsigned int)(h)) * 0xF4243) ^ (unsigned char)(c))

    /* this improves dispersion */
    h += ~(h << 9);
    h ^=  (h >> 14);
    h +=  (h << 4);
    h ^=  (h >> 10);

    return (int)h;
}

/**
 * Check for equality of non-null reference x and possibly-null y.
 */
static int eq(const char* x, const char* y, size_t xlen, size_t ylen)
{
    size_t i;
    if(x == y) return 1;
    if(xlen != ylen) return 0;
    for(i = 0; i < xlen; i++) {
        if(x[i] != y[i]) return 0;
    }
    return 1;
}


/**
 * Returns index for hash code h.
 */
static int index_for(int h, int length)
{
    return ( h & (length-1) );
}

/**
 * Rehashes the contents of this map into a new array with a
 * larger capacity. This method is called automatically when the
 * number of keys in this map reaches its threshold.
 *
 * The new_capacity the new capacity, MUST be a power of two;
 * it must be greater than current capacity.
 *
 * If current capacity is HT_MAXIMUM_CAPACITY, this method does not
 * resize the map, but but sets threshold to 0x7FFFFFFF (max int value).
 * This has the effect of preventing future calls.
 */
static void resize( _hash* ht, int new_capacity )
{
    ht_entry** new_table;
    int old_capacity = ht->table_length;
    int j, n;

    if( old_capacity == HT_MAXIMUM_CAPACITY ) {
        ht->threshold = 0x7FFFFFFF;
        return;
    }

    new_table = (ht_entry**) malloc( new_capacity * sizeof(ht_entry*) );
    memset(new_table, 0, new_capacity * sizeof(ht_entry*)); //清零



    /* Transfer all entries from current table to new_table. */
    n = ht->table_length;
    for( j = 0; j < n; j++ ) {
        ht_entry* e = ht->table[j];
        if( e != NULL ) {
            do {
                ht_entry* next = e->next;
                int i = index_for( e->hash, new_capacity );
                e->next = new_table[i];
                new_table[i] = e;
                e = next;
            } while( e != NULL );
        }
    }

    free( ht->table );
    ht->table = new_table;
    ht->table_length = new_capacity;
    ht->threshold = (int)( new_capacity * ht->load_factor );
}

/**
 * Add a new entry with the specified key, value and hash code to
 * the specified bucket. It is the responsibility of this
 * method to resize the table if appropriate.
 *
 * Subclass overrides this to alter the behavior of put method.
 */
static void add_entry( _hash* ht,
                       int hash,
                       const char* key,
                       void* value,
                       size_t keylen,
                       int index )
{
    ht_entry* e = (ht_entry*)malloc( sizeof(ht_entry) );
    e->hash = hash;
    e->key = key;
    e->keylen = keylen;
    e->value = value;
    e->next = ht->table[index];
    ht->table[index] = e;

    if( ht->size++ >= ht->threshold )
        resize( ht, 2 * ht->table_length );
}



/* hash table iterator */
typedef struct {
    const _hash* ht;
    int ht_idx;
    const ht_entry* _ht_entry;
} _htiter;

static void _hti_movenext( _htiter* hti );

/**
 * Creates a hash table read-only iterator that points to the first entry
 * from the hash table
 **/
htiter create_htiter( se_hash ht )
{
    _htiter* hti;

    ASSERT(ht);

    hti = (_htiter*) malloc( sizeof(_htiter) );
    memset(hti, 0, sizeof(*hti)); //清零

    hti->ht = (const _hash*) ht;
    hti->ht_idx = -1;
    /* move to the first ht entry (if any) */
    _hti_movenext( hti );

    return (se_hash) hti;
}

/**
 * Destroys a hash table read-only iterator
 **/
void destroy_htiter( htiter _hti )
{
    _htiter* hti = (_htiter*) _hti;

    if (hti)
        free(hti);
}

/**
 * Returns 1 if the hash table iterator has more entries, or 0 otherwise
 **/
int hti_hasnext( htiter _hti )
{
    _htiter* hti = (_htiter*) _hti;

    return hti->ht_idx < hti->ht->table_length && hti->_ht_entry;
}

/**
 * Output the current entry from the hash table and move the iterator forward
 **/
void hti_next( htiter _hti,
                   const char** key, /* hash table entry's key, can be NULL */
                   void** value )   /* hash table entry's value, can be NULL */
{
    _htiter* hti = (_htiter*) _hti;

    ASSERT( hti_hasnext(hti) );

    /* output current entry */
    if( key )
        *key = hti->_ht_entry->key;
    if( value )
        *value = hti->_ht_entry->value;

    /* move next the iterator */
    _hti_movenext( hti );
}

static void _hti_movenext( _htiter* hti )
{
    ASSERT(hti);
    ASSERT(hti->ht_entry || hti->ht_idx == -1);

    if( hti->_ht_entry && hti->_ht_entry->next )
        hti->_ht_entry = hti->_ht_entry->next;
    else {
        for( ++hti->ht_idx;
             hti->ht_idx < hti->ht->table_length;
             ++hti->ht_idx ) {
            if( hti->ht->table[hti->ht_idx] ) {
                hti->_ht_entry = hti->ht->table[hti->ht_idx];
                return;
            }
        }
        hti->_ht_entry = NULL;
    }
}
