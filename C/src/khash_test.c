//
// Created by Hood on 5/17/2019.
//

// To run this program: `echo a bc a cd bc|./this_prog`
#include <stdio.h>
#include <string.h>
#include "khash.h"
KHASH_SET_INIT_STR(str)



int main(int argc, char *argv[])
{
    char s[4096]; // max string length: 4095 characters
    khash_t(str) *h;
    khint_t k;
    h = kh_init(str);
    while (scanf("%s", s) > 0) {
        int absent;
        k = kh_put(str, h, s, &absent);
        if (absent) kh_key(h, k) = strdup(s);
        // else, the key is not touched; we do nothing
    }
    printf("# of distinct words: %d\n", kh_size(h));
    // IMPORTANT: free memory allocated by strdup() above
    for (k = 0; k < kh_end(h); ++k)
        if (kh_exist(h, k))
            free((char*)kh_key(h, k));
    kh_destroy(str, h);
    return 0;
}

#include "khash.h"
KHASH_MAP_INIT_INT(int_to_char_map, char)      // instantiate structs and methods
int main() {
    int absent, is_missing;
    khint_t k;
    khash_t(int_to_char_map) *h = kh_init(int_to_char_map);  // allocate a hash table
    k = kh_put(int_to_char_map, h, 5, &absent);  // insert a key to the hash table
    if (!ret) kh_del(int_to_char_map, h, k);
    kh_value(h, k) = 10;             // set the value
    k = kh_get(int_to_char_map, h, 10);          // query the hash table
    is_missing = (k == kh_end(h));   // test if the key is present
    k = kh_get(int_to_char_map, h, 5);
    kh_del(int_to_char_map, h, k);               // remove a key-value pair
    for (k = kh_begin(h); k != kh_end(h); ++k) {  // traverse
        if (kh_exist(h, k)) { // test if a bucket contains data
            kh_value(h, k) = 1;
        }
    }
    kh_destroy(int_to_char_map, h);              // deallocate the hash table
    return 0;
}