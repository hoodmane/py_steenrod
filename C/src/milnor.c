//
// Created by Hood on 5/8/2019.
//

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

#include "khash.h"

#include "combinatorics.h"
#include "FpVector.h"
#include "milnor.h"
#include "milnor_private.h"

// Milnor Algebra
// A MilnorAlgebra object has a prime and is either generic or not generic.
// It also contains a profile which specifies some subHopfalgebra of the Steenrod algebra.
// If the Profile parameter is null, then the full Steenrod algebra is used.
//
// The MilnorAlgebra allocates tables to represent the basis. In order to use, first you must call
// MilnorAlgebra_generateBasis(max_degree). Then the methods _getBasis, getDimension, and multiply 
// implement the algebra structure.
MilnorAlgebra *MilnorAlgebra__initializeFields(MilnorAlgebraInternal *algebra, uint p, bool generic, Profile *profile);
MilnorAlgebra *MilnorAlgebra_construct(uint p, bool generic, Profile *profile){
    uint profile_ppart_length = 0;
    if(profile != NULL){
        profile_ppart_length = profile->p_part_length + 1;
    }
    uint num_products;
    if(generic){
        num_products = 2;
    } else {
        num_products = 3;
    }
    size_t algebra_size =         
        sizeof(MilnorAlgebraInternal) 
        + profile_ppart_length * sizeof(uint) // Store the P list part of the profile here
        + sizeof(FiltrationOneProductList) // These are for listing out the products.
            + 2 * num_products * sizeof(uint); // more product space
    MilnorAlgebraInternal *algebra = malloc(algebra_size);

    // This next assignment gets written over in initializer so we'll need to do it again...
    algebra->public_algebra.profile.p_part = profile_ppart_length > 0 ? (uint*)(algebra + 1) : NULL;
    Algebra *inner_algebra = &algebra->public_algebra.algebra; 
    inner_algebra->product_list = (FiltrationOneProductList*)((uint*)(algebra + 1) + profile_ppart_length);
    inner_algebra->product_list->degrees = (int*)(inner_algebra->product_list + 1);
    inner_algebra->product_list->indices = (uint*)inner_algebra->product_list->degrees + num_products;

    assert((char*)(inner_algebra->product_list->indices + num_products) == ((char *)algebra) + algebra_size);
    MilnorAlgebra__initializeFields(algebra, p, generic, profile);
    return (MilnorAlgebra*)algebra;
}

MilnorAlgebra *MilnorAlgebra__initializeFields(MilnorAlgebraInternal *algebra, uint p, bool generic, Profile *profile){
    initializePrime(p);
    algebra->public_algebra.algebra.p = p;
    algebra->public_algebra.generic = generic;
    // Initialize Profile
    if(profile == NULL){
        algebra->public_algebra.profile.restricted = false;
        algebra->public_algebra.profile.truncated = false;
        algebra->public_algebra.profile.p_part_length = 0;
        algebra->public_algebra.profile.q_part = -1;
        algebra->public_algebra.profile.p_part_length = 0;
        algebra->public_algebra.profile.p_part = NULL;
    } else {
        algebra->public_algebra.profile = *profile; // Oh no, this writes over p_part.
        algebra->public_algebra.profile.p_part = (uint*)(algebra + 1); // Fix it.
        memcpy(algebra->public_algebra.profile.p_part, profile->p_part, profile->p_part_length * sizeof(uint));
    }
    algebra->public_algebra.profile.generic = generic;
    // Fill in the Algebra function pointers.
    algebra->public_algebra.algebra.computeBasis = MilnorAlgebra_generateBasis;
    algebra->public_algebra.algebra.getDimension = MilnorAlgebra_getDimension;
    algebra->public_algebra.algebra.multiplyBasisElements = MilnorAlgebra_multiply;
    algebra->public_algebra.sort_order = NULL;

    // Basis tables
    algebra->public_algebra.algebra.max_degree = 0;

    algebra->P_table = NULL;
    algebra->P_table_by_P_length = NULL;
    algebra->P_table_max_degree = 0;

    algebra->Q_table = NULL;
    algebra->Q_table_max_tau = 0;

    algebra->basis_table = NULL;
    algebra->basis_element_to_index_map = NULL;

    // Products
    // Length field has to match with amount of space we decided to allocate for this
    // in constructor function directly before this.
    FiltrationOneProductList *product_list = algebra->public_algebra.algebra.product_list;
    if(generic){
        product_list->length = 2;
        product_list->degrees[0] = 1; // beta
        product_list->indices[0] = 0;        
        product_list->degrees[1] = 2 * (p - 1); // P1
        product_list->indices[1] = 0;
    } else {
        product_list->length = 3;
        product_list->degrees[0] = 1; // Sq1
        product_list->indices[0] = 0; 
        product_list->degrees[1] = 2; // Sq2
        product_list->indices[1] = 0; 
        product_list->degrees[2] = 4; // Sq4
        product_list->indices[2] = 0;  
    }
    MilnorAlgebra__generateName((MilnorAlgebra*)algebra);
    return (MilnorAlgebra*)algebra;
}

void MilnorAlgebra_free(MilnorAlgebra *algebra){
    MilnorAlgebra_freeBasis(algebra);
    free(algebra);
}

uint MilnorAlgebra_getDimension(Algebra *this, int degree, int excess __attribute__((unused))) {
    MilnorAlgebraInternal *algebra = (MilnorAlgebraInternal *) this;
    assert(algebra->basis_table != NULL);
    assert(degree < this->max_degree);
    if(degree < 0){
        return 0;
    }
    return algebra->basis_table[degree].length;
}

MilnorBasisElement_list MilnorAlgebra_getBasis(MilnorAlgebra *public_algebra, int degree){
    MilnorAlgebraInternal *algebra = (MilnorAlgebraInternal *) public_algebra;
    assert(algebra->basis_table != NULL && degree < public_algebra->algebra.max_degree);
    return algebra->basis_table[degree];
}




// TODO: Break this up into several methods.
// Generate the PTable through degree new_max_degree.
// This is only needed in order to help generate the basis table, so it can in principle be deallocated immediately afterwards.
// However, we keep it around so that if we need to extend the table we don't have to redo our work.
void MilnorAlgebra__generatePpartTable(MilnorAlgebraInternal *algebra, int new_max_degree) {
    uint p = algebra->public_algebra.algebra.p;
    uint *xi_degrees = getXiDegrees(p);

    uint profile_list[MAX_XI_TAU];
    for(uint idx = 0;  idx < MAX_XI_TAU; idx ++){
        profile_list[idx] = Profile__getExponent(algebra->public_algebra.profile, p, idx) - 1;
    }
    uint old_max_degree = algebra -> P_table_max_degree;
    if(new_max_degree < old_max_degree){
        return;
    }
    algebra -> P_table_max_degree = new_max_degree;
    algebra -> P_table = realloc(algebra -> P_table , new_max_degree * sizeof(P_part_list));
    algebra -> P_table_by_P_length = realloc(algebra -> P_table_by_P_length, new_max_degree * sizeof(P_part_list[MAX_XI_TAU]));
    P_part_list *degree_table = algebra->P_table;
    P_part_list **tau_table = algebra->P_table_by_P_length;

    // If old_max_degree is 0 we need to set up the base case for our recursion.
    // In degree 0 there's a single partition of length 0.
    if(old_max_degree == 0){
        old_max_degree++;
        uint *p_part = calloc(1, sizeof(uint));
        // In order to deallocate correctly, the zero dimensional part has to be the same as the rest.
        // degree_table and tau_table share P_parts, each P_part appears once in degree_table and once in tau_table/
        // degree_table and tau_table DO NOT share P_part_lists, so even though the same list appears in each,
        // we need to make sure to allocate two separate ones to avoid double free on teardown.
        degree_table[0].length = 1;
        degree_table[0].list = malloc(sizeof(P_part));
        degree_table[0].list->degree = 0;
        degree_table[0].list->length = 0;
        degree_table[0].list->p_part = p_part;

        tau_table[0] = (P_part_list*) calloc(10, sizeof(P_part_list));
        tau_table[0][0].length = 1;
        tau_table[0][0].list = malloc(sizeof(P_part));
        tau_table[0][0].list->degree = 0;
        tau_table[0][0].list->length = 0;
        tau_table[0][0].list->p_part = p_part;
    }

    // Populate the new table.
    // We iterate over each degree and for each degree n iterate over each xi_i until we hit one of degree larger
    // than n and look for monomials of degree n - |xi_i|. Copy each such monomial, increment the power of xi_i
    // and add it to the degree n component.

    P_part degree_list_buffer[MAX_DIMENSION];
    uint degree_list_length;
    P_part tau_list_buffer[MAX_DIMENSION];
    uint tau_list_length;
    for(int current_degree = old_max_degree; current_degree < new_max_degree; current_degree++){
        degree_list_length = 0;
        tau_table[current_degree] = (P_part_list*) calloc(MAX_XI_TAU, sizeof(P_part_list));

        // Now populate the tables.
        // Loop over the highest degree xi_i in the new monomial
        for(uint xi = 0; xi < MAX_XI_TAU; xi ++){
            uint xi_degree = xi_degrees[xi];
            if(xi_degree > current_degree){
                break;
            }
            int remaining_degree = current_degree - xi_degree;
            tau_list_length = 0;
            // Iterate over the highest xi_i in the old monomial
            for(uint xi_old = 0; xi_old <= xi; xi_old ++) {
                P_part_list remaining_degree_list = tau_table[remaining_degree][xi_old];
                // Iterate over the monomials of degree n - |xi| with highest xi xi_old.
                for(uint i = 0; i < remaining_degree_list.length; i++) {
                    P_part old_p_part = remaining_degree_list.list[i];
                    if( profile_list[xi] == 0 
                        || (xi_old == xi && old_p_part.p_part[xi] == profile_list[xi])){
                        continue;
                    }
                    // Copy the old p_part to a new p_part with length xi + 1.
                    P_part new_p_part;
                    new_p_part.degree = current_degree;
                    new_p_part.length = xi + 1;
                    new_p_part.p_part = (uint *) calloc(new_p_part.length, sizeof(uint));
                    memcpy(new_p_part.p_part, old_p_part.p_part, old_p_part.length * sizeof(uint));
                    // Increment the exponent of xi_i
                    new_p_part.p_part[xi]++;
                    // Add our new_p_part to the tables.
                    degree_list_buffer[degree_list_length] = new_p_part;
                    tau_list_buffer[tau_list_length] = new_p_part;
                    degree_list_length++;
                    tau_list_length++;
                }
            }
            tau_table[current_degree][xi].length = tau_list_length;
            tau_table[current_degree][xi].list = malloc(tau_list_length * sizeof(P_part));
            memcpy(tau_table[current_degree][xi].list, tau_list_buffer, tau_list_length * sizeof(P_part));
        }
        degree_table[current_degree].length = degree_list_length;
        degree_table[current_degree].list = malloc(degree_list_length * sizeof(P_part));
        memcpy(degree_table[current_degree].list, degree_list_buffer, degree_list_length * sizeof(P_part));
    }
}

// Free the PpartTable.
void MilnorAlgebra__freePPartTable(MilnorAlgebraInternal *algebra){
    if(algebra->P_table == NULL){
        return;
    }
    uint max_degree = algebra -> P_table_max_degree;
    P_part_list *degree_table = algebra->P_table;
    P_part_list **tau_table = algebra->P_table_by_P_length;
    for(int current_degree = 0; current_degree < max_degree; current_degree++){
        // Now populate the tables.
        // Loop over the highest degree xi_i in the new monomial
        for(uint xi = 0; xi < MAX_XI_TAU; xi ++){
            free(tau_table[current_degree][xi].list);
        }
        for(uint idx = 0; idx < degree_table[current_degree].length; idx++) {
            free(degree_table[current_degree].list[idx].p_part);
        }
        free(degree_table[current_degree].list);
        free(tau_table[current_degree]);
    }
    free(algebra->P_table);
    free(algebra->P_table_by_P_length);
    algebra->P_table = NULL;
    algebra->P_table_by_P_length = NULL;
    algebra -> P_table_max_degree = -1;
}


// Generate the MilnorBasisQpartTable through degree n.
// TODO: Break this into smaller parts?
// Generate the QTable through degree new_max_degree.
// This is only needed in order to help generate the basis table, so it could in principle be deallocated immediately afterwards.
// However, we keep it around so that if we need to extend the table we don't have to redo our work.
void MilnorAlgebra__generateQpartTable(MilnorAlgebraInternal *algebra, int new_max_degree){
    uint p = algebra->public_algebra.algebra.p;
    uint q = 2*(p-1);
    uint *tau_degrees = getTauDegrees(p);
    uint profile = algebra->public_algebra.profile.q_part;

    uint old_max_tau = algebra->Q_table_max_tau;
    // New max tau is number of tau_i's less than max_degree.
    uint new_max_tau = 0;
    for( ; tau_degrees[new_max_tau] < new_max_degree; new_max_tau++){}

    if(new_max_tau <= old_max_tau){
        return;
    }
    algebra->Q_table_max_tau = new_max_tau;

    // Initialize map residue --> Q_part_list
    if(old_max_tau == 0) {
        algebra->Q_table = (Q_part_list *) calloc(q, sizeof(Q_part_list));
    }
    uint available_taus = 0;
    for(uint i = 0; i < new_max_tau; i++){
        available_taus += (profile >> i) & 1;
    }

    profile = ~profile;

    // Calculate how many elements are going to land in each bin. There are available_taus choose i monomials in tau_0, ..., tau_(n-1)
    // of length i and these land in residue class i % q.
    // We could replace this with a buffer...
    uint residue_bin_length[q];
    memset(residue_bin_length, 0, q * sizeof(uint));
    Q_part residue_bin_buffer[q][MAX_DIMENSION/q];

    if(old_max_tau == 0) {
        // degree 0 is an edge case because n = 0 is the only integer such that there is not a unique nonnegative
        // exponent i with 2^i <= n < 2^(i+1). It goes in residue 0.
        residue_bin_buffer[0][0] = (Q_part){0, 0};
        residue_bin_length[0] ++;
    }

    uint total = 0;
    //The total starts out as tau_degrees[old_max_tau] but we update by tau_degrees[old_max_tau] - Sum(smaller tau_i's).
    //So initialize total = Sum(smaller tau_i's)
    for(uint i = 0; i < old_max_tau; i++){
        total += tau_degrees[i];
    }
    uint bit_string_min = 1 << old_max_tau;
    uint bit_string_max = 1 << new_max_tau;
    //The residue starts out as 1, but we update by 1 - # of trailing 0's.
    //On the first pass, # of trailing 0's is old_max_tau, so initialize residue = old_max_tau
    int residue = old_max_tau % q;
    for(uint bit_string = bit_string_min;  bit_string < bit_string_max; bit_string++) {
        // Iterate over trailing zero bits
        uint v = (bit_string ^ (bit_string - 1)) >> 1;
        uint c = 0; // c is counting the number of zero bits at the end of bit_string.
        for (; v != 0; c++) {
            v >>= 1;
            total -= tau_degrees[c];
        }
        total += tau_degrees[c];
        // residue is the number of 1's in bit_string mod q.
        // c bits were unset when we incremented and one new bit was set so we have to add 1 - c.
        residue += 1 - c;
        if((bit_string & profile) != 0){
            continue;
        }
        residue = ModPositive(residue, q);
        residue_bin_buffer[residue][residue_bin_length[residue]] = (Q_part){total, bit_string};
        residue_bin_length[residue]++;
    }

    Q_part_list *table = algebra->Q_table;
    for (uint residue = 0; residue < q; residue++) {
        table[residue].list = (Q_part*) realloc(table[residue].list, (table[residue].length + residue_bin_length[residue]) * sizeof(Q_part));
        memcpy(
                table[residue].list + table[residue].length,
                residue_bin_buffer[residue],
                residue_bin_length[residue] * sizeof(Q_part)
        );
        table[residue].length = residue_bin_length[residue];
    }
}

// Free the QpartTable.
// Thankfully a lot less elaborate than freeing PpartTable.
void MilnorAlgebra__freeQPartTable(MilnorAlgebraInternal *algebra){
    uint p = algebra->public_algebra.algebra.p;
    uint q = 2*p - 2;
    Q_part_list *table = algebra->Q_table;
    for(uint residue = 0; residue < q; residue++) {
        free(table[residue].list);
    }
    free(algebra->Q_table);
    algebra->Q_table_max_tau = 0;
    algebra->Q_table = NULL;
}

void MilnorAlgebra__generateBasisElementToIndexMap(MilnorAlgebraInternal *algebra, int old_max_degree, int max_degree);
// Generate the basis for MilnorAlgebra. The main work is already contained in the methods
// that generate the P part and Q part. This is mostly just putting stuff in the right spots.
// We need to make the basis_table which is a list of the basis elements in each degree
// and basis_element_to_index_map which is a khash map from MilnorBasisElements to indices.
// Here we allocate the memory for these structures and dispatch to the generic and 2 cases.
void MilnorAlgebra_generateBasis(Algebra *this, int max_degree) {
    MilnorAlgebraInternal *algebra = (MilnorAlgebraInternal*) this;
    uint p = this->p;
    initializePrime(p);
    uint old_max_degree = this->max_degree;
    algebra->basis_table = (MilnorBasisElement_list *) realloc(
            algebra->basis_table,
            max_degree * sizeof(MilnorBasisElement_list)
    );
    algebra->public_algebra.algebra.max_degree = max_degree;
    algebra->basis_element_to_index_map = (khash_t(monomial_index_map) **)realloc(
            algebra->basis_element_to_index_map,
            max_degree * sizeof(khash_t(monomial_index_map) *)
        );
    khash_t(monomial_index_map) **name_table = algebra->basis_element_to_index_map;
    for(uint k = old_max_degree; k < max_degree; k++){
        name_table[k] = kh_init(monomial_index_map);
    }
    if(algebra->public_algebra.generic){
        MilnorAlgebra__generateBasisGeneric(algebra, old_max_degree, max_degree);
    } else {
        MilnorAlgebra__generateBasis2(algebra, old_max_degree, max_degree);
    }
    MilnorAlgebra__generateBasisElementToIndexMap(algebra, old_max_degree, max_degree);
}

void MilnorAlgebra_freeBasis(MilnorAlgebra *public_algebra){
    MilnorAlgebraInternal *algebra = (MilnorAlgebraInternal*) public_algebra;
    MilnorBasisElement_list *table = algebra->basis_table;
    khash_t(monomial_index_map) **name_table = algebra->basis_element_to_index_map;
    for(int degree = 0; degree < algebra->public_algebra.algebra.max_degree; degree++) {
        free(table[degree].list);
        khint_t bin;
        for (bin = 0; bin < kh_end(name_table[degree]); ++bin) {
            if (kh_exist(name_table[degree], bin)){
                // Free the key objects (the cast to char* casts away a const).
                free((char *) kh_key(name_table[degree], bin));
            }
        }
    }

    MilnorAlgebra__freePPartTable(algebra);
    if(algebra->public_algebra.generic){
        MilnorAlgebra__freeQPartTable(algebra);
    }

    for(uint k = 0; k < algebra->public_algebra.algebra.max_degree; k++){
        kh_destroy(monomial_index_map, algebra->basis_element_to_index_map[k]);
    }
    free(algebra->basis_element_to_index_map);
    free(algebra->basis_table);
}

// In the nongeneric case, we just copy the P_part table into basis_table and
// populate the basis element to index map with the inverse.
void MilnorAlgebra__generateBasis2(MilnorAlgebraInternal *algebra, int old_max_degree, int new_max_degree){
    MilnorAlgebra__generatePpartTable(algebra, new_max_degree);

    MilnorBasisElement_list *table = algebra->basis_table;
    MilnorBasisElement_list current_degree_list;
    for(int degree = old_max_degree; degree < new_max_degree; degree++){
        P_part_list p_parts = algebra->P_table[degree];
        current_degree_list.length = 0;
        current_degree_list.list = (MilnorBasisElement*) malloc(p_parts.length * sizeof(MilnorBasisElement));
        for(uint i = 0; i < p_parts.length; i++){
            P_part x = p_parts.list[i];
            current_degree_list.list[current_degree_list.length] = (MilnorBasisElement){0, 0, degree, x.length, x.p_part};
            current_degree_list.length ++;
        }
        if(algebra->public_algebra.sort_order != NULL){
            printf("sort? %llx\n", (uint64)algebra->public_algebra.sort_order);
            qsort(current_degree_list.list, current_degree_list.length, sizeof(MilnorBasisElement), algebra->public_algebra.sort_order);
        }
        table[degree] = current_degree_list;
    }
}

// Get the basis in degree n for the generic steenrod algebra at the prime p.
// We just put together the "even part" of the basis and the "Q part".
// First we select our Q_part because there is at most one possible Q_part in each degree.
// Our Q_part can vary over any Q_part that has the right residue class mod 2p-2 and
// has degree bounding by the degree we're looking for. Then we take the left over degree
// and add all possible P parts.
void MilnorAlgebra__generateBasisGeneric(MilnorAlgebraInternal *algebra, int old_max_degree, int new_max_degree){
    uint p = algebra->public_algebra.algebra.p;
    uint q = 2 * (p - 1);
    // Generate P-part and Q-part table. Notice that P-part table has degree 
    // divided by the factor of 2p-2.
    MilnorAlgebra__generatePpartTable(algebra, (new_max_degree + q - 1) / q);
    MilnorAlgebra__generateQpartTable(algebra, new_max_degree);
    MilnorBasisElement_list *table = algebra->basis_table;
    for(int degree = old_max_degree; degree < new_max_degree; degree++){
        // Pick Q_part with appropriate residue.
        Q_part_list q_list = algebra->Q_table[degree % q];
        uint degree_list_length = 0;
        // There will eventually be a buffer overflow in a large enough dimension.
        // TODO: check for this?
        MilnorBasisElement degree_list_buffer[MAX_DIMENSION];
        for(uint i = 0; i < q_list.length; i++) {
            Q_part q_part = q_list.list[i];
            uint q_deg = q_part.degree;
            // The q_parts in a given residue class are ordered increasing in degree
            // so if this one is too big so are all later ones.
            if(q_deg > degree){
                break;
            }
            uint q_bit_string = q_part.bit_string;
            uint p_deg = (degree - q_deg);
            P_part_list p_list = algebra->P_table[p_deg / q];
            // Use up the leftover degree in a P part. There are multiple P parts in a 
            // given degree.
            for(uint j = 0; j < p_list.length; j++) {
                P_part p_part = p_list.list[j];
                // Now populate table (index --> basis element)
                // and map (basis element --> index)
                MilnorBasisElement m = (MilnorBasisElement){q_part.degree, q_bit_string, p_deg, p_part.length, p_part.p_part};
                degree_list_buffer[degree_list_length] = m;
                degree_list_length ++;
            }
        }
        if(algebra->public_algebra.sort_order != NULL){
            qsort(degree_list_buffer, degree_list_length, sizeof(MilnorBasisElement), algebra->public_algebra.sort_order);
        }        
        table[degree].length = degree_list_length;
        table[degree].list = malloc(degree_list_length * sizeof(MilnorBasisElement));
        // Copy from stack to heap.
        memcpy(table[degree].list, degree_list_buffer, degree_list_length * sizeof(MilnorBasisElement));
    }
}

void MilnorAlgebra__generateBasisElementToIndexMap(MilnorAlgebraInternal *algebra, int old_max_degree, int max_degree){
    khash_t(monomial_index_map) **name_table = algebra->basis_element_to_index_map;
    MilnorBasisElement_list *basis_table = algebra->basis_table;
    for(int n = old_max_degree; n < max_degree; n++){
        // We somewhat arbitrarily choose a load factor of 0.5 for the table
        // khash will resize the map if the load factor rises above 0.77.
        // At some point maybe we should profile the hash table and find out what the 
        // best load factor is. We could also try using perfect hash...
        kh_resize(monomial_index_map, name_table[n], basis_table[n].length * 2);
        for(uint i=0; i<basis_table[n].length; i++){
            char key[200];
            MilnorAlgebra_basisElement_toKey(key, &basis_table[n].list[i]);
            int absent;
            khint_t bin = kh_put(monomial_index_map, name_table[n], key, &absent);
            assert(absent);
            kh_key(name_table[n], bin) = strdup(key);
            kh_val(name_table[n], bin) = i;
        }
    }
}

MilnorBasisElement MilnorAlgebra_basisElement_fromIndex(MilnorAlgebra *public_algebra, int degree, uint idx) {
    MilnorAlgebraInternal *algebra = (MilnorAlgebraInternal*) public_algebra;
    return algebra->basis_table[degree].list[idx];
}

// Precondition: b is a valid MBE for the algebra -- it actually appears as a key in 
uint MilnorAlgebra_basisElement_toIndex(MilnorAlgebra *public_algebra, MilnorBasisElement b){
    MilnorAlgebraInternal *algebra = (MilnorAlgebraInternal*) public_algebra;
    int degree = b.q_degree + b.p_degree;
    assert(degree < public_algebra->algebra.max_degree);
    kh_monomial_index_map_t *map = algebra->basis_element_to_index_map[degree];
    char key[200];
    MilnorAlgebra_basisElement_toKey(key, &b);
    khint_t bin = kh_get(monomial_index_map, map, key);
    if(bin == kh_end(map)){
        MilnorAlgebra_basisElement_toString(key, public_algebra, &b);
        printf("Uh-oh, not here. degree: %d, elt: '%s'\n", degree, key);
        assert(false);
        return 0;
    }
    uint result = kh_val(map, bin);
    assert(result < algebra->basis_table[degree].length);
    assert(algebra->basis_table[degree].list[result].q_degree == b.q_degree);
    assert(algebra->basis_table[degree].list[result].q_part == b.q_part);    
    assert(algebra->basis_table[degree].list[result].p_degree == b.p_degree);
    assert(algebra->basis_table[degree].list[result].p_length == b.p_length);
    return result;
}

// 
// Multiplication
// 

// Initializes an len(r)+1 by len(s)+1 matrix
// Puts r along the first column and s along the first row and zeroes everywhere else.
void milnor_matrix_initialize(uint M[MAX_XI_TAU][MAX_XI_TAU], P_part r, P_part s) {
    M[0][0] = 0; // shouldn't really matter
    memcpy(M[0] + 1, s.p_part, s.length * sizeof(uint));
    for(uint i = 0; i < r.length; i++){
        M[i+1][0] = r.p_part[i];
        memset(M[i+1] + 1, 0, (s.length)* sizeof(uint));
    }
}

// This seems to move an i x j block of M back to the first row and column.
// To be honest, I don't really know what the point is, but the milnor_matrices
// function was a little long and this seemed like a decent chunk to extract.
// At least it contains all of the steps that modify M so that seems like a good thing.
void milnor_matrix_step_update_matrix(uint M[MAX_XI_TAU][MAX_XI_TAU], P_part r, P_part s, uint i, uint j, uint x)  {
    uint cols = s.length + 1;
    for(uint row = 1; row < i; row ++){
        M[row][0] = r.p_part[row-1];
        for(uint col = 1; col < cols; col++){
            M[0][col] += M[row][col];
            M[row][col] = 0;
        }
    }
    for(uint col = 1; col < j; col++){
        M[0][col] += M[i][col];
        M[i][col] = 0;
    }
    M[0][j] --;
    M[i][j] ++;
    M[i][0] = x;
}



// Generator for Milnor matrices. milnor_product_even iterates over this.
// Uses the same algorithm Monks does in his Maple package to iterate through
// the possible matrices: see
// https://monks.scranton.edu/files/software/Steenrod/steen.html
bool milnor_matrix_step(uint  p, uint M[MAX_XI_TAU][MAX_XI_TAU], P_part r, P_part s) {
    uint rows = r.length + 1;
    uint cols = s.length + 1;
    for(uint i = 1; i < rows; i++){
        uint total = M[i][0];
        uint p_to_the_j = 1;
        for(uint j = 1; j < cols; j++){
            p_to_the_j *= p;
            if(total < p_to_the_j){
                // We don't have enough weight left in the entries above this one in the column to increment this cell.
                // Add the weight from this cell to the total, we can use it to increment a cell lower down.
                total += M[i][j] * p_to_the_j;
                continue;
            }

            // Check if any entry in column j above row i is nonzero. I'm still not sure why tbh.
            for(uint k = 0; k < i; k++){
                if(M[k][j] != 0){
                    // If so, we found our next matrix.
                    milnor_matrix_step_update_matrix(M, r, s, i, j, total - p_to_the_j);
                    return true;
                }
            }
            // All the cells above this one are zero so we didn't find our next matrix.
            // Add the weight from this cell to the total, we can use it to increment a cell lower down.
            total += M[i][j] * p_to_the_j;
        }
    }
    return false;
}

int max(int a, int b){
    if(a > b){
        return a;
    }
    return b;
}

int min(int a, int b){
    if(a < b){
        return a;
    }
    return b;
}



// Handles the multiplication in the even subalgebra of the Steenrod algebra P.
// When p = 2, this is isomorphic to the whole Steenrod algebra so this method does everything.
void MilnorAlgebra__multiplyEven(MilnorAlgebraInternal *algebra, uint *result, MilnorBasisElement r_elt, MilnorBasisElement s_elt){   
    uint p = algebra->public_algebra.algebra.p;
    P_part r, s;
    r.degree = r_elt.p_degree;
    r.length = r_elt.p_length;
    r.p_part = r_elt.p_part;
    s.degree = s_elt.p_degree;
    s.length = s_elt.p_length;
    s.p_part = s_elt.p_part;
    int output_degree = (r_elt.p_degree + s_elt.p_degree);
    uint output_dimension = algebra->basis_table[output_degree].length;
    memset(result, 0, output_dimension *sizeof(uint));
    uint rows = r.length + 1;
    uint cols = s.length + 1;
    uint number_of_diagonals = r.length + s.length;

    uint M[MAX_XI_TAU][MAX_XI_TAU];
    milnor_matrix_initialize(M, r, s);

    do {
        uint coeff = 1;
        uint diagonal_sums[MAX_XI_TAU];
        uint max_nonzero_diagonal = 0;
        for(uint diagonal_index = 1; diagonal_index <= number_of_diagonals; diagonal_index++){
            // We're going to iterate along the diagonal and copy the diagonal into nth_diagonal
            // and put the sum into nth_diagonal_sum.
            uint i_min = max(0, diagonal_index - cols + 1);
            uint i_max = min(1 + diagonal_index, rows);
            uint diagonal_length = i_max - i_min;
            uint diagonal[diagonal_length];
            uint diagonal_sum = 0;
            uint index = 0;
            // Iterate along diagonal
            for(uint i = i_min; i < i_max;  i++){
                diagonal[index] = M[i][diagonal_index - i];
                diagonal_sum += diagonal[index];
                index ++;
            }
            diagonal_sums[diagonal_index - 1] = diagonal_sum;
            if(diagonal_sum == 0) {
                // If we dropped the continue, the only thing underneath that affects correctness is max_nonzero_diagonal = diagonal.
                // The multinomial will be 1 so the coefficient will not change and everything else will do nothing.
                continue;
            }
            max_nonzero_diagonal = diagonal_index;
            coeff *= Multinomial(p, diagonal_length, diagonal);
            coeff = coeff % p;
            if(coeff == 0){
                break;
            }
        }
        if(coeff != 0) {
            MilnorBasisElement m = (MilnorBasisElement){0, 0, output_degree, max_nonzero_diagonal, diagonal_sums};
            uint idx = MilnorAlgebra_basisElement_toIndex((MilnorAlgebra*)algebra, m);
            result[idx] = (result[idx] + coeff) % p;
        }
    } while(milnor_matrix_step(p, M, r, s));
}



// Reduce m1 * f = (Q_e0 Q_e1 ... P(r1, r2, ...)) * (Q_f0 Q_f1 ...) into the form Sum of Q's * P's
// Result is represented as dictionary of pairs of tuples.
void MilnorAlgebra__multiplyQpart(MilnorAlgebraInternal *algebra, uint *output, MilnorBasisElement m1, uint f) {
    uint p = algebra->public_algebra.algebra.p;
    uint q;
    q = 2*p-2;
    uint *tau_degrees = getTauDegrees(p);
    uint *xi_degrees = getXiDegrees(p);
    int p_degree = m1.p_degree;
    int q_degree = m1.q_degree;
    uint idx = MilnorAlgebra_basisElement_toIndex((MilnorAlgebra*)algebra, m1);

    int result_degree = p_degree + q_degree;
    uint result_dimension = algebra->basis_table[result_degree].length;
    
    // We move one Q over at a time and generate a series of intermediate results.
    // Thus, we need a result and an old_result.
    // result and old_result have these two "memory" arrays as backing, we swap them
    // after commuting each Q.
    uint memory1[MAX_DIMENSION], memory2[MAX_DIMENSION];
    uint old_result_dimension, old_result_degree;

    uint *result = memory1;
    uint *old_result = memory2;

    memset(result, 0, result_dimension * sizeof(uint));
    result[idx] = 1;
    uint p_to_the_k = 0;
    for(uint k = 0; ( f & ~((1 << k) - 1) ) != 0; k++){
        p_to_the_k *= p;
        if(k == 0){
            p_to_the_k = 1;
        }
        if((f & (1 << k)) == 0){
            continue;
        }
        q_degree += tau_degrees[k];
        // swap result and old result
        uint *swap_temp = old_result;
        old_result = result;
        old_result_degree = result_degree;
        old_result_dimension = result_dimension;
        result = swap_temp;
        result_degree = p_degree + q_degree;
        result_dimension = algebra->basis_table[result_degree].length;
        // Zero result
        memset(result, 0, result_dimension * sizeof(uint));
        for(uint idx = 0; idx < old_result_dimension; idx++){
            if(old_result[idx] == 0){
                continue;
            }
            MilnorBasisElement mono = MilnorAlgebra_basisElement_fromIndex((MilnorAlgebra*)algebra, old_result_degree, idx);
            for(uint i = 0; i < mono.p_length + 1; i++){
                // There is already a Q of this index on the other side, so moving another one over will give 0.
                if((mono.q_part & (1 << (k+i))) != 0){
                    continue;
                }
                // Make sure mono.p_part[i - 1] is large enough to deduct p^k from it
                if(i > 0 && mono.p_part[i - 1] < p_to_the_k ){
                    continue;
                }

                uint q_mono = mono.q_part;
                int q_degree = mono.q_degree;
                uint *p_mono = mono.p_part;
                uint p_mono_length = mono.p_length;
                int p_degree = mono.p_degree;

                uint new_p_mono[MAX_XI_TAU];
                if(i > 0){
                    p_mono[i - 1] -= p_to_the_k;
                    uint new_p_mono_length = p_mono_length;
                    for( ; new_p_mono_length > 0 && p_mono[new_p_mono_length - 1] == 0; new_p_mono_length -- ){}
                    memcpy(new_p_mono, p_mono, new_p_mono_length * sizeof(uint));
                    p_mono[i - 1] += p_to_the_k;
                    p_degree -= q * p_to_the_k * xi_degrees[i-1];
                    p_mono = new_p_mono;
                    p_mono_length = new_p_mono_length;
                }
                // Add the new Q.
                q_mono |= 1 << (k + i);
                q_degree += tau_degrees[k+i];

                uint larger_Qs = 0;
                uint v = q_mono >> (k + i + 1);
                for( ; v != 0; v >>= 1){
                    larger_Qs += v & 1;
                }
                uint coeff = MinusOneToTheN(p, larger_Qs) * old_result[idx];
                MilnorBasisElement m = (MilnorBasisElement){q_degree, q_mono, p_degree, p_mono_length, p_mono};
                uint out_idx = MilnorAlgebra_basisElement_toIndex((MilnorAlgebra*)algebra, m);
                result[out_idx] += coeff;
                result[out_idx] = modPLookup(p, result[out_idx]);
            }
        }
    }
    memcpy(output, result, result_dimension*sizeof(uint));
}


// Product of Milnor basis elements defined by m1 and m2 at the prime p.
//
// INPUT:
//
// - m1 - pair of tuples (e,r), where e is an increasing tuple of
//   non-negative integers and r is a tuple of non-negative integers
// - m2 - pair of tuples (f,s), same format as m1
// - p -- odd prime number
//
// OUTPUT:
//
// Dictionary of terms of the form (tuple: coeff), where 'tuple' is
// a pair of tuples, as for r and s, and 'coeff' is an integer mod p.
//
// This computes the product of the Milnor basis elements
// $Q_{e_1} Q_{e_2} ... P(r_1, r_2, ...)$ and
// $Q_{f_1} Q_{f_2} ... P(s_1, s_2, ...)$.
void MilnorAlgebra__multiplyFull(MilnorAlgebraInternal *algebra, uint *result, MilnorBasisElement m1, MilnorBasisElement m2) {
    uint p = algebra->public_algebra.algebra.p;
    MilnorBasisElement s = (MilnorBasisElement) {0, 0, m2.p_degree, m2.p_length, m2.p_part};
    int m1_times_f_degree = m1.q_degree + m1.p_degree + m2.q_degree;
    uint m1_times_f_dimension = algebra->basis_table[m1_times_f_degree].length;
    uint m1_times_f[m1_times_f_dimension];
    memset(m1_times_f, 0, m1_times_f_dimension * sizeof(uint));
    // uint out_degree = m1_times_f_degree + m2.p_degree;
    // uint out_dimension = algebra->basis_table[out_degree].length;
    MilnorAlgebra__multiplyQpart(algebra, m1_times_f, m1, m2.q_part);

    // Now for the Milnor matrices.  For each entry '(e,r): coeff' in answer,
    // multiply r with s.  Record coefficient for matrix and multiply by coeff.
    // Store in 'result'.
    if (m2.p_length == 0) {
        memcpy(result, m1_times_f, m1_times_f_dimension * sizeof(uint));
        return;
    }
    for (uint i = 0; i < m1_times_f_dimension; i++) {
        if (m1_times_f[i] == 0) {
            continue;
        }
        uint coeff = m1_times_f[i];
        MilnorBasisElement er_mono = MilnorAlgebra_basisElement_fromIndex((MilnorAlgebra *) algebra, m1_times_f_degree, i);
        MilnorBasisElement r = (MilnorBasisElement) {0, 0, er_mono.p_degree, er_mono.p_length, er_mono.p_part};
        int prod_degree = r.p_degree + s.p_degree;
        uint prod_dimension = algebra->basis_table[prod_degree].length;
        uint prod[prod_dimension];
        MilnorAlgebra__multiplyEven(algebra, prod, r, s);
        // Vector *prod_test = Vector_construct(p, prod_dimension, 0);
        // Vector_pack(prod_test, prod);
        // MilnorElement_print("prod: %s\n", algebra, prod_degree, prod_test);
        // Vector_free(prod_test);
        for (uint j = 0; j < prod_dimension; j++) {
            // Watch out: many of the entries in prod correspond to basis vectors that
            // have nonzero Q part. The result of MilnorAlgebra__multiplyEven is guaranteed to have
            // 0's there. We could definitely make this slicker.
            if(prod[j] == 0){
                continue;
            }
            MilnorBasisElement m = MilnorAlgebra_basisElement_fromIndex((MilnorAlgebra *) algebra, prod_degree, j);
            m.q_degree = er_mono.q_degree;
            m.q_part = er_mono.q_part;
            uint out_idx = MilnorAlgebra_basisElement_toIndex((MilnorAlgebra *) algebra, m);
            result[out_idx] += coeff * prod[j];
            result[out_idx] = result[out_idx] % p;
        }
        // Vector_pack(result_test, result);
        // MilnorElement_print("result: %s\n\n", algebra, out_degree, result_test);
    }
    // Vector_free(result_test);
}


// Multiply the basis elements corresponding to r_index and s_index in the Milnor algebra determined by algebra.
// Note that since profile functions determine subalgebras, the product
// doesn't need to care about the profile function at all.
void MilnorAlgebra_multiply(Algebra *this, Vector *result, uint coeff, int r_degree, uint r_index, int s_degree, uint s_index, int excess __attribute__((unused)))  {
    MilnorAlgebraInternal *algebra = (MilnorAlgebraInternal*) this;
    MilnorBasisElement r, s;
    assert(r_index < MilnorAlgebra_getDimension(this, r_degree, excess));
    assert(s_index < MilnorAlgebra_getDimension(this, s_degree, excess));

    r = MilnorAlgebra_basisElement_fromIndex((MilnorAlgebra*)algebra, r_degree, r_index);
    s = MilnorAlgebra_basisElement_fromIndex((MilnorAlgebra*)algebra, s_degree, s_index);
    int output_degree = r_degree + s_degree;
    uint output_dimension = algebra->basis_table[output_degree].length;
    assert(output_dimension == result->dimension);
    assert(output_degree < algebra->public_algebra.algebra.max_degree);
    uint product_array[output_dimension];
    memset(product_array, 0, output_dimension * sizeof(uint));

    if(algebra->public_algebra.generic){
        MilnorAlgebra__multiplyFull(algebra, product_array, r, s);
    } else {
        MilnorAlgebra__multiplyEven(algebra, product_array, r, s);
    }
    Vector_addArray(result, product_array, coeff);
}
/**/




/**
int main() {
    char buffer[10000];

//    MilnorBasisElement m, n;
//    m.q_degree = 1;
//    m.q_part = 1;
//    m.p_degree = 8;
//    m.p_length = 1;
//    m.p_part = (uint[1]){1};
//
//    n.q_degree = 1;
//    n.q_part = 1;
//    n.p_degree = 8;
//    n.p_length = 1;
//    n.p_part = (uint[1]){1};
    uint p = 2;
    MilnorAlgebra *algebra = MilnorAlgebra_construct(p, p!=2, NULL);

//    generateMilnorBasisQpartTable(algebra, 20);
//    MilnorAlgebra__freeQPartTable(algebra);
//    generateMilnorBasisQpartTable(algebra, 20);
//    generateMilnorBasisQpartTable(algebra, 76);
//
//    Q_part_list *Q_table = algebra->Q_table;
//    for(uint int residue = 0; residue < 2*3 - 2; ++residue) {
//        Q_part_list v = Q_table[residue];
//        printf("residue: %ld\n", residue);
//        for(int i = 0; i < v.length ; i++){
//            Q_part q_part = v.list[i];
//            printf(" %ld , %ld\n", q_part.bit_string, q_part.degree);
//        }
//        printf("\n");
//    }
//    MilnorAlgebra__freeQPartTable(algebra);

//    MilnorAlgebra__generatePpartTable(algebra, 10);
//    P_part_list *P_table = algebra->P_table;
//    for(int n = 0; n < 10; n++){
//        printf("deg: %d \n", n);
//        for(int i = 0; i < P_table[n].length; i++){
//            P_part v = P_table[n].list[i];
//            printf("   [");
//            for(int j = 0; j < v.length; j++){
//                printf("%ld, ", v.p_part[j]);
//            }
//            printf("]\n");
//        }
//        printf("\n");
//    }
//    MilnorAlgebra__freePPartTable(algebra);


    //GenerateMilnorBasis(algebra, 76);
    //MilnorAlgebra__freeQPartTable(algebra);
    MilnorAlgebra_generateBasis((Algebra*)algebra, 101);

//    for(int degree = 0; degree < algebra->max_degree; degree++){
//        MilnorBasisElement_list basis_list = GetMilnorAlgebraBasis(algebra, degree);
//        printf("degree: %d \n", degree);
//        for(int i = 0; i < basis_list.length; i++){
//            MilnorBasisElement v = basis_list.list[i];
//            milnor_basis_element_to_string(buffer, &v);
//            printf("  %d: '%s'\n", i , buffer);
//        }
//        printf("\n");
//    }

    MilnorBasisElement A, B;
    //(A5.Q(1) * A5.P(2), A5.P(1, 1))
    A = MilnorAlgebra_basisElement_fromString(algebra, "P(10, 4)");
    B = MilnorAlgebra_basisElement_fromString(algebra, "P(10, 13, 2, 1)");
    uint A_idx = MilnorAlgebra_basisElement_toIndex(algebra, A);
    uint B_idx = MilnorAlgebra_basisElement_toIndex(algebra, B);
    uint A_deg = A.p_degree + A.q_degree;
    uint B_deg = B.p_degree + B.q_degree;

    uint output_degree = A_deg + B_deg;
    uint output_dimension = MilnorAlgebra_getDimension((Algebra*)algebra, output_degree, 0);
    Vector *result = Vector_construct(algebra->algebra.p, output_dimension, 0);
    MilnorAlgebra_multiply((Algebra*)algebra, result, 1, A_deg, A_idx, B_deg, B_idx, 0);
    char *str1 = buffer;
    int len1 = MilnorAlgebra_element_toString(str1, algebra, output_degree, result);
    char *str2 = str1 + len1 + 1;
    int len2 = MilnorAlgebra_basisElement_toString(str2, algebra, &A);
    char *str3 = str2 + len2 + 1;
    MilnorAlgebra_basisElement_toString(str3, algebra,  &B);
    printf("%s * %s = %s\n", str2, str3, str1);
    free(A.p_part);
    free(B.p_part);
    Vector_free(result);
    MilnorAlgebra_free(algebra);
}
//*/