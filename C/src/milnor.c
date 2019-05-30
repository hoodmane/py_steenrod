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

#define MAX_DIMENSION 10000


void printMilnorBasisElement(MilnorBasisElement *b){
    char buffer[1000];
    milnor_basis_element_to_string(buffer, b);
}


MilnorAlgebra *constructMilnorAlgebra(uint p, bool generic, Profile *profile){
    MilnorAlgebraInternal *algebra = malloc(sizeof(MilnorAlgebraInternal));
    initializePrime(p);
    algebra->public_algebra.p = p;
    algebra->public_algebra.algebra.p = p;
    algebra->public_algebra.algebra.vectorInterface = (p == 2) ? Vector2Interface : VectorGenericInterface;
    algebra->public_algebra.generic = generic;
    if(profile == NULL){
        algebra->public_algebra.profile.truncated = false;
        algebra->public_algebra.profile.p_part_length = 0;
        algebra->public_algebra.profile.q_part = -1;
    } else {
        algebra->public_algebra.profile = *profile;
    }

    algebra->public_algebra.algebra.computeBasis = GenerateMilnorBasis;
    algebra->public_algebra.algebra.getDimension = GetMilnorAlgebraDimension;
    algebra->public_algebra.algebra.multiplyBasisElements = MilnorProduct;

    algebra->P_table = NULL;
    algebra->P_table_by_P_length = NULL;
    algebra->P_table_max_degree = 0;

    algebra->Q_table = NULL;
    algebra->Q_table_max_tau = 0;

    algebra->basis_table = NULL;
    algebra->public_algebra.max_degree = -1;
    algebra->basis_element_to_index_map = NULL;

    return (MilnorAlgebra*)algebra;
}

void freeMilnorAlgebra(MilnorAlgebra *algebra){
    freeMilnorBasis(algebra);
    free(algebra);
}

uint GetMilnorAlgebraDimension(Algebra *public_algebra, uint degree) {
    MilnorAlgebraInternal *algebra = (MilnorAlgebraInternal *) public_algebra;
    return algebra->basis_table[degree].length;
}

MilnorBasisElement_list GetMilnorAlgebraBasis(MilnorAlgebra *public_algebra, uint degree){
    MilnorAlgebraInternal *algebra = (MilnorAlgebraInternal *) public_algebra;
    return algebra->basis_table[degree];
}



// TODO: Break this up into several methods.
// Generate the PTable through degree new_max_degree.
// This is only needed in order to help generate the basis table, so it can in principle be deallocated immediately afterwards.
// However, we keep it around so that if we need to extend the table we don't have to redo our work.
void generateMilnorBasisPpartTable(MilnorAlgebraInternal *algebra, uint new_max_degree) {
    uint p = algebra->public_algebra.p;
    uint *xi_degrees = getXiDegrees(p);

    uint profile_list[MAX_XI_TAU];
    for(uint idx = 0;  idx < MAX_XI_TAU; idx ++){
        profile_list[idx] = getProfileExponent(algebra->public_algebra.profile, p, idx) - 1;
    }

    uint old_max_degree = algebra -> P_table_max_degree;
    if(new_max_degree < old_max_degree){
        return;
    }
    algebra -> P_table_max_degree = new_max_degree;
    algebra -> P_table = realloc(algebra -> P_table , (new_max_degree + 1) * sizeof(P_part_list));
    algebra -> P_table_by_P_length = realloc(algebra -> P_table_by_P_length, (new_max_degree + 1) * sizeof(P_part_list[MAX_XI_TAU]));
    P_part_list *degree_table = algebra->P_table;
    P_part_list **tau_table = algebra->P_table_by_P_length;

    // If old_max_degree is 0 we need to set up the base case for our recursion.
    // In degree 0 there's a single partition of length 0.
    if(old_max_degree == 0){
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
    for(uint current_degree = old_max_degree + 1; current_degree <= new_max_degree; current_degree++){
        degree_list_length = 0;
        tau_table[current_degree] = (P_part_list*) calloc(MAX_XI_TAU, sizeof(P_part_list));

        // Now populate the tables.
        // Loop over the highest degree xi_i in the new monomial
        for(uint xi = 0; xi < MAX_XI_TAU; xi ++){
            uint xi_degree = xi_degrees[xi];
            if(xi_degree > current_degree){
                break;
            }
            uint remaining_degree = current_degree - xi_degree;
            tau_list_length = 0;
            // Iterate over the highest xi_i in the old monomial
            for(uint xi_old = 0; xi_old <= xi; xi_old ++) {
                P_part_list remaining_degree_list = tau_table[remaining_degree][xi_old];
                // Iterate over the monomials of degree n - |xi| with highest xi xi_old.
                for(uint i = 0; i < remaining_degree_list.length; i++) {
                    P_part old_p_part = remaining_degree_list.list[i];
                    if(xi_old == xi && old_p_part.p_part[xi] == profile_list[xi]){
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
void freeMilnorBasisPpartTable(MilnorAlgebraInternal *algebra){
    if(algebra->P_table == NULL){
        return;
    }
    uint max_degree = algebra -> P_table_max_degree;
    P_part_list *degree_table = algebra->P_table;
    P_part_list **tau_table = algebra->P_table_by_P_length;
    for(uint current_degree = 0; current_degree <= max_degree; current_degree++){
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
void generateMilnorBasisQpartTable(MilnorAlgebraInternal *algebra, uint new_max_degree){
    uint p = algebra->public_algebra.p;
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
    Q_part residue_bin_buffer[q][MAX_DIMENSION];

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
    uint residue = old_max_tau;
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
    for (int residue = 0; residue < q; residue++) {
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
void freeMilnorBasisQPartTable(MilnorAlgebraInternal *algebra){
    uint p = algebra->public_algebra.p;
    uint q = 2*p - 2;
    Q_part_list *table = algebra->Q_table;
    for (int residue = 0; residue < q; residue++) {
        free(table[residue].list);
    }
    free(algebra->Q_table);
    algebra->Q_table_max_tau = 0;
    algebra->Q_table = NULL;
}

bool GenerateMilnorBasis(Algebra *public_algebra, uint max_degree) {
    MilnorAlgebraInternal *algebra = (MilnorAlgebraInternal*) public_algebra;
    uint p = algebra->public_algebra.p;
    initializePrime(p);
    uint old_max_degree = algebra->public_algebra.max_degree;
    algebra->basis_table = (MilnorBasisElement_list *) realloc(
            algebra->basis_table,
            (max_degree + 1) * sizeof(MilnorBasisElement_list)
    );
    algebra->public_algebra.max_degree = max_degree;
    algebra->basis_element_to_index_map = (khash_t(monomial_index_map) **)realloc(
            algebra->basis_element_to_index_map,
            (max_degree + 1) * sizeof(khash_t(monomial_index_map) *)
        );
    khash_t(monomial_index_map) **name_table = algebra->basis_element_to_index_map;
    for(uint k = old_max_degree + 1; k <= max_degree; k++){
        name_table[k] = kh_init(monomial_index_map);
    }
    if(algebra->public_algebra.generic){
        GenerateMilnorBasisGeneric(algebra, old_max_degree, max_degree);
    } else {
        GenerateMilnorBasis2(algebra, old_max_degree, max_degree);
    }
    return true;
}

void freeMilnorBasis(MilnorAlgebra *public_algebra){
    MilnorAlgebraInternal *algebra = (MilnorAlgebraInternal*) public_algebra;
    MilnorBasisElement_list *table = algebra->basis_table;
    khash_t(monomial_index_map) **name_table = algebra->basis_element_to_index_map;
    for(uint degree = 0; degree <= algebra->public_algebra.max_degree; degree++) {
        free(table[degree].list);
        khint_t bin;
        for (bin = 0; bin < kh_end(name_table[degree]); ++bin) {
            if (kh_exist(name_table[degree], bin)){
                free((char *) kh_key(name_table[degree], bin));
            }
        }
    }

    freeMilnorBasisPpartTable(algebra);
    if(algebra->public_algebra.generic){
        freeMilnorBasisQPartTable(algebra);
    }

    for(uint k = 0; k <= algebra->public_algebra.max_degree; k++){
        kh_destroy(monomial_index_map, algebra->basis_element_to_index_map[k]);
    }
    free(algebra->basis_element_to_index_map);
    free(algebra->basis_table);
}

void GenerateMilnorBasis2(MilnorAlgebraInternal *algebra, uint old_max_degree, uint new_max_degree){
    generateMilnorBasisPpartTable(algebra, new_max_degree);

    MilnorBasisElement_list *table = algebra->basis_table;
    khash_t(monomial_index_map) **name_table = algebra->basis_element_to_index_map;

    MilnorBasisElement_list current_degree_list;
    for(uint degree = old_max_degree + 1; degree <= new_max_degree; degree++){
        P_part_list p_parts = algebra->P_table[degree];
        current_degree_list.length = 0;
        current_degree_list.list = (MilnorBasisElement*) malloc(p_parts.length * sizeof(MilnorBasisElement));
        for(uint i = 0; i < p_parts.length; i++){
            P_part x = p_parts.list[i];
            MilnorBasisElement m = (MilnorBasisElement){0, 0, degree, x.length, x.p_part};
            char key[200];
            milnor_basis_element_to_key(key, &m);
            int absent;
            khint_t bin = kh_put(monomial_index_map, name_table[degree], key, &absent);
            if(absent){
                kh_key(name_table[degree], bin) = strdup(key);
            }
            kh_val(name_table[degree], bin) = current_degree_list.length;
            current_degree_list.list[current_degree_list.length] = m;
            current_degree_list.length ++;
        }
        table[degree] = current_degree_list;
    }
}

// Get the basis in degree n for the generic steenrod algebra at the prime p.
// We just put together the "even part" of the basis and the "Q part".
void GenerateMilnorBasisGeneric(MilnorAlgebraInternal *algebra, uint old_max_degree, uint new_max_degree){
    uint p = algebra->public_algebra.p;
    uint q = 2 * (p - 1);
    generateMilnorBasisPpartTable(algebra, new_max_degree / q);
    generateMilnorBasisQpartTable(algebra, new_max_degree);
    MilnorBasisElement_list *table = algebra->basis_table;
    khash_t(monomial_index_map) **name_table = algebra->basis_element_to_index_map;

    for(uint degree = old_max_degree + 1; degree <= new_max_degree; degree++){
        // p_deg records the desired degree of the P part of the basis element.
        // Since p-parts are always divisible by 2p-2, we divide by this first.
        // pow(p, -1) returns 1, so min_q_deg is 0 if q divides n evenly.
        Q_part_list q_list = algebra->Q_table[degree % q];
        uint degree_list_length = 0;
        MilnorBasisElement degree_list_buffer[MAX_DIMENSION];
        for(uint i = 0; i < q_list.length; i++) {
            Q_part q_part = q_list.list[i];
            uint q_deg = q_part.degree;
            if(q_deg > degree){
                break;
            }
            uint q_bit_string = q_part.bit_string;
            uint p_deg = (degree - q_deg);
            P_part_list p_list = algebra->P_table[p_deg / q];
            for(uint j = 0; j < p_list.length; j++) {
                P_part p_part = p_list.list[j];
                MilnorBasisElement m = (MilnorBasisElement){q_part.degree, q_bit_string, p_deg, p_part.length, p_part.p_part};
                char key[200];
                milnor_basis_element_to_key(key, &m);
                int absent;
                khint_t bin = kh_put(monomial_index_map, name_table[degree], key, &absent);
                if(absent){
                    kh_key(name_table[degree], bin) = strdup(key);
                }
                kh_value(name_table[degree], bin) = degree_list_length;

                degree_list_buffer[degree_list_length] = m;
                degree_list_length ++;
            }
        }
        table[degree].length = degree_list_length;
        table[degree].list = malloc(degree_list_length * sizeof(MilnorBasisElement));
        memcpy(table[degree].list, degree_list_buffer, degree_list_length * sizeof(MilnorBasisElement));
    }
}



MilnorBasisElement GetMilnorBasisElementFromIndex(MilnorAlgebra *public_algebra, uint degree, uint idx) {
    MilnorAlgebraInternal *algebra = (MilnorAlgebraInternal*) public_algebra;
    return algebra->basis_table[degree].list[idx];
}

uint GetIndexFromMilnorBasisElement(MilnorAlgebra *public_algebra, MilnorBasisElement b){
    MilnorAlgebraInternal *algebra = (MilnorAlgebraInternal*) public_algebra;
    kh_monomial_index_map_t *map = algebra->basis_element_to_index_map[b.q_degree + b.p_degree];
    char key[200];
    milnor_basis_element_to_key(key, &b);
    khint_t bin = kh_get(monomial_index_map, map, key);
    if(bin == kh_end(map)){
        milnor_basis_element_to_string(key, &b);
        printf("Uh-oh, not here. degree: %d, elt: '%s'\n", b.q_degree + b.p_degree, key);
        return -1;
    }
    return kh_val(map, bin);
}

// Initializes an len(r)+1 by len(s)+1 matrix
// Puts r along the first column and s along the first row and zeroes everywhere else.
void initialize_milnor_matrix(uint M[MAX_XI_TAU][MAX_XI_TAU], P_part r, P_part s) {
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
void step_milnor_matrix_update_matrix(uint M[MAX_XI_TAU][MAX_XI_TAU], P_part r, P_part s, int i, int j, int x)  {
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
bool step_milnor_matrix(uint  p, uint M[MAX_XI_TAU][MAX_XI_TAU], P_part r, P_part s) {
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
                    step_milnor_matrix_update_matrix(M, r, s, i, j, total - p_to_the_j);
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

//Remove trailing zeroes from the list l.
//func remove_trailing_zeroes(l []int) []int {
//for i := len(l) - 1; i >= 0; i-- {
//if l[i] != 0 {
//return l[:i+1]
//}
//}
//return l[:0]
//}

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
void MilnorProductEven(MilnorAlgebraInternal *algebra, uint *result, MilnorBasisElement r_elt, MilnorBasisElement s_elt){   
    uint p = algebra->public_algebra.p;
    P_part r, s;
    r.degree = r_elt.p_degree;
    r.length = r_elt.p_length;
    r.p_part = r_elt.p_part;
    s.degree = s_elt.p_degree;
    s.length = s_elt.p_length;
    s.p_part = s_elt.p_part;
    uint output_degree = (r_elt.p_degree + s_elt.p_degree);
    uint output_dimension = algebra->basis_table[output_degree].length;
    memset(result, 0, output_dimension *sizeof(uint));
    uint rows = r.length + 1;
    uint cols = s.length + 1;
    uint number_of_diagonals = r.length + s.length;

    uint M[MAX_XI_TAU][MAX_XI_TAU];
    initialize_milnor_matrix(M, r, s);

    do {
        uint coeff = 1;
        uint diagonal_sums[MAX_XI_TAU];
        uint max_nonzero_diagonal = 0;
        for(int diagonal_index = 1; diagonal_index <= number_of_diagonals; diagonal_index++){
            // We're going to iterate along the diagonal and copy the diagonal into nth_diagonal
            // and put the sum into nth_diagonal_sum.
            int i_min = max(0, diagonal_index - cols + 1);
            int i_max = min(1 + diagonal_index, rows);
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
            uint idx = GetIndexFromMilnorBasisElement((MilnorAlgebra*)algebra, m);
            result[idx] = coeff;
        }
    } while(step_milnor_matrix(p, M, r, s));
}



// Reduce m1 * f = (Q_e0 Q_e1 ... P(r1, r2, ...)) * (Q_f0 Q_f1 ...) into the form Sum of Q's * P's
// Result is represented as dictionary of pairs of tuples.
void MilnorProductFullQpart(MilnorAlgebraInternal *algebra, uint *output, MilnorBasisElement m1, uint f) {
    uint p = algebra->public_algebra.p;
    uint q;
    q = 2*p-2;
    uint *tau_degrees = getTauDegrees(p);
    uint *xi_degrees = getXiDegrees(p);
    uint p_degree = m1.p_degree;
    uint q_degree = m1.q_degree;
    uint idx = GetIndexFromMilnorBasisElement((MilnorAlgebra*)algebra, m1);

    uint result_degree = p_degree + q_degree;
    uint result_dimension = algebra->basis_table[result_degree].length;
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
        uint *swap_temp = old_result;
        old_result = result;
        old_result_degree = result_degree;
        old_result_dimension = result_dimension;
        result = swap_temp;
        result_degree = p_degree + q_degree;
        result_dimension = algebra->basis_table[result_degree].length;
        memset(result, 0, result_dimension * sizeof(uint));
        for(uint idx = 0; idx < old_result_dimension; idx++){
            if(old_result[idx] == 0){
                continue;
            }
            MilnorBasisElement mono = GetMilnorBasisElementFromIndex((MilnorAlgebra*)algebra, old_result_degree, idx);
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
                uint q_degree = mono.q_degree;
                uint *p_mono = mono.p_part;
                uint p_mono_length = mono.p_length;
                uint p_degree = mono.p_degree;

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
                uint out_idx = GetIndexFromMilnorBasisElement((MilnorAlgebra*)algebra, m);
                result[out_idx] += coeff;
                result[out_idx] = modPLookup(p, result[out_idx]);
                char buffer[1000];
                array_to_string(buffer, result, result_dimension);
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
void MilnorProductFull(MilnorAlgebraInternal *algebra, uint *result, MilnorBasisElement m1, MilnorBasisElement m2) {
    uint p = algebra->public_algebra.p;
    MilnorBasisElement s = (MilnorBasisElement) {0, 0, m2.p_degree, m2.p_length, m2.p_part};
    uint m1_times_f_degree = m1.q_degree + m1.p_degree + m2.q_degree;
    uint m1_times_f_dimension = algebra->basis_table[m1_times_f_degree].length;
    uint m1_times_f[m1_times_f_dimension];
    memset(m1_times_f, 0, m1_times_f_dimension * sizeof(uint));
    MilnorProductFullQpart(algebra, m1_times_f, m1, m2.q_part);

    printMilnorBasisElement(&m2);
    // char buffer[200];
    // milnor_element_to_string(buffer, (MilnorAlgebra *) algebra, &m1_times_f);

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
        MilnorBasisElement er_mono = GetMilnorBasisElementFromIndex((MilnorAlgebra *) algebra, m1_times_f_degree, i);
        MilnorBasisElement r = (MilnorBasisElement) {0, 0, er_mono.p_degree, er_mono.p_length, er_mono.p_part};
        uint prod_degree = r.p_degree + s.p_degree;
        uint prod_dimension = algebra->basis_table[prod_degree].length;
        uint prod[prod_dimension];
        MilnorProductEven(algebra, prod, r, s);
        for (uint j = 0; j < prod_dimension; j++) {
            uint c = prod[j];
            MilnorBasisElement m = GetMilnorBasisElementFromIndex((MilnorAlgebra *) algebra, prod_degree, j);
            m.q_degree = er_mono.q_degree;
            m.q_part = er_mono.q_part;
            uint out_idx = GetIndexFromMilnorBasisElement((MilnorAlgebra *) algebra, m);
            result[out_idx] = (coeff * c) % p;
        }
    }
}


// Multiply the basis elements corresponding to r_index and s_index in the Milnor algebra determined by algebra.
// Note that since profile functions determine subalgebras, the product
// doesn't need to care about the profile function at all.
void MilnorProduct(Algebra *public_algebra, Vector *result, uint coeff, uint r_degree, uint r_index, uint s_degree, uint s_index)  {
    MilnorAlgebraInternal *algebra = (MilnorAlgebraInternal*) public_algebra;
    MilnorBasisElement r, s;
    r = GetMilnorBasisElementFromIndex((MilnorAlgebra*)algebra, r_degree, r_index);
    s = GetMilnorBasisElementFromIndex((MilnorAlgebra*)algebra, s_degree, s_index);
    // if(r_degree + s_degree != result->degree){
    //     printf("Result has degree %ld but should have degree %ld\n", result->degree, r.q_degree + r.p_degree + s.q_degree + s.p_degree);
    //     return;
    // }
    uint output_degree = r_degree + s_degree;
    uint output_dimension = algebra->basis_table[output_degree].length;
    assert(output_dimension == result->dimension);
    uint product_array[output_dimension];
    memset(product_array, 0, output_dimension * sizeof(uint));
    if(algebra->public_algebra.generic){
        MilnorProductFull(algebra, product_array, r, s);
    } else {
        MilnorProductEven(algebra, product_array, r, s);
    }
    public_algebra->vectorInterface.addArray(result, product_array, coeff);
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

    MilnorAlgebra *algebra = constructMilnorAlgebra(3, true, NULL);
    VectorInterface vectorInteterface = algebra->algebra.vectorInterface;

//    generateMilnorBasisQpartTable(algebra, 20);
//    freeMilnorBasisQPartTable(algebra);
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
//    freeMilnorBasisQPartTable(algebra);

//    generateMilnorBasisPpartTable(algebra, 10);
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
//    freeMilnorBasisPpartTable(algebra);


    //GenerateMilnorBasis(algebra, 76);
    //freeMilnorBasisQPartTable(algebra);
    GenerateMilnorBasis((Algebra*)algebra, 100);

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

    MilnorBasisElement A,B;
    //(A5.Q(1) * A5.P(2), A5.P(1, 1))
    A = milnor_basis_element_from_string(algebra, " Q(0)");
    B = milnor_basis_element_from_string(algebra, "Q(0, 1) * P(3)");
    uint A_idx = GetIndexFromMilnorBasisElement(algebra, A);
    uint B_idx = GetIndexFromMilnorBasisElement(algebra, B);
    uint A_deg = A.p_degree + A.q_degree;
    uint B_deg = B.p_degree + B.q_degree;

    uint output_degree = A_deg + B_deg;
    uint output_dimension = GetMilnorAlgebraDimension((Algebra*)algebra, output_degree);
    Vector *result = vectorInteterface.construct(algebra->p, output_dimension);
    MilnorProduct((Algebra*)algebra, result, 1, A_deg, A_idx, B_deg, B_idx);
    string str1 = buffer;
    int len1 = milnor_element_to_string(str1, algebra, output_degree, result);
    string str2 = str1 + len1 + 1;
    int len2 = milnor_basis_element_to_string(str2, &A);
    string str3 = str2 + len2 + 1;
    milnor_basis_element_to_string(str3, &B);
    printf("%s * %s = %s\n", str2, str3, str1);
    free(A.p_part);
    free(B.p_part);
    freeVector(result);
    freeMilnorAlgebra(algebra);
}
**/