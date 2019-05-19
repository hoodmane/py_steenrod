//
// Created by Hood on 5/8/2019.
//

#include "milnor_datatypes.h"
#include "milnor.h"
#include "combinatorics.h"

#include "khash.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define MAX_ALGEBRA_DIMENSION 500

int main() {
    char buffer[100];

    initializePrime(3);
    MilnorAlgebra * algebra;
    MilnorAlgebra a;
    algebra = & a;
    algebra->p = 3;
    algebra->generic = true;
    constructMilnorAlgebra(algebra);
//    algebra->p = 3;
//    algebra->generic = true;
//    algebra->profile.restricted = true;
//    algebra->profile.q_part = (~0);
//    algebra->profile.p_part_length = 1;
//    algebra->profile.p_part = (unsigned long*)malloc(sizeof(unsigned long));
//    algebra->profile.p_part[0] = 5;
//    algebra->profile.truncated = false;
//    generateMilnorBasisQpartTable(algebra, 100);
//    Q_part_list* Q_table = algebra->Q_table;
//    for(unsigned long int residue = 0; residue < 2*3 - 2; ++residue) {
//        Q_part_list v = Q_table[residue];
//        printf("residue: %ld\n", residue);
//        for(int i = 0; i < v.length ; i++){
//            Q_part q_part = v.list[i];
//            printf(" %ld , %ld\n", q_part.bit_string, q_part.degree);
//        }
//        printf("\n");
//    }
//
//    generateMilnorBasisPpartTable(algebra, 10);
//    P_part_list * P_table = algebra->P_table;
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

    GenerateMilnorBasis(algebra, 76);
    MilnorBasisElement m;
    m.p_degree = 0;
    m.p_length = 0;
    m.q_degree = 5;
    m.q_part = 2;
    printf("test:\n");
    MonomialIndex idx1 = GetIndexFromMilnorBasisElement(algebra, m);
    GetMilnorBasisElementFromIndex(algebra, idx1);


//    MilnorBasisElement_list * basis_table = algebra->basis_table;
//    for(int n = 0; n < algebra->basis_max_degree; n++){
//        printf("deg: %d \n", n);
//        for(int i = 0; i < basis_table[n].length; i++){
//            MilnorBasisElement v = basis_table[n].list[i];
//            milnor_basis_element_to_string(buffer, &v);
//            printf("  %d: '%s'\n", i , buffer);
//        }
//        printf("\n");
//    }

//    P_part r, s;
//    r.length = 3;
//    unsigned long r_list[3] = {2, 5, 1};
//    r.p_part = r_list;
//    s.length = 3;
//    unsigned long s_list[3] = {2, 2, 2};
//    s.p_part = s_list;
//    auto M = initialize_milnor_matrix(r, s);
//    do {
//        cout << "[\n";
//        for(int row = 0; row <= r.length; row++) {
//            cout << "  [";
//            for (int col = 0; col <= s.length; col++) {
//                cout << M[row][col] << ", ";
//            }
//            cout << "]\n";
//        }
//        cout << "]\n\n";
//    } while(step_milnor_matrix(2, M, r, s));
    m.p_degree = 0;
    m.p_length = 0;
    m.q_degree = 5;
    m.q_part = 2;
    MonomialIndex idx = GetIndexFromMilnorBasisElement(algebra, m);

    MilnorBasisElement A,B;
    A = milnor_basis_element_from_string(algebra, "P(2)");
    B = milnor_basis_element_from_string(algebra, "Q(1)");
    MilnorElement * result = allocateMilnorElement(algebra, A.p_degree + A.q_degree + B.p_degree + B.q_degree);
    MilnorProduct(algebra, result, A, B);
    string str1 = buffer;
    int len1 = milnor_element_to_string(str1, algebra, result);
    string str2 = str1 + len1 + 1;
    int len2 = milnor_basis_element_to_string(str2, &A);
    string str3 = str2 + len2 + 1;
    milnor_basis_element_to_string(str3, &B);
    printf("%s * %s = %s\n", str2, str3, str1);
}

/**/



//string to_string(MilnorBasisElement m) {
//    return fmt.Sprintf("{%v %v}", m.q_part, m.p_part)
//}



// TODO: Break this up into several methods.
void generateMilnorBasisPpartTable(MilnorAlgebra * algebra, unsigned long n) {
    unsigned long p = algebra->p;
    unsigned long * xi_degrees = getXiDegrees(p);

    unsigned long profile_list[MAX_XI_TAU];
    for(long idx = 0;  idx < MAX_XI_TAU; idx ++){
        profile_list[idx] = getProfileExponent(algebra->profile, p, idx) - 1;
    }

    unsigned long old_max_degree = algebra -> P_table_max_degree;
    unsigned long new_max_degree = n;
    printf("old_max: %ld, new_max: %ld\n", old_max_degree, new_max_degree);
    if(new_max_degree < old_max_degree){
        return;
    }
    algebra -> P_table_max_degree = new_max_degree;
    algebra -> P_table = (P_part_list*) realloc(algebra -> P_table , (new_max_degree + 1) * sizeof(P_part_list));
    algebra -> P_table_by_P_length = (P_part_list**) realloc(algebra -> P_table_by_P_length, (new_max_degree + 1) * sizeof(P_part_list[MAX_XI_TAU]));
    P_part_list* degree_table = algebra->P_table;
    P_part_list** tau_table = algebra->P_table_by_P_length;

    // If old_max_degree is 0 we need to set up the base case for our recursion.
    // In degree 0 there's a single partition of length 0.
    P_part_list degree_list, tau_list;
    if(old_max_degree == 0){
        degree_list.length = 1;
        degree_list.list = (P_part*)malloc(sizeof(P_part));
        degree_list.list->degree = 0;
        degree_list.list->length = 0;
        degree_list.list->p_part = (unsigned long*)calloc(1, sizeof(unsigned long));
        degree_table[0] = degree_list;
        tau_table[0] = (P_part_list*) calloc(10, sizeof(P_part_list));
        tau_table[0][0] = degree_list;
    }

    // Populate the new table.
    // We iterate over each degree and for each degree n iterate over each xi_i until we hit one of degree larger
    // than n and look for monomials of degree n - |xi_i|. Copy each such monomial, increment the power of xi_i
    // and add it to the degree n component.


    for(unsigned long current_degree = old_max_degree + 1; current_degree <= new_max_degree; current_degree++){
        // On the first pass we compute the size of the degree list.
        // The next pile of loops will have the same structure as this one but do more.
        unsigned long degree_list_capacity = 0;

        // Loop over the highest degree xi_i in the new monomial
        for(int xi = 0; xi < MAX_XI_TAU; xi ++){
            unsigned long xi_degree = xi_degrees[xi];
            if (xi_degree > current_degree) {
                break;
            }
            unsigned long remaining_degree = current_degree - xi_degree;
            // For every monomial m of degree n - |xi_i| such that the max xi appearing in m is <= i
            // we'll get a monomial xi_i * m of degree n.
            // Loop over the highest degree xi_i in the old monomial
            for(int xi_old = 0; xi_old <= xi; xi_old ++) {
                degree_list_capacity += tau_table[remaining_degree][xi_old].length;
            }
        }

        tau_table[current_degree] = (P_part_list*) calloc(10, sizeof(P_part_list));

        degree_list.list = (P_part*)malloc(degree_list_capacity*sizeof(P_part));
        degree_list.length = 0;

        // Now populate the tables.
        // Loop over the highest degree xi_i in the new monomial
        for(unsigned long xi = 0; xi < MAX_XI_TAU; xi ++){
            unsigned long xi_degree = xi_degrees[xi];
            if(xi_degree > current_degree){
                break;
            }
            unsigned long remaining_degree = current_degree - xi_degree;
            unsigned long tau_list_capacity = 0;
            // xi_old is the max xi present in the remaining degree monomial.
            // We do another degree calculation for tau_list
            for(int xi_old = 0; xi_old <= xi; xi_old ++) {
                tau_list_capacity += tau_table[remaining_degree][xi_old].length;
            }
            tau_list.list = (P_part*)malloc(tau_list_capacity*sizeof(P_part));
            tau_list.length = 0;
            // Now everything is allocated so we can do the real work.
            // Iterate over the highest xi_i in the old monomial
            for(unsigned long xi_old = 0; xi_old <= xi; xi_old ++) {
                P_part_list remaining_degree_list = tau_table[remaining_degree][xi_old];
                // Iterate over the monomials of degree n - |xi| with highest xi xi_old.
                for (unsigned long i = 0; i < remaining_degree_list.length; i++) {
                    P_part old_p_part = remaining_degree_list.list[i];
                    if(xi_old == xi && old_p_part.p_part[xi] == profile_list[xi]){
                        continue;
                    }
                    // Copy the old p_part to a new p_part with length xi + 1.
                    P_part new_p_part;
                    new_p_part.degree = current_degree;
                    new_p_part.length = xi + 1;
                    new_p_part.p_part = (unsigned long *) calloc(new_p_part.length, sizeof(unsigned long));
                    memcpy(new_p_part.p_part, old_p_part.p_part, old_p_part.length * sizeof(unsigned long));
                    // Increment the exponent of xi_i
                    new_p_part.p_part[xi]++;
                    // Add our new_p_part to the tables.
                    degree_list.list[degree_list.length] = new_p_part;
                    tau_list.list[tau_list.length] = new_p_part;
                    degree_list.length++;
                    tau_list.length++;
                }
            }
            tau_table[current_degree][xi] = tau_list;
        }
        degree_table[current_degree] = degree_list;
    }
}

// Generate the MilnorBasisQpartTable through degree n.
void generateMilnorBasisQpartTable(MilnorAlgebra * algebra, unsigned long max_degree){
    unsigned long p = algebra->p;
    unsigned long q = 2*(p-1);
    unsigned long * tau_degrees = getTauDegrees(p);
    unsigned long profile = algebra->profile.q_part;

    unsigned long old_max_tau = algebra->Q_table_max_tau;
    // New max tau is number of tau_i's less than max_degree.
    long new_max_tau = 0;
    for( ; tau_degrees[new_max_tau] < max_degree; new_max_tau++){}

    if(new_max_tau <= old_max_tau){
        return;
    }
    algebra->Q_table_max_tau = new_max_tau;

    // Initialize map residue --> Q_part_list
    if(old_max_tau == 0) {
        algebra->Q_table = (Q_part_list *) calloc(q, sizeof(Q_part_list));
    }
    Q_part_list* table = algebra->Q_table;

    long available_taus = 0;
    for(long i = 0; i < new_max_tau; i++){
        available_taus += (profile >> i) & 1;
    }

    profile = ~profile;

    // Calculate how many elements are going to land in each bin. There are available_taus choose i monomials in tau_0, ..., tau_(n-1)
    // of length i and these land in residue class i % q.
    long capacities[q];
    memset(capacities, 0, q * sizeof(long));
    long n_choose_i = 1;
    capacities[0] += n_choose_i;
    for(int i = 1; i <= available_taus; i++){
        n_choose_i *= (available_taus + 1 - i);
        n_choose_i /= i;
        capacities[i % q] += n_choose_i;
    }
    // Allocate more space for the lists, but don't change the length yet because we need to know where to insert our
    // new elements. So for now capacities list knows the capacity, the length attribute stores the actual number of
    // elements present in the list.
    for (int residue = 0; residue < q; residue++) {
        table[residue].list = (Q_part*) realloc(table[residue].list, capacities[residue]* sizeof(Q_part));
    }

    if(old_max_tau == 0) {
        // degree 0 is an edge case because n = 0 is the only integer such that there is not a unique nonnegative
        // exponent i with 2^i <= n < 2^(i+1). It goes in residue 0.
        table[0].list[0] = (Q_part){0, 0};
        table[0].length++;
    }

    unsigned long total = 0;
    //The total starts out as tau_degrees[old_max_tau] but we update by tau_degrees[old_max_tau] - Sum(smaller tau_i's).
    //So initialize total = Sum(smaller tau_i's)
    for(unsigned long i = 0; i < old_max_tau; i++){
        total += tau_degrees[i];
    }
    unsigned long bit_string_min = 1 << old_max_tau;
    unsigned long bit_string_max = 1 << new_max_tau;
    //The residue starts out as 1, but we update by 1 - # of trailing 0's.
    //On the first pass, # of trailing 0's is old_max_tau, so initialize residue = old_max_tau
    unsigned long residue = old_max_tau;
    for(unsigned long bit_string = bit_string_min;  bit_string < bit_string_max; bit_string++) {
        // Iterate over trailing zero bits
        unsigned long v = (bit_string ^ (bit_string - 1)) >> 1;
        unsigned long c = 0; // c is counting the number of zero bits at the end of bit_string.
        for (; v != 0; c++) {
            v >>= 1;
            total -= tau_degrees[c];
        }
        total += tau_degrees[c];
        // residue is the number of 1's in bit_string mod q.
        // c bits were unset when we incremented and one new bit was set so we have to add 1 - c.
        residue += 1 - c;
        if(bit_string & profile != 0){
            continue;
        }
        residue = ModPositive(residue, q);
        table[residue].list[table[residue].length] = (Q_part){total, bit_string};
        table[residue].length++;
    }
}

void freeMilnorAlgebraPQTables(MilnorAlgebra * algebra){

}

void freeMilnorAlgebra(MilnorAlgebra * algebra){
    //algebra ->
}

void GenerateMilnorBasis(MilnorAlgebra * algebra, unsigned long max_degree) {
    printf("Start MilnorBasis\n");
    printf("algebra: %lx, %lx, %lx, %lx\n", (long)algebra->basis_table, (long)algebra->P_table_max_degree, (long)algebra->P_table, (long)algebra->Q_table);
    unsigned long p = algebra->p;
    printf("p: %ld\n", p);
    initializePrime(p);
    unsigned long old_max_degree = algebra->basis_max_degree;
    algebra->basis_table = (MilnorBasisElement_list *) realloc(
            algebra->basis_table,
            (max_degree + 1) * sizeof(MilnorBasisElement_list)
    );
    algebra->basis_max_degree = max_degree;
    algebra->basis_name_to_index_map = (khash_t(monomial_index_map) **)realloc(
            algebra->basis_name_to_index_map,
            (max_degree + 1)* sizeof(khash_t(monomial_index_map) *)
        );
    MilnorBasisElement_list * table = algebra->basis_table;
    khash_t(monomial_index_map) ** name_table = algebra->basis_name_to_index_map;
    for(long k = old_max_degree + 1; k <= max_degree; k++){
        name_table[k] = kh_init(monomial_index_map);
    }
    if(algebra->generic){
        GenerateMilnorBasisGeneric(algebra, old_max_degree, max_degree);
    } else {
        GenerateMilnorBasis2(algebra, old_max_degree, max_degree);
    }
}

void GenerateMilnorBasis2(MilnorAlgebra * algebra, unsigned long old_n, unsigned long n){
    generateMilnorBasisPpartTable(algebra, n);
    unsigned long p = algebra->p;
    Profile profile = algebra->profile;

    MilnorBasisElement_list * table = algebra->basis_table;
    khash_t(monomial_index_map) ** name_table = algebra->basis_name_to_index_map;

    MilnorBasisElement_list current_degree_list;
    for(unsigned long degree = old_n + 1; degree <= n; degree++){
        P_part_list p_parts = algebra->P_table[degree];
        current_degree_list.length = 0;
        current_degree_list.list = (MilnorBasisElement*) malloc(p_parts.length * sizeof(MilnorBasisElement));
        for(unsigned long i = 0; i < p_parts.length; i++){
            P_part x = p_parts.list[i];
            MilnorBasisElement m = (MilnorBasisElement){0, 0, degree, x.length, x.p_part};
            current_degree_list.list[current_degree_list.length] = m;
            current_degree_list.length ++;
            char key[200];
            milnor_basis_element_to_string(key, &m);
            int absent;
            khint_t k = kh_put(monomial_index_map, name_table[degree], key, &absent);
            if(absent){
                kh_key(name_table[degree], k) = strdup(key);
            }
            kh_val(name_table[degree], k) = (MonomialIndex){degree, current_degree_list.length};
        }
        table[degree] = current_degree_list;
    }
}


// Get the basis in degree n for the generic steenrod algebra at the prime p.
// We just put together the "even part" of the basis and the "Q part".
void GenerateMilnorBasisGeneric(MilnorAlgebra * algebra, unsigned long old_max_degree, unsigned long new_max_degree){
    unsigned long p = algebra->p;
    unsigned long q = 2*(p-1);
    generateMilnorBasisPpartTable(algebra, new_max_degree / q);
    generateMilnorBasisQpartTable(algebra, new_max_degree);
    MilnorBasisElement_list* table = algebra->basis_table;
    khash_t(monomial_index_map) ** name_table = algebra->basis_name_to_index_map;

    MilnorBasisElement_list current_degree_list;
    for(unsigned long degree = old_max_degree + 1; degree <= new_max_degree; degree++){
        // p_deg records the desired degree of the P part of the basis element.
        // Since p-parts are always divisible by 2p-2, we divide by this first.
        // pow(p, -1) returns 1, so min_q_deg is 0 if q divides n evenly.
        Q_part_list q_list = algebra->Q_table[degree % q];
        // Let's count the length first.
        unsigned long table_size = 0;
        for(long i = 0; i < q_list.length; i++) {
            Q_part q_part = q_list.list[i];
            unsigned long q_deg = q_part.degree;
            if(q_deg > degree){
                break;
            }
            unsigned long p_deg = (degree - q_deg) ;
            P_part_list p_list = algebra->P_table[p_deg / q];
            table_size += p_list.length;
        }

        current_degree_list.length = 0;
        current_degree_list.list = (MilnorBasisElement*) malloc(table_size * sizeof(MilnorBasisElement));
        for(long i = 0; i < q_list.length; i++) {
            Q_part q_part = q_list.list[i];
            unsigned long q_deg = q_part.degree;
            if(q_deg > degree){
                break;
            }
            unsigned long q_bit_string = q_part.bit_string;
            unsigned long p_deg = (degree - q_deg);
            P_part_list p_list = algebra->P_table[p_deg / q];
            for(long j = 0; j < p_list.length; j++) {
                P_part p_part = p_list.list[j];
                MilnorBasisElement m = (MilnorBasisElement){q_part.degree, q_bit_string, p_deg, p_part.length, p_part.p_part};
                char key[200];
                milnor_basis_element_to_string(key, &m);
                int absent;
                khint_t k = kh_put(monomial_index_map, name_table[degree], key, &absent);
                if(absent){
                    kh_key(name_table[degree], k) = strdup(key);
                }
                kh_value(name_table[degree], k) = (MonomialIndex){degree, current_degree_list.length};

                current_degree_list.list[current_degree_list.length] = m;
                current_degree_list.length ++;
            }
        }
        table[degree] = current_degree_list;
    }
    MilnorBasisElement m;
    m.p_degree = 0;
    m.p_length = 0;
    m.q_degree = 5;
    m.q_part = 2;
    GetIndexFromMilnorBasisElement(algebra, m);
}


MilnorBasisElement GetMilnorBasisElementFromIndex(MilnorAlgebra * algebra, MonomialIndex idx) {
    return algebra->basis_table[idx.degree].list[idx.index];
}

MonomialIndex GetIndexFromMilnorBasisElement(MilnorAlgebra * algebra, MilnorBasisElement b){
    kh_monomial_index_map_t * map = algebra->basis_name_to_index_map[b.q_degree + b.p_degree];
    char buffer[200];
    milnor_basis_element_to_string(buffer, &b);
    khint_t k = kh_get(monomial_index_map, map, buffer);
    if(k == kh_end(map)){
        printf("Uh-oh, not here. degree: %ld, elt: '%s'\n", b.q_degree + b.p_degree, buffer);
        return (MonomialIndex){-1, -1};
    }
    return kh_val(map, k);
}


unsigned long** allocate_milnor_matrix(unsigned long rows, unsigned long cols)  {
    unsigned long** M = (unsigned long**) calloc(1, rows * sizeof(unsigned long*) + rows * cols * sizeof(unsigned long));
    for(long row = 0; row < rows; row++){
        M[row] = (unsigned long*)(M + rows) + row * cols;
    }
    return M;
}

// Initializes an len(r)+1 by len(s)+1 matrix
// Puts r along the first column and s along the first row and zeroes everywhere else.
unsigned long** initialize_milnor_matrix(P_part r, P_part s) {
    unsigned long rows = r.length + 1;
    unsigned long cols = s.length + 1;
    unsigned long** M = allocate_milnor_matrix(rows, cols);
    memcpy(M[0] + 1, s.p_part, s.length * sizeof(unsigned long));
    for(long i = 0; i < r.length; i++){
        M[i+1][0] = r.p_part[i];
    }
    return M;
}

// This seems to move an i x j block of M back to the first row and column.
// To be honest, I don't really know what the point is, but the milnor_matrices
// function was a little long and this seemed like a decent chunk to extract.
// At least it contains all of the steps that modify M so that seems like a good thing.
void step_milnor_matrix_update_matrix(unsigned long **M, P_part r, P_part s, int i, int j, int x)  {
    unsigned long rows = r.length + 1;
    unsigned long cols = s.length + 1;
    for(long row = 1; row < i; row ++){
        M[row][0] = r.p_part[row-1];
        for(long col = 1; col < cols; col++){
            M[0][col] += M[row][col];
            M[row][col] = 0;
        }
    }
    for(long col = 1; col < j; col++){
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
bool step_milnor_matrix(unsigned long  p, unsigned long ** M, P_part r, P_part s) {
    long rows = r.length + 1;
    long cols = s.length + 1;
    for(long i = 1; i < rows; i++){
        unsigned long total = M[i][0];
        unsigned long p_to_the_j = 1;
        for(long j = 1; j < cols; j++){
            p_to_the_j *= p;
            if(total < p_to_the_j){
                // We don't have enough weight left in the entries above this one in the column to increment this cell.
                // Add the weight from this cell to the total, we can use it to increment a cell lower down.
                total += M[i][j] * p_to_the_j;
                continue;
            }

            // Check if any entry in column j above row i is nonzero. I'm still not sure why tbh.
            for(long k = 0; k < i; k++){
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

long max(long a, long b){
    if(a > b){
        return a;
    }
    return b;
}

long min(long a, long b){
    if(a < b){
        return a;
    }
    return b;
}



// Handles the multiplication in the even subalgebra of the Steenrod algebra P.
// When p = 2, this is isomorphic to the whole Steenrod algebra so this method does everything.
void MilnorProductEven(MilnorAlgebra * algebra, MilnorElement * result, MilnorBasisElement r_elt, MilnorBasisElement s_elt){
    unsigned long p = algebra->p;
    P_part r, s;
    r.degree = r_elt.p_degree;
    r.length = r_elt.p_length;
    r.p_part = r_elt.p_part;
    s.degree = s_elt.p_degree;
    s.length = s_elt.p_length;
    s.p_part = s_elt.p_part;
    unsigned long output_degree = (r_elt.p_degree + s_elt.p_degree);
    unsigned long rows = r.length + 1;
    unsigned long cols = s.length + 1;
    unsigned long number_of_diagonals = r.length + s.length;
    unsigned long** M = initialize_milnor_matrix(r, s);

    if(result->degree != output_degree){
        printf("Result has degree %ld but output has degree %ld. ", result->degree, output_degree);
        printf("I don't know how we handle errors though.\n");
        return;
    }
    memset(result->vector, 0, result->algebra_dimension * sizeof(long));
    do {
        // check diagonals
        long coeff = 1; // I think this needs to be signed.
        unsigned long *diagonal_sums = (unsigned long*)malloc(number_of_diagonals * sizeof(unsigned long));
        unsigned long max_nonzero_diagonal = 0;
        for(unsigned long diagonal_index = 1; diagonal_index <= number_of_diagonals; diagonal_index++){
            // We're going to iterate along the diagonal and copy the diagonal into nth_diagonal
            // and put the sum into nth_diagonal_sum.
            unsigned long i_min = max(0, diagonal_index - cols + 1);
            unsigned long i_max = min(1 + diagonal_index, rows);
            unsigned long diagonal_length = i_max - i_min;
            unsigned long diagonal[diagonal_length];
            unsigned long diagonal_sum = 0;
            long index = 0;
            // Iterate along diagonal
            for(long i = i_min; i < i_max;  i++){
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
        if(coeff == 0) {
            free(diagonal_sums);
        } else {
            diagonal_sums = (unsigned long*)realloc(diagonal_sums, max_nonzero_diagonal * sizeof(unsigned long));
            MilnorBasisElement m = (MilnorBasisElement){0, 0, output_degree, max_nonzero_diagonal, diagonal_sums};
            MonomialIndex idx = GetIndexFromMilnorBasisElement(algebra, m);
            addBasisElementToMilnorElement(result, idx, coeff);
        }
    } while(step_milnor_matrix(p, M, r, s));
}



// Reduce m1 * f = (Q_e0 Q_e1 ... P(r1, r2, ...)) * (Q_f0 Q_f1 ...) into the form Sum of Q's * P's
// Result is represented as dictionary of pairs of tuples.
void MilnorProductFullQpart(MilnorAlgebra * algebra, MilnorElement * output, MilnorBasisElement m1, unsigned long f) {
    unsigned long p = algebra->p;
    unsigned long q;
    if(algebra->generic){
        q = 2*p-2;
    } else {
        q = 1;
    }
    unsigned long * tau_degrees = getTauDegrees(p);
    unsigned long * xi_degrees = getXiDegrees(p);
    unsigned long p_degree = m1.p_degree;
    unsigned long q_degree = m1.q_degree;
    MonomialIndex idx = GetIndexFromMilnorBasisElement(algebra, m1);
    MilnorElement result, old_result;
    result.algebra = algebra;
    result.degree = p_degree + q_degree;
    result.algebra_dimension = result.algebra->basis_table[result.degree].length;
    long result_vector1[MAX_ALGEBRA_DIMENSION];
    long result_vector2[MAX_ALGEBRA_DIMENSION];
    memset(result_vector1, 0, result.algebra_dimension * sizeof(long));
    result.vector = result_vector1;
    result.vector[idx.index] = 1;
    old_result.vector = result_vector2;
//    char buffer[200]; // For test print statements
    unsigned long p_to_the_k = 0;
    for(unsigned long k = 0; ( f & ~((1 << k) - 1) ) != 0; k++){
        p_to_the_k *= p;
        if(k == 0){
            p_to_the_k = 1;
        }
        if((f & (1 << k)) == 0){
            continue;
        }
        q_degree += tau_degrees[k];
        long * swap_temp = old_result.vector;
        old_result = result;
        result.vector = swap_temp;
        result.degree = p_degree + q_degree;
        result.algebra_dimension = result.algebra->basis_table[result.degree].length;
        memset(result.vector, 0, result.algebra_dimension * sizeof(long));
//        milnor_element_to_string(buffer, algebra, &old_result);
//        printf("old_result: %s\n", buffer);
        for(long idx = 0; idx < old_result.algebra_dimension; idx++){
            if(old_result.vector[idx] == 0){
                continue;
            }
            MonomialIndex monomialIndex;
            monomialIndex.degree = old_result.degree;
            monomialIndex.index = idx;
            MilnorBasisElement mono = GetMilnorBasisElementFromIndex(algebra, monomialIndex);
//            milnor_basis_element_to_string(buffer, &mono);
//            printf("    MONO: %s,   q_deg: %ld\n", buffer, mono.q_degree);
            unsigned long coeff = old_result.vector[idx];
            for(unsigned long i = 0; i < mono.p_length + 1; i++){
                if((mono.q_part & (1 << k+i)) != 0){
                    continue;
                }
                // Make sure mono.p_part[i - 1] is large enough to deduct p^k from it
                if(i > 0 && mono.p_part[i - 1] < p_to_the_k ){
                    continue;
                }
//                milnor_basis_element_to_string(buffer, &mono);
//                printf("        mono: %s, i: %ld\n",buffer, i);

                unsigned long q_mono = mono.q_part;
                unsigned long q_degree = mono.q_degree;
                unsigned long *p_mono = mono.p_part;
                unsigned long p_mono_length = mono.p_length;
                unsigned long p_degree = mono.p_degree;

                if(i > 0){
                    p_mono[i - 1] -= p_to_the_k;
                    unsigned long new_p_mono_length = p_mono_length;
                    for( ; p_mono[new_p_mono_length - 1] == 0; new_p_mono_length -- ){}
                    unsigned long * new_p_mono = (unsigned long*)malloc(new_p_mono_length * sizeof(unsigned long));
                    memcpy(new_p_mono, p_mono, new_p_mono_length * sizeof(unsigned long));
                    p_mono[i - 1] += p_to_the_k;
                    p_degree -= q * p_to_the_k * xi_degrees[i-1];
                    p_mono = new_p_mono;
                    p_mono_length = new_p_mono_length;
                }
                // insert(q_mono, len(q_mono) - qs_gt_k_plus_i, k+i)
                q_mono += (1 << k + i);
                q_degree += tau_degrees[k+i];

                unsigned long larger_Qs = 0;
                unsigned long v = q_mono >> k + i + 1;
                for( ; v != 0; v >>= 1){
                    larger_Qs += v & 1;
                }
                coeff *= MinusOneToTheN(larger_Qs);
                MilnorBasisElement m = (MilnorBasisElement){q_degree, q_mono, p_degree, p_mono_length, p_mono};

//                milnor_basis_element_to_string(buffer, &m);
//                printf("        MBE: %s\n", buffer);
//                printf("        q_deg: %ld, tau_deg(%ld) = %ld\n", q_degree, k+i, tau_degrees[k+i]);
                MonomialIndex out_idx = GetIndexFromMilnorBasisElement(algebra, m);
                addBasisElementToMilnorElement(&result, out_idx, coeff);
//                milnor_element_to_string(buffer, algebra, &result);
//                printf("        result: %s\n", buffer);
//                printf("  new vector[%ld]: %ld\n", out_idx.index, result.vector[out_idx.index]);
            }
        }
    }
//    milnor_element_to_string(buffer, algebra, &result);
//    printf("result.degree: %ld, result: %s\n", result.degree, buffer);
    assignMilnorElement(output, &result);
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
void MilnorProductFull(MilnorAlgebra * algebra, MilnorElement * result, MilnorBasisElement m1, MilnorBasisElement m2) {
    unsigned long p = algebra-> p;
    MilnorBasisElement s = (MilnorBasisElement){0, 0, m2.p_degree, m2.p_length, m2.p_part};
    MilnorElement m1_times_f;
    m1_times_f.algebra = algebra;
    m1_times_f.degree = m1.q_degree + m1.p_degree + m2.q_degree;
    m1_times_f.algebra_dimension = algebra->basis_table[m1_times_f.degree].length;
    long m1_times_f_vector[m1_times_f.algebra_dimension];
    memset(m1_times_f_vector, 0, m1_times_f.algebra_dimension * sizeof(long));
    m1_times_f.vector = m1_times_f_vector;
    MilnorProductFullQpart(algebra, &m1_times_f, m1, m2.q_part);
    // Now for the Milnor matrices.  For each entry '(e,r): coeff' in answer,
    // multiply r with s.  Record coefficient for matrix and multiply by coeff.
    // Store in 'result'.
//    char buffer[200];
//    milnor_element_to_string(buffer, algebra, &m1_times_f);
//    printf("m1_times_f: %s\n", buffer);
    if(m2.p_length == 0){
        assignMilnorElement(result, &m1_times_f);
        return;
    }
    unsigned long output_degree = m1.p_degree + m1.q_degree + m2.q_degree + m2.q_degree;
    for(long i = 0; i < m1_times_f.algebra_dimension; i++){
        if(m1_times_f.vector[i] == 0){
            continue;
        }
        MonomialIndex er_idx = (MonomialIndex){m1_times_f.degree, i};
        long coeff = m1_times_f.vector[i];
        MilnorBasisElement er_mono = GetMilnorBasisElementFromIndex(algebra, er_idx);
        MilnorBasisElement r = (MilnorBasisElement){0, 0, er_mono.p_degree, er_mono.p_length, er_mono.p_part};
        MilnorElement prod;
        prod.algebra = algebra;
        prod.degree = r.p_degree + s.p_degree;
        prod.algebra_dimension = algebra->basis_table[prod.degree].length;
        long prod_vector[prod.algebra_dimension];
        memset(prod_vector, 0, prod.algebra_dimension * sizeof(long));
        prod.vector = prod_vector;
        MilnorProductEven(algebra, &prod, r, s);
        for(long j = 0; j < prod.algebra_dimension; j++){
            MonomialIndex idx = (MonomialIndex){prod.degree, j};
            long c = prod.vector[j];
            MilnorBasisElement m = GetMilnorBasisElementFromIndex(algebra, idx);
            m.q_degree = er_mono.q_degree;
            m.q_part = er_mono.q_part;
            MonomialIndex out_idx = GetIndexFromMilnorBasisElement(algebra, m);
            addBasisElementToMilnorElement(result, out_idx, coeff*c);
        }
    }
}

// Multiplication of Milnor basis elements in the non generic case.
void MilnorProduct2(MilnorAlgebra * algebra, MilnorElement * result, MilnorBasisElement r, MilnorBasisElement s){
    MilnorProductEven(algebra, result, r, s);
}

void MilnorProductGeneric(MilnorAlgebra * algebra, MilnorElement * result, MilnorBasisElement r, MilnorBasisElement s ) {
    MilnorProductFull(algebra, result, r, s);
}


// Multiply r and s in the Milnor algebra determined by algebra.
// Note that since profile functions determine subalgebras, the product
// doesn't need to care about the profile function at all.
void  MilnorProduct(MilnorAlgebra * algebra, MilnorElement * result, MilnorBasisElement r, MilnorBasisElement s)  {
//    if(r.q_degree + r.p_degree + s.q_degree + s.p_degree != result->degree){
//        printf("Result has degree %ld but should have degree %ld\n", result->degree, r.q_degree + r.p_degree + s.q_degree + s.p_degree);
//        return;
//    }
    if(algebra->generic){
        MilnorProductFull(algebra, result, r, s);
    } else {
        MilnorProductEven(algebra, result, r, s);
    }
}
/**/
