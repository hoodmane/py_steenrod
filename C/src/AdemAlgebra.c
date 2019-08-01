#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "khash.h"
#include "combinatorics.h"
#include "AdemAlgebra.h"


typedef struct {
    uint length;
    uint *list;
}  bockstein_list;

static bockstein_list bockstein_table[MAX_XI_TAU] = {0};
static uint bockstein_table_inner[(1<<MAX_XI_TAU) - 1];

static void initializeBocksteinTable(){
    if(bockstein_table[0].list != NULL){
        return;
    }    
    uint n_choose_k = 1;
    uint *bockstein_table_ptr = bockstein_table_inner;
    for(uint k=1; k<=MAX_XI_TAU; k++){
        bockstein_table[k-1].length = 0;
        bockstein_table[k-1].list = bockstein_table_ptr;
        bockstein_table_ptr += n_choose_k;
        n_choose_k *= (MAX_XI_TAU + 1 - k);
        n_choose_k /= k; 
    }

    for(uint i=0; i<((1<<MAX_XI_TAU) - 1); i++){
        uint bits_set = __builtin_popcount(i);
        bockstein_table[bits_set].list[bockstein_table[bits_set].length] = i;
        bockstein_table[bits_set].length ++;
    }
}

int AdemAlgebra_excessSortOrder(const void *a, const void *b){
    AdemBasisElement *x = (AdemBasisElement*) a;
    AdemBasisElement *y = (AdemBasisElement*) b;
    return x->excess - y->excess; 
}

// We need this for generic basis generation.
int AdemAlgebra_lengthSortOrder(const void *a, const void *b){
    AdemBasisElement *x = (AdemBasisElement*) a;
    AdemBasisElement *y = (AdemBasisElement*) b;
    return -x->P_length + y->P_length; 
}

KHASH_MAP_INIT_STR(monomial_index_map, uint)

typedef struct {
    AdemAlgebra public_algebra;
    uint basis_max_degree;
    AdemBasisElement_list *even_basis_table;
    AdemBasisElement_list *basis_table;
    khash_t(monomial_index_map) **basis_element_to_index_map; // degree -> admissible sequence -> index
    Vector ****multiplication_table;// degree -> first square -> admissibile sequence idx -> result vector
    uint **excess_table;
    AdemBasisElement_list filtrationOneProduct_basisElements;
} AdemAlgebraInternal;

static void AdemAlgebra__initializeFields(AdemAlgebraInternal *algebra, uint p, bool generic, bool unstable);
uint AdemAlgebra__generateName(AdemAlgebra *algebra); // defined in adem_io
AdemAlgebra *AdemAlgebra_construct(uint p, bool generic, bool unstable){
    uint num_products;
    if(generic){
        num_products = 2;
    } else {
        num_products = 3;
    }
    size_t algebra_size =         
        sizeof(AdemAlgebraInternal) 
        + num_products * sizeof(FiltrationOneProduct); // product space

    AdemAlgebraInternal *algebra = malloc(algebra_size);
    Algebra *inner_algebra = &algebra->public_algebra.algebra;
    inner_algebra->product_list.length = num_products;
    inner_algebra->product_list.list = (FiltrationOneProduct*)(algebra + 1);    
    assert((char*)(inner_algebra->product_list.list + num_products) == ((char *)algebra) + algebra_size);

    AdemAlgebra__initializeFields(algebra, p, generic, unstable);
    AdemAlgebra__generateName((AdemAlgebra*)algebra);
    return (AdemAlgebra*)algebra;    
}

static void AdemAlgebra__initializeFields(AdemAlgebraInternal *algebra, uint p, bool generic, bool unstable){
    initializePrime(p);
    algebra->public_algebra.algebra.p = p;
    algebra->public_algebra.generic = generic;
    algebra->public_algebra.unstable = unstable;
    algebra->public_algebra.sort_order = NULL;
    // Fill in the Algebra function pointers.
    algebra->public_algebra.algebra.computeBasis = AdemAlgebra_generateBasis;
    algebra->public_algebra.algebra.multiplyBasisElements = AdemAlgebra_multiply;
    algebra->public_algebra.algebra.basisElementToString = AdemAlgebra_basisElementIndex_toString;
    if(unstable){
        algebra->public_algebra.sort_order = AdemAlgebra_excessSortOrder;
        algebra->public_algebra.algebra.getDimension = AdemAlgebra_getDimension_unstable;
    } else {
        algebra->public_algebra.algebra.getDimension = AdemAlgebra_getDimension;
    }

    algebra->public_algebra.algebra.max_degree = 0;

    algebra->even_basis_table = NULL;
    algebra->basis_table = NULL;
    algebra->basis_element_to_index_map = NULL;
    algebra->multiplication_table = NULL;
    algebra->excess_table = NULL;

    // Products
    // Length field has to match with amount of space we decided to allocate for this
    // in constructor function directly before this.
    FiltrationOneProduct *product_list = algebra->public_algebra.algebra.product_list.list;
    uint num_basis_elements;
    uint total_p_length;
    if(generic){
        num_basis_elements = 2;
        total_p_length = 1;
    } else {
        num_basis_elements = 3;
        total_p_length = 3;
    }
    algebra->filtrationOneProduct_basisElements.length = num_basis_elements;
    AdemBasisElement **product_basisElement_table = malloc(
        num_basis_elements*(sizeof(AdemBasisElement*) + sizeof(AdemBasisElement)) + total_p_length*sizeof(uint)
    );

    AdemBasisElement *basis_elt_ptr = (AdemBasisElement*)(product_basisElement_table + num_basis_elements);
    for(uint i = 0; i < num_basis_elements; i++){
        product_basisElement_table[i] = basis_elt_ptr;
        basis_elt_ptr ++;
    }
    uint *Ps_ptr = (uint*)basis_elt_ptr;

    algebra->filtrationOneProduct_basisElements.list = product_basisElement_table;
    if(generic){
        AdemBasisElement *b = *product_basisElement_table;
        // beta
        b->degree = 1;
        b->bocksteins = 1;
        b->P_length = 0;
        b->Ps = NULL;
        b++;
        // P1
        b->degree = 2*(p-1);
        b->bocksteins = 0;
        b->P_length = 1;
        b->Ps = Ps_ptr;
        b->Ps[0] = 1;
        b++;
    } else {
        AdemBasisElement *b = *product_basisElement_table;
        for(int hi = 0; hi < 3; hi++){
            uint degree = 1 << hi; // degree is 2^hi 
            b->degree = degree;
            b->bocksteins = 0;
            b->P_length = 1;
            b->Ps = Ps_ptr;
            Ps_ptr += b->P_length;
            b->Ps[0] = degree;
            b++;
        }
    }
    if(generic){
        product_list[0].type = "a0";
        product_list[0].degree = 1;
        product_list[0].index = 0;
        product_list[1].type = "h0";
        product_list[1].degree = 2 * (p - 1); // P1
        product_list[1].index = 0;
    } else {
        product_list[0].type = "h0";
        product_list[0].degree = 1; // Sq1
        product_list[0].index = 0; 
        product_list[1].type = "h1";
        product_list[1].degree = 2; // Sq2
        product_list[1].index = 0; 
        product_list[2].type = "h2";
        product_list[2].degree = 4; // Sq4
        product_list[2].index = 0;  
    }
}

void AdemAlgebra_free(AdemAlgebra * algebra){
    AdemAlgebra_freeBasis(algebra);
    free(algebra);
}

static void AdemAlgebra__generateBasisGeneric(AdemAlgebraInternal *algebra, int old_max_degree, int max_degree);
static void AdemAlgebra__generateBasis2(AdemAlgebraInternal *algebra, int old_max_degree, int max_degree);
static void AdemAlgebra__generateBasisEven(AdemAlgebraInternal *algebra, int old_max_degree, int max_degree);
static void AdemAlgebra__generateBasisElementToIndexMap(AdemAlgebraInternal *algebra, int old_max_degree, int max_degree);
static void AdemAlgebra__generateMultiplicationTable(AdemAlgebraInternal *algebra, int old_max_degree, int max_degree);
static void AdemAlgebra__generateExcessTable(AdemAlgebraInternal *algebra, int old_max_degree, int max_degree);

void AdemAlgebra_generateBasis(Algebra *this, int max_degree){
    AdemAlgebraInternal *algebra = (AdemAlgebraInternal*) this;
    if(max_degree <= this->max_degree){
        return;
    }
    initializePrime(this->p);
    initializeBocksteinTable();
    uint old_max_degree = this->max_degree;
    if(algebra->public_algebra.generic){
        // generateMultiplcationTableGeneric sometimes goes over by one due to its bockstein logic.
        // rather than testing for this, we take the lazy way out and calculate everything else out one extra step.
        max_degree ++;
        if(old_max_degree > 0){
            old_max_degree ++; // If we've done work before, we also did that one extra step.
        }
    }
    algebra->even_basis_table = realloc(
            algebra->even_basis_table,
            max_degree * sizeof(AdemBasisElement_list)
    );
    this->max_degree = max_degree;
    algebra->basis_element_to_index_map = realloc(
            algebra->basis_element_to_index_map,
            max_degree * sizeof(khash_t(monomial_index_map) *)
        );
    khash_t(monomial_index_map) **name_table = algebra->basis_element_to_index_map;
    for(int k = old_max_degree; k < max_degree; k++){
        name_table[k] = kh_init(monomial_index_map);
    }
    if(algebra->public_algebra.generic){
        AdemAlgebra__generateBasisGeneric(algebra, old_max_degree, max_degree);
    } else {
        AdemAlgebra__generateBasis2(algebra, old_max_degree, max_degree);
    }
    AdemAlgebra__generateBasisElementToIndexMap(algebra, old_max_degree, max_degree);
    if(algebra->public_algebra.generic){
        // AdemAlgebra__generateMultiplicationTable consumes the one extra degree we computed in the generic case
        max_degree --;
        if(old_max_degree > 0){
            old_max_degree --;
        }
    }
    AdemAlgebra__generateMultiplicationTable(algebra, old_max_degree, max_degree);
    if(algebra->public_algebra.generic){
        this->max_degree--;
    }
    if(algebra->public_algebra.unstable){
        AdemAlgebra__generateExcessTable(algebra, old_max_degree, max_degree);
    }
    // Make sure product_list reflects sort order.
    for(uint i = 0; i < algebra->filtrationOneProduct_basisElements.length; i++){
        AdemBasisElement *b = algebra->filtrationOneProduct_basisElements.list[i];
        if(b->degree < max_degree){
            this->product_list.list[i].index = AdemAlgebra_basisElement_toIndex(&algebra->public_algebra, b);
        }
    }
}

static void AdemAlgebra__generateBasis2(AdemAlgebraInternal *algebra, int old_max_degree, int max_degree){
    AdemAlgebra__generateBasisEven(algebra, old_max_degree, max_degree);
    if(algebra->public_algebra.sort_order == NULL){
        algebra->basis_table = algebra->even_basis_table;
        return;
    }
    algebra->basis_table = realloc(algebra->basis_table, max_degree * sizeof(AdemBasisElement_list));
    for(int n=old_max_degree; n<max_degree; n++){
        uint cur_basis_len = algebra->even_basis_table[n].length;
        AdemBasisElement basisElementBuffer[cur_basis_len];
        uint cur_Ps_len = 0;
        for(uint i=0; i<cur_basis_len; i++){
            cur_Ps_len += algebra->even_basis_table[n].list[i]->P_length;
            basisElementBuffer[i] = *algebra->even_basis_table[n].list[i];
        }
        qsort(basisElementBuffer, cur_basis_len, sizeof(AdemBasisElement), algebra->public_algebra.sort_order);

        AdemBasisElement **basis_memory = malloc(
            cur_basis_len * sizeof(AdemBasisElement*)
            + cur_basis_len * sizeof(AdemBasisElement)
            + cur_Ps_len * sizeof(uint)
        );
        char *basis_memory_ptr = (char *)(basis_memory + cur_basis_len);
        for(uint i = 0; i < cur_basis_len; i++){
            AdemBasisElement *basis_elt_target = (AdemBasisElement*)basis_memory_ptr;
            basis_memory[i] = basis_elt_target;
            basis_memory_ptr += sizeof(AdemBasisElement);
            *basis_elt_target = basisElementBuffer[i];
            basis_elt_target->Ps = (uint*)basis_memory_ptr;
            memcpy(basis_memory_ptr, basisElementBuffer[i].Ps, basisElementBuffer[i].P_length * sizeof(uint));
            basis_memory_ptr += basisElementBuffer[i].P_length * sizeof(uint);
        }
        algebra->basis_table[n].list = basis_memory;
        algebra->basis_table[n].length = cur_basis_len;        
    }
}

// This algorithm is sort of clever, but it might be worse than the straightforward way.
// Currently our approach is to pick the bocksteins and the P's separately and merge.
// This approach is what we did for the Milnor basis but it feels worse here.
// 
// The other approach is to use as base case:
// If n % q == 0 add P^(n/q). If n % q == 1 add bP^{n/q} and P^{n/q}b. If n%q == 2, add bP^{n/q}b.
// Then let last go from n/(p+1) - epsilon to 1, find lower elements and append P and Pb to the end.
static void AdemAlgebra__generateBasisGeneric_degreen(AdemAlgebraInternal *algebra, int n);
static void AdemAlgebra__generateBasisGeneric(AdemAlgebraInternal *algebra, int old_max_degree, int max_degree){
    AdemAlgebra__generateBasisEven(algebra, old_max_degree, max_degree);    
    algebra->basis_table = (AdemBasisElement_list *) realloc(
            algebra->basis_table,
            max_degree * sizeof(AdemBasisElement_list)
    );
    // Handle the bocksteins
    for(int n = old_max_degree; n < max_degree; n++){
        AdemAlgebra__generateBasisGeneric_degreen(algebra, n);
    }
}


// Now handle the bocksteins.
// We have our Ps in even_basis_table and they contain in their bockstein field
// a bit flag that indicates where bocksteins are allowed to go.
static void AdemAlgebra__generateBasisGeneric_degreen(AdemAlgebraInternal *algebra, int n){
    uint p = algebra->public_algebra.algebra.p;
    uint q = 2*(p-1);        
    uint cur_basis_len = 0;
    uint cur_Ps_len = 0;
    AdemBasisElement basisElementBuffer[MAX_DIMENSION];
    AdemBasisElement *basis_elt_ptr = basisElementBuffer;
    uint PsBuffer[MAX_DIMENSION]; //* MAX_XI_TAU];
    uint *ps_buffer_ptr = PsBuffer;
    uint residue = n % q;
    // First we need to know how many bocksteins we'll use so we know how much degree
    // to assign to the Ps. The Ps all have degree divisible by q=2p-2, so num_bs needs to
    // be congruent to degree mod q.
    for(uint num_bs = residue; num_bs < MAX_XI_TAU && num_bs <= (uint)n; num_bs += q){
        uint P_deg = (n - num_bs)/q;
        AdemBasisElement_list P_list = algebra->even_basis_table[P_deg];
        for(int i = P_list.length-1; i>=0; i--){
            AdemBasisElement *P = P_list.list[i];
            // We pick our P first.
            if(P->P_length + 1 < num_bs){ // Not enough space to fit the bs.
                continue; // Ps ordered in descending length, so none of the later ones will have space either
            }
            uint bflags_length = bockstein_table[num_bs].length;
            uint *bflags_list = bockstein_table[num_bs].list;
            for(uint j=0; j<bflags_length; j++){
                char bocksteins = bflags_list[j];
                // clz is count leading zeros. We can't have 1's past 1 + P->P_length.
                // clz has weird behavior when the input is zero, so bit or with 1 to prevent that.
                if((uint)(32 - __builtin_clz(bocksteins | 1)) > P->P_length + 1){
                    // Too large of a b. We sorted the Ps in descending length order so we can break now.
                    break;
                }
                // P->bocksteins contains 1 in locations where the sequence is "just barely admissible" and so 
                // adding a bockstein would make it inadmissible.
                if((bocksteins & P->bocksteins) != 0){
                    continue;
                }
                // Okay, everything's good with this bocksteins, P pair so let's add it to our basis.
                // Write new basis element to basis element buffer
                basis_elt_ptr->degree = n;
                basis_elt_ptr->excess = 2*P->excess; // Ps contribute 2 to excess
                basis_elt_ptr->excess += (bocksteins & 1); // leading bockstein increases excess by 1
                basis_elt_ptr->excess -= __builtin_popcount(bocksteins & ((1<<P->P_length) & (~1))); // remaining bocksteins reduce excess by 1
                basis_elt_ptr->bocksteins = bocksteins;
                basis_elt_ptr->P_length = P->P_length;
                basis_elt_ptr->Ps = ps_buffer_ptr;
                memcpy(basis_elt_ptr->Ps, P->Ps, P->P_length * sizeof(uint));
                
                ps_buffer_ptr += basis_elt_ptr->P_length;
                cur_Ps_len += basis_elt_ptr->P_length;
                basis_elt_ptr ++;
                cur_basis_len ++;   
            }
        }
    }

    if(algebra->public_algebra.sort_order != NULL){
        qsort(basisElementBuffer, cur_basis_len, sizeof(AdemBasisElement), algebra->public_algebra.sort_order);
    }
    // Commit to heap
    AdemBasisElement **basis_memory = malloc(
        cur_basis_len * sizeof(AdemBasisElement*)
        + cur_basis_len * sizeof(AdemBasisElement)
        + cur_Ps_len * sizeof(uint)
    );
    char *basis_memory_ptr = (char *)(basis_memory + cur_basis_len);
    for(uint i = 0; i < cur_basis_len; i++){
        AdemBasisElement *basis_elt_target = (AdemBasisElement*)basis_memory_ptr;
        basis_memory[i] = basis_elt_target;
        basis_memory_ptr += sizeof(AdemBasisElement);
        *basis_elt_target = basisElementBuffer[i];
        basis_elt_target->Ps = (uint*)basis_memory_ptr;
        memcpy(basis_memory_ptr, basisElementBuffer[i].Ps, basisElementBuffer[i].P_length * sizeof(uint));
        basis_memory_ptr += basisElementBuffer[i].P_length * sizeof(uint);
    }
    algebra->basis_table[n].list = basis_memory;
    algebra->basis_table[n].length = cur_basis_len;
}

static void AdemAlgebra__generateBasisEven_degreen(AdemAlgebraInternal *algebra, int n);
static void AdemAlgebra__generateBasisEven(AdemAlgebraInternal *algebra, int old_max_degree, int max_degree){
    if(old_max_degree == 0){
        algebra->even_basis_table[0].length = 1;
        algebra->even_basis_table[0].list = malloc(sizeof(AdemBasisElement*) + sizeof(AdemBasisElement));
        algebra->even_basis_table[0].list[0] = (AdemBasisElement*)(algebra->even_basis_table[0].list + 1);
        algebra->even_basis_table[0].list[0]->degree = 0;
        algebra->even_basis_table[0].list[0]->excess = 0;
        algebra->even_basis_table[0].list[0]->bocksteins = 0;
        algebra->even_basis_table[0].list[0]->P_length = 0;
        algebra->even_basis_table[0].list[0]->Ps = NULL;
        old_max_degree++;
    }

    for(int n = old_max_degree; n < max_degree; n ++){
        AdemAlgebra__generateBasisEven_degreen(algebra, n);
    }
}

static void AdemAlgebra__generateBasisEven_degreen(AdemAlgebraInternal *algebra, int n){
    uint p = algebra->public_algebra.algebra.p;
    uint cur_basis_len = 0;
    uint cur_Ps_len = 0;
    
    AdemBasisElement basisElementBuffer[MAX_DIMENSION];
    AdemBasisElement *basis_elt_ptr = basisElementBuffer;
    uint PsBuffer[MAX_DIMENSION];// * MAX_XI_TAU];
    uint *ps_buffer_ptr = PsBuffer;
    // Put Sqn into the list.
    basis_elt_ptr->degree = n;
    basis_elt_ptr->excess = n;
     // bocksteins is a bitmask that has a 1 where bocksteins are illegal.
     // For a single square, the bottom two bocksteins are legal and all others illegal.
    basis_elt_ptr->bocksteins = algebra->public_algebra.generic ? -1 : 0;
    basis_elt_ptr->bocksteins <<= 2; 
    basis_elt_ptr->P_length = 1;
    basis_elt_ptr->Ps = ps_buffer_ptr;
    ps_buffer_ptr[0] = n;
    ps_buffer_ptr += basis_elt_ptr->P_length;
    cur_Ps_len += basis_elt_ptr->P_length;
    basis_elt_ptr ++;
    cur_basis_len ++;

    // last = last term. We append (last,) to the end of
    // elements of degree n - last whose own last square is
    // at least p * last.
    // In order for this to be possible, this means that p last <= n - last, 
    // or (p+1) * last <= n or last <= n/(p+1). We order the squares in decreasing
    // order of their last element so that as we walk over the previous basis
    // when we find a square whose end is too small, we can break.
    for(uint last = n/(p+1); last >= 1; last --){
        uint previous_basis_length = algebra->even_basis_table[n-last].length;
        AdemBasisElement **previous_basis = algebra->even_basis_table[n-last].list;
        for(uint i = 0; i < previous_basis_length; i++){
            AdemBasisElement *prev_elt = previous_basis[i];
            uint old_last_sq = prev_elt->Ps[previous_basis[i]->P_length - 1];
            if(old_last_sq < p * last){
                break;
            }
            // Write new basis element to basis element buffer
            basis_elt_ptr->degree = prev_elt->degree + last;
            basis_elt_ptr->excess = prev_elt->excess - (p-1)*last;
            // We're using bocksteins as a bit mask:
            // A bit in bocksteins shall be set if it's illegal for a bockstein to occur there.
            if(algebra->public_algebra.generic){
                basis_elt_ptr->bocksteins = prev_elt->bocksteins; 
                basis_elt_ptr->bocksteins |= ((old_last_sq == p*last) << prev_elt->P_length);
                basis_elt_ptr->bocksteins &= ~(1 << (prev_elt->P_length +1));
            } else {
                basis_elt_ptr->bocksteins = 0;
            }
            basis_elt_ptr->P_length = prev_elt->P_length + 1;
            basis_elt_ptr->Ps = ps_buffer_ptr;
            memcpy(basis_elt_ptr->Ps, prev_elt->Ps, prev_elt->P_length * sizeof(uint));
            basis_elt_ptr->Ps[basis_elt_ptr->P_length - 1] = last;

            ps_buffer_ptr += basis_elt_ptr->P_length;
            cur_Ps_len += basis_elt_ptr->P_length;
            basis_elt_ptr ++;
            cur_basis_len ++;
        }
    }
    // Okay now we need to move our basis from the stack to the heap.
    AdemBasisElement **basis_memory = malloc(
        cur_basis_len * sizeof(AdemBasisElement*)
        + cur_basis_len * sizeof(AdemBasisElement)
        + cur_Ps_len * sizeof(uint)
    );
    char *basis_memory_ptr = (char *)(basis_memory + cur_basis_len);
    for(uint i = 0; i < cur_basis_len; i++){
        AdemBasisElement *basis_elt_target = (AdemBasisElement*)basis_memory_ptr;
        basis_memory[i] = basis_elt_target;
        basis_memory_ptr += sizeof(AdemBasisElement);
        *basis_elt_target = basisElementBuffer[i];
        basis_elt_target->Ps = (uint*)basis_memory_ptr;
        memcpy(basis_memory_ptr, basisElementBuffer[i].Ps, basisElementBuffer[i].P_length * sizeof(uint));
        basis_memory_ptr += basisElementBuffer[i].P_length * sizeof(uint);
    }
    algebra->even_basis_table[n].list = basis_memory;
    algebra->even_basis_table[n].length = cur_basis_len;
}

static void AdemAlgebra__generateBasisElementToIndexMap(AdemAlgebraInternal *algebra, int old_max_degree, int max_degree){
    khash_t(monomial_index_map) **name_table = algebra->basis_element_to_index_map;
    AdemBasisElement_list *basis_table = algebra->basis_table;
    for(int n = old_max_degree; n < max_degree; n++){
        // We somewhat arbitrarily choose a load factor of 0.5 for the table
        // khash will resize the map if the load factor rises above 0.77.
        // At some point maybe we should profile the hash table and find out what the 
        // best load factor is. We could also try using perfect hash...
        kh_resize(monomial_index_map, name_table[n], basis_table[n].length * 2);
        for(uint i=0; i<basis_table[n].length; i++){
            char key[200];
            AdemAlgebra_basisElement_toKey(key, basis_table[n].list[i]);
            int absent;
            khint_t bin = kh_put(monomial_index_map, name_table[n], key, &absent);
            assert(absent);
            kh_key(name_table[n], bin) = strdup(key);
            kh_val(name_table[n], bin) = i;
        }
    }
}

static void AdemAlgebra__generateMultiplicationTableGeneric(AdemAlgebraInternal *algebra, int old_max_degree, int max_degree);
static void AdemAlgebra__generateMultiplicationTableGeneric_step(AdemAlgebraInternal *algebra, int n, uint x, uint idx);
static void AdemAlgebra__generateMultiplicationTable2(AdemAlgebraInternal *algebra, int old_max_degree, int max_degree);
static void AdemAlgebra__generateMultiplicationTable2_step(AdemAlgebraInternal *algebra, int n, uint x, uint idx);

static void AdemAlgebra__generateMultiplicationTable(AdemAlgebraInternal *algebra, int old_max_degree, int max_degree){
    if(algebra->public_algebra.generic){
        AdemAlgebra__generateMultiplicationTableGeneric(algebra, old_max_degree, max_degree);
    } else {
        AdemAlgebra__generateMultiplicationTable2(algebra, old_max_degree, max_degree);
    }
}

static void AdemAlgebra__generateMultiplicationTable2(AdemAlgebraInternal *algebra, int old_max_degree, int max_degree){
    // degree -> first_square -> admissibile sequence idx -> result vector
    algebra->multiplication_table = realloc(
            algebra->multiplication_table,
            max_degree * sizeof(Vector ***)
        );
    uint p = algebra->public_algebra.algebra.p;
    if(old_max_degree == 0){
        old_max_degree++;
    }
    for(int n=old_max_degree; n < max_degree; n++){
        uint output_dimension = algebra->basis_table[n].length;
        size_t vector_size = Vector_getSize(p, output_dimension, 0);
        uint total_outputs = 0;
        for(int i=1; i<=n; i++){
            total_outputs += AdemAlgebra_getDimension((Algebra*)algebra, n-i, -1);
        }
        size_t table_size = (n+1) * sizeof(Vector **);
        table_size += total_outputs * sizeof(Vector*);
        table_size += total_outputs * vector_size;
        Vector ***top_of_table = (Vector***)malloc(table_size);
        Vector ***current_ptr_1 = top_of_table;
        Vector **current_ptr_2 = (Vector **)(current_ptr_1 + n + 1);
        char *current_ptr_3 = (char*)(current_ptr_2 + total_outputs);
        current_ptr_1 ++; // Skip first_square = 0
        for(int i = 1; i <= n; i++){
            *current_ptr_1 = current_ptr_2;
            uint dimension = AdemAlgebra_getDimension((Algebra*)algebra, n-i, -1);
            for(uint j=0; j < dimension; j++){
                *current_ptr_2 = Vector_initialize(p, &current_ptr_3, output_dimension, 0);
                current_ptr_2++;
            }
            current_ptr_1++;
        }
        assert(current_ptr_1 == top_of_table + n + 1);
        assert(current_ptr_2 == (Vector**)current_ptr_1 + total_outputs);
        assert((char*)current_ptr_3 == (char*)top_of_table + table_size); 
        algebra->multiplication_table[n] = top_of_table;
    }
    for(int n=old_max_degree; n < max_degree; n++){     
        for(uint x = n; x > 0; x--){
            for(uint idx = 0; idx < AdemAlgebra_getDimension((Algebra*)algebra, n-x, -1); idx++){
                AdemAlgebra__generateMultiplicationTable2_step(algebra, n, x, idx);
            }
        }         
    }
}

static void AdemAlgebra__generateMultiplicationTableGeneric(AdemAlgebraInternal *algebra, int old_max_degree, int max_degree){
    // degree -> first_square -> admissibile sequence idx -> result vector
    algebra->multiplication_table = realloc(
            algebra->multiplication_table,
            (max_degree + 1) * sizeof(Vector ***)
        );
    uint p = algebra->public_algebra.algebra.p;
    uint q = 2*p-2;
    if(old_max_degree==0){
        old_max_degree++;
    }
    for(int n=old_max_degree; n < max_degree + 1; n++){
        uint output_dimension = algebra->basis_table[n].length;
        size_t vector_size = Vector_getSize(p, output_dimension, 0);
        uint total_outputs = 0;
        for(uint i=0; i<=n/q; i++){
            for(uint b = 0; b <= 1; b++){
                total_outputs += AdemAlgebra_getDimension((Algebra*)algebra, n-q*i-b, -1);
            }
        }
        uint num_entries = 2*(n/q + 1);
        size_t table_size = num_entries * sizeof(Vector **);
        table_size += total_outputs * sizeof(Vector*);
        table_size += total_outputs * vector_size;
        Vector ***top_of_table = (Vector***)malloc(table_size);
        Vector ***current_ptr_1 = top_of_table;
        Vector **current_ptr_2 = (Vector **)(current_ptr_1 + num_entries);
        char *current_ptr_3 = (char*)(current_ptr_2 + total_outputs);
        for(uint i = 0; i <= n/q; i++){
            for(uint b=0; b<=1; b++){
                *current_ptr_1 = current_ptr_2;
                uint dimension = AdemAlgebra_getDimension((Algebra*)algebra, n - q*i - b, -1);
                for(uint j=0; j < dimension; j++){
                    *current_ptr_2 = Vector_initialize(p, &current_ptr_3, output_dimension, 0);
                    current_ptr_2++;
                }
                current_ptr_1++;
            }
        }
        assert(current_ptr_1 == top_of_table + num_entries);
        assert(current_ptr_2 == (Vector**)current_ptr_1 + total_outputs);
        assert((char*)current_ptr_3 == (char*)top_of_table + table_size); 
        algebra->multiplication_table[n] = top_of_table;
    }
    for(int n=old_max_degree; n < max_degree; n++){    
        for(int x = n/q; x >= 0; x--){
            for(uint idx = 0; idx < AdemAlgebra_getDimension((Algebra*)algebra, n - q*x, -1); idx++){
                AdemAlgebra__generateMultiplicationTableGeneric_step(algebra, n, x, idx);
            }
        }         
    }
}

static void AdemAlgebra__generateMultiplicationTable2_step(AdemAlgebraInternal *algebra, int n, uint x, uint idx){
    Vector *result = algebra->multiplication_table[n][x][idx];
    AdemBasisElement working_elt;
    uint working_elt_Ps[MAX_XI_TAU + 1];
    working_elt.Ps = working_elt_Ps;
    working_elt.degree = n;    
    AdemBasisElement *cur_basis_elt = AdemAlgebra_basisElement_fromIndex((AdemAlgebra*)algebra, n-x, idx);
    // Enough space to fit Sq^x * (rest of Sqs)
    working_elt.bocksteins = 0;
    working_elt.P_length = cur_basis_elt->P_length + 1;
    memcpy(working_elt.Ps + 1, cur_basis_elt->Ps, cur_basis_elt->P_length * sizeof(uint));
    working_elt.Ps[0] = x;      
    // Be careful to deal with the case that cur_basis_elt has length 0            
    // If the length is 0 or the sequence is already admissible, we can just write a 1 in the answer
    // and continue.
    if(cur_basis_elt->P_length == 0 || x >= 2*cur_basis_elt->Ps[0]){
        uint out_idx = AdemAlgebra_basisElement_toIndex((AdemAlgebra*)algebra, &working_elt);
        Vector_addBasisElement(result, out_idx, 1);
        return;
    }
    uint y = working_elt.Ps[1];              
    // We only needed the extra first entry to perform the lookup if our element
    // happened to be admissible. Otherwise, take the rest of the list and forget about it.
    working_elt.P_length --;
    working_elt.Ps ++;
    for(uint j=0; j < 1 + x/2; j++){
        if(Binomial2(y - j - 1, x - 2*j) == 0){
            continue;
        }
        if(j==0){
            working_elt.Ps[0] = x + y;
            // In this case the result is guaranteed to be admissible so we can immediately add it to result
            uint out_idx = AdemAlgebra_basisElement_toIndex((AdemAlgebra*)algebra, &working_elt);
            Vector_addBasisElement(result, out_idx, 1);
            continue;
        }
        AdemBasisElement rest_of_term = working_elt;
        rest_of_term.degree = n - x - y;
        // Take rest of list
        rest_of_term.Ps += 1; 
        rest_of_term.P_length -= 1; 
        // Now we need to reduce Sqj * (rest of Sqs)
        // The answer to this is in the table we're currently making.
        uint rest_of_term_idx = AdemAlgebra_basisElement_toIndex((AdemAlgebra*)algebra, &rest_of_term);
        // total degree -> first sq -> idx of rest of squares
        Vector *rest_reduced = algebra->multiplication_table[j + rest_of_term.degree][j][rest_of_term_idx];
        for(
            VectorIterator it = Vector_getIterator(rest_reduced);
            it.has_more;
            it = Vector_stepIterator(it)
        ){
            if(it.value == 0){
                continue;
            }
            // Reduce Sq^{x+y-j} * whatever square using the table in the same degree, larger index
            // Since we're doing the first squares in decreasing order and x + y - j > x, 
            // we already calculated this.
            Vector_add(result, algebra->multiplication_table[n][x+y-j][it.index], 1);
        }
    }
}

static void AdemAlgebra__generateMultiplicationTableGeneric_step(AdemAlgebraInternal *algebra, int n, uint x, uint idx){
    uint p = algebra->public_algebra.algebra.p;
    uint q = 2*p-2;
    uint x_index = (x<<1);
    Vector *result = algebra->multiplication_table[n][x_index][idx];
    Vector *beta_result = algebra->multiplication_table[n + 1][x_index + 1][idx];
    AdemBasisElement working_elt;
    uint working_elt_Ps[MAX_XI_TAU + 1];
    working_elt.Ps = working_elt_Ps;
    working_elt.degree = n;
    AdemBasisElement *cur_basis_elt = AdemAlgebra_basisElement_fromIndex((AdemAlgebra*)algebra, n-q*x, idx);  
    // Enough space to fit Sq^x * (rest of Sqs)
    uint x_len = x>0 ? 1 : 0;
    working_elt.P_length = cur_basis_elt->P_length + x_len;
    working_elt.bocksteins = cur_basis_elt->bocksteins << x_len;
    memcpy(working_elt.Ps + x_len, cur_basis_elt->Ps, cur_basis_elt->P_length * sizeof(uint));
    if(x>0){
        working_elt.Ps[0] = x;
    }
    // Be careful to deal with the case that cur_basis_elt has length 0            
    // If the length is 0 or the sequence is already admissible, we can just write a 1 in the answer
    // and continue.
    uint b = (cur_basis_elt->bocksteins) & 1;
    if(cur_basis_elt->P_length == 0 || x==0 || x >= p*cur_basis_elt->Ps[0] + b){
        uint out_idx = AdemAlgebra_basisElement_toIndex((AdemAlgebra*)algebra, &working_elt);
        Vector_addBasisElement(result, out_idx, 1);
        if(working_elt.bocksteins & 1){
            // Two bocksteins run into each other (only possible when x=0)
            return;
        }
        working_elt.bocksteins |= 1;
        working_elt.degree ++;
        out_idx = AdemAlgebra_basisElement_toIndex((AdemAlgebra*)algebra, &working_elt);
        Vector_addBasisElement(beta_result, out_idx, 1);
        return;
    }
    uint y = cur_basis_elt->Ps[0];     
    // We only needed the extra first entry to perform the lookup if our element
    // happened to be admissible. Otherwise, take the rest of the list and forget about it.
    working_elt.P_length --;
    working_elt.Ps ++;
    working_elt.bocksteins >>= 1;
    // Adem relation
    for(uint e1 = 0;  e1 <= b; e1++){
        uint e2 = b - e1;
        // e1 and e2 determine where a bockstein shows up.
        // e1 determines if a bockstein shows up in front 
        // e2 determines if a bockstein shows up in middle
        // So our output term looks like b^{e1} P^{x+y-j} b^{e2} P^{j}
        for(uint j=0; j <= x/p; j++){
            uint c = BinomialOdd(p, (y-j) * (p-1) + e1 - 1, x - p*j - e2);
            if(c == 0){
                continue;
            }
            c *= MinusOneToTheN(p, x + j + e2);
            c = c % p;            
            if(j==0){
                AdemBasisElement temp_elt = working_elt;
                if(e2 & (temp_elt.bocksteins >> 1)){
                    // Two bocksteins run into each other:
                    // P^x b P^y b --> P^{x+y} b b = 0
                    continue;
                }                
                temp_elt.Ps[0] = x + y;
                // Mask out bottom bit of original bocksteins.
                temp_elt.bocksteins &= ~1;
                // Now either the front bit or idx + 1 might need to be set depending on e1 and e2.
                temp_elt.bocksteins |= e1;
                temp_elt.bocksteins |= e2 << 1;                
                // In this case the result is guaranteed to be admissible so we can immediately add it to result
                uint out_idx = AdemAlgebra_basisElement_toIndex((AdemAlgebra*)algebra, &temp_elt);
                Vector_addBasisElement(result, out_idx, c);
                if(e1==0){
                    temp_elt.bocksteins |= 1;
                    temp_elt.degree ++;
                    out_idx = AdemAlgebra_basisElement_toIndex((AdemAlgebra*)algebra, &temp_elt);
                    Vector_addBasisElement(beta_result, out_idx, c);
                }
                continue;
            }
            AdemBasisElement rest_of_term = working_elt;
            rest_of_term.degree = n - q*(x + y) - b;
            // Take rest of list
            rest_of_term.Ps ++; 
            rest_of_term.P_length --; 
            rest_of_term.bocksteins >>= 1;
            // Now we need to reduce b^{e2} P^j * (rest of term)
            // The answer to this is in the table we're currently making.
            uint rest_of_term_idx = AdemAlgebra_basisElement_toIndex((AdemAlgebra*)algebra, &rest_of_term);
            uint bj_idx = (j<<1) + e2;
            // (rest of term) has degree n - q*(x + y) - b, 
            // b^{e2} P^j has degree q*j + e2, so the degree of the product is the sum of these two quantities.
            int bj_degree = q*j + e2;
            int bPj_rest_degree =  rest_of_term.degree + bj_degree;
            // total degree ==> b^eP^j ==> rest of term idx ==> Vector
            Vector *rest_of_term_reduced = algebra->multiplication_table[bPj_rest_degree][bj_idx][rest_of_term_idx];
            for(
                VectorIterator it = Vector_getIterator(rest_of_term_reduced);
                it.has_more;
                it = Vector_stepIterator(it)
            ){
                if(it.value == 0){
                    continue;
                }
                // Reduce P^{x+y-j} * whatever square using the table in the same degree, larger index
                // Since we're doing the first squares in decreasing order and x + y - j > x, 
                // we already calculated this.
                bj_idx = ((x+y-j) << 1) + e1;
                Vector *output_vector = algebra->multiplication_table[n][bj_idx][it.index];
                Vector_add(result, output_vector, (c*it.value)%p);
                for(
                    VectorIterator it_output = Vector_getIterator(output_vector);
                    it_output.has_more;
                    it_output = Vector_stepIterator(it_output)
                ){
                    if(it_output.value == 0){
                        continue;
                    }
                    AdemBasisElement *z = AdemAlgebra_basisElement_fromIndex((AdemAlgebra *)algebra, n, it_output.index);
                    if((z->bocksteins & 1) == 0){
                        AdemBasisElement bz = *z;
                        bz.bocksteins |= 1;
                        bz.degree ++;
                        uint idx = AdemAlgebra_basisElement_toIndex((AdemAlgebra *)algebra, &bz);
                        Vector_addBasisElement(beta_result, idx, (it_output.value * c*it.value)%p);
                    }
                }
            }
        }
    }
}

void AdemAlgebra__generateExcessTable(AdemAlgebraInternal *algebra, int old_max_degree, int max_degree){
    algebra->excess_table = realloc(algebra->excess_table, sizeof(uint*)*max_degree);
    for(int n=old_max_degree; n<max_degree; n++){
        uint dim = AdemAlgebra_getDimension((Algebra*)algebra, n, -1);
        algebra->excess_table[n] = malloc(n * sizeof(uint));
        uint cur_excess = 0;
        for(uint i=0; i < dim; i++){
            AdemBasisElement *elt = AdemAlgebra_basisElement_fromIndex((AdemAlgebra*)algebra, n, i);
            for(int j=cur_excess; j<elt->excess; j++){
                algebra->excess_table[n][j] = i;
            }
            cur_excess = elt->excess;
        }
        for(int j=cur_excess; j<n; j++){
            algebra->excess_table[n][j] = dim;
        }        
    }
}

void AdemAlgebra_freeBasis(AdemAlgebra *public_algebra){
    AdemAlgebraInternal *algebra = (AdemAlgebraInternal *) public_algebra;
    for(int n = 0; n < public_algebra->algebra.max_degree; n++){
        free(algebra->basis_table[n].list);
    }
}

uint AdemAlgebra_getDimension(Algebra *this, int degree, int excess __attribute__((unused))){
    assert(degree < (int)this->max_degree);
    if(degree < 0){
        return 0;
    }
    AdemAlgebraInternal *algebra = (AdemAlgebraInternal*) this;
    return algebra->basis_table[degree].length;
}

uint AdemAlgebra_getDimension_unstable(Algebra *this, int degree, int excess){
    assert(degree < this->max_degree);
    if(degree < 0){
        return 0;
    }
    AdemAlgebraInternal *algebra = (AdemAlgebraInternal*) this;
    if(excess >= degree){
        return algebra->basis_table[degree].length;
    }
    return algebra->excess_table[degree][excess];
}

AdemBasisElement_list AdemAlgebra_getBasis(AdemAlgebra *algebra, int degree){
    assert(degree < algebra->algebra.max_degree);
    return ((AdemAlgebraInternal*)algebra)->basis_table[degree];
}

AdemBasisElement *AdemAlgebra_basisElement_fromIndex(AdemAlgebra *public_algebra, int degree, uint index){
    AdemAlgebraInternal *algebra = (AdemAlgebraInternal *)public_algebra;
    assert(degree < public_algebra->algebra.max_degree);
    assert(index < algebra->basis_table[degree].length);
    AdemBasisElement *result = algebra->basis_table[degree].list[index];
    assert(AdemAlgebra_basisElement_toIndex(public_algebra, result) == index);
    return result;
}

uint AdemAlgebra_basisElement_toIndex(AdemAlgebra *public_algebra,  AdemBasisElement *b){
    assert(b->degree < public_algebra->algebra.max_degree);
    AdemAlgebraInternal *algebra = (AdemAlgebraInternal *)public_algebra;
    kh_monomial_index_map_t *map = algebra->basis_element_to_index_map[b->degree];
    char key[200];
    AdemAlgebra_basisElement_toKey(key, b);
    khint_t bin = kh_get(monomial_index_map, map, key);
    if(bin == kh_end(map)){
        AdemAlgebra_basisElement_toString(key, public_algebra, b);
        printf("Uh-oh, not here. degree: %d, elt: '%s'\n", b->degree, key);
        assert(false);
        return 0;
    }
    uint result = kh_val(map, bin);
    return result;
}


static void AdemAlgebra__makeMonoAdmissible2(AdemAlgebra *algebra, Vector *result, AdemBasisElement *monomial, int idx, int leading_degree, int excess, bool stop_early);
static void AdemAlgebra__makeMonoAdmissibleGeneric(AdemAlgebra *algebra, Vector *result, uint coeff, AdemBasisElement *monomial, int idx, int leading_degree, int excess, bool stop_early);

void AdemAlgebra_multiply(Algebra *this, Vector *result, uint coeff, 
                        int r_degree, uint r_index, 
                        int s_degree, uint s_index, int excess)
{
    if(coeff == 0){
        return;
    }
    AdemAlgebra *algebra = (AdemAlgebra*)this;
    assert(r_index < Algebra_getDimension(this, r_degree, excess + s_degree));
    assert(s_index < Algebra_getDimension(this, s_degree, excess));

    if(s_degree == 0){
        // If s is of length 0 then max_idx "r->P_length" is off the edge of the list and it segfaults.
        // Avoid this by returning early in this case.
        Vector_addBasisElement(result, r_index, coeff);
        return;
    }
    AdemBasisElement *r = AdemAlgebra_basisElement_fromIndex(algebra, r_degree, r_index);
    AdemBasisElement *s = AdemAlgebra_basisElement_fromIndex(algebra, s_degree, s_index);
    AdemBasisElement monomial;
    monomial.P_length = r->P_length + s->P_length;
    monomial.degree = r->degree + s->degree;

    monomial.bocksteins = 0;
    if(algebra->generic
        && (r->bocksteins >> (r->P_length)) & (s->bocksteins & 1)){
        // If there is a bockstein at the end of r and one at the beginning of s, these run into each other
        // and the output is 0.
        return;
    } else if(algebra->generic){
        monomial.bocksteins = r->bocksteins;
        monomial.bocksteins |= s->bocksteins << (r->P_length);
    }

    // if(algebra->generic && s->P_length == 0){ 
    //     // In this case s is just a bockstein. This causes the same trouble as the 
    //     // s is the identity case we already covered (because s->P_length == 0), 
    //     // so we just handle it separately.
    //     monomial.Ps = r->Ps;
    //     uint idx = AdemAlgebra_basisElement_toIndex(algebra, &monomial);
    //     Vector_addBasisElement(result, idx, coeff);
    //     return;
    // }
    uint memory[monomial.P_length];
    memcpy(memory, r->Ps, r->P_length * sizeof(uint));
    memcpy(memory + r->P_length, s->Ps, s->P_length * sizeof(uint));
    monomial.Ps = memory;
    if(algebra->generic){
        // If r ends in a bockstein, we need to move it over because we consider
        // the monomial from right to left in chunks like bP^i. The b from the end of r gets donated
        // to the P from the beginning of s.
        int leading_degree = r->degree - ((r->bocksteins >> r->P_length) & 1);
        AdemAlgebra__makeMonoAdmissibleGeneric(algebra, result, coeff, &monomial, r->P_length - 1, leading_degree, excess, true);
    } else {
        AdemAlgebra__makeMonoAdmissible2(algebra, result, &monomial, r->P_length - 1, r->degree, excess, true);
    }
}

void AdemAlgebra_makeMonoAdmissible(AdemAlgebra *algebra, Vector *result, uint coeff, AdemBasisElement *monomial, int excess){
    assert(monomial->degree < algebra->algebra.max_degree);
    uint q = algebra->generic ? 2*algebra->algebra.p-2 : 1;
    int leading_degree = monomial->degree - q * monomial->Ps[monomial->P_length - 1];
    int idx = monomial->P_length - 2;    
    if(algebra->generic){
        leading_degree -= (monomial->bocksteins >> (monomial->P_length - 1)) & 1;
        AdemAlgebra__makeMonoAdmissibleGeneric(algebra, result, coeff, monomial, idx, leading_degree, excess, false);
    } else {
        AdemAlgebra__makeMonoAdmissible2(algebra, result, monomial, idx, leading_degree, excess, false);
    }
}

#define MAX(a, b) ((a)>(b)) ? (a) : (b)
#define MIN(a, b) ((a)<(b)) ? (a) : (b)
/**
 * Reduce a Steenrod monomial at the prime 2.
 * Arguments:
 *    algebra -- an Adem algebra. This would be a method of class AdemAlgebra.
 *    result  -- Where we put the result
 *    monomial -- a not necessarily admissible Steenrod monomial which we will reduce. 
 *                We destroy monomial->Ps.
 *    idx -- the only index to check for inadmissibility in the input (we assume that we've gotten
 *           our input as a product of two admissible sequences.)
 *    leading_degree -- the degree of the squares between 0 and idx (so of length idx + 1)
 */
static void AdemAlgebra__makeMonoAdmissible2(
    AdemAlgebra *public_algebra, Vector *result, AdemBasisElement *monomial,
    int idx, int leading_degree, int excess, bool stop_early
){
    AdemAlgebraInternal *algebra = (AdemAlgebraInternal*)public_algebra;
    // Check for admissibility;
    while(idx < 0 || (uint)idx == monomial->P_length - 1 || monomial->Ps[idx] >= 2*monomial->Ps[idx + 1]){
        if(idx < 0 || stop_early){
            // Admissible so write monomial to result.
            uint idx = AdemAlgebra_basisElement_toIndex(public_algebra, monomial);
            // If excess is too large, quit. It's faster to check this by comparing idx to dimension
            // than to use fromIndex because fromIndex  dereferences a hash map.
            if(public_algebra->unstable && idx >= Algebra_getDimension((Algebra*)algebra, monomial->degree, excess)){
                return;
            }
            Vector_addBasisElement(result, idx, 1);
            return;
        }
        leading_degree -= monomial->Ps[idx];
        idx --;
    }
    AdemBasisElement tail_of_monomial = *monomial;
    tail_of_monomial.P_length -= idx + 1;
    tail_of_monomial.Ps += idx + 1;
    tail_of_monomial.degree -= leading_degree;
    uint x = monomial->Ps[idx];
    uint adm_idx = AdemAlgebra_basisElement_toIndex(public_algebra, &tail_of_monomial);
    int tail_degree = tail_of_monomial.degree + x;
    Vector *reduced_tail = algebra->multiplication_table[tail_degree][x][adm_idx];
    uint dim = Algebra_getDimension((Algebra*)algebra, tail_degree, excess);
    for(
        VectorIterator it = Vector_getIterator(reduced_tail);
        it.index < dim;
        it = Vector_stepIterator(it)
    ){
        if(it.value == 0){
            continue;
        }
        AdemBasisElement *cur_tail_basis_elt = AdemAlgebra_basisElement_fromIndex(public_algebra, tail_degree, it.index);
        tail_of_monomial.degree = monomial->degree;
        tail_of_monomial.P_length = idx + cur_tail_basis_elt->P_length;
        uint temp_Ps[tail_of_monomial.P_length];
        tail_of_monomial.Ps = temp_Ps;
        memcpy(tail_of_monomial.Ps, monomial->Ps, idx * sizeof(*tail_of_monomial.Ps));
        memcpy(tail_of_monomial.Ps + idx, cur_tail_basis_elt->Ps, cur_tail_basis_elt->P_length * sizeof(*cur_tail_basis_elt->Ps));        
        AdemAlgebra__makeMonoAdmissible2(public_algebra, result, &tail_of_monomial, idx - 1, leading_degree - x, excess, stop_early);
    }
}


static void AdemAlgebra__makeMonoAdmissibleGeneric(
    AdemAlgebra *public_algebra, Vector *result, uint coeff, AdemBasisElement *monomial,
    int idx, int leading_degree, int excess, bool stop_early
    ){
    AdemAlgebraInternal *algebra = (AdemAlgebraInternal*)public_algebra;    
    // array_print("mono: (%s,", monomial->Ps, monomial->P_length);
    // printf("%x) coeff: %d\n", monomial->bocksteins, coeff);
    uint p = public_algebra->algebra.p;
    uint q = 2*p-2;
    // Check for admissibility    
    uint b1 = (monomial->bocksteins >> idx) & 1;
    uint b2 = (monomial->bocksteins >> (idx+1)) & 1;
    while(idx < 0 || (uint)idx == monomial->P_length - 1 || monomial->Ps[idx] >= p*monomial->Ps[idx + 1] + b2){
        if(idx < 0 || stop_early){
            // Admissible so write monomial to result.
            uint idx = AdemAlgebra_basisElement_toIndex(public_algebra, monomial);
            if(public_algebra->unstable && idx >= Algebra_getDimension((Algebra*)algebra, monomial->degree, excess)){
                return;
            }        
            Vector_addBasisElement(result, idx, coeff);
            return;
        }
        leading_degree -= q * monomial->Ps[idx];
        leading_degree -= (monomial->bocksteins >> idx) & 1;
        // monomial->degree += monomial->Ps[idx];
        // monomial->degree += (monomial->bocksteins >> idx) & 1;
        idx --;
    }    
    AdemBasisElement tail_of_monomial = *monomial;
    tail_of_monomial.P_length -= idx + 1;
    tail_of_monomial.Ps += idx + 1;
    tail_of_monomial.bocksteins >>= idx + 1;
    tail_of_monomial.degree -= leading_degree;
    // Notice how much we avoid bockstein twiddling here. It's all hidden in multiplication_table =)
    uint x = monomial->Ps[idx];
    uint bx = (x << 1) + b1;
    uint adm_idx = AdemAlgebra_basisElement_toIndex(public_algebra, &tail_of_monomial);
    int tail_degree = tail_of_monomial.degree + q*x + b1;
    Vector *reduced_tail = algebra->multiplication_table[tail_degree][bx][adm_idx];
    uint dim = Algebra_getDimension((Algebra*)algebra, tail_degree, excess);    
    for(
        VectorIterator it = Vector_getIterator(reduced_tail);
        it.index < dim;
        it = Vector_stepIterator(it)
    ){
        if(it.value == 0){
            continue;
        }
        AdemBasisElement *cur_tail_basis_elt = AdemAlgebra_basisElement_fromIndex(public_algebra, tail_degree, it.index);
        tail_of_monomial.degree = monomial->degree;
        tail_of_monomial.P_length = idx + cur_tail_basis_elt->P_length;
        uint temp_Ps[tail_of_monomial.P_length];
        tail_of_monomial.Ps = temp_Ps;
        memcpy(tail_of_monomial.Ps, monomial->Ps, idx * sizeof(*tail_of_monomial.Ps));
        memcpy(tail_of_monomial.Ps + idx, cur_tail_basis_elt->Ps, cur_tail_basis_elt->P_length * sizeof(*cur_tail_basis_elt->Ps));        
        tail_of_monomial.bocksteins = monomial->bocksteins & ((1<<idx)-1);
        tail_of_monomial.bocksteins |= cur_tail_basis_elt->bocksteins << idx;
        int new_leading_degree = leading_degree - q*x - b1;
        AdemAlgebra__makeMonoAdmissibleGeneric(public_algebra, result, (coeff * it.value)%p, &tail_of_monomial, idx - 1, new_leading_degree, excess, stop_early);
    }
}

// This is for the javascript bindings to save us the trouble of computing offsets into
// AdemBasisElements.
AdemBasisElement *AdemAlgebra_basisElement_construct(uint degree, uint P_length, uint *Ps, uint bocksteins){
    AdemBasisElement *result = malloc(sizeof(AdemBasisElement));
    result->degree = degree;
    result->P_length = P_length;
    result->Ps = Ps;
    result->bocksteins = bocksteins;
    return result;
}

void AdemAlgebra_basisElement_free(AdemBasisElement *b){
    free(b);
}


uint AdemAlgebra_basisElement_getPlength(AdemBasisElement *b){
    return b->P_length;
}

uint *AdemAlgebra_basisElement_getPs(AdemBasisElement *b){
    return b->Ps;
}

uint AdemAlgebra_basisElement_getBocksteins(AdemBasisElement *b){
    return b->bocksteins;
}


// for e1 in range(1 + bockstein):
//     e2 = bockstein - e1
//     for j in range(1 + (A - e1)//p):
//         coeff = combinatorics.binomial_odd((B-j) * (p-1) - 1 + e1, A - p*j - e2, p)
//         coeff *= (-1)**(A+j + e2)
//         coeff = coeff % p
//         if coeff != 0 and j == 0:
//             result[(e1, A+B, e2)] = coeff
//         elif coeff != 0 and j != 0:
//             result[(e1, A+B-j, e2, j, 0)] = coeff 

/**
int main(){
    char buffer[5000];
    uint p = 2;
    bool generic = p != 2;
    generic = true;
    uint max_degree = 40;    
    AdemAlgebra *algebra = AdemAlgebra_construct(p, generic, false);
    algebra->sort_order = AdemAlgebra_excessSortOrder;
    AdemAlgebra_generateBasis((Algebra*)algebra, max_degree);
    AdemAlgebra_generateBasis((Algebra*)algebra, 45);
    // for(uint i=0; i < max_degree; i++){
    //     AdemBasisElement_list basis = AdemAlgebra_getBasis(algebra, i);
    //     printf("degree %d:\n", i);
    //     for(uint j=0; j<basis.length; j++){
    //         AdemAlgebra_basisElement_print("   %s", algebra, basis.list[j]);
    //         printf(" excess: %d\n", basis.list[j]->excess);
    //     }
    //     printf("\n\n");
    // }
    // return 0;
    
    AdemBasisElement *A, *B;
    A = AdemAlgebra_basisElement_fromString(algebra, "P2");
    // A = AdemAlgebra_basisElement_fromString(algebra, "P13");
    // AdemBasisElement unit = (AdemBasisElement){0};
    // A = &unit;
    B = AdemAlgebra_basisElement_fromString(algebra, "P18");
    // B = AdemAlgebra_basisElement_fromString(algebra, "P6 b P1");
    uint A_idx = AdemAlgebra_basisElement_toIndex(algebra, A);
    uint B_idx = AdemAlgebra_basisElement_toIndex(algebra, B);
    

    uint output_degree = A->degree + B->degree;
    uint output_dimension = AdemAlgebra_getDimension((Algebra*)algebra, output_degree, 0);
    Vector *result = Vector_construct(algebra->algebra.p, output_dimension, 0);
    AdemAlgebra_multiply((Algebra*)algebra, result, 1, A->degree, A_idx, B->degree, B_idx, 0);
    char *str1 = buffer;
    int len1 = AdemAlgebra_element_toString(str1, algebra, output_degree, result);
    char *str2 = str1 + len1 + 1;
    int len2 = AdemAlgebra_basisElement_toString(str2, algebra,  A);
    char *str3 = str2 + len2 + 1;
    AdemAlgebra_basisElement_toString(str3, algebra, B);
    printf("%s * %s = %s\n", str2, str3, str1);
    Vector_free(result);
    return 0;
}
//*/