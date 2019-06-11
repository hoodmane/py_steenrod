
#ifndef CSTEENROD_ADEM_H
#define CSTEENROD_ADEM_H

#include "algebra.h"


typedef struct {
    Algebra algebra;
    bool generic;
    bool unstable;
    // This will be passed to q_sort and determines our monomial ordering.
    // return <0 if a should go before b, >0 if a should come after b, and 
    // 0 if you don't care. It's probably a good idea to impose a total order
    // for increased consistency.
    int (*sort_order)(const void *a, const void *b); 
} AdemAlgebra;

typedef struct {
    uint degree;
    uint excess;
    uint bocksteins;
    uint P_length;
    uint *Ps;
} AdemBasisElement;

typedef struct {
    uint length;
    AdemBasisElement **list;
} AdemBasisElement_list;

uint AdemAlgebra_element_toString(char *buffer, AdemAlgebra *algebra, uint degree, Vector *m);
void AdemAlgebra_element_print(char *fmt_string, AdemAlgebra *algebra, uint degree, Vector *m);

uint AdemAlgebra_basisElement_toKey(char *buffer, AdemBasisElement *b);
uint AdemAlgebra_basisElement_toString(char *buffer, AdemAlgebra *algebra, AdemBasisElement *b);
void AdemAlgebra_basisElement_print(char *fmt_string, AdemAlgebra *algebra, AdemBasisElement *b);
AdemBasisElement *AdemAlgebra_basisElement_fromString(AdemAlgebra *algebra, char *elt_string);

AdemAlgebra *AdemAlgebra_construct(uint p, bool generic, bool unstable);
void AdemAlgebra_free(AdemAlgebra *algebra);

void AdemAlgebra_generateBasis(Algebra *this, uint max_degree);
void AdemAlgebra_freeBasis(AdemAlgebra *algebra);

uint AdemAlgebra_getDimension(Algebra *this, int degree, uint excess);
AdemBasisElement_list AdemAlgebra_getBasis(AdemAlgebra *algebra, uint degree);
AdemBasisElement *AdemAlgebra_basisElement_fromIndex(AdemAlgebra *algebra, uint degree, uint index);
uint AdemAlgebra_basisElement_toIndex(AdemAlgebra *algebra,  AdemBasisElement *b);

void AdemAlgebra_multiply(Algebra *this, Vector *result, uint coeff, uint r_degree, uint r_index, uint s_degree, uint s_index, uint excess);


#endif //CSTEENROD_ADEM_H
