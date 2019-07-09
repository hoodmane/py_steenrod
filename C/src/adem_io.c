#include <stdio.h>
#include <string.h>
#include "adem.h"


uint AdemAlgebra__generateName(AdemAlgebra *algebra){
    char buffer[1000];
    uint len = 0;
    len += sprintf(buffer + len, "AdemAlgebra(p=%d", algebra->algebra.p);
    if(algebra->generic != (algebra->algebra.p != 2)){
        len += sprintf(buffer + len, ", generic=%s", algebra->generic ? "true" : "false");
    }
    len += sprintf(buffer + len, ")");
    char *result = malloc((len + 1) * sizeof(char));
    memcpy(result, buffer, (len + 1) * sizeof(char));
    algebra->algebra.name = result;
    return len;
}

uint AdemAlgebra_basisElement_toKey(char *buffer, AdemBasisElement *b){
    // Note that we ignore the degree and excess fields.
    // These are calculated fields, and we don't want to force calculating code to
    // always recalculate them.
    uint len = 0;
    memcpy(buffer, &b->bocksteins, sizeof(uint));
    len += sizeof(uint);
    memcpy(buffer + len, b->Ps, b->P_length * sizeof(uint));
    len += b->P_length * sizeof(uint);
    // Ensure that none of the chars in our key are 0.
    // A zero would terminate our key early.
    for(uint i = 0; i < len; i++){
        buffer[i] = (buffer[i] << 1) | 1;
    }
    buffer[len] = '\0';
    return len;
}

uint AdemAlgebra_basisElementIndex_toString(Algebra *this, char *buffer, int degree, uint idx){
    AdemAlgebra *algebra = (AdemAlgebra*)this;
    AdemBasisElement *b = AdemAlgebra_basisElement_fromIndex(algebra, degree, idx);
    return AdemAlgebra_basisElement_toString(buffer, algebra, b);
}

uint AdemAlgebra_basisElement_toString(char *buffer, AdemAlgebra *algebra, AdemBasisElement *b){
    printf("AdemAlgebra_basisElement_toString\n");
    uint len = 0;
    bool generic = algebra->generic;
    char P_or_Sq[4];
    if(generic){
        sprintf(P_or_Sq, "P");
    } else {
        sprintf(P_or_Sq, "Sq");
    }
    uint bockstein;
    for(uint i = 0; i < b->P_length; i++){
        bockstein = (b->bocksteins >> i) & 1;
        if(bockstein && bockstein){
            len += sprintf(buffer + len, "b");
        }
        len += sprintf(buffer + len, "%s%d ", P_or_Sq, b->Ps[i]);
    }
    bockstein = (b->bocksteins >> b->P_length) & 1;
    if(bockstein && bockstein){
        len += sprintf(buffer + len, "b");
    } else if(b->P_length == 0){
        buffer[0] = '1';
        buffer[1] = '\0';
        len = 1;
    } else { // delete trailing space.
        len --;
        buffer[len] = '\0';
    }
    return len;
}

void AdemAlgebra_basisElement_print(char *fmt_string, AdemAlgebra *algebra, AdemBasisElement *b){
    char buffer[2000];
    AdemAlgebra_basisElement_toString(buffer, algebra, b);
    printf(fmt_string, buffer);
}

// Parses squares like: "P3 b P1" or "P4 P2 P1".
AdemBasisElement *AdemAlgebra_basisElement_fromString(AdemAlgebra *algebra, char *elt_string){
    char *string, *tofree;
    tofree = string = strdup(elt_string);
    char * Op;
    uint q = 1;
    if(algebra->generic){
        q = 2*algebra->algebra.p - 2;
    }    
    uint P_part[20] = {0};
    uint bocksteins = 0;
    uint len = 0;
    uint degree = 0;
    while ((Op = strsep(&string, " ")) != NULL){
        uint num_position = 1;
        if(Op[0] == 'b'){
            bocksteins |= 1 << len;
            num_position ++;
            degree ++;
        }
        if(Op[1] == '\0'){ // In this case we ended with a b.
            continue; 
        }
        if(Op[0]=='S'){//S for Sq. We need to skip the q too.
            num_position++;
        }
        char *end;
        P_part[len] = strtoul(Op + num_position, &end, 10);
        degree += q*P_part[len];
        len++;
    }
    free(tofree);
    tofree = NULL;
    string = NULL;

    // TODO: instead of mallocing a rogue element, use the copy in the master table.
    AdemBasisElement ABE;
    ABE.degree = degree;
    ABE.bocksteins = bocksteins;
    ABE.P_length = len;
    ABE.Ps = P_part;
    uint idx = AdemAlgebra_basisElement_toIndex(algebra, &ABE);
    AdemBasisElement *result = AdemAlgebra_basisElement_fromIndex(algebra, degree, idx);
    return result;
}


uint AdemAlgebra_element_toString(Algebra *this, char *buffer, int degree, Vector *m){
    AdemAlgebra *algebra = (AdemAlgebra *)this;
    uint len = 0;
    for(
        VectorIterator it = Vector_getIterator(m); 
        it.has_more; 
        it = Vector_stepIterator(it)
    ){
        if(it.value == 0){
            continue;
        }
        if(it.value != 1) {
            len += sprintf(buffer + len, "%d * ", it.value);
        }
        AdemBasisElement *b = AdemAlgebra_basisElement_fromIndex(algebra, degree, it.index);
        len += AdemAlgebra_basisElement_toString(buffer + len, algebra, b);
        len += sprintf(buffer + len, " + ");
    }
    if(len == 0){
        len += sprintf(buffer + len, "0");
    } else {
        // Remove trailing " + "
        len -= 3;
        buffer[len] = '\0';
    }
    return len;
}

void AdemAlgebra_element_print(char *fmt_string, AdemAlgebra *algebra, int degree, Vector *m){
    char buffer[1000];
    AdemAlgebra_element_toString((Algebra *)algebra, buffer, degree, m);
    printf(fmt_string, buffer);
}