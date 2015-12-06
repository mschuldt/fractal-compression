#ifndef COUNT_OPS_H
#define COUNT_OPS_H

#define COUNT_OPS 0

#include <cstdio>
#include <cstdlib>

extern long long int op_counter;

void inc_op(int n);

#if COUNT_OPS
#define INC_OP(n) inc_op(n)
#else
#define INC_OP(n)
#endif


#endif
