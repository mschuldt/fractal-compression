#ifndef COUNT_OPS_H
#define COUNT_OPS_H

#define COUNT_OPS false

#include <cstdio>
#include <cstdlib>


#define downsample_i 1
#define execute_i 2
#define getscalefactore_i 3
#define geterror_i 4
#define getaveragepixel_i 5
#define last_index 6

extern long long int op_counter[last_index];

void inc_op(int n, int i);
long long int total_ops(void);

#if COUNT_OPS
#define INC_OP(n) inc_op(n, 0)
#define INC_OP2(n, i) inc_op(n, i)
#else
#define INC_OP(n)
#define INC_OP2(n, i)
#endif


#endif
