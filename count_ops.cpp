#include "count_ops.h"

long long int op_counter = 0;

void inc_op(int n){
  op_counter += n;
  if (op_counter < 0){
    printf("op_counter overflow\n");
    exit(0);
  }
}

