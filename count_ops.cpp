#include "count_ops.h"

long long int op_counter[last_index] = {0};

void inc_op(int n){
  op_counter[0] += n;
  if (op_counter[0] < 0){
    printf("op_counter overflow\n");
    exit(0);
  }
}

void inc_op(int n, int i){
  if (i >= last_index){
    printf("op_counter(%d, %d): invalid index\n", n, i);
    exit(0);
  }
  op_counter[i] += n;
  if (op_counter[i] < 0){
    printf("op_counter[%d] overflow\n", i);
    exit(0);
  }
}

long long int total_ops(){
  long long int sum = 0;
  for (int i=0; i < last_index;i++){
    sum += op_counter[i];
  }
  return sum;
}

