#include <sys/time.h>
#include <time.h>
#include <cstdio> 

int main(int argc, char *argv[]){
  struct timeval tv;
  gettimeofday(&tv, 0);
  double time = tv.tv_sec +1e-6 * tv.tv_usec;
  printf("Current time = %f seconds. \n", time);

  int sum = 0;
  for(int i = 0; i <10; i++){
  	sum+=i;
  }
  gettimeofday(&tv, 0);
  double end = tv.tv_sec +1e-6 * tv.tv_usec;
  printf("END time = %f seconds. \n", end);


  printf("total execution time = %f seconds. \n", end-time);


}