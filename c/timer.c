#include <time.h>
#include <stdio.h>

void start_timing(clock_t *time)
{
	*time=clock();
}

void stop_timing(clock_t *time, const char*string)
{
  	printf("Time taken for %s =%lg s\n", string, (double) (clock()-*time)/((double) CLOCKS_PER_SEC)); 	
}
