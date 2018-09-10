#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include "Seed.h"
void seed(void)
{
	int fd, buf;
	 
	if((fd = open("/dev/random", O_RDONLY)) < 0) {
	        perror("/dev/random");
	        //return(1);
	 }
	 read(fd, &buf, sizeof(buf));
	 close(fd);
	 srand(buf);
}
