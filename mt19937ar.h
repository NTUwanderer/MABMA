/********************************
  Author: Tian-Li Yu
  E-mail: tianliyu@nextumi.com
********************************/
#ifndef _MT199937AR_H
#define _MT199937AR_H

extern void init_by_array (unsigned long init_key[], int key_length);
                                 //[0,1]
extern double genrand_real1 (void);
                                 //[0,1)
extern double genrand_real2 (void);
                                 //(0,1)
extern double genrand_real3 (void);
                                 //[0,1)
extern double genrand_res53 (void);

#endif
