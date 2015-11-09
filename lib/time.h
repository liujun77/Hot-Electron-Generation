#ifndef jga_time          // ************************************************
#define jga_time 1        // ***  jga/time.h                     20-v-99  ***
                          // ***                                          ***
#include "jga.h"          // ***                                          ***
                          // ***                                          ***
// ***  ----------------  // ***  --------------------------------------  ***
// ***                                                                    ***
// ***  Author: F. J. Garcia de Abajo, CSIC, Spain                        ***
// ***                                                                    ***
// ***  ----------------------------------------------------------------  ***
// ***                                                                    ***
// ***  Time tracking                                                     ***
// ***                                                                    ***
// ***  This permits to print the total CPU time since the program        ***
// ***  was invoked, and the CPU since the last time it was printed.      ***
// ***  Notice that time_update() has to be called regularly to avoid     ***
// ***  loss of time track, since clock cannot count indefinitely.        ***
// ***                                                                    ***
// ***  ----------------------------------------------------------------  ***
// ***                                                                    ***
// ***  Routines:  time_print(fout);  print (update) time to file fout    ***
// ***             time_print();      print time to standard output       ***
// ***             time_update();     update time counter to prevent      ***
// ***                                loss of track                       ***
// ***                                                                    ***
// **************************************************************************


#define time_seconds  0.000001  // =0.000001 for SUN C++;  =0.01 for GNU g++
numero  time_initial=clock(),   // last output of clock()
        time_partial=0,         // CPU time since last print out of time
        time_total=0;           // total CPU time

// --------------------------------------------------------------------------

void time_update(void)
{
  numero time_temp=clock();
  time_partial+=time_temp-time_initial;
  time_initial=time_temp;
}

// --------------------------------------------------------------------------

void time_print(FILE *fout)
{
  time_update();  time_total+=time_partial;

  fprintf(fout, "--- Time:  total=%4.2f sec, partial= %4.2f sec.\n",
          time_total*time_seconds,
          time_partial*time_seconds);

  time_partial=0;
}

// --------------------------------------------------------------------------

void time_print(void) {time_print(foutput);}

#endif  // ******************************************************************
