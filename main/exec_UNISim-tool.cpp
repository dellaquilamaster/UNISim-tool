/*
 * UNISim-tool Framework Main Program
 * Created by: Daniele Dell'Aquila (daniele.dellaquila@irb.hr)
 * 10 June 2019
 * 
 * 
 */

#include <stdio.h>
#include <stdlib.h>

#include <UNISFramework.h>
#include <UNISLogo.h>

int main (int argc, char ** argv)
{ 
  //
  //Printing logo
  PrintLogo();
  //
  
  //
  //Creation of the framework
  UNISFramework TheFramework;
  //
  
  //
  //Reading input
  if(TheFramework.ReadInput(argc,argv)<1) {
    printf("Error: invalid input!\nPlease specify number of events.\nAborting.\n\n");
    exit(2);
  }
  //
  
  //
  //Configuring the framework
  if(TheFramework.ConfigureFramework()<=0) {
    printf("Error: failed to configure the framework. Check configuration file. Aborting\n\n");
    exit(1);
  }
  //
  
  //
  //Printing configuration
  TheFramework.PrintConfiguration();
  //

  //
  //Output Tree Initialization
  TheFramework.InitTree();
  //
  
  //
  //Running the framework
  TheFramework.ProcessIterations();
  //
  
  //
  //Terminating the process
  TheFramework.EndProcess();
  //
  
  return 0; 
}
