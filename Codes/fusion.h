#include<iostream>
#include"PARAMETER.h"
#include"cuLayer.h"
void usage(char *command);
int parseParameters(char *fname, CuLayer *psensor,PARAMETER *par);
void Re_fusion(CuLayer *psensor,PARAMETER *par);