#include<stdio.h>
#include"fusion.h"
void usage(char *command) {
  printf("Usage: %s <IN_PARAMETER_FILE>\n", command);
  printf("<IN_PARAMETER_FILE> is a text input file which contains:\n");
  printf("STNLFFM_PARAMETER_START\n");
  printf("\tNUM_IN_PAIRS = \n");
 // printf("\tThe_pf_band_of_Landsat_for_calculating = \n");
//  printf("\tThe_pc_band_of_MODIS_for_calculating = \n");
  printf("\tIN_PAIR_MODIS_FNAME = \n");
  printf("\tIN_PAIR_LANDSAT_FNAME = \n");
  printf("\tIN_PDAY_MODIS_FNAME = \n");
  printf("\tOUT_PDAY_LANDSAT_FNAME = \n");
  printf("\tis_limit_a_CalcuRela\n");
  printf("\tThe_width_of_searching_window=\n");
   printf("\tFilter_parameters=\n");
  printf("\tG_Type =(optional)\n");
 // printf("\t =\n");
 // printf("\tLandsat_sensor_error= \n");
 // printf("\tMODIS_sensor_error = \n"); 
  printf("STNLFFM_PARAMETER_END\n");
}
int parseParameters(char *fname, CuLayer *psensor,PARAMETER *par)
{
	  int   i, k, total_pairs=0;
  char  buffer[4000] = "\0";
  char  *label =NULL;
  par->NUM_PREDICTIONS=0;
  char  *tokenptr =NULL;
  char readpath[4000];
  string argName;
  string t;
  char  *separator = "= ,";
  FILE  *in;  
	//char  buffer[1000] = "\0";
	if((in=fopen(fname,"r"))==NULL) 
	{
    printf("Can't open input %s\n", fname);
    return -1;
    }
	fscanf(in, "%s", buffer);
	if(strcmp(buffer, "STNLFFM_PARAMETER_START") != 0)
	{
    printf("This is not a valid input file\n");
    return -1;
    }
	int nn=0;
	while(1) 
	{
		nn++;
		if(nn>1000)
		{
			cerr<<"配置文件应该用STNLFFM_PARAMETER_END结尾"<<endl;
				break;
		}
		memset(buffer,0,1000);
		if(fgets(buffer, 1000, in)==NULL)
			continue;
		if(strcmp(buffer, "STNLFFM_PARAMETER_END") == 0) break;
		tokenptr = strtok(buffer, separator);
		label=tokenptr;
		if(strcmp(label,"#") == 0) continue;
		while(tokenptr != NULL) 
		{
			tokenptr = strtok(NULL, separator);
			if(strcmp(label, "NUM_IN_PAIRS") == 0) 
			{
				par->NUM_PAIRS = atoi(tokenptr);
				if(par->NUM_PAIRS<=0)
				{
					cerr<<"参考影像至少一对"<<endl;
					return -1;
				}
			}
			else if(strcmp(label, "IN_PAIR_LANDSAT_FNAME") == 0)
				for(i=0; i<par->NUM_PAIRS; i++) 
				{
					sscanf(tokenptr, "%s",readpath);
					psensor[i].Read(readpath);
					tokenptr = strtok(NULL, separator);
				}
			else if(strcmp(label, "IN_PAIR_MODIS_FNAME") == 0)
				for(i=par->NUM_PAIRS; i<2*par->NUM_PAIRS; i++) 
				{
					sscanf(tokenptr, "%s", readpath);
					psensor[i].Read(readpath);
					tokenptr = strtok(NULL, separator);
				}
			else if(strcmp(label, "IN_PDAY_MODIS_FNAME") == 0)
			{
				 k = 0;
				do 
				{
				 sscanf(tokenptr, "%s",  psensor[2*(par->NUM_PAIRS+k)].outpath);
				 tokenptr = strtok(NULL, separator);
				// psensor[2*(par->NUM_PAIRS+k)].Read(readpath);
				// sscanf(tokenptr, "%s",  psensor[2*(par->NUM_PAIRS+k)+1].outpath);
				 k++;
				}while(tokenptr != NULL);
				if(par->NUM_PREDICTIONS == 0)
					par->NUM_PREDICTIONS = k;
				else if(k != par->NUM_PREDICTIONS) 
				{
					printf("\nnumber of IN_PDAY_MODIS_MASK does not match IN_PDAY_MODIS_FNAME\n");
					return -1;
				}
			}
			else if(strcmp(label, "OUT_PDAY_LANDSAT_FNAME") == 0) 
			{
				k = 0;
				do
				{
				sscanf(tokenptr, "%s",  psensor[2*(par->NUM_PAIRS+k)+1].outpath);
				//psensor[2*par->NUM_PAIRS+1].resize(psensor[0].getWidth(),psensor[0].getHeight());
				 tokenptr = strtok(NULL, separator);
				k++;
				}while(tokenptr != NULL);
				if(par->NUM_PREDICTIONS == 0)
					par->NUM_PREDICTIONS = k;
				else if(k != par->NUM_PREDICTIONS) 
				{
					printf("\nnumber of IN_PDAY_MODIS_MASK does not match IN_PDAY_MODIS_FNAME\n");
					return -1;
				}
			}
			else if(strcmp(label, "The_width_of_searching_window") == 0)
				par->WIN_SIZE=atoi(tokenptr);
			else if(strcmp(label, "is_limit_a_CalcuRela") == 0)
				par->r=atoi(tokenptr);
			else if(strcmp(label, "Filter_parameters") == 0)
				par->h=atof(tokenptr);
			else if(strcmp(label, "L_ERR") == 0)
				par->L_ERR=atof(tokenptr);
			else if(strcmp(label, "M_ERR") == 0)
				par->M_ERR=atof(tokenptr);
			else if(strcmp(label, "patchSize") == 0)
				par->pathSize=atoi(tokenptr);
			else if(strcmp(label, "d") == 0)
				par->d=atof(tokenptr);
			else if(strcmp(label, "p") == 0)
				par->p=atoi(tokenptr);
			else if(strcmp(label, "G_Type") == 0)
				sscanf(tokenptr, "%s",par->G_Type);
			label = tokenptr;
		}
		
	}
	for(i=1; i<2*par->NUM_PAIRS; i++) 
    {
		if(psensor[0].getHeight()!=psensor[i].getHeight()||psensor[0].getWidth()!=psensor[i].getWidth())
		{
			cerr<<"输入影像图幅范围不匹配"<<endl;
			return -1;
		}
	}
	fclose(in);
  return 0;
}