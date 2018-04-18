//#include"cuLayer.h"
typedef struct { 
 // char outpath[1000];
 // int A;
	int p;      //产品类型
 float gamma;
 float d;    //光谱经验参数
  char G_Type[50];
 // int pf;
 // int pc;
  float L_ERR;
  float M_ERR;
  int pathSize;
  int NUM_PREDICTIONS;
 // float uncertain;
  float h;                       /*滤波参数*/
  int r;   
  float m;
  float l;
  int NUM_PAIRS;             /* number of input data pair */ 
  int WIN_SIZE;              /* window size determined by the MAX_SEARCH_DIS */
} PARAMETER;
