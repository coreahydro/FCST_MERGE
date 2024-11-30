#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include "read_bin_data.h"

SITE_HEADER h_r;

float* read_bin_(char *filename,int *nmax, int *nt, double data[*nt][*nmax])
{
    gzFile fp = NULL;
    int skip = 0;
    short int temp,fst_time;
    short int **data_temp;    
    int head_size,data_size,i,xy,t,z,n;
    int xydim,tdim;
    short int spval=-32768;
    
    for(i=0;i<256;i++){if(filename[i]==' '){filename[i]='\0';break;} }


    if( (fp = gzopen(filename, "r")) == NULL ){return NULL;}

    //head_size = sizeof(short int) * 1024;
    xydim = *nmax;
    tdim = *nt;

    head_size = sizeof(short int) * 1024;
//    gzread(fp,&h_r, sizeof(short int) * 1024);  // Read header information
/*    printf("%d \n",sizeof(SITE_HEADER));
    gzread(fp,&h_r,sizeof(h_r));  // Read header information        
    printf("%d %d \n",h_r.xdim,h_r.ydim);
    printf("%f %f \n",h_r.lamc_slat1,h_r.lamc_slat2);    
    printf("%f %f \n",h_r.lamc_olon,h_r.lamc_olat);
    printf("%f %f \n",h_r.lamc_xo,h_r.lamc_yo);
    printf("%f \n",h_r.lamc_grid);
    printf("%d \n",h_r.lead_time);
    printf("%d \n",h_r.write_dt);
    printf("%d \n",h_r.num_dataset);
    printf("%d \n",head_size);      
    printf("%d \n",data_size);
    printf("%d %d \n",xydim,tdim);*/
    skip = head_size;
    gzseek(fp, skip, SEEK_CUR);
    data_temp = (short int **)malloc(sizeof(short int *)*tdim); 
   	n=0;
    for(t=0;t<=tdim;t++) //Be careful!! "t<=tdim"
            {
           	gzread(fp, &fst_time, sizeof(short int) * 1);
         		//if(t == 10) {printf("%d \n",fst_time);}
            data_temp[t] = (short int *)malloc(sizeof(short int)*xydim); 
         		gzread(fp, data_temp[t], sizeof(short int)*xydim);
            for(xy=0;xy<xydim;xy++) //xdim
            {
              temp=data_temp[t][xy];
            	if(temp < 0 || temp > 10000 )
            	 data[t][xy]=-99.9;
            	else
            	 data[t][xy]=temp/100.; //100.
            	//if(t == 0) {printf("%d %d %10.2f \n",n,data_temp[t],data[t][xy]);}
          	}
//            printf("%d \n",t);
//          printf("\n");
        		}    
    gzclose(fp);
    //printf("%d %d \n",*nmax,*nt);

    //return imsi_data;
       }

