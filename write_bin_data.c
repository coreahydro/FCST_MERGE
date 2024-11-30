#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
//#include "read_bin_data.h"
#include "write_bin_data.h"

MGRF_HEADER m_h;

float write_bin_(char *filename,char *fcst_info,char *nwp_info,int *nmax, int *nx, int *ny, int *nt, double data[*nt][*nmax],int *nt_out)
{
    gzFile fp = NULL;
    int skip = 0;
    double temp;
    short int fst_time;
    short int **data_temp;    
    int head_size,data_size,i,ix,iy,xy,t,z,n;
    int xydim,xdim,ydim,tdim;
    short int spval=-32768;

    for(i=0;i<256;i++){if(filename[i]==' '){filename[i]='\0';break;} }

    //head_size = sizeof(short int) * 1024;
    xydim = *nmax;
    xdim = *nx;
    ydim = *ny;    
    tdim = *nt_out;
    //tdim = *nt;
		//strcpy(m_h.info1,fcst_info);
		strncpy(m_h.info1,fcst_info,18);
		//strcpy(m_h.info2,nwp_info);
		strncpy(m_h.info2,nwp_info,18);

		/*printf("%s \n",fcst_info);
		printf("%s \n",nwp_info);		
		printf("%s \n",m_h.info1);
		printf("%s \n",m_h.info2);
		printf("%d \n",xdim);
		printf("%d \n",ydim);*/

    fp = gzopen(filename,"w");
    if (fp == NULL)
	 {
        fprintf(stderr,"ERROR : cannot open file : %s\n",filename);
        return -12345;
     	 }
    else
    	{
    gzwrite(fp,&m_h,sizeof(MGRF_HEADER));

    data_temp = (short int **)malloc(sizeof(short int *)*tdim);
    for(t=0;t<=tdim;t++) 
		{

      data_temp[t] = (short int *)malloc(sizeof(short int )*xydim);
    //printf("%d \n",sizeof(data_temp[t]));

      for(iy=0;iy<ydim;iy++) //xdim
      {
	      for(ix=0;ix<xdim;ix++) //xdim
  	    {
      		//xy=ix+iy*xdim;            // Write data with lower left cell as (1, 1)
      		xy=ix+(ydim-iy-1)*xdim;     // Write data with upper left cell as (1, 1)     		
		      temp=data[t][xy];
         	if(temp < 0 || temp > 10000 )
         	 data_temp[t][xy]=-999;
         	else
       	 data_temp[t][xy]=temp*100; //*100; //100.
         //if(t == 0) {printf("%d %10.4f \n",data_temp[t][xy],data[t][xy]);}

      		gzwrite(fp, &data_temp[t][xy], sizeof(short int));
	  	  }
			}    
    // Write data with lower left cell as (1, 1)
    /*  for(xy=0;xy<xydim;xy++) //xdim
      {
				      temp=data[t][xy];
            	if(temp < 0 || temp > 10000 )
            	 data_temp[t][xy]=-999;
            	else
            	 data_temp[t][xy]=temp*10; //*100; //100.
 	            //if(t == 0) {printf("%d %10.4f \n",data_temp[t][xy],data[t][xy]);}

      		gzwrite(fp, &data_temp[t][xy], sizeof(short int));
	    }*/

 		//gzwrite(fp, &data[t], sizeof(double)*xydim);
 		//gzwrite(fp, &data_temp[t], sizeof(short int)*xydim);
		}
       gzclose(fp);

        }
}



