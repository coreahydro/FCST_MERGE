typedef struct
{
    short int time;
    short int day;
    short int month;
    short int year;
    short int nx;
    short int ny;
    short int xydim;
    short int lat_d;
    short int lat_m;
    short int lat_s;
    short int lon_d;
    short int lon_m;
    short int lon_s;
    short int height;
    short int lead_time;
    short int write_dt;
    short int num_dataset;
    short int xmargin;
    short int ymargin;    
}FCST_NEW_HEAD;

typedef struct
{
    int xdim; 
    int ydim; 
    float lamc_slat1;
    float lamc_slat2;
    float lamc_olon;
    float lamc_olat;
    float lamc_xo;
    float lamc_yo;
    float lamc_grid;
    short int lead_time;
    short int write_dt;
    short int num_dataset;     
}SITE_HEADER;

extern SITE_HEADER h_r;

