#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <stdlib.h>
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tag.h"
#include "view.h"

// define constants
#define MAX_DIGITS 8

// define physical constants and window size
double drop_dia; // m
double SIGMA = 72.8e-3; //N/m
double box_length; //m
double start_height; //m

// Define relevant dimensionless numbers
double Re;
double We;

// create scalar field f0[] to hold the fluid interface
scalar f0[];

// define max refinement level for grid
int MAXLEVEL = 8;

// define error tolerance
double uemax = 0.1;

// define initial velocity
double u0; // m/s

// define end time for the simulation
double t_end;// = 0.03; 

// boundary conditions
u.t[bottom] = dirichlet(0.);
u.n[bottom] = dirichlet(0.);
f[bottom] = 0;

// Create a name for the directory to save the image files
char path_str[20];

// define a name for the logfile and a pointer to the logfile
char logfile[30];
FILE *logfptr;

// define a name for the output directory which will be created
char dirname[30];

int main(int argc, char *argv[]){
    init_grid(64);
   
    // take drop dia and u0 from command line arguments 
    drop_dia = atof(argv[1]);
    u0 = atof(argv[2]); 

    // define the box_length and start_height based on droplet diameter
    box_length = 20*drop_dia;
    start_height = 5*drop_dia;

    origin(-box_length/2, 0);
    
    // define the end time of the simulation
    t_end = 1.3 * start_height / u0;

    size(box_length);
    
    // define the densities and viscosities of the two fluids
        // fluid 1: Water at room temp
        // fluid 2: Air at room temp
    rho1 = 997.; // kg/m^3
    rho2 = 1.293; // kg/m^3
    mu1 = 0.89e-3; // Pa*s
    mu2 = 1.8e-5; // Pa*s
     
    // create a folder for saving the outputs
    sprintf(dirname, "v=%g__D=%g", u0, drop_dia);
    mkdir(dirname, 0755);

    // Create the filename and create the file for logging the simulation progress
    sprintf(logfile, "%s/log.log", dirname);
    logfptr = fopen(logfile, "w");
    fclose(logfptr);

    // define the surface tension coefficient
    f.sigma = SIGMA;
    
    // Calculate Reynolds number for the impact droplet
    Re = rho1 * u0 * drop_dia / mu1;
    We = rho1 * drop_dia *sq(u0) / f.sigma;

    run();
}

event init(t=0){
    // Print the Reynold's number and Weber number for this test case
    printf("Re: %d, We:%d\n", (int) round(Re), (int) round(We));

    // refine grid within a circular region with radius sqrt(2) times the radius of the droplet
    refine((sq(x)+sq(y-start_height) < 2*sq(drop_dia/2)) && (level<MAXLEVEL));

    // set up the interface between the fluids.
        // f0 takes on the value of 1 inside the circle and zero outside
    fraction(f0, sq(drop_dia/2) - sq(x)-sq(y-start_height));
    f0.refine = f0.prolongation = fraction_refine;
    restriction({f0});

    // set the value of the field f[] to f0[] at each point in the domain
    foreach(){
        f[] = f0[];
        u.y[] = -u0*f[];
    }
    boundary(all);
}

event log_status(i++){
    // print out the iteration number and time 
    char logstr[20];
    sprintf(logstr, "%d %g\n", i, t);

    logfptr = fopen(logfile, "a");
    fputs(logstr, logfptr);
    fclose(logfptr);
}

event end(t=t_end){
    // Print the Reynold's number and Weber number for this test case
    printf("Re: %d, We:%d\n", (int) round(Re), (int) round(We));
    char donestr[20];
    sprintf(donestr, "Re:%d, We:%d\n", (int) round(Re), (int) round(We));
 
    logfptr = fopen(logfile, "a");
    fputs(donestr, logfptr);
    fclose(logfptr);
}

#if 0
event initial_graphics_display(i=0){
    // set the view parameters for the display 
    view(tx=0, ty=-0.5, width=800, height=800);
    clear();

    // draw the interface between the two fluids
    draw_vof("f");

    //create numbered grid lines
    box();

    // save the resu
    save("init_img.png");
}
#endif

event movie(i+=5){
    view(tx=0, ty=-0.5, width=800, height=800);
    clear();
    draw_vof("f");
    squares("u.y", linear=true, spread=10);
    box();
        
    char img_index[MAX_DIGITS] = "00000000";
    
    char filename[50];
    if(i != 0){
       
        int num_digits = floor(log10(i)) + 1;
        printf("%d\n", num_digits);
        char num_as_str[num_digits];
        sprintf(num_as_str, "%d", i);
    
        int start_ind = MAX_DIGITS - num_digits;
        for(int j=0; j<=num_digits-1;j++){
            img_index[start_ind + j] = num_as_str[j];
        }
    }
    //sprintf(filename, "output_vids/Re=%d_We=%d.mp4", (int) round(Re), (int) round(We));
    char png_save_path[50];
    sprintf(png_save_path, "%s/%s_t=%g.png", dirname, img_index, t);
    save(png_save_path);
}

event adapt(i++){
    adapt_wavelet({f, p, u}, (double[]){0.01,0.01,0.01,0.01,0.01}, maxlevel=MAXLEVEL);
}

