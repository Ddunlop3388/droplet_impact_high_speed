#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <stdlib.h>
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "contact.h"
#include "tension.h"
#include "view.h"
#include "reduced.h"
#include "tag.h"
#include "navier-stokes/conserving.h"

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
double t_dl; // dimensionless time

// create scalar field f0[] to hold the fluid interface
scalar f0[];

// define max refinement level for grid
int MAXLEVEL;

// define error tolerance
double uemax = 0.1;

// define initial velocity
double u0; // m/s

// define end time for the simulation
double t_end;

// set the equilibrium contact angle
vector h[];
double theta0 = 90.;
h.t[bottom] = contact_angle(theta0 * pi/180.);
h.r[bottom] = contact_angle(theta0 * pi/180.);

// boundary conditions --> no slip
u.t[bottom] = dirichlet(0.);
u.r[bottom] = dirichlet(0.);
u.n[bottom] = dirichlet(0.);
f[bottom] = 0;

// create a name for the big folder holding all simulation results
char save_location_name[30];

// Create a name for the directory to save the image files
char path_str[20];

// define a name for the logfile and a pointer to the logfile
char logfile[120];
FILE *logfptr;

// define a name for the output directory which will be created
char dirname[110];

int main(int argc, char *argv[]){
    init_grid(64);
     
    // define the acceleration due to gravity
    G.y = -9.81;

    // take drop dia and u0 from command line arguments
    drop_dia = atof(argv[1]);
    u0 = atof(argv[2]);
    sprintf(save_location_name, "refinement_level%d", atoi(argv[3]));

    MAXLEVEL = atoi(argv[3]);
    // define the box_length and start_height based on droplet diameter
    box_length = 6.0*drop_dia;
    start_height = 0.5 * drop_dia;

    origin(-box_length/2, 0, -box_length/2);

    // define the end time of the simulation
    t_end = 10.0*start_height / u0;

    size(box_length);

    // define the densities and viscosities of the two fluids
        // fluid 1: Water at room temp
        // fluid 2: Air at room temp
    rho1 = 998.; // kg/m^3
    rho2 = 1.2; // kg/m^3
    mu1 = 1.002e-3; // Pa*s
    mu2 = 1.820e-5; // Pa*s

    // create the save location name for the big folder for all simulatinos
    //sprintf(save_location_name, "refinement_level%d", MAXLEVEL);
    //mkdir(save_location_name, 0755);

    // create a folder for saving the outputs
    sprintf(dirname, "%s/v=%g__D=%g", save_location_name, u0, drop_dia);
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

double max_in_array(double arr[], int n){
    double max = arr[0];

    for(int i=1; i<n; i++){
        if (arr[i] > max){
            max = arr[i];
        }
    }
    return max;
}

event init(t=0){
    // Print the Reynold's number and Weber number for this test case
    printf("Re: %d, We:%d\n", (int) round(Re), (int) round(We));

    // refine grid within a spherical region with radius sqrt(2) times the radius of the droplet
    refine((sq(x)+sq(y-start_height) + sq(z) < 2*sq(drop_dia/2)) && (level<MAXLEVEL));

    // set up the interface between the fluids.
        // f0 takes on the value of 1 inside the circle and zero outside
    fraction(f0, sq(drop_dia/2) - sq(x)-sq(y-start_height) - sq(z));
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
    
    t_dl = u0 * t / drop_dia; 
    char logstr[50];
    sprintf(logstr, "%d %g %g\n", i, t, t_dl);

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

event movie_front(i+=5){
    view(tx=0, ty=-0.5, width=1100, height=1100, theta=pi);
  
    clear();
    draw_vof("f");
    //squares("u.y", linear=true, spread=10);
    box();

    char img_index[MAX_DIGITS] = "00000000";

    if(i != 0){

        int num_digits = floor(log10(i)) + 1;
        char num_as_str[num_digits];
        sprintf(num_as_str, "%d", i);

        int start_ind = MAX_DIGITS - num_digits;
        for(int j=0; j<=num_digits-1;j++){
            img_index[start_ind + j] = num_as_str[j];
        }
    }
    //sprintf(filename, "output_vids/Re=%d_We=%d.mp4", (int) round(Re), (int) round(We));
    char png_save_path[170];
    sprintf(png_save_path, "%s/%s_t=%g_tau=%g_front.png", dirname, img_index, t, t_dl);
    save(png_save_path);
}

event movie_top(i+=5){
    view(width=1100, height=1100, camera = "top");
    
    clear();
    draw_vof("f");
    
    char img_index[MAX_DIGITS] = "00000000";

    if (i != 0){
        int num_digits = floor(log10(i)) + 1;
        char num_as_str[num_digits];
        sprintf(num_as_str, "%d", i);

        int start_ind = MAX_DIGITS - num_digits;

        for(int j=0; j <= num_digits - 1; j++){
            img_index[start_ind + j] = num_as_str[j];
        }
    }
    char png_save_path[170];
    sprintf(png_save_path, "%s/%s_t=%g_tau=%g_top.png", dirname, img_index, t, t_dl);
    save(png_save_path);
}

event adapt(i++){
    adapt_wavelet({f, p, u}, (double[]){0.01,0.01,0.01,0.01,0.01}, maxlevel=MAXLEVEL);

    boundary(all);
}

event update_contact_angle(i++){
    scalar m[];
    foreach()
        m[] = f[] > 1e-3;
    int num_droplets = tag(m);
    
    //printf("%d\n", num_droplets); 
    double v[num_droplets];
    coord b[num_droplets];

    for (int j=0; j< num_droplets; j++){
        v[j] = b[j].x = b[j].y = b[j].z = 0.;
    }    

    foreach(serial){
        if (m[] > 0){
            int j = m[] - 1;
            v[j] += dv() * f[];
            coord p = {x, y, z};
            foreach_dimension()
                b[j].x += dv() * f[] * p.x;
        }
    }
    
    printf("%g\n", max_in_array(v, num_droplets));
    #if _MPI
        MPI_Allreduce(MPI_IN_PLACE, v, num_droplets, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, b, 3*num_droplets, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    #endif
    #if 0
    for (int j=0; j < num_droplets; j++){
        fprintf(fout, "%d %g %d %g %g %g %g\n", i, t, j, v[j], b[j].x/v[j], b[j].y/v[j], b[j].z/v[j]);
    }   
    #endif
}

