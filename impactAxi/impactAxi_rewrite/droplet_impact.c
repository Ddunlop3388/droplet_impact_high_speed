#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <stdlib.h>
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tag.h"
#include "view.h"

#define MAX_DIGITS 8

double drop_dia;
double SIGMA = 72.8e-3;
double box_length;
double start_height;

double Re;
double We;

scalar f0[];

int MAXLEVEL;

double uemax = 0.1;

double u0;

double t_end;

u.t[left] = dirichlet(0.);
u.n[left] = dirichlet(0.);
f[left] = 0.;

char save_location_name[30];

char path_str[20];

char logfile[90];

FILE *logfptr;

char dirname[80];

int main(int argc, char *argv[]){
    init_grid(64);

    drop_dia = atof(argv[1]);
    u0 = atof(argv[2]);
    sprintf(save_location_name, "refinement_level%d", atoi(argv[3]));
    
    MAXLEVEL = atoi(argv[3]);

    box_length = 10. * drop_dia;
    start_height = 5. * drop_dia;
    
    t_end = 1.3 * start_height / u0;

    size(box_length);

    rho1 = 997.;
    rho2 = 1.293;
    mu1 = 0.89e-3;
    mu2 = 1.8e-5;

    sprintf(dirname, "%s/v=%g__D=%g", save_location_name, u0, drop_dia);
    mkdir(dirname, 0755);

    sprintf(logfile, "%s/log.log", dirname);
    logfptr = fopen(logfile, "w");
    fclose(logfptr);

    f.sigma = SIGMA;

    Re = rho1 * u0 * drop_dia / mu1;
    We = rho1 * drop_dia * sq(u0) / f.sigma;

    run();

    
}

event init(i=0){
    printf("Re: %d, We: %d\n", (int) round(Re), (int) round(We));

    refine((sq(y) + sq(x-start_height) < 2*sq(drop_dia/2)) && (level<MAXLEVEL));

    fraction(f0, sq(drop_dia/2) - sq(y) - sq(x - start_height));

    f0.refine = f0.prolongation = fraction_refine;

    restriction({f0});

    foreach(){
        f[] = f0[];
        u.x[] = -u0*f[];
    }
    boundary(all);
}

#if 0
event init_display(i=0){
    printf("display_section\n");
    view(tx=-0.5, ty=-0.5, width=800, height=800);
    clear();

    draw_vof("f");
    squares("u.x", linear=true, spread=10);

    box();

    save("init_img.png");
}
#endif

event log_status(i++){
    char logstr[20];
    sprintf(logstr, "%d %g\n", i, t);

    logfptr = fopen(logfile, "a");
    fputs(logstr, logfptr);
    fclose(logfptr);
}

event end(t=t_end){
    printf("Re: %d, We: %d\n", (int) round(Re), (int) round(We));
    char donestr[20];
    sprintf(donestr, "Re: %d, We: %d\n", (int) round(Re), (int) round(We));

    logfptr = fopen(logfile, "a");
    fputs(donestr, logfptr);
    fclose(logfptr);
}

event movie(i+=5){
    view(tx=-0.5, ty=-0.5, width=800, height=800);
    
    clear();
    draw_vof("f");
    //squares("u.x", linear=true, spread=10);
    box();

    char img_index[MAX_DIGITS] = "00000000";

    if (i != 0){
        int num_digits = floor(log10(i)) + 1;
        printf("%d\n", num_digits);
        char num_as_str[num_digits];
        sprintf(num_as_str, "%d", i);

        int start_ind = MAX_DIGITS - num_digits;

        for(int j=0; j<=num_digits-1;j++){
            img_index[start_ind + j] = num_as_str[j];
        }

    }
    char png_save_path[100];
    sprintf(png_save_path, "%s/%s_t=%g.png", dirname, img_index, t);
    save(png_save_path);
}

event adapt(i++){
    adapt_wavelet({f, p, u}, (double[]) {0.01, 0.01, 0.01, 0.01, 0.01}, maxlevel=MAXLEVEL);
    boundary(all);
}
