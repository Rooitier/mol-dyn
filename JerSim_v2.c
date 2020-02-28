#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* Defintions */

#define PartNo 50 // Number of particles
#define Dim 2 // Box dimension
#define TwoPi 6.283185307 // 2Pi
#define Filenamelength 100 // Max characters for a file name eg. test1
#define sig 0.6 
#define eps 0.001
#define shift 0.065

/*-----------------*/


/* Box-Muller transformation is used to produce random samples for the Maxwell distribution which will be used to provide the intial values of the speed of the particles */

double sampgauss(double temp){

	double u1, u2, z0;

	u1 = drand48();
	u2 = drand48();

	z0 = sqrt(-2.0 * log(u1)) * cos(TwoPi * u2);

	return (sqrt(temp) * z0);


}
/* L-J potential force */

void lj_force (double pos[PartNo][Dim], double f_net[PartNo][Dim]){
    double f_lj, r;
    int i, j, d;
    for (i=0; i<PartNo; i++){
		
	for (d=0; d < Dim; d++){
		f_net[i][d] = 0.0;
            for (j = i+1; j < PartNo; j++){
		

		r = sqrt(pow(pos[i][0]-pos[j][0],2) + pow(pos[i][1]-pos[j][1],2)); // Distance calculation

        if (r <= 2.5*sig) { // Lower truncation
            f_lj = (-4*eps*(1/r)*(6*pow(sig,6)/pow(r,6)-12*pow(sig,12)/pow(r,12))+shift*0.0163*eps) * (pos[i][d]-pos[j][d])/(r);
        }
		else if (r > 2.5*sig){
		 	f_lj = 0.0;
		} 
                else { // Otherwise just compute normal LJ_force on each particle
                    f_lj = -4*eps*(1/r)*(12*pow(sig,6)/pow(r,6)-12*pow(sig,12)/pow(r,12)) * (pos[i][d]-pos[j][d])/r;
                }
                f_net[i][d] += f_lj;
                f_net[j][d] -= f_lj;
            }
        }                 
    }
}
void lj_pot (double pos[PartNo][Dim], double v_pot[PartNo][Dim]){
    int i, j, d;
    double v_temp;
    for (i=0; i<PartNo; i++){
        for (d=0; d < Dim; d++){
            v_pot[PartNo][Dim] = 0.0;
            for (j=i+1; j < PartNo; j++){
              r = sqrt(pow(pos[i][0]-pos[j][0],2) + pow(pos[i][1]-pos[j][1],2));  
            if (r <= 2.5*sig){
                v_temp = 4*eps*((pow(sig,12)/pow(r,12))-(pow(sig,6)/pow(r,6))) + 0.0163*eps
            } 
            else if (r > 2.5*sig){
                v_temp = 0.0;
            }
            else
            {
                v_temp = 4*eps*((pow(sig,12)/pow(r,12))-(pow(sig,6)/pow(r,6)));
            }
            }
        }
    v_pot[i][d] += v_temp; 
    }
}

/* L-J potential force in y */

/*void dljy (double posi[2], double posj[2] epsilon, sigma){
  r = pow(posi[1]-posj[1],2) + pow(posi[2]-posj[2],2)
    return -4*ep*(12*pow(sig,11)/pow(r,7)-24*pow(sig,12)/pow(r,23))*(posj[2]-posi[2])/r

}


/* Initialisation of the position and velocities of the particles. (Changes the entries in the array that calls it??) ./a.out 0.0001 0.5 100 1.1 1.2 1.0 datafile will run each of the elements seperated by commas*/
/* Keeping the centre of velocity constant so the box doesn't drift */

void initpart (double pos[PartNo][Dim], double vel[PartNo][Dim], double temperature, double boxdims[Dim]){

	int i, d;
    double VCoM[Dim];
	srand48((long)time(NULL));

	for (i=0; i<PartNo; i++) {

		for (d=0; d<Dim; d++) {

			pos[i][d] = drand48()*boxdims[d]; // Places the particle at a random position within 				the boxes length
			vel[i][d] = sampgauss(temperature); // Sample velocity from the Maxwell distribution 				using Box-Muller transform. *This can be hard-coded to test the interaction potentials 				later on. 
		
            VCoM[d] += vel[i][d];
        }
	}
    for ( d = 0; d < Dim; d++){
        VCoM[d] /= PartNo;
    }
    for ( i = 0; i < PartNo; i++){
        for ( d = 0; d < Dim; d++){
            vel[i][d] -= VCoM[d];
        }
    }
}

/* We now wish to determine the force acting on each particle by Newton II. This will be done by using the LJ 6-12 potential and velocity Verlet algorithm. For now the force on the particles is zero, thus the acceleration is set to zero. */

/* void calaccel (double pos[PartNo][Dim], double acc[PartNo][Dim]){
    double f_net;
	int i, j, d; // Particles pairs (i,j) and Cartesian component index (d).
	acc[i][j][d] = 0.0;
	for (i=0; i<PartNo; i++){
        for(j=i+1; j<PartNo; j++){
		    for (d=0; d<Dim; d++){
			    acc[i][d] += lj_force(pos);
                acc[j][d] -= lj_force(pos);
            
            }
		}
	}
}




/* Integration step. Implementing the Verlet algorithm for updating positions, velocities - with provided acc. This section of codes also introduces a way for the boundary of the box to be established - the particles density remaining constant in this way */

void move_part(double pos[PartNo][Dim], double vel[PartNo][Dim], double acc[PartNo][Dim], double dt, double *t, double *KinE, double boxdims[Dim]){

	int i, d;
	double out_pos, newvel, newpos, newacc, f_net[PartNo][Dim];

	*KinE = 0.0;
    lj_force(pos, f_net);
	for (i=0; i<PartNo; i++){
		for (d=0; d<Dim; d++){
            out_pos = pos[i][d];
            

			while (out_pos > boxdims[d]){
			
				out_pos -=  boxdims[d]; // Check the positions of the particles in the loop and updates it. 

            }
            while (out_pos < 0)
            {
                out_pos += boxdims[d];
            }

            newpos = out_pos + vel[i][d] * dt;
            vel[i][d] += 0.5 * f_net[i][d] * dt;
            
		pos[i][d] = newpos;
		f_net[i][d] = 0.0;

		
		}
	}

    for (i=0; i<PartNo; i++){
        for (d=0; d<Dim; d++){
            newvel += 0.5 * dt * f_net[i][d];
        
        *KinE += 0.5*newvel*newvel;
        vel[i][d] = newvel;
        
        }
    }
    *t += dt;
} 


/* void verlet_update (double pos[PartNo][Dim], double vel[PartNo][Dim], double acc[PartNo][Dim], double dt, double)
    
    int i, d;
    double f_net[PartNo][Dim], dt, t;
    move_part(pos, vel, acc, dt, &t, &KinE, boxdims); */







void writepositions(double pos[PartNo][Dim], char filename[Filenamelength], double boxdims[Dim]){
    
    int i;
    char filenameupdated[Filenamelength];
    
    FILE *posfile;
    
    strcpy(filenameupdated, filename);
    strcat(filenameupdated, "_mathematica.nb"); // data written to filename that is formated and designated a mathematica notebook
    
    posfile = fopen(filenameupdated, "a");

    fprintf(posfile, "ListPlot[{");
    
    for (i=0; i<(PartNo-1); i++) {
        fprintf(posfile, "{ %f, %f },\n", pos[i][0], pos[i][1]);
    }
    
    fprintf(posfile, "{ %f, %f }\n", pos[PartNo-1][0], pos[PartNo-1][1]);
    
    fprintf(posfile, "}, PlotStyle -> PointSize[Large], \n PlotRange -> {{0, %f}, {0, %f}}, AspectRatio -> 1],\n", boxdims[0], boxdims[1]);
    fclose(posfile);
    
}

/* The energies are written to a simple text file: timestamp and kinetic energy. */

void writeenergies(double KinE, double t, char filename[Filenamelength]){
    
    char filenameupdated[Filenamelength];
    
    FILE *KineFile;
    
    strcpy(filenameupdated, filename);
    strcat(filenameupdated, "_KE_list.txt"); // data written to filename that is pure text
    
    KineFile = fopen(filenameupdated, "a");

    fprintf(KineFile, "%f %f\n", t, KinE);
    
    fclose(KineFile);
    
}



/* Creating a mathematica to animate the motions of the particles */





void tidyupmathematicafile(double pos[PartNo][Dim], char filename[Filenamelength], double boxdims[Dim]){
    
    int i;
    char filenameupdated[Filenamelength];
    
    FILE *posfile;
    
    strcpy(filenameupdated, filename);
    strcat(filenameupdated, "_mathematica.nb");
    
    posfile = fopen(filenameupdated, "a");

    fprintf(posfile, "ListPlot[{");
    
    for (i=0; i<(PartNo-1); i++) {
        fprintf(posfile, "{ %f, %f },\n", pos[i][0], pos[i][1]);
    }
    
    fprintf(posfile, "{ %f, %f }\n", pos[PartNo-1][0], pos[PartNo-1][1]);
    
    fprintf(posfile, "}, PlotStyle -> PointSize[Large], \n PlotRange -> {{0, %f}, {0, %f}}, AspectRatio -> 1]}]\n", boxdims[0], boxdims[1]);

    fclose(posfile);
    
}

void startmathematica(char filename[Filenamelength]){
    
    char filenameupdated[Filenamelength];
    
    FILE *posfile;
    
    strcpy(filenameupdated, filename);
    strcat(filenameupdated, "_mathematica.nb"); 
    
    posfile = fopen(filenameupdated, "w");

    fprintf(posfile, "ListAnimate[{");
    
    fclose(posfile);
    
}


/*Main function */

int main(int argc, const char *argv[]) {

    double boxdims[Dim]; // Array for x and y dimensions of the box
    double initi_temp; // Temperature for distribution sampling
    double curr_time=0.0, deltat, endtime; // Time initialisation
    double VCoM =0.0; // Centre-of-velocity
    double Kinetic=0.0; // Temporary storage of kinetic energy
    double positions[PartNo][Dim]; // 2D array for particle positions
    double velocities[PartNo][Dim]; // 2D array for particle vel
    double accelerations[PartNo][Dim]; // 2D array for particle accelerations
    char outputfilename[Filenamelength];
    int outputinterval, stepssinceoutput=0;


    deltat = atof(argv[1]);
    endtime = atof(argv[2]);
    outputinterval = atof(argv[3]);
    boxdims[0] = atof(argv[4]);
    boxdims[1] = atof(argv[5]);
    initi_temp = atof(argv[6]);
    strcpy(outputfilename, argv[7]);


    printf("code has read in delta t= %.4e end time=%.4e interval for data=%d\n box x=%.4e box y=%.4e initial t=%.4e file=%s\n",deltat, endtime, outputinterval, boxdims[0], boxdims[1], initi_temp, outputfilename);
    
    
    initpart(positions, velocities, initi_temp, boxdims);
    
    startmathematica(outputfilename);
    
    while (curr_time < endtime) {
        
        // calaccel(positions, accelerations);
        
        move_part(positions, velocities, accelerations, deltat, &curr_time, &Kinetic, boxdims);
        
        stepssinceoutput += 1;
        
        if (stepssinceoutput == outputinterval) {
            stepssinceoutput=0;
            
            writepositions (positions, outputfilename, boxdims);
            writeenergies (Kinetic, curr_time, outputfilename);
            
        }
        
    }
    
    tidyupmathematicafile(positions, outputfilename, boxdims);
    
    
}




















