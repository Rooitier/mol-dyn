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
#define sig 1
#define eps 1
#define r_shift 2.5*sig


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
void lj_force (double pos[PartNo][Dim], double acc[PartNo][Dim]){
    double f_lj, r;
    int i, j, d;
    for (i=0; i<PartNo; i++){
	  	for (d=0; d < Dim; d++){ // Removed zero'd acceleration
        for (j = i+1; j < PartNo; j++){
		    		r = sqrt(pow(pos[i][0]-pos[j][0],2) + pow(pos[i][1]-pos[j][1],2)); // Distance calculation

        		if (r <= 2.5*sig) { // Lower truncation
            		f_lj = (-4*eps*(1/(r+r_shift))*(6*pow(sig,6)/pow(r+r_shift,6)-12*pow(sig,12)/pow(r+r_shift,12))) * (pos[i][d]-pos[j][d])/r;
        		}
						else if (r > 2.5*sig){
		 						f_lj = 0.0;
						}
        		acc[i][d] += f_lj;
        		acc[j][d] -= f_lj;
						}
				}
		}
}

/* L-J potential */
void lj_pot (double pos[PartNo][Dim], double *v_pot){
    int i, j;
    double r;
    *v_pot = 0.0;
    for (i=0; i<PartNo; i++){
        for (j=i+1; j < PartNo; j++){
            r = sqrt(pow(pos[i][0]-pos[j][0],2) + pow(pos[i][1]-pos[j][1],2));
            if (r <= 2.5*sig){
            		*v_pot += 4*eps*((pow(sig,12)/pow(r+r_shift,12))-(pow(sig,6)/pow(r+r_shift,6))) + 0.0163*eps;
            }
            else if (r > 2.5*sig){
            		*v_pot += 0.0;
            }
        }
    }
}

/* Periodic interactions: creating periodic images of particles */
int per_int (double pos[PartNo][Dim], double boxdims[Dim], double  per_arr[3*PartNo][Dim]){

    int i, k=0;

		for (i = 0; i < PartNo; i++)
		{ 
        // Corners
    	if (pos[i][0] < 2.5*sig && pos[i][1] < 2.5*sig) // Bottom left
    	{
        	per_arr[k][0] = pos[i][0] + boxdims[0];
        	per_arr[k][1] = pos[i][1] + boxdims[1];
        	k++;
    	}
        if (pos[i][0] > boxdims[0] - 2.5*sig && pos[i][1] < 2.5*sig) // Bottom right
    	{
        	per_arr[k][0] = pos[i][0] - boxdims[0];
        	per_arr[k][1] = pos[i][1] + boxdims[1];
        	k++;
    	}
    	if (pos[i][0] < 2.5*sig && pos[i][1] > boxdims[1] - 2.5*sig) // Top left
    	{
        	per_arr[k][0] = pos[i][0] + boxdims[0];
        	per_arr[k][1] = pos[i][1] - boxdims[1];
        	k++;
    	}
    	if (pos[i][0] > boxdims[0] - 2.5*sig && pos[i][1] > boxdims[1] - 2.5*sig) // Top right
    	{
        	per_arr[k][0] = pos[i][0] - boxdims[0];
        	per_arr[k][1] = pos[i][1] - boxdims[1];
        	k++;
    	} 
        // Edges
        if (pos[i][0] < boxdims[0] && pos [i][1] < 2.5*sig) // Bottom
    	{
        	per_arr[k][0] = pos[i][0];
        	per_arr[k][1] = pos[i][1] + boxdims[1];
        	k++;
    	}
    	if (pos[i][0] < 2.5*sig && pos[i][1] < boxdims[1]) // Left
    	{
        	per_arr[k][0] = pos[i][0] + boxdims[0];
        	per_arr[k][1] = pos[i][1];
        	k++;
    	}
    	if (pos[i][0] < boxdims[0] && pos[i][1] > boxdims[1] - 2.5*sig) // Top
    	{
        	per_arr[k][0] = pos[i][0];
        	per_arr[k][1] = pos[i][1] - boxdims[1];
        	k++;
    	}
    	if (pos[i][0] > boxdims[0] - 2.5*sig && pos[i][1] > boxdims[1]) // Right
    	{
        	per_arr[k][0] = pos[i][0] - boxdims[0];
        	per_arr[k][1] = pos[i][1];
        	k++;
    	}
		}
		return k;
}

/* Periodic force */
void per_forces (double pos[PartNo][Dim], double acc[PartNo][Dim], double per_arr[3*PartNo][Dim], double boxdims[Dim]){

int i, k, l, d;
double r, f_lj;
l = per_int(pos, boxdims, per_arr);
for (k = 0; k < l ; k++){
    for (i = 0; i < PartNo; i++) {
        	r = sqrt(pow(pos[i][0]-per_arr[k][0],2) + pow(pos[i][1]-per_arr[k][1],2));
        	for (d = 0; d < Dim; d++) {
        		if (r <= 2.5*sig) { // Lower truncation
            	f_lj = (-4*eps*(1/(r+r_shift))*(6*pow(sig,6)/pow(r+r_shift,6)-12*pow(sig,12)/pow(r+r_shift,12))) * (pos[i][d]-per_arr[k][d])/r;// segfault was here you had pos[i][d]-pos[k][d] instead of pos[i][d]-per_arr[k][d]
        		}
						else if (r > 2.5*sig){
		 					f_lj = 0.0;
						}
        		acc[i][d] += f_lj;
        	}
        }
    }
}

void lj_per_pot (double pos[PartNo][Dim], double *v_pot, double per_arr[3*PartNo][Dim], double boxdims[Dim]){
    int i, l, k;
    double r;
	l = per_int(pos, boxdims, per_arr);
    for (k=0; k<l; k++){
        for (i=0; i < PartNo; i++){
            r = sqrt(pow(pos[i][0]-per_arr[k][0],2) + pow(pos[i][1]-per_arr[k][1],2));
            if (r <= 2.5*sig){
            		*v_pot += 0.5*(4*eps*((pow(sig,12)/pow(r+r_shift,12))-(pow(sig,6)/pow(r+r_shift,6))) + 0.0163*eps);
            }
            else if (r > 2.5*sig){
            		*v_pot += 0.0;
            }
        }
    }
}
/* Initialisation of the position and velocities of the particles. (Changes the entries in the array that calls it??) ./a.out 0.0001 0.5 100 1.1 1.2 1.0 datafile will run each of the elements seperated by commas*/
/* Keeping the centre of velocity constant so the box doesn't drift */
void initpart (double pos[PartNo][Dim], double vel[PartNo][Dim], double temperature, double boxdims[Dim]){

		int i, d;
    double VCoM[Dim];

		srand48((long)time(NULL));
		for (i=0; i<PartNo; i++) {

				for (d=0; d<Dim; d++) {

						pos[i][d] = drand48()*boxdims[d]; // Places the particle at a random position within the boxes length
						vel[i][d] = sampgauss(temperature); // Sample velocity from the Maxwell distribution using Box-Muller transform. *This can be hard-coded to test the interaction potentials 				later on.

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

/* Integration algorithm */
void move_part(double pos[PartNo][Dim], double vel[PartNo][Dim], double acc[PartNo][Dim], double dt, double *t, double *KinE, double boxdims[Dim], double per_arr[3*PartNo][Dim]){

	int i, d;
	double newvel, newpos, newacc;

	*KinE = 0.0;

		for (i=0; i<PartNo; i++){
				for (d=0; d<Dim; d++){

            	newpos = pos[i][d] + vel[i][d] * dt;
            	vel[i][d] += 0.5 * acc[i][d] * dt;

				while (newpos > boxdims[d]){

									newpos -=  boxdims[d]; // Check the positions of the particles in the loop and updates it.

            	}
            	while (newpos < 0)
            	{
                	newpos += boxdims[d];
            	}

    	pos[i][d] = newpos;
        acc[i][d] = 0.0;

			    }
		}
	per_forces(pos, acc, per_arr, boxdims);
    lj_force(pos,acc);

    for (i=0; i<PartNo; i++){
        for (d=0; d<Dim; d++){
            vel[i][d] += 0.5 * dt * acc[i][d];
        }
        *KinE += 0.5*(vel[i][0]*vel[i][0] + vel[i][1]*vel[i][1]);
    }
    *t += dt;
}


void box_p_v (double pos[PartNo][Dim], double acc[PartNo][Dim], double temperature, double pressure, double n_step){

int i, n, j;
double temp_pv;
temp_pv = pressure; 
    for (i = 0; i < PartNo; i++) // Calculating pressure for each time step
    {
        temp_pv += PartNo*temperature - (1/Dim) * (pos[i][0]*acc[i][0] + pos[i][1]*acc[i][1]); 
    }
    pressure = (1/n_step) * temp_pv; // Value of the average pressure
    printf("pressure = %f\n", pressure);
}




/* The energies are written to a simple text file: timestamp, kinetic energy and potetial energy. */
void writeenergies(double KinE, double t, char filename[Filenamelength], double v_pot){

    char filenameupdated[Filenamelength];

    FILE *KineFile;

    strcpy(filenameupdated, filename);
    strcat(filenameupdated, "_E_list.txt"); // data written to filename that is pure text space seperated into 4 columns

    KineFile = fopen(filenameupdated, "a");

    fprintf(KineFile, "%f %f %f %f\n", t, KinE, v_pot, KinE + v_pot);

    fclose(KineFile);

}



/* Creating a mathematica to animate the motions of the particles
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
/* void tidyupmathematicafile(double pos[PartNo][Dim], char filename[Filenamelength], double boxdims[Dim]){

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
*/


/*Main function */

int main(int argc, const char *argv[]) {

    double boxdims[Dim]; // Array for x and y dimensions of the box
    double initi_temp; // Temperature for distribution sampling
    double curr_time=0.0, deltat, endtime, steps; // Time initialisation
    double VCoM[Dim] = {0}; // Centre-of-velocity
    double Kinetic=0.0; // Temporary storage of kinetic energy
    double Potential=0.0; // Storage of potential energy
    double pressure = 0.0; // Storage of pressure
    double positions[PartNo][Dim]; // 2D array for particle positions
    double velocities[PartNo][Dim]; // 2D array for particle vel
    double accelerations[PartNo][Dim] = {0}; // 2D array for particle accelerations
    double per_arr[3*PartNo][Dim]; // 2D array for periodic ghost particles
    char outputfilename[Filenamelength];
    int outputinterval, stepssinceoutput=0;

    steps = atof(argv[2])/atof(argv[1]);
    deltat = atof(argv[1]);
    endtime = atof(argv[2]);
    outputinterval = atof(argv[3]);
    boxdims[0] = atof(argv[4]);
    boxdims[1] = atof(argv[5]);
    initi_temp = atof(argv[6]);
    strcpy(outputfilename, argv[7]);

    printf("code has read in delta t= %.4e end time=%.4e interval for data=%d\n box x=%.4e box y=%.4e initial t=%.4e file=%s\n",deltat, endtime, outputinterval, boxdims[0], boxdims[1], initi_temp, outputfilename);
    printf("n = %f\n", steps);

    initpart(positions, velocities, initi_temp, boxdims);
    lj_force(positions, accelerations);
    per_forces(positions, accelerations, per_arr, boxdims);        
    box_p_v(positions, accelerations, initi_temp, pressure, steps);
    //startmathematica(outputfilename);
    while (curr_time < endtime) {

        lj_force(positions, accelerations);
        per_forces(positions, accelerations, per_arr, boxdims);
        move_part(positions, velocities, accelerations, deltat, &curr_time, &Kinetic, boxdims, per_arr);

        stepssinceoutput += 1;

        if (stepssinceoutput == outputinterval) {
            stepssinceoutput=0;
						lj_pot(positions, &Potential);
						lj_per_pot(positions, &Potential, per_arr, boxdims);

                        // writepositions (positions, outputfilename, boxdims);
                        writeenergies (Kinetic, curr_time, outputfilename, Potential);

        }

    }

    //tidyupmathematicafile(positions, outputfilename, boxdims);


}
