//-----------------------------------------------------------------------
// Particle Filters Simulation    ---- PARALLEL WITH OPEN MP 
//-----------------------------------------------------------------------
//  Written by: Javier Pastorino
//  Updated in Dec-2016
//-----------------------------------------------------------------------
#include <iostream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
using namespace std;

//-----------------------------------------------------------------------
//   Data Structures and Constants
//-----------------------------------------------------------------------

#define MAX_PART_TO_PRINT 20

struct Particle {
	int x,y;
	double weight;
	bool choosen;
};

struct Robot {
	int x,y;
};

Robot theRobot;  /*Variable that register the real position of the robot. Used to simulate the probabilities of the particles and compare if the algorithm succeeds.*/
int dx,dy;
bool showOutput=false;

/********************************************************************/
bool GetUserInput(int argc, char *argv[],int& spaceDimention, long& particleQuantity){
	/*Gets the parameters from input*/
	bool isOK = true;
	int robotX,robotY;

	if(argc < 4) 
	{
		cout << "Arguments:<spaceDimention> <particleQuantity> <showOutput>" << endl;
		cout << "spaceDimention : Space Matrix size [ spaceDimention X spaceDimention]" << endl;
		cout << "particleQuantity : number of particles to create" << endl;
		cout << "showOutput : 0|1 if show iteration output. Summary will always be shown." << endl;
		isOK = false;
	}
	else 
	{
		//get spaceDimention
		spaceDimention = atoi(argv[1]);
		if (spaceDimention <=0) 
		{
			cout << "Space size must be larger than 0" <<endl;
			isOK = false;
		}
		//get particleQuantity
		particleQuantity = atol(argv[2]);
		if (particleQuantity <= 0) 
		{
			cout << "Particles must be more than 0" <<endl;
			isOK = false;
		}

		//get showOutput
		int SO = atoi(argv[3]);
		if (SO ==0)	showOutput=false;
		else 
			if (SO ==1)	showOutput=true;
			else {
				cout << "showOutput should be 0 or 1" <<endl;
				isOK = false;
			}


		if (isOK){
			/*Read Robot Initial Location*/
			cout<<"Select robot position X (0.."<<spaceDimention-1<<"): ";  cin>>robotX;
			cout<<"Select robot position Y (0.."<<spaceDimention-1<<"): ";  cin>>robotY;
			
			theRobot.x=robotX;	theRobot.y=robotY;
		}
	}
	return isOK;
}



void updateRobotMovements(){
  	srand (time(NULL));/* initialize random seed: */

	if ( (rand()%100 ) > 50) {
		dx = 1;
		dy=0; 
	}   
	else {
		dx = 0;
		dy=1;
	}
}

/********************************************************************/
double estimateParticleWeight(Particle aParticle){
	double distance = sqrt( pow( theRobot.x - aParticle.x ,2) + pow( theRobot.y - aParticle.y  ,2) );

	if (distance == 0)
		return 1;
	else
		return 1/distance;
}



/********************************************************************/
void drawFirstParticleSet(Particle* particleSpace, int spaceDimention, long particleQuantity){
	/*Draw the first <particleQuantity> particles inside the matrix [spaceDimention X spaceDimention]*/
	int x,y;

  	srand (time(NULL));/* initialize random seed: */

	for ( int i=0; i<particleQuantity; i++){
 		x = rand() % spaceDimention;/* generate secret number between 1 and spaceDimention: */
 		y = rand() % spaceDimention;/* generate secret number between 1 and spaceDimention: */
 		particleSpace[i].x=x;
 		particleSpace[i].y=y;
 		particleSpace[i].weight=0;
 		particleSpace[i].choosen=false;
	}
}


/********************************************************************/
void printMatrixParticles (Particle* particleSpace, long particleQuantity, int spaceDimention){
	
	/*Prints the current particles to screen*/
	if (spaceDimention <= MAX_PART_TO_PRINT && showOutput){
		long particlesUnderRobot=0;

		for (int i=0; i < spaceDimention; i++){
			printf("ROW [%2i] ", i );
			for (int j=0; j < spaceDimention; j++){
				long count=0;

				for ( int k = 0; k< particleQuantity; k++ ){
					if ( particleSpace[k].x == i && particleSpace[k].y == j)	
						count++;
				}

				if ( theRobot.x == i && theRobot.y == j){
					cout << "R.";
					particlesUnderRobot=count;
				}
				else{
					if (count >0) printf("%2i ", count); 
					else cout <<"--";
				}

			}
			cout << endl;
		}
		printf("Particles Under Robot: %2i ", particlesUnderRobot); 
	}
}



/********************************************************************/
void printParticles (Particle* particleSpace, long particleQuantity, int spaceDimention){
	/*Prints the current particles to screen*/
	if (spaceDimention <= MAX_PART_TO_PRINT && showOutput){
		for ( int i = 0; i < particleQuantity; i++ )
		    cout<<"particle "<< i << " X:"<<particleSpace[i].x << " Y:" << particleSpace[i].y << " Weight:"<< particleSpace[i].weight<<" Choosen:"<<particleSpace[i].choosen<<endl;
	}
}


/********************************************************************/
void displayInitialConfiguration(int spaceDimention, long particleQuantity){
	system ("clear");
	cout << "Parallel Implementation w/OpenMP using "<< omp_get_num_procs() << " threads/CPUs" << endl;
	cout << "Simulation Configuration:" <<endl << "Space Dimention:" <<spaceDimention<<endl <<"Number of Particles:"<<particleQuantity<<endl;
	cout << "Robot initial position (x,y) = ("<<theRobot.x<<","<<theRobot.y<<")"<<endl; 
	cout <<"---------------------------------------------------"<<endl;
	cout <<"Press any key to start...\n";
	std::cin.ignore();
}




/********************************************************************/
void normalizeWeights(Particle* particleSpace, long particleQuantity, double normWeight){
	/*Normalizes the weight of the particles so Sum(Wi)==1.  If normWeight == 0, we calculate the normalization factor as the sum of the weight. */

	long i,k;
	#pragma omp parallel shared(particleQuantity,particleSpace,normWeight) firstprivate(i,k) 
	{
		#pragma omp single 
		{
			if (normWeight == 0){
				for (i=0; i<particleQuantity; i++){
					normWeight += particleSpace[i].weight;
				}
			}
		}

		#pragma omp for schedule(static)
		for (k=0; k<particleQuantity; k++){
			particleSpace[k].weight = (1 / normWeight) * particleSpace[k].weight;
		}
	}
}



/********************************************************************/
void printParticleProbability(Particle* particleSpace, long particleQuantity , int spaceDimention){
	/* Prints a summary of the particles and the probabylity for each one*/

	long      summaryQty = 0;										//Stores the numeber of particles.
	double   *particleProbability = new double[particleQuantity]; 	//Stores the probability
	long     *particleNumber = new long[particleQuantity]; 			//Stores the number of particles in an specific cell
	Particle *particleSummary = new Particle[particleQuantity];  	//Stores the particles 

	for (long i=0; i < particleQuantity; i++){
		particleProbability[i]=0;
		particleNumber[i]=0;
		particleSummary[i].x = -1;
		particleSummary[i].y = -1;
	}


	for (long i=0; i<particleQuantity; i++){
		bool found=false;
		long index=0;
		long j=0;
		while (j<summaryQty and !found ){  /*Search for the particle in the Summary array*/
			if ( particleSpace[i].x == particleSummary[j].x 
				 && 
				 particleSpace[i].y == particleSummary[j].y    ){
				found=true;
				index=j;
			}
			j++;
		}

		if (! found){ /*must add it*/
			index=summaryQty;
			particleSummary[index].x = particleSpace[i].x;
			particleSummary[index].y = particleSpace[i].y;
			summaryQty++;
		}
		particleProbability[index] += particleSpace[i].weight;
		particleNumber[index] ++;
	}

	system ("clear");
	cout << "Simulation Configuration:" <<endl << "Space Dimention:" <<spaceDimention<<endl <<"Number of Particles:"<<particleQuantity<<endl;
	cout << "Robot Final position (x,y) = ("<<theRobot.x<<","<<theRobot.y<<")"<<endl; 
	cout <<"---------------------------------------------------"<<endl;

	cout<<"Particle Summary:\n";
	for (long j=0; j < summaryQty; j++){
		printf("Position (%3i) (x,y)=(%5i,%5i) #Particles: %7i Probability of robot here: %6f%% \n", j,particleSummary[j].x,particleSummary[j].y,particleNumber[j],((double)particleProbability[j]*100) );
	}

	delete[] particleProbability;
	delete[] particleNumber;
	delete[] particleSummary;
}



/********************************************************************/
void estimateParticlesWeight(Particle* particleSpace, long  particleQuantity){
	double normWeight=0;

	for (long i=0;i<particleQuantity;i++){
		particleSpace[i].weight = estimateParticleWeight(particleSpace[i]);
		normWeight += particleSpace[i].weight;
	}

	normalizeWeights(particleSpace, particleQuantity, normWeight);

}



/********************************************************************/
void prefixCalculation(Particle *particles, double *prefix, long N)
{
	double *threadPrefix;
 
 	#pragma omp parallel shared(particles,prefix,threadPrefix,N)
    {
		const int ithread  = omp_get_thread_num();
        const int nthreads = omp_get_num_threads();
  
        #pragma omp single
        {                
        	/*
        	Stores the max of the particles.weight section calculated by each thread to be sum up to the subsequents parts of the prefix.
        	For K threads will store K+1 values. First one (i==0 )has value 0. Then one per thread. Value for thread K I wont need it but keep it for alignment */
        	threadPrefix = new double[nthreads+1];   
         	threadPrefix[0] = 0;
        }

        /*  Calculates a subarray prefix.  */
        double sum=0;
        #pragma omp for schedule(static)
        for (long i=0; i<N; i++) { 
            sum += particles[i].weight;
            prefix[i] = sum;
        }

        /* Each array will record the last value on threadPrefix at threadId +1 */
        threadPrefix[ithread+1]=sum;


        #pragma omp barrier  //Thread Syncronization.

        /* Make a prefix of the Thread Sums */
        #pragma omp single
        {
        	//for (int k=0;k<nthreads;k++){cout<< threadPrefix[k] << '|';}
        	//cout<<endl;

        	for (long k=1; k<nthreads; k++){ //Starting at position == 1
        		threadPrefix[k] += threadPrefix[k-1];
        	}

        	//for (int k=0;k<nthreads;k++){cout<< threadPrefix[k] << '|';}
        	//cout<<endl;

        }
		
        /* Now sumup the subarray prefixes with the subtotal calculated and get the final total prefix. */
		#pragma omp for schedule(static)
        for (long i=0; i<N; i++) { 
            /* Last value of partialPrefix will never be used as the last threadid (K) wrote in position K+1 */
            prefix[i] += threadPrefix[ithread];

        }
 	}
	
}

/********************************************************************/
void applyParticleFilters(Particle* particleSpace, int spaceDimention, long particleQuantity){


	long j,k;
	long currentNewParticle, myCurrentNewParticle;
	double normWeight=0;
	double randomProbability;
	double *cumulativeWeight     = new double[particleQuantity];	// Stores the prefix of the particle weight

	
	srand (time(NULL));/* initialize random seed: */

	randomProbability = (double)(((rand() << 15) + rand()) & ((1 << 24) - 1)) / (1 << 24);
    randomProbability = randomProbability * ( (double)1 / (double) particleQuantity );

    /*Calculating the cumulative weight -- Prefix  -- Parallel Algorithm*/
    prefixCalculation(particleSpace,cumulativeWeight,particleQuantity);

	
	j=0;
	currentNewParticle=0;
	myCurrentNewParticle=0;

	///#pragma omp parallel for
	#pragma omp parallel shared(particleQuantity,randomProbability,cumulativeWeight,particleSpace,currentNewParticle,normWeight) firstprivate(j,k) private(myCurrentNewParticle)  
	{
		#pragma omp for schedule(static) 
		for (k=0; k < particleQuantity; k++){ /*For each particle*/
			
			myCurrentNewParticle=k;

			double uk = (double) randomProbability +  ( (double) k / particleQuantity ) ;

			while (uk > cumulativeWeight[j]) {
				j ++;
			}

			//APPLY RANDOM MOVEMENT.

			particleSpace[myCurrentNewParticle].x = particleSpace[j].x + dx;
			particleSpace[myCurrentNewParticle].y = particleSpace[j].y + dy;
			particleSpace[myCurrentNewParticle].choosen = true;

			//Boundary control.
			if (particleSpace[myCurrentNewParticle].x < 0 || 
				particleSpace[myCurrentNewParticle].x > spaceDimention-1)	{
				particleSpace[myCurrentNewParticle].x = particleSpace[myCurrentNewParticle].x + (dx * - 2);	
			}
			if (particleSpace[myCurrentNewParticle].y < 0 || 
				particleSpace[myCurrentNewParticle].y > spaceDimention-1)	{
				particleSpace[myCurrentNewParticle].y = particleSpace[myCurrentNewParticle].y + (dy * - 2);	
			}

			particleSpace[myCurrentNewParticle].weight = estimateParticleWeight(particleSpace[myCurrentNewParticle]);
		}
	}
	normWeight=0; //Calculation in sequence is better than having a critical secction for update in parallel

	normalizeWeights(particleSpace,particleQuantity, normWeight);
 
	delete[] cumulativeWeight;
}





//********************************************************************
// Main Program
//********************************************************************
int main(int argc, char *argv[])
{
	/*********************/
	/**** Variables      */
	int spaceDimention;
	long particleQuantity;
	Particle *particleSpace;
	float runtime;

	if ( GetUserInput(argc,argv,spaceDimention,particleQuantity) == false ) return 1;  
	displayInitialConfiguration(spaceDimention,particleQuantity);  /*Prints initial configuration.*/


	runtime = omp_get_wtime();


	updateRobotMovements();

	particleSpace = new Particle[particleQuantity];  /*Allocates memory for the particles.*/ 

	drawFirstParticleSet(particleSpace,spaceDimention,particleQuantity); /*draw the first set of particles sparsed on the grid.*/

	estimateParticlesWeight(particleSpace,particleQuantity);

	printMatrixParticles(particleSpace,particleQuantity,spaceDimention); /*prints the inital particles (depending on size)*/
	printParticles(particleSpace,particleQuantity,spaceDimention); /*prints the inital particles (depending on size)*/


	/*ITERATE!!!*/

	int iterationsToDo=0;

	//for (int i=0; i <300; i++){  ///////////TODO!!!! Resolve when to stop.
	while (iterationsToDo < (spaceDimention*0.5) ){

		//cout <<iterationsToDo<<endl;

		if (iterationsToDo % (spaceDimention/2) == 0 )
			updateRobotMovements();

		if (showOutput){	cout <<"Iteration No."<<iterationsToDo<<" Press any key to continue...\n";	std::cin.ignore();	}

		theRobot.x = theRobot.x + dx;	theRobot.y = theRobot.y + dy;	/*Robot Moves.*/

		//Boundary control.
		if (theRobot.x < 0 || theRobot.x > spaceDimention-1)	{theRobot.x = theRobot.x + (dx * -2);	dx = dx * -1;}
		if (theRobot.y < 0 || theRobot.y > spaceDimention-1)	{theRobot.y = theRobot.y + (dy * -2);	dy = dy * -1;}
		
		applyParticleFilters(particleSpace, spaceDimention,particleQuantity);
		
		printMatrixParticles(particleSpace,particleQuantity,spaceDimention); /*prints the inital particles (depending on size)*/

		iterationsToDo++;

	}

	cout<<endl;

	printMatrixParticles(particleSpace,particleQuantity,spaceDimention); /*prints the inital particles (depending on size)*/

	cout<<endl<<endl;

	printParticleProbability(particleSpace, particleQuantity, spaceDimention );

	cout<<endl<<endl;
	runtime = omp_get_wtime() - runtime;
	cout<< "Program runs in " << setiosflags(ios::fixed) << setprecision(2) << runtime << " seconds\n"; 

	delete[] particleSpace;

	cout <<"---------------------------------------------------"<<endl;
	cout <<"-----            Simulation Ended              ----"<<endl;
	cout <<"---------------------------------------------------"<<endl;
	return 0;
}