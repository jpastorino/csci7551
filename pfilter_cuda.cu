//-----------------------------------------------------------------------
// Particle Filters Simulation   ---  SEQUENTIAL ALGORITHM
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
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>

using namespace std;

//-----------------------------------------------------------------------
//   Data Structures and Constants
//-----------------------------------------------------------------------

#define THREADS_PER_BLOCK 1024

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


/********************************************************************/
void updateRobotMovements(int &dx, int &dy){
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


/**************************** Run in CPU ****************************************/
double cpuEstimateParticleWeight(Particle aParticle){
	/* Estimates the weigh of a particle being in the robots position. */
	double distance = sqrt( pow( theRobot.x - aParticle.x ,2) + pow( theRobot.y - aParticle.y  ,2) );

	if (distance == 0)
		return 1;
	else
		return 1/distance;
}

/*****************************  Run in CPU ***************************************/
void estimateParticlesWeight(Particle* particleSpace, long  particleQuantity){
	/*For first initialization calculates the particle weight.  Could be improved in parallel*/

	double normWeight=0;

	for (long i=0;i<particleQuantity;i++){
		particleSpace[i].weight = cpuEstimateParticleWeight(particleSpace[i]);
		normWeight += particleSpace[i].weight;
	}

	
	for (long i=0; i<particleQuantity; i++){
		normWeight += particleSpace[i].weight;
	}

	for (long i=0; i<particleQuantity; i++){
		particleSpace[i].weight = (1 / normWeight) * particleSpace[i].weight;
	}

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
double calculateRandomProbability(long particleQuantity){
	double randomProbability=0;
	srand (time(NULL));/* initialize random seed: */
	randomProbability = (double)(((rand() << 15) + rand()) & ((1 << 24) - 1)) / (1 << 24);
	randomProbability = randomProbability * ( (double)1 / (double) particleQuantity );
	return randomProbability;
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
	cout << "Simulation Configuration:" <<endl << "Space Dimention:" <<spaceDimention<<endl <<"Number of Particles:"<<particleQuantity<<endl;
	cout << "Robot initial position (x,y) = ("<<theRobot.x<<","<<theRobot.y<<")"<<endl; 
	cout <<"---------------------------------------------------"<<endl;
	cout <<"Press any key to start...\n";
	std::cin.ignore();
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


/******************************** RUNS IN GPU ************************************/
__global__ void normalizeWeights(Particle* particleSpace, long particleQuantity){  
	/*Normalizes the weight of the particles so Sum(Wi)==1.  we calculate the normalization factor as the sum of the weight. */
	unsigned int threadId = ( blockIdx.x * blockDim.x ) + threadIdx.x;

	__shared__ double normWeight;

	if (threadId == 0){
		//do SERIAL Reduction
		normWeight=0;
		for (long i=0; i<particleQuantity; i++)
			normWeight+=particleSpace[i].weight;
	}
	__syncthreads();

	if (threadId<particleQuantity)
		particleSpace[threadId].weight = (1 / normWeight) * particleSpace[threadId].weight;
	
}


/******************************** RUNS IN GPU ************************************/
__global__ void prefixCalculation(Particle *particles, double *prefix, long N){
	unsigned int threadId = ( blockIdx.x * blockDim.x ) + threadIdx.x;

	if (threadId == 0){  //Dummy in serial to check everything.
		//*Calculating the cumulative weight -- Prefix ****
		prefix[0]=particles[0].weight;
		for (long k=1; k<N; k++){	
			prefix[k] = prefix[k-1]+particles[k].weight;
		}
	} 
}




/**************************** Run in GPU ****************************************/
__device__ double estimateParticleWeight(Particle aParticle, Robot theRobot){
	/* Estimates the weigh of a particle being in the robots position. */
	double distance = sqrt( (double)(pow( (double)(theRobot.x - aParticle.x) ,2) + pow( (double)(theRobot.y - aParticle.y)  ,2)) );

	if (distance == 0)
		return 1;
	else
		return 1/distance;
}

/******************************** RUNS IN GPU ************************************/
__global__ void applyParticleFilters(Particle* particleSpace, double* cumulativeWeight, int spaceDimention, long particleQuantity, int dx, int dy, Robot theRobot, double randomProbability){

	unsigned int threadId = ( blockIdx.x * blockDim.x ) + threadIdx.x;
	long j=0;
	 
	if (threadId < particleQuantity)   //Control access in particle array boundaries
	{
 		double uk = (double) randomProbability +  ( (double) threadId / particleQuantity ) ;

		while (uk > cumulativeWeight[j] && j<particleQuantity-1) {	j ++;	}

		//APPLY RANDOM MOVEMENT.
		particleSpace[threadId].x = particleSpace[j].x + dx;
		particleSpace[threadId].y = particleSpace[j].y + dy;
		particleSpace[threadId].choosen = true;

		//Boundary control.
		if (particleSpace[threadId].x < 0 || 
			particleSpace[threadId].x > spaceDimention-1)	{
			particleSpace[threadId].x = particleSpace[threadId].x + (dx * - 2);	
		}
		if (particleSpace[threadId].y < 0 || 
			particleSpace[threadId].y > spaceDimention-1)	{
			particleSpace[threadId].y = particleSpace[threadId].y + (dy * - 2);	
		}

		particleSpace[threadId].weight = estimateParticleWeight(particleSpace[threadId], theRobot);
	}
}





//********************************************************************
// Main Program
//********************************************************************
int main(int argc, char *argv[])
{
	/*********************/
	/**** Variables      */
	int dx,dy;						//Robot movements
	int spaceDimention;
	long particleQuantity;
	Particle *particleSpace;
	Particle *d_particleSpace;
		double *cumulativeWeight;  

	double *d_cumulativeWeight;  
	double randomProbability;


	float runtime;

	if ( GetUserInput(argc,argv,spaceDimention,particleQuantity) == false ) return 1;  

	//Configure GPU Thread distribution.
	int numOfBlocks = particleQuantity / THREADS_PER_BLOCK + ((particleQuantity%THREADS_PER_BLOCK)?1:0);

	displayInitialConfiguration(spaceDimention,particleQuantity);  /*Prints initial configuration.*/
	cout <<"RUNNIN CUDA WITH BLOCK:"<<numOfBlocks<<"  AND THREADS:"<<THREADS_PER_BLOCK<<endl;

	runtime = clock()/(float)CLOCKS_PER_SEC;
	
	//Initialize the robot movements.
	updateRobotMovements(dx,dy);

	particleSpace = new Particle[particleQuantity];  /*Allocates memory for the particles.*/ 
	cumulativeWeight=new double[particleQuantity];
	
	//Allocate memory on device for the particles and prefix.
	cudaMalloc((void**)&d_particleSpace,    particleQuantity*sizeof(Particle));
	cudaMalloc((void**)&d_cumulativeWeight, particleQuantity*sizeof(double));  //Prefix Calculation.
 
	drawFirstParticleSet(particleSpace,spaceDimention,particleQuantity); /*draw the first set of particles sparsed on the grid.*/

	estimateParticlesWeight(particleSpace,particleQuantity);//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	printMatrixParticles(particleSpace,particleQuantity,spaceDimention); /*prints the inital particles (depending on size)*/
	printParticles(particleSpace,particleQuantity,spaceDimention);       /*prints the inital particles (depending on size)*/


	//Copy Particles to the GPU. After that all procedure will occur on the GPU until finished.
	cudaMemcpy(d_particleSpace, particleSpace, particleQuantity*sizeof(Particle), cudaMemcpyHostToDevice);

	/*ITERATE!!!*/
	int iterationsToDo=0;

	while (iterationsToDo < (spaceDimention*0.5) ){

		if (iterationsToDo % (spaceDimention/2) == 0 )
			updateRobotMovements(dx,dy);

		if (showOutput){	cout <<"Iteration No."<<iterationsToDo<<" Press any key to continue...\n";	std::cin.ignore();	}

		theRobot.x = theRobot.x + dx;	theRobot.y = theRobot.y + dy;	/*Robot Moves.*/

		//Boundary control.
		if (theRobot.x < 0 || theRobot.x > spaceDimention-1)	{theRobot.x = theRobot.x + (dx * -2);	dx = dx * -1;}
		if (theRobot.y < 0 || theRobot.y > spaceDimention-1)	{theRobot.y = theRobot.y + (dy * -2);	dy = dy * -1;}
		
		prefixCalculation<<<numOfBlocks, THREADS_PER_BLOCK>>>(d_particleSpace, d_cumulativeWeight, particleQuantity); 	//Lets Calculate Prefix, 
		cudaThreadSynchronize();  //SyncThreads to continue with the prefix calculated.

		randomProbability = calculateRandomProbability(particleQuantity);
		applyParticleFilters<<<numOfBlocks,THREADS_PER_BLOCK>>> (d_particleSpace, d_cumulativeWeight, spaceDimention, particleQuantity,dx,dy, theRobot, randomProbability);		//Lets apply the particle filter
		cudaThreadSynchronize();  //Sync threads to continue next step
		
		normalizeWeights<<<numOfBlocks,THREADS_PER_BLOCK>>> (d_particleSpace, particleQuantity);//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		cudaThreadSynchronize();  //Sync threads to continue to next iteration.

		//Could not print as data is in GPU,otherwise had to copyback:	printMatrixParticles(particleSpace,particleQuantity,spaceDimention); /*prints the inital particles (depending on size)*/

		iterationsToDo++;
	}
	//Copy back the data from GPU to CPU.
	cudaMemcpy(particleSpace, d_particleSpace, particleQuantity*sizeof(Particle), cudaMemcpyDeviceToHost);

	
	//////// Display Information.
	cout<<endl;
	printMatrixParticles(particleSpace,particleQuantity,spaceDimention); 
	
	cout<<endl<<endl;
	printParticleProbability(particleSpace, particleQuantity, spaceDimention );

	cout<<endl<<endl;
	runtime = clock()/(float)CLOCKS_PER_SEC - runtime;
	cout<< "Program runs in " << setiosflags(ios::fixed) << setprecision(2) << runtime << " seconds\n"; 

	cudaFree(d_particleSpace);
	cudaFree(d_cumulativeWeight);
	delete[] particleSpace;

	cout <<"---------------------------------------------------"<<endl;
	cout <<"-----            Simulation Ended              ----"<<endl;
	cout <<"---------------------------------------------------"<<endl;
	return 0;
}