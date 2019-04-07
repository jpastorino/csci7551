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


/********************************************************************/
void updateRobotMovements(){
	/* Update the way the robot moves, so the particles.*/
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
	/* Simulation. estimates the probability of the particle being the robots based on what it senses.
	   In this simulation, that probability is just the inverse of the euclidean distance to the robot. 
	   Other Simulation methods colud be incorporated, such as the use of beacons. */
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
	cout << "Simulation Configuration:" <<endl << "Space Dimention:" <<spaceDimention<<endl <<"Number of Particles:"<<particleQuantity<<endl;
	cout << "Robot initial position (x,y) = ("<<theRobot.x<<","<<theRobot.y<<")"<<endl; 
	cout <<"---------------------------------------------------"<<endl;
	cout <<"Press any key to start...\n";
	std::cin.ignore();
}




/********************************************************************/
void normalizeWeights(Particle* particleSpace, long particleQuantity, double normWeight){
	/*Normalizes the weight of the particles so Sum(Wi)==1.  If normWeight == 0, we calculate the normalization factor as the sum of the weight. */
	if (normWeight == 0){
		for (long i=0; i<particleQuantity; i++){
			normWeight += particleSpace[i].weight;
		}
	}

	for (long i=0; i<particleQuantity; i++){
		particleSpace[i].weight = (1 / normWeight) * particleSpace[i].weight;
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
	/*For the first run this function just estimates all particles weights.*/
	double normWeight=0;

	for (long i=0;i<particleQuantity;i++){
		particleSpace[i].weight = estimateParticleWeight(particleSpace[i]);
		normWeight += particleSpace[i].weight;
	}

	normalizeWeights(particleSpace, particleQuantity, normWeight);
}




/********************************************************************/
void applyParticleFilters(Particle* particleSpace, int spaceDimention, long particleQuantity){
	/* This is the Particle filter algorithm.
	   Includes the prefix calculation for cumulative weights.
	   The resampligm.
	   and the weight normalization.*/

	long j;
	long currentNewParticle;
	double normWeight=0;
	double randomProbability;
	double *cumulativeWeight     = new double[particleQuantity];	// Stores the prefix of the particle weight

	
	srand (time(NULL));/* initialize random seed: */

	randomProbability = (double)(((rand() << 15) + rand()) & ((1 << 24) - 1)) / (1 << 24);
    randomProbability = randomProbability * ( (double)1 / (double) particleQuantity );



	 
	/*Calculating the cumulative weight -- Prefix */
	cumulativeWeight[0]=particleSpace[0].weight;
	for (long k=1; k<particleQuantity; k++){	
		cumulativeWeight[k] = cumulativeWeight[k-1]+particleSpace[k].weight;
	}

	j=0;
	currentNewParticle=0;


	double uk =0;
	for (long k=0; k < particleQuantity; k++){ /*For each particle*/

		uk = (double) randomProbability +  (( (double) k )/ particleQuantity ) ;

		while (uk > cumulativeWeight[j]){
			j ++;
		}

		//APPLY RANDOM MOVEMENT.
		particleSpace[currentNewParticle].x = particleSpace[j].x + dx;
		particleSpace[currentNewParticle].y = particleSpace[j].y + dy;
		particleSpace[currentNewParticle].choosen = true;

		//Boundary control.
		if (particleSpace[currentNewParticle].x < 0 || particleSpace[currentNewParticle].x > spaceDimention-1)	{
			particleSpace[currentNewParticle].x = particleSpace[currentNewParticle].x + (dx * - 2);	
		}
		if (particleSpace[currentNewParticle].y < 0 || particleSpace[currentNewParticle].y > spaceDimention-1)	{
			particleSpace[currentNewParticle].y = particleSpace[currentNewParticle].y + (dy * - 2);	
		}

		particleSpace[currentNewParticle].weight = estimateParticleWeight(particleSpace[currentNewParticle]);
		normWeight += particleSpace[currentNewParticle].weight;

		currentNewParticle++;
	}

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


	runtime = clock()/(float)CLOCKS_PER_SEC;


	updateRobotMovements();

	particleSpace = new Particle[particleQuantity];  						/*Allocates memory for the particles.*/ 

	drawFirstParticleSet(particleSpace,spaceDimention,particleQuantity); 	/*draw the first set of particles sparsed on the grid.*/

	estimateParticlesWeight(particleSpace,particleQuantity);

	printMatrixParticles(particleSpace,particleQuantity,spaceDimention); 	/*prints the inital particles (depending on size)*/
	printParticles(particleSpace,particleQuantity,spaceDimention); 			/*prints the inital particles (depending on size)*/


	/*ITERATE!!!*/

	int iterationsToDo=0;

	while (iterationsToDo < (spaceDimention*0.5) ){

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
	runtime = clock()/(float)CLOCKS_PER_SEC - runtime;
	cout<< "Program runs in " << setiosflags(ios::fixed) << setprecision(2) << runtime << " seconds\n"; 

	delete[] particleSpace;

	cout <<"---------------------------------------------------"<<endl;
	cout <<"-----            Simulation Ended              ----"<<endl;
	cout <<"---------------------------------------------------"<<endl;
	return 0;
}