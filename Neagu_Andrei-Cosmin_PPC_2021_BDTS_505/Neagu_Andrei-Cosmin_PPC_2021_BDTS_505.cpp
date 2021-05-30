#include <iostream>
#include <iomanip>
#include <ctime>
#include <string>
#include "mpi.h"
#include "stdafx.h"

using namespace std;

/**
 * @file Neagu_Andrei-Cosmin_PPC_2021_BDTS_505.cpp
 *
 * @brief Circuit Satisfiability using MPI // mpiexec -n [numberOfProcessors] Neagu_Andrei-Cosmin_PPC_2021_BDTS_505.exe
 *
 * @author Neagu Andrei-Cosmin - Baze de date si tehnologii software - Anul II, Grupa 505
 *
 */

class Circuit {
public:
	Circuit();
	virtual void createBinaryArray(int integerValue, int binaryArray[], int arraySize);
	virtual bool getCircuitOutput(int binaryArray[], int numberOfLogicalVariables);
};

class TimeUtil {
public:
	TimeUtil();
	virtual string getElapsedTime(clock_t startTime, int processId);
};

Circuit::Circuit() {}

/**
 * Creates a binary array by converting an integer value.
 *
 * @param integerValue - the integer value.
 * @param binaryArray -  the binary array.
 * @param arraySize - the size of the binary array.
 */
void Circuit::createBinaryArray(int integerValue, int binaryArray[], int arraySize) {
	for (int i = arraySize - 1; i >= 0; i--) {
		binaryArray[i] = integerValue % 2;
		integerValue /= 2;
	}
}

/**
 * Obtains the output of the circuit, given an array of binary values by performing logical operations on those.
 * The method has three predefined circuit logical operations.
 *
 * @param binaryArray - the binary array containing our logical variables.
 * @param numberOfLogicalVariables - the array size / number of logical variables.
 * @return a boolean value - the output of the circuit, either 1 (true) or 0 (false).
 */
bool Circuit::getCircuitOutput(int circuit[], int numberOfLogicalVariables) {
	const bool circuitOneOutput =
		((((circuit[0] || circuit[1]) && (!circuit[1] || !circuit[3])) && ((circuit[2] || circuit[3]) && (!circuit[3] || !circuit[4])))
		&&
		(((circuit[4] || !circuit[5]) && (circuit[5] || circuit[6])) && ((circuit[5] || !circuit[6]) && (circuit[7] || !circuit[8]))))
		&&
		((((circuit[8] || circuit[9]) && (circuit[8] || !circuit[9])) && ((!circuit[9] || !circuit[10]) && (circuit[10] || circuit[11])))
		&&
		(((circuit[9] || circuit[11]) && (circuit[12] || circuit[13])) && ((!circuit[7] || !circuit[13]) && ((circuit[14] || circuit[15]) && (circuit[6] || !circuit[15])))));

	const bool circuitTwoOutput =
		(((((circuit[0] && circuit[1]) && (circuit[2] && circuit[3])) && ((circuit[4] && circuit[5]) && (circuit[6] && circuit[7]))))
		&&
		((((circuit[8] && circuit[9]) && (circuit[10] && circuit[11])) && ((circuit[12] && circuit[13]) && (circuit[14] && circuit[13])))));

	const bool circuitThreeOutput =
		((((circuit[0] && circuit[1]) && (!(circuit[2] || circuit[3]) && (circuit[4] && circuit[5]))) ^ !circuit[5])
		&&
		((circuit[6] && (!circuit[6] && circuit[5])) ^ !(circuit[7] || (circuit[8] ^ circuit[9]))));

	switch (numberOfLogicalVariables) {
	case 16:
		return circuitOneOutput;
	case 15:
		return circuitTwoOutput;
	case 10:
		return circuitThreeOutput;
	default:
		throw exception("The inserted values are not supported by this method!");
	}
}

TimeUtil::TimeUtil() {};

/**
 * Calculates the elapsed time from a specific start time until the moment the method is being called for a process ID.
 *
 * @param startTime - the start time of the process.
 * @param processId - the ID of the process.
 * @return a string - the message including the time elapsed measured in seconds for the specified process.
 */
string TimeUtil::getElapsedTime(clock_t startTime, int processId) {
	double elapsedTime = double(clock() - startTime) / CLOCKS_PER_SEC;
	string message = "Process " + to_string(processId) + " finished after " + to_string(elapsedTime) + " seconds.\n";
	return message;
}

void main(int argc, char* argv[]) {
	const int logicalVariablesSize = 16;
	int binaryArray[logicalVariablesSize];

	Circuit* circuit = new Circuit();
	TimeUtil* timeUtil = new TimeUtil();

	MPI_Init(&argc, &argv);

	int processId;
	MPI_Comm_rank(MPI_COMM_WORLD, &processId);

	int numberOfProcesses;
	MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcesses);

	clock_t startTimeOfCurrentProcess = clock();

	int initialLowerBoundary = 0;
	int initialUpperBoundary = pow(2, logicalVariablesSize);

	if (processId == 0) {
		cout << "\n* Number of processes: " << numberOfProcesses << " *\n";
		cout << "\n* Number of logical variables: " << logicalVariablesSize << " *\n";;
		cout << "\n* Number of binary arrays: " << initialUpperBoundary << " *\n";
	}

	int lowerBoundary = ((numberOfProcesses - processId) * initialLowerBoundary + (processId) * initialUpperBoundary) / (numberOfProcesses);
	int upperBoundary = ((numberOfProcesses - processId - 1) * initialLowerBoundary + (processId + 1) * initialUpperBoundary) / (numberOfProcesses);
	int numberOfSolutionsForCurrentProcess = 0;

	cout << "\n* Processor with ID " << processId << " starts now and will cover values between: [" << lowerBoundary << "," << upperBoundary << ")\n";

	for (int currentInteger = lowerBoundary; currentInteger < upperBoundary; currentInteger++) {
		circuit->createBinaryArray(currentInteger, binaryArray, logicalVariablesSize);
		bool circuitOutput = circuit->getCircuitOutput(binaryArray, logicalVariablesSize);
		if (circuitOutput) {
			numberOfSolutionsForCurrentProcess++;
			cout << setw(10) << numberOfSolutionsForCurrentProcess << setw(10) << currentInteger << setw(10);
			for (int currentBinaryValue = 0; currentBinaryValue < logicalVariablesSize; currentBinaryValue++) {
				cout << binaryArray[currentBinaryValue];
			}
			cout << "\n";
		}		
	}

	cout << timeUtil->getElapsedTime(startTimeOfCurrentProcess, processId);

	int numberOfSolutionsForAllProcesses;
	MPI_Reduce(&numberOfSolutionsForCurrentProcess, &numberOfSolutionsForAllProcesses, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	if (processId == 0) {
		cout << "\n * The program has found: " << numberOfSolutionsForAllProcesses << " satisfiable solutions. *\n\n";
	}

	MPI_Finalize();

	return;
}
