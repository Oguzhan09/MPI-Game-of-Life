/**
 * Student Name: Oğuzhan Kırlar
 * Student Number: 2014400195
 * Compile Status: Compiling
 * Program Status: Working
 */


#include <mpi.h>
#include <fstream>
#include <algorithm>
#include <iostream>

using namespace std;

int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    // Find out rank, size
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    //the number of slave processes.
    int N = world_size-1;
    int s = 360/N;
    ifstream inputFile(argv[1]);        // Input file stream object

    if(world_rank == 0){
        //input array.
        int X[360][360];
        //reading the file and assigning corresponding array element.
        int current_number = 0;
        for (int i = 0; i < 360; i++) {
            for (int j = 0; j < 360; j++) {
                inputFile >> current_number;
                X[i][j] = (current_number);
            }
        }

        int sending_array [s][360];
        //broadcasting to slave processes.
        for(int k = 1 ; k <= N ; ++k){
            for (int i = 0; i < s; i++) {
                for (int j = 0; j < 360; ++j) {
                    sending_array[i][j] = X [(k-1)*(s)+i][j];
                }
            }
            MPI_Send(sending_array, 360*s, MPI_INT, k, 0, MPI_COMM_WORLD);
        }
        // Close the file.
        inputFile.close();
        // output array after processing.
        int outputArr[360][360];
        //recieves subarrays and join them according to their process ids.
        for (int k = 1; k <= N; ++k) {
            int* subarr = NULL;
            subarr = (int *)malloc(sizeof(int) * 360*360/N);
            MPI_Recv(subarr, 360*360/N, MPI_INT, k, k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i = 0; i < 360 / N; ++i) {
                for (int j = 0; j < 360; ++j) {
                    outputArr[(k-1)*(360/N)+i][j] = subarr[i*360+j];
                }
            }
        }
        ofstream outputfile;
        outputfile.open(argv[2]);
        // writing to output file.
        for (int n = 0; n < 360; ++n) {
            for (int i = 0; i < 360; ++i) {
                outputfile << outputArr[n][i] << " ";
            }
            outputfile << endl;
        }
        outputfile.close();
    }else{
        int sub_arr [s][360];
        //gets the corresponding subarray from master process.
        MPI_Recv(sub_arr, 360*s, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        int generations = atoi(argv[3]);
        int Z[s][360];
        copy(&sub_arr[0][0], &sub_arr[0][0] + s * 360, &Z[0][0]);


        int to_down[360];	int to_up[360]; int from_down[360]; int from_up[360]; //arrays to send and to receive
        for (int k = 1; k <=generations; ++k) {
            if (world_rank!=world_size-1) {// all except for last send down
                for (int j=0; j<360; j++) to_down[j]=sub_arr[s-1][j];
                MPI_Send(&to_down, 360, MPI_INT, world_rank+1, 1, MPI_COMM_WORLD);
            } else {
                for (int k=0; k<360; k++) from_down[k]=0;  // last one generates empty stripe "from down"
            }

            if (world_rank != 1) { // all except for first receive from up
                MPI_Recv(&from_up, 360, MPI_INT, world_rank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            } else {
                for (int k=0; k<360; k++) from_up[k]=0; // first one generates empty line "from up"
            }

            if (world_rank != 1) {// all except for first send up
                for (int j=0; j<360; j++) to_up[j]=sub_arr[0][j];
                MPI_Send(&to_up, 360, MPI_INT, world_rank-1, 1, MPI_COMM_WORLD);
            }

            if (world_rank != world_size-1) {// all except for last receive from down
                MPI_Recv(&from_down, 360, MPI_INT, world_rank+1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            //COUNTING NEIGHBORS
            int sum=0; // sum of neighbours
            for (int x=0; x<s; x++) {
                for (int y=0; y<360; y++) {
                    if (x==0 && y==0) //upper-left cell
                        sum = sub_arr[x+1][y]+sub_arr[x+1][y+1]+sub_arr[0][y+1]+from_up[0]+from_up[1];
                    else if (x==0 && y==359) //upper-right cell
                        sum = sub_arr[x][y-1]+sub_arr[x+1][y-1]+sub_arr[x+1][y]+from_up[359]+from_up[358];
                    else if (x==s-1 && y==0) //lower-left cell
                        sum = sub_arr[x][y+1]+sub_arr[x-1][y+1]+sub_arr[x-1][y]+from_down[0]+from_down[1];
                    else if (x==s-1 && y==359) //lower-right cell
                        sum = sub_arr[x-1][y]+sub_arr[x-1][y-1]+sub_arr[x][y-1]+from_down[359]+from_down[358];
                    else {
                        if (y==0) // leftmost line, not corner
                            sum=sub_arr[x-1][y]+sub_arr[x-1][y+1]+sub_arr[x][y+1]+sub_arr[x+1][y+1]+sub_arr[x+1][y];
                        else if (y==359) //rightmost line, not corner
                            sum=sub_arr[x-1][y]+sub_arr[x-1][y-1]+sub_arr[x][y-1]+sub_arr[x+1][y-1]+sub_arr[x+1][y];
                        else if (x==0) //uppermost line, not corner
                            sum=sub_arr[x][y-1]+sub_arr[x+1][y-1]+sub_arr[x+1][y]+sub_arr[x+1][y+1]+sub_arr[x][y+1]+from_up[y-1]+from_up[y]+from_up[y+1];
                        else if (x==s-1) //lowermost line, not corner
                            sum=sub_arr[x-1][y-1]+sub_arr[x-1][y]+sub_arr[x-1][y+1]+sub_arr[x][y+1]+sub_arr[x][y-1]+from_down[y-1]+from_down[y]+from_down[y+1];
                        else //general case, any cell within
                            sum=sub_arr[x-1][y-1]+sub_arr[x-1][y]+sub_arr[x-1][y+1]+sub_arr[x][y+1]+sub_arr[x+1][y+1]+sub_arr[x+1][y]+sub_arr[x+1][y-1]+sub_arr[x][y-1];
                    }
                    //PUT THE NEW VALUE OF A CELL
                    if (sub_arr[x][y]==0 && sum==3) Z[x][y]=1;
                    else if (sub_arr[x][y]==1 && sum>3) Z[x][y]=0;
                    else if (sub_arr[x][y]==1 && sum<2) Z[x][y]=0;
                }
            }
            // copy
            for (int x=0; x<s; x++)
                for (int y=0; y<360; y++)
                    sub_arr[x][y]=Z[x][y];
        }
        MPI_Send(sub_arr, 360*s, MPI_INT, 0, world_rank, MPI_COMM_WORLD);
    }
    MPI_Finalize();
}