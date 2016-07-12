#include "mpi.h" 
#include <stdio.h> 
#include <string.h>
#include <stdlib.h> 

#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include <sys/time.h>

typedef struct { 
   double* values;   //Array of node values in one row
   long int* col;    //Array of col_inds in one row
   long int vcount;  //Number of nodes in one row
   int rowId;
} rowData;

typedef struct {
   rowData* rowary;  //Stores data for each row
   int rcount;       //Number of rows
} row;

struct pNode {
   int val;
   struct pNode *next;
   struct pNode *last;
};

struct cNode {
   int val;
   int numofReq;
   struct pNode *reqLink;
   struct cNode *next;
};

typedef struct {
   double val;
   int rowId;
}  mValue;

int main(int argc, char **argv) { 
   int numprocs; 
   int myid; 
   MPI_Status stat;
   int c;
   int i = 0;

   int iterationCount = 0;
   struct timeval start, end;

   MPI_Init(&argc,&argv); 
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs); 
   MPI_Comm_rank(MPI_COMM_WORLD,&myid); 

   int *multiplyindex; 
   int indexsize;
   int numofRows = 0;
   double epsilon = 0;

   double *recvMatrixData;
   int recvSize;

   char finish = 1;

   if(myid == 0) { 
      FILE *pFile = fopen(argv[1], "r");

      if (pFile==NULL) {
         printf("File error");
         fputs ("File error", stderr); 
         exit(1);
      }

      int size;
      fseek (pFile, 0, SEEK_END);
      size = ftell(pFile);
      rewind (pFile);

      char st[20];
      int readRowsPtrs = 0;
      long int* rownum = (long int *) malloc(size * sizeof(long int));

      printf("READING LAST LINE\n");
      while (fscanf(pFile, "%20s", st) != EOF) {
         if (readRowsPtrs == 1) {
            long int v;
            v = strtol(st, NULL, 0);
            rownum[i] = v;
            i++;
         }

         if (strcmp(st, "row_ptrs:") == 0) {
            readRowsPtrs = 1;
         }
      }
      rewind(pFile);
      printf("LAST LINE READ\n");

      row matrixdata;
      matrixdata.rowary = malloc(i * sizeof(rowData));
      matrixdata.rcount = i-1;

      int *outGoingEdgeCount = (int *) malloc(matrixdata.rcount * sizeof(int));

      //Initialize matrix 
      for (c=0; c < matrixdata.rcount; c++) { //i is the total amount of rows in the matrix
         long int rowSize;
         rowSize = rownum[c+1] - rownum[c];

         matrixdata.rowary[c].values = malloc(rowSize * sizeof(double));
         matrixdata.rowary[c].col = malloc(rowSize * sizeof(long int));
         matrixdata.rowary[c].vcount = rowSize;

         outGoingEdgeCount[c] = 0;
      }
      free(rownum);
      
      i = 0;
      c = 0;      //Matrix row incrementor
      char str[20];
      int newline = 0;
      char ch;
      rewind(pFile);

      printf("READING col_inds LINE\n");
      readRowsPtrs = 0;
      while (fscanf(pFile, "%20s", st) != EOF) {
         if (strcmp(st, "row_ptrs:") == 0) {
            break;
         }

         if (readRowsPtrs == 1) {
            int v;
            v = strtol(st, NULL, 0);
            outGoingEdgeCount[v] = outGoingEdgeCount[v] + 1;
         }

         if (strcmp(st, "col_inds:") == 0) {
            readRowsPtrs = 1;
         }
      }
      rewind(pFile);
      printf("col_inds LINE READ\n");

      //Read matrix file and populate matrix data structure
      printf("POPULATING MATRIX\n");
      while (fscanf(pFile, "%s%c", str, &ch) != EOF && newline < 3) {
         if(isdigit(str[0])) {
            if(i == matrixdata.rowary[c].vcount) {
               i = 0;
               c++;
            }
            if(c > matrixdata.rcount) {
               printf("Number of rows is more than allocated rows\n");
               exit(1);
            }

            //Determine what the current values represent
            if(newline == 1) {
               //Convert String values to double strtod 
               double v;
               v = strtod(str, NULL);
               matrixdata.rowary[c].values[i] = v;
            } else if (newline == 2) { //if(strcmp(header, "col_inds:") == 0) {
               //Convert String values to int strtol
               long int v;
               v = strtol(str, NULL, 0);
               matrixdata.rowary[c].col[i] = v;
               matrixdata.rowary[c].values[i] = matrixdata.rowary[c].values[i]/outGoingEdgeCount[v];
            } 
            i++; //Counts the total number of nodes
         } else {
            c=0;
            i=0;
            newline++;
         }
         
         if (ch == '\n') { //END OF LINE
            c = 0;
            i = 0;
         }
      }
      printf("END OF FILE\n");
      fclose(pFile);
      printf("FILE CLOSED\n");
      free(outGoingEdgeCount);

      //Read partition file
      char *dot = strrchr(argv[2], '.');
      if(!dot || dot == argv[2]) dot = "";
      else dot = dot + 1;
      
      int partitionCount = 0;
      if (isdigit(dot[0])) {
         partitionCount = atoi(dot);
      } else {
         printf("INVALID PARTITION FILE\n");
         exit(1);
      }

      struct pNode *partitionArray = malloc(partitionCount * sizeof(struct pNode));
      //Initialize linked list
      for (i=0;i<partitionCount;i++){
         partitionArray[i].next = 0;
         partitionArray[i].last = 0;
         partitionArray[i].val = -1;
      }

      FILE *partFile = fopen(argv[2], "r");

      if (partFile==NULL) {
         printf("File error");
         fputs ("File error", stderr); 
         exit(1);
      }

      i = 0;
      int in;

      //Read Partition file and store where rows should be borken up into
      printf("READING PARTITION FILE\n");
      while (fscanf(partFile, "%d", &in) != EOF) {
         //in is the processor number
         struct pNode *newNode = &partitionArray[in];
         struct pNode *firstNode = &partitionArray[in];

         if (firstNode->last != 0) {
            newNode = firstNode->last; //Go to last node
         } 
         if (newNode->val != -1) {
            if (newNode != 0) {
               while (newNode->next != 0){
                  newNode = newNode->next;
               }
            }

            newNode->next = malloc(sizeof(struct pNode));
            newNode = newNode->next;
            if (newNode == 0){
               printf("NO MEMORY\n");
               exit(1);
            }
         }
         
         newNode->next = 0;
         newNode->val = i; //i is the line number
         firstNode->last = newNode;
         //sizeofEachPartition[c] = sizeofEachPartition[c] + 1;
         i++;
      }
      printf("FINISHED READING PARTITION FILE\n");

      //multiplyindex will contain information about which node is on which processor. Used for vector
      multiplyindex = (int *) malloc(i * sizeof(int));

      i = 0;
      rewind(partFile);
      while (fscanf(partFile, "%d", &in) != EOF) {
         multiplyindex[i] = in;
         i++;
      }
      indexsize = i;
      fclose(partFile);

      //Partition the Matrix and send to other processors
      //Data is combined into one array. 
      //|rowTotal|rowNum|size|val|col|val|col|...|rowNum|size|val|col|...
      printf("PARTITION MATRIX\n");
      for (i=0;i<partitionCount;i++){
         struct pNode *newNode = &partitionArray[i];
         size = 1;
         int rc = 0;
         while (newNode != NULL) {
            size = size+2; //Includes the row number and row size
            size = size+(matrixdata.rowary[newNode->val].vcount*2);
            newNode = newNode->next;
            rc++;
         }
         double *partitionedData = (double *) malloc(size * sizeof(double));
         newNode = &partitionArray[i];
         c = 0;
         partitionedData[c] = rc;
         c++;
         while (newNode != NULL) {
            partitionedData[c] = newNode->val;
            c++;
            partitionedData[c] = matrixdata.rowary[newNode->val].vcount*2;
            c++;
            int m;
            for (m=0; m < matrixdata.rowary[newNode->val].vcount; m++) {
               partitionedData[c] = matrixdata.rowary[newNode->val].values[m];
               c++;
               partitionedData[c] = matrixdata.rowary[newNode->val].col[m];
               c++;
            }
            newNode = newNode->next;
         }

         //sending partitioned data
         if (i == 0) {
            recvMatrixData = (double *) malloc(size * sizeof(double));
            recvSize = size;
            recvMatrixData = partitionedData;
         } else {
            MPI_Send(partitionedData, size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);   
         }
      }

      //printf("BCAST multiplyindex\n");
      numofRows = matrixdata.rcount;
      MPI_Bcast(&numofRows, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(multiplyindex, indexsize, MPI_INT, 0, MPI_COMM_WORLD);
 
      for (c=0; c < matrixdata.rcount; c++) { 
         free(matrixdata.rowary[c].col);
         free(matrixdata.rowary[c].values);
      }
      free(partitionArray);
      free(matrixdata.rowary);
   } 
   
   if(myid != 0) {
      //Recv partition data
      MPI_Probe(0, 0, MPI_COMM_WORLD, &stat);
      MPI_Get_count(&stat, MPI_DOUBLE, &recvSize);
      recvMatrixData = (double *) malloc(recvSize * sizeof(double));
      MPI_Recv(recvMatrixData, recvSize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &stat);
      
      //Recv multiplyIndex
      MPI_Bcast(&numofRows, 1, MPI_INT, 0, MPI_COMM_WORLD);
      multiplyindex = (int *) malloc(numofRows * sizeof(int));
      MPI_Bcast(multiplyindex, numofRows, MPI_INT, 0, MPI_COMM_WORLD);
   } 

   row pmatrixdata;
   int rowCount = 0;
   int mc = 0;
   i = 0;

   //Columns which are needed for multiplication. 
   struct pNode *colLink = malloc(sizeof(struct pNode)); 
   struct pNode *newColNode;
   int colCount = 0; //Number of columns needed to multiply

   int *allColumnList = (int*) malloc(numofRows * sizeof(int));
   for(c=0; c<numofRows; c++) {
      allColumnList[c] = 0;
   }

   colLink->next = 0;
   colLink->last = 0;
   colLink->val  = -1;

   newColNode = colLink;

   //Have list of columns used. As I read in the data, store unique columns 
   //printf("ID(%d) COMPUTE COLs NEEDED\n", myid);
   for(c=0; c<recvSize; c++){
      if (c == 0) {
         pmatrixdata.rowary = malloc(recvMatrixData[0] * sizeof(rowData));
         pmatrixdata.rcount = recvMatrixData[0];   
      }else if (c - rowCount == 1) {
         if (c > 1) 
            mc++;
         
         i = 0;
         pmatrixdata.rowary[mc].rowId = (int)recvMatrixData[c];
      } else if (c - rowCount==2){
         rowCount = recvMatrixData[c];
         pmatrixdata.rowary[mc].values = malloc(rowCount/2 * sizeof(double));
         pmatrixdata.rowary[mc].col = malloc(rowCount/2 * sizeof(long int));
         pmatrixdata.rowary[mc].vcount = rowCount/2;
         rowCount = rowCount + c;
      } else if (c % 2) {
         pmatrixdata.rowary[mc].values[i] = recvMatrixData[c];
         pmatrixdata.rowary[mc].col[i] = (int)recvMatrixData[c+1];
         if (allColumnList[(int)recvMatrixData[c+1]] == 0) {
            allColumnList[(int)recvMatrixData[c+1]] = 1;
            if (newColNode->val != -1) {
               if (newColNode != 0) {
                  while (newColNode->next != 0){
                     newColNode = newColNode->next;
                  }
               }
               newColNode->next = malloc(sizeof(struct pNode));
               newColNode = newColNode->next;
               if (newColNode == 0){
                  printf("NO MEMORY\n");
                  exit(1);
               }
            }
            newColNode->next = 0;
            newColNode->val = (int)recvMatrixData[c+1];
            colCount++;
         }
         c++;
         i++;
      }
   }
   free(allColumnList);

   //Create multiply vector. 
   double *resultArray = (double *) malloc(numofRows * sizeof(double)); //malloc(pmatrixdata.rcount * sizeof(double));
   double *multiplyArray = malloc(numofRows * sizeof(double));
   //printf("ID(%d) INITIALIZE MULTIPLY ARRAY\n", myid);
   for(c=0; c<numofRows; c++) {
      multiplyArray[c] = 1/(double)numofRows;
      resultArray[c] = 0;
   }

   bool pageRankComputed = false;
   int calcCount = 0;
   double normalizedValue = 0;
   
   struct pNode *recColArray = malloc(numprocs * sizeof(struct pNode)); //array of linked lists, each array element is a processor
   int *recColCount = malloc(numprocs * sizeof(int));
   for (i=0; i<numprocs; i++) {
      recColArray[i].next = 0;
      recColArray[i].last = 0;
      recColArray[i].val  = -1;
      recColCount[i] = 0;
   }

   newColNode = colLink; //iterate through columns which are needed for multiplication. 

   //loop through required column linked list, newColNode
   //printf("POPULATE LIST OF PROC TO REQ COLS\n");
   while (newColNode != NULL) {
      int pmultiply = multiplyindex[newColNode->val];
      if (pmultiply != myid) {
         struct pNode *newRecColNode = &recColArray[pmultiply];
         struct pNode *firstNode = &recColArray[pmultiply];
         if (firstNode->last != 0) {
            newRecColNode = firstNode->last;
         }

         if (newRecColNode->val != -1){
            if (newRecColNode != 0) {
               while (newRecColNode->next != 0) {
                  newRecColNode = newRecColNode->next;
               }
            }
            newRecColNode->next = malloc(sizeof(struct pNode));
            newRecColNode = newRecColNode->next;
            if (newRecColNode == 0) {
               printf("NO MEMORY\n");
               exit(1);
            }
         }
         newRecColNode->next = 0;
         newRecColNode->val = newColNode->val;
         firstNode->last = newRecColNode;

         //Increment number if columns needed for each processor
         recColCount[pmultiply] = recColCount[pmultiply] + 1; 
      }
      newColNode = newColNode->next;
   }

   MPI_Barrier(MPI_COMM_WORLD);
   //printf("ID(%d) PERFORM COMPUTATION\n", myid);
   if (myid == 0)
      gettimeofday(&start, NULL);
   
   while(calcCount == 0 || finish != 0) {
      double sumValue = 0;
      MPI_Barrier(MPI_COMM_WORLD);
      if (calcCount != 0) {
         for(c=1; c<numprocs; c++) {
            int requestFromProcessorNumber = (myid-c+numprocs)%numprocs;
            int requestToProcessorNumber = (myid+c)%numprocs;
            int recColSize;
            struct pNode *subNode = &recColArray[requestToProcessorNumber];
            
            if (subNode != NULL) {
               int *neededColValues = (int *)malloc(recColCount[requestToProcessorNumber] * sizeof(int));
               if (subNode->val != -1) {
                  i=0;
                  while (subNode != NULL) {
                     neededColValues[i] = subNode->val;
                     subNode = subNode->next;
                     i++;
                  }
                  MPI_Barrier(MPI_COMM_WORLD);
                  
                  MPI_Request request[2];

                  MPI_Isend(neededColValues, recColCount[requestToProcessorNumber], MPI_INT, requestToProcessorNumber, 4, MPI_COMM_WORLD, &request[0]);
                  //MPI_Send(neededColValues, recColCount[requestToProcessorNumber], MPI_INT, requestToProcessorNumber, 4, MPI_COMM_WORLD);

                  MPI_Probe(requestFromProcessorNumber, 4, MPI_COMM_WORLD, &stat);
                  MPI_Get_count(&stat, MPI_INT, &recColSize);
                  int *requestFromProcessorColArray = (int *) malloc(recColSize * sizeof(int));
                  MPI_Irecv(requestFromProcessorColArray, recColSize, MPI_INT, requestFromProcessorNumber, 4, MPI_COMM_WORLD, &request[1]);
                  //MPI_Recv(requestFromProcessorColArray, recColSize, MPI_INT, requestFromProcessorNumber, 4, MPI_COMM_WORLD, &stat);
                  
                  MPI_Waitall(2, request, NULL);

                  double *sendRecColArray = (double *) malloc(recColSize * sizeof(double));
                  for (i=0;i<recColSize;i++) {
                     sendRecColArray[i] = multiplyArray[requestFromProcessorColArray[i]];
                  } 
                  MPI_Isend(sendRecColArray, recColSize, MPI_DOUBLE, requestFromProcessorNumber, 5, MPI_COMM_WORLD, &request[0]);
                  //MPI_Send(sendRecColArray, recColSize, MPI_DOUBLE, requestFromProcessorNumber, 5, MPI_COMM_WORLD);

                  double *recvRecColArray = (double *) malloc(recColCount[requestToProcessorNumber] * sizeof(double));
                  MPI_Irecv(recvRecColArray, recColCount[requestToProcessorNumber], MPI_DOUBLE, requestToProcessorNumber, 5, MPI_COMM_WORLD, &request[1]);
                  //MPI_Recv(recvRecColArray, recColCount[requestToProcessorNumber], MPI_DOUBLE, requestToProcessorNumber, 5, MPI_COMM_WORLD, &stat);

                  MPI_Waitall(2, request, NULL);

                  subNode = &recColArray[requestToProcessorNumber];
                  for (i=0; i<recColCount[requestToProcessorNumber]; i++) {
                     multiplyArray[neededColValues[i]] = recvRecColArray[i];
                  }

                  free(requestFromProcessorColArray);
                  free(sendRecColArray);
                  free(recvRecColArray);
               }
               free(neededColValues);
            } 

            MPI_Barrier(MPI_COMM_WORLD);
         }
      }
      
      //Loop through each row
      for (c=0; c < pmatrixdata.rcount; c++) {
         double rowSum = 0;
         //Loop through each column
         for (i=0; i<pmatrixdata.rowary[c].vcount; i++) {
            rowSum = rowSum + (pmatrixdata.rowary[c].values[i] * multiplyArray[pmatrixdata.rowary[c].col[i]]);
         }
         resultArray[pmatrixdata.rowary[c].rowId] = rowSum;
         sumValue = sumValue + pow(rowSum, 2);
      }
      
      if(myid == 0) {
         for(i=1;i<numprocs;i++) {
            double receivedValue = 0;
            MPI_Recv(&receivedValue, 1, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &stat); 
            sumValue = sumValue + receivedValue;
         }
         normalizedValue = sqrt(sumValue);
         
         iterationCount++;
         if (epsilon != 0 && ((epsilon - normalizedValue < 10E-5) && (epsilon - normalizedValue > -10E-5))) {
            finish = 0; //End Calculations
         } 
         epsilon = normalizedValue;

         MPI_Bcast(&normalizedValue, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
         MPI_Bcast(&finish, 1, MPI_CHAR, 0, MPI_COMM_WORLD);
      } else {
         MPI_Send(&sumValue, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
         MPI_Bcast(&normalizedValue, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
         MPI_Bcast(&finish, 1, MPI_CHAR, 0, MPI_COMM_WORLD);
         //MPI_Recv(&normalizedValue, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, &stat);
      }
      
      for (c=0; c < pmatrixdata.rcount; c++) {
         multiplyArray[pmatrixdata.rowary[c].rowId] = resultArray[pmatrixdata.rowary[c].rowId];
         resultArray[pmatrixdata.rowary[c].rowId] = 0;
      }
      calcCount++;
   } 

   //END OF COMPUTATION
   MPI_Barrier(MPI_COMM_WORLD);
   if (myid == 0)
      gettimeofday(&end, NULL);

   //Combine multiply array and write to file
   if (myid != 0) {
      double *finalMultiplyArray = malloc((pmatrixdata.rcount*2) * sizeof(double));
      i=0;
      for (c=0; c < pmatrixdata.rcount; c++) {
         finalMultiplyArray[i] = pmatrixdata.rowary[c].rowId;
         i++;
         finalMultiplyArray[i] = multiplyArray[pmatrixdata.rowary[c].rowId];
         i++;
      }
      MPI_Send(finalMultiplyArray, pmatrixdata.rcount*2, MPI_DOUBLE, 0, 7, MPI_COMM_WORLD);
      free(finalMultiplyArray);
   } else {
      printf("id(0), normalizedValue = %.12f\n", normalizedValue);  
      printf("iterations: %d\n", iterationCount);
      for (c=0; c<numprocs; c++) {
         if (c!=myid) {
            MPI_Probe(MPI_ANY_SOURCE, 7, MPI_COMM_WORLD, &stat);
            int recColSize;

            MPI_Get_count(&stat, MPI_DOUBLE, &recColSize);

            double *recvFinalMultiplyArray = (double *) malloc(recColSize * sizeof(double));
            MPI_Recv(recvFinalMultiplyArray, recColSize, MPI_DOUBLE, stat.MPI_SOURCE, 7, MPI_COMM_WORLD, NULL);

            for (i=0; i<recColSize; i++) {
               if ((i % 2) == 0) {
                  multiplyArray[(int)recvFinalMultiplyArray[i]] = recvFinalMultiplyArray[i+1];
               }
            }
            free(recvFinalMultiplyArray);
         }
      }
      long runTime = (end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec);
      printf("TIME: %ld\n", runTime);
      FILE *fp;

      fp = fopen("pagerank.result", "w+");
      fprintf(fp, "time: %ld\n", runTime);
      fprintf(fp, "node_Id pagerank \n");
      for (c=0; c < numofRows; c++) { 
         fprintf(fp, "%d %.12f\n", c, multiplyArray[c]);
      }
      fclose(fp);
   }

   free(recColArray);
   free(recColCount);
   free(resultArray);
   free(colLink);
   free(multiplyArray);
   for (c=0; c < pmatrixdata.rcount; c++) { 
      free(pmatrixdata.rowary[c].col);
      free(pmatrixdata.rowary[c].values);
   }
   free(pmatrixdata.rowary);

   MPI_Finalize(); 
}