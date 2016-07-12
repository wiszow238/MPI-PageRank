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
} rowData;

typedef struct {
   rowData* rowary;  //Stores data for each row
   int rcount;       //Number of rows
} row;

int main(int argc, char **argv) { 
   struct timeval start, end;
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

   int i = 0;
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
   
   double *multiplyArray = (double *) malloc(matrixdata.rcount * sizeof(double));
   int *outGoingEdgeCount = (int *) malloc(matrixdata.rcount * sizeof(int));

   int c;
   //Initialize the matrix and the multiply array
   for (c=0; c < matrixdata.rcount; c++) { //i is the total amount of rows in the matrix
      long int rowSize;
      rowSize = rownum[c+1] - rownum[c];

      multiplyArray[c] = 1/(double)matrixdata.rcount;

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

   printf("CALCULATING NORMS\n");
   double *normalValue = (double *) malloc(2 * sizeof(double));
   double *resultArray = (double *) malloc(matrixdata.rcount * sizeof(double));
   normalValue[1] = 0;
   double epsilon = 0;
   gettimeofday(&start, NULL);

   int a = 0;
   do {
      double norm = 0;
      double sumValue = 0;
      //Loop through each row
      for (c=0; c < matrixdata.rcount; c++) {
         double rowSum = 0;
         //Loop through each column
         for (i=0; i < matrixdata.rowary[c].vcount; i++) {
            rowSum = rowSum + (matrixdata.rowary[c].values[i] * multiplyArray[matrixdata.rowary[c].col[i]]);
         }
         resultArray[c] = rowSum;
         sumValue = sumValue + pow(rowSum, 2);
      }

      norm = sqrt(sumValue);
      for (c=0; c < matrixdata.rcount; c++) {
         multiplyArray[c] = resultArray[c];
         resultArray[c] = 0;
      }
      normalValue[0] = normalValue[1];
      normalValue[1] = norm;
      a++;
      //printf("norm: %.12f\n", norm);
   } while(normalValue[0] == 0 || ((normalValue[0] - normalValue[1]) > 10E-5) || ((normalValue[1] - normalValue[0]) > 10E-5));

   printf("iterations: %d\n", a);
   printf("epsilon - normalizedValue = %f\n", normalValue[0] - normalValue[1]);
   printf("normalizedValue = %f\n", normalValue[1]);
   gettimeofday(&end, NULL);
   long runTime = (end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec);
   printf("TIME: %ld\n", runTime);
   FILE *fp;

   fp = fopen("pagerank.result", "w+");
   fprintf(fp, "time: %ld\n", runTime);
   fprintf(fp, "node_Id pagerank \n");
   for (c=0; c < matrixdata.rcount; c++) { 
      fprintf(fp, "%d %.12f\n", c, multiplyArray[c]);
      free(matrixdata.rowary[c].col);
      free(matrixdata.rowary[c].values);
   }
   
   fclose(fp);
   free(normalValue);
   free(resultArray);
   free(matrixdata.rowary);
   free(multiplyArray);
}