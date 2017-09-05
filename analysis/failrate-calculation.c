#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

int main() {

FILE *f;


f = fopen("trim10.csv","r");

int lines;
float c;

lines = 0;
while (!feof(f)) {
	c = fgetc(f);
	if (c == '\n') {
		lines++;
	}
}

printf("Reading number of lines... \n");
printf("%i\n",lines);


//float array[lines];

float *array;

array = (float *)malloc(sizeof(float)*lines);

rewind(f);
// char buffer[100];
 
// while (fgets(buffer, sizeof buffer, f) != NULL)
//   fputs(buffer, stdout);

// setvbuf(stdout, NULL, _IOLBF, 0);
printf("Rewind to sof... \n");
int j;
float num;
j=0;

while (!feof(f)) {
	fscanf(f,"%f",&array[j]);
	// printf("%f\n",array[j]);
	j++;
}

// for (j = 0; j < lines; j++) {
// 	fscanf(f,"%1f \n", &array[j]);
// 	printf("%1f \n", array[j]);
// }
fclose(f);


printf("Successfully read file into array.\n");

// printf("%i\n",lines);



int g;

int len = lines;

while (true) {
	


	int i = 0;
	lines -= g;
	int m = 1;

	while (i < lines  && m != (double)floor( (double) array[i])) {
		i++;
	} 

	printf("Finding imax... \n");
	printf("First imax: %d\n",i);
	float k = array[i];
	printf("Value at imax: %f\n",k);

	//array slicing format that maps to y = x[n:m]
	//memcpy(y, x + n, (m - n) * sizeof(*y))

	bool zerosfinder;
	g = i ;
	
	for (g; g<lines; g++) {

		zerosfinder = (floor(array[g]) <= 1);
		if ((zerosfinder)==false) {
			printf("Index where the zeroes stop: %d\n",g);
			break;
		}
	}

	g +=1;

	float *tmp;

	tmp = (float *)malloc(sizeof(float)*(lines-g));

	memcpy(tmp,array+g,(lines-g) * sizeof(*tmp));

	memcpy(array,tmp,sizeof(*tmp));



	int v = 0;
	int n = 5;
	while (v < lines  && n != (double)ceil( (double) array[v])) {
		v++;
	} 

	float *slc;

	slc = (float *)malloc(sizeof(float)*(g-v));

	memcpy(slc,array+v,(g-v) * sizeof(*array));

	printf("Slice array test %f\n", slc[0]);
	

	int u = 0;
	while (u < lines  && 10 != (double)floor( (double) slc[u])) {
		u++;
	} 
	printf("back above n at index: %d \n",u);

	break;
}


// if (i < lines) {
//    printf("Number found at the location = %d", i + 1);
// } else {
//    printf("Number not found");
// }




return 0;
}

