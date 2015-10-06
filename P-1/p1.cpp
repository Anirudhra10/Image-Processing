/*
	B490/B659 Project 1 Skeleton Code    (1/2015)
	
	Be sure to read over the project document and this code (should you choose to use it) before 
	starting to program. 
	

	Compiling:
		A simple console command 'make' will excute the compilation instructions
		found in the Makefile bundled with this code. It will result in an executable
		named p1.

	Running:
		The executable p1 takes commands of the form:
			./p1 problem_ID input_File ouput_File additional_Arguments
	
		Some examples:

			./p1 2.1 input.png out.png 
			
			This runs the 'averageGrayscale' function defined in this file and described
			in the project doc Part 2, problem 1 on the input. The output is saved as out.png. 
		
			----

			./p1 4.2b input.jpg output.png 0.5

			This runs the Gaussian filter function (gaussianFilter) with sigma = 0.5 on input.jpg.
			This is problem 2b from Part 4 in the project documentation. 

*/


//Link to the header file
#include "CImg.h"
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <cstdlib>
#include <cmath>
#define PI 3.14159265359

//Use the cimg namespace to access the functions easily
using namespace cimg_library;
using namespace std;


//Part 2 - Basic image operations
CImg<double> averageGrayscale(CImg<double> input);
CImg<double> simpleBW(CImg<double> input);
CImg<double> advancedBW(CImg<double> input);

//Part 3 - Adding noise
CImg<double> uniformNoise(CImg<double> input);
CImg<double> gaussianNoise(CImg<double> input, double sigma);
CImg<double> saltAndPepperNoise(CImg<double> input);

//Part 4 - Filtering
CImg<double> filter(CImg<double> input, CImg<double> filter,int filterSize);
CImg<double> meanFilter(CImg<double> input, int filterSize);
CImg<double> gaussianFilter(CImg<double> input, double sigma);
CImg<double> medianFilter(CImg<double> input, int size);



int main(int argc, char **argv){


	if(argc < 4){
		cout << "Insufficent number of arguments. Please see documentation" << endl;
		cout << "p1 problemID inputfile outputfile" << endl;
		return -1;
	}
	
	
	char* inputFile = argv[2];
	char* outputFile = argv[3];
	cout << "In: " << inputFile <<"  Out: " << outputFile << endl;
	
	CImg<double> input(inputFile);
	CImg<double> output; 

	if(!strcmp(argv[1], "2.1")){	
		cout << "# Problem 2.1 - Average Grayscale" << endl;
		if(input.spectrum() != 3){	cout << "INPUT ERROR: Input image is not a color image!" << endl;return -1;}
		output = averageGrayscale(input);
		output.save(outputFile);
	}
	else if(!strcmp(argv[1], "2.2a")){
		cout << "# Problem 2.1a - Simple Threshold Black and White" << endl;
		if(input.spectrum() != 1){	cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;return -1;}
		output = simpleBW(input);
		output.save(outputFile);
	}
	else if(!strcmp(argv[1], "2.2b")){
		cout << "# Problem 2.2b -typedef Advanced Threshold Black and White" << endl;
		if(input.spectrum() != 1){	cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;return -1;}
		output = advancedBW(input);
		output.save(outputFile);
	}
	else if(!strcmp(argv[1], "3.1")){
		cout << "# Problem 3.1 - Uniform Noise" << endl;
		if(input.spectrum() != 1){	cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;return -1;}
		output = uniformNoise(input);
		output.save(outputFile);
	}
	else if(!strcmp(argv[1], "3.2")){
		cout << "# Problem 3.2 - Gaussian Noise" << endl;
		if(input.spectrum() != 1){	cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;return -1;}
		if(argc != 5){	cout << "INPUT ERROR: Provide sigma as additional argument!" << endl;return -1;}
		output = gaussianNoise(input, atof(argv[4]));
		output.save(outputFile);
	}
	else if(!strcmp(argv[1], "3.3")){
		cout << "# Problem 3.3 - Salt & Pepper Noise" << endl;
		if(input.spectrum() != 1){	cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;return -1;}
		output = saltAndPepperNoise(input);
		output.save(outputFile);
	}
	else if(!strcmp(argv[1], "4.2a")){
		cout << "# Problem 4.2a - Mean Filter Noise" << endl;
		if(input.spectrum() != 1){	cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;return -1;}
		if(argc != 5){	cout << "INPUT ERROR: Provide filter size as additional argument!" << endl;return -1;}
		output = meanFilter(input, atoi(argv[4]));
		output.save(outputFile);
	}
	else if(!strcmp(argv[1], "4.2b")){
		cout << "# Problem 4.2b - Gaussian Filter Noise" << endl;
		if(input.spectrum() != 1){	cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;return -1;}
		if(argc != 5){	cout << "INPUT ERROR: Provide sigma as additional argument!" << endl;return -1;}
		output = gaussianFilter(input, atof(argv[4]));
		output.save(outputFile);
	}
	else if(!strcmp(argv[1], "4.3")){
		cout << "# Problem 4.3 - Median Noise" << endl;
		if(input.spectrum() != 1){	cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;return -1;}
		if(argc != 5){	cout << "INPUT ERROR: Provide filter size as additional argument!" << endl;return -1;}
		output = medianFilter(input, atoi(argv[4]));
		output.save(outputFile);
	}
	else if(!strcmp(argv[1], "5")){
		cout << "# Problem 5 - Noise Removal Analysis" << endl;
		//FIXME You will need to implement this section yourself
	}
	else if(!strcmp(argv[1], "6.1")){
		cout << "# Problem 6.1 - Separable Kernel Convolutions" << endl;
		
		std::clock_t start = std::clock();
     	//FIXME You will need to implement this section yourself
    	std::cout << "Time Elapsed: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
	
	}
	else if(!strcmp(argv[1], "6.2")){
		cout << "# Problem 6.2 - Dynamic Box Filter" << endl;
		
		std::clock_t start = std::clock();
     	//FIXME You will need to implement this section yourself
    	std::cout << "Time Elapsed: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;

	}
	else if(!strcmp(argv[1], "6.3")){
		cout << "# Problem 6.3 - Fast Gaussian Smoothing" << endl;
		
		std::clock_t start = std::clock();
     	//FIXME You will need to implement this section yourself
    	std::cout << "Time Elapsed: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;

	}
	else{
		cout << "Unknown input command" << endl; 
	}
	return 0;
}


//Part 2 - Basic image operations
CImg<double> averageGrayscale(CImg<double> input){
	//Creates a new grayscale image with same size as the input initialized to all 0s (black)
	CImg<double> output(input.width(), input.height(), 1, 1); 

	//FIXME fill in new grayscale image
	for(int x=0;x<output.width();x++){
		for(int y=0;y<output.height();y++){
				output(x,y,0,0)=0.3*input(x,y,0,0)+0.6*input(x,y,0,1)+0.1*input(x,y,0,2);
}
}

	
	return output;
}

CImg<double> simpleBW(CImg<double> input){
	//Creates a new grayscale image with same size as the input initialized to all 0s (black)
	CImg<double> output(input, "xyzc", 0);  
/*	for(int x=0;x<output.width();x++){
                for(int y=0;y<input.height();y++){
				if(input(x,y,0,0)>127||input(x,y,0,0)==127){
					output(x,y,0,0)=255;}
					else
					output(x,y,0,0)=0;
}
}*/
	//FIXME Do stuff
       for(int i=0;i<output.width();i++)
	{
	for(int j=0;j<output.height();j++)
	{
		if(input(i,j,0)>127)
			output(i,j,0)=255;
		else
			output(i,j,0)=0;
	}
}	
	return output;
}

CImg<double> advancedBW(CImg<double> input){
	double e;
	//Creates a new grayscale image with same size as the input initialized to all 0s (black)
	        CImg<double> output(input.width(), input.height(), 1, 1);
	
	//                //FIXME Do stuff
	for(int x=0;x<output.width();x++){
                for(int y=0;y<input.height();y++){
			
			if(input(x,y,0)>127)
                               output(x,y,0)=255;
                        else
                               output(x,y,0)=0;
			
		e=input(x,y,0)-output(x,y,0);   //compute the error value

			if(x+1<input.width())                             //pixel to the right
				input(x+1,y,0)=input(x+1,y,0)+0.4*e;

			if(y+1<input.height())                            //pixel at the bottom
				input(x,y+1,0)=input(x,y+1,0)+0.3*e;

			if(x+1<input.width() && y+1<input.height())       //pixel to the right and below
				input(x+1,y+1,0)=input(x+1,y+1,0)+0.1*e;

			if(x>0 && y+1<input.height())	                  //pixel to the left and below
				input(x-1,y+1,0)=input(x-1,y+1,0)+0.2*e;
}
}	
	return output;
}

//Part 3 - Adding noise
CImg<double> uniformNoise(CImg<double> input){
	double noise;
	//Creates a new image with same size as the input initialized to all 0s (black)
	CImg<double> output(input, "xyzc", 0); 
//		for(int x=0;x<output.width();x++){
//                    for(int y=0;y<input.height();y++){
int random_num,sum;				
	//FIXME Add noise
for (int i=0;i<input.width();i++)
{
	for(int j=0;j<input.height();j++)
	{
		random_num= (rand()%40)-20;
		sum = input(i,j,0)+random_num;
		if(sum>255)
			sum =255;
		if(sum<0)
			sum=0;
	
		output(i,j,0) = sum;	
	}
}

	return output;
}


CImg<double> gaussianNoise(CImg<double> input, double sigma){

	//Creates a new image with same size as the input initialized to all 0s (black)
	CImg<double> output(input, "xyzc", 0); 

	//FIXME Add noise
       
	return output;
}

CImg<double> saltAndPepperNoise(CImg<double> input){

	//Creates a new image with same size as the input initialized to all 0s (black)
	CImg<double> output(input, "xyzc", 0); 

	//FIXME Add noise
int random_pixel,black_white;
for (int i=0;i<input.width();i++)
{
        for(int j=0;j<input.height();j++)
        {
		output(i,j,0) = input(i,j,0);

                random_pixel = (rand()%10)+1;       //generate a random number between 1-10 to achieve a probability of (0.1)-> Low probability
		
		if(random_pixel == 5)  //5 has a probability of 0.1 (any number between 1-10 has a probability of 0.1)-> This is to choose a pixel with P=0.1 
		{
			black_white = (rand()%2)+1;  //to set the value of pixel to 0 or 255
			
			if(black_white == 1)        //has P=0.5
				output(i,j,0) = 0;
			else
				output(i,j,0) = 255; 
		}
	}
}

	return output;
}

CImg<double> reflect_edge(CImg<double> input,int filterSize){

double sum = 0;
int fsize;
int width,height;

fsize = (filterSize-1)/2;
width = input.width();
height = input.height();

CImg<double> inter(input.width()+2*fsize,input.height()+2*fsize);

for(int i=0;i<input.width();i++)        //copy the input image to an intermediate img of a bigger dimension
        for(int j=0;j<input.height();j++)
                inter(i+fsize,j+fsize,0)=input(i,j,0);

for(int i=0; i<fsize;i++)
{
        for(int j=0; j<height;j++)
        {       if(j<fsize)   //the corner squares
                {
                        inter(fsize-i-1,fsize-j-1,0) = input(i,j,0); //top left square
                        inter(i+fsize+width,fsize-j-1,0)  = input(width-i-1,j,0);//top right square
                }

                inter(i,j+fsize,0)=input(fsize-i-1,j,0); //left column
                inter(i+fsize+width,j+fsize,0)=input(width-i-1,j,0); //right column
        }
	
	for(int j=0; j<width;j++)
        {
                if(j<fsize)
                {
                        inter(fsize-i-1,j+fsize+height,0) = input(i,height-j-1,0);//bottom left square
                        inter(i+fsize+width,j+fsize+height,0) = input(width-i-1,height-j-1,0);//bottom right square
                }
                inter(j+fsize,i,0)=input(j,fsize-i-1,0); //top row
                inter(j+fsize,i+fsize+height,0)=input(j,height-i-1,0);//bottom row
        }
}
	
return inter;
}
//Part 4 - Filtering
CImg<double> filter(CImg<double> input, CImg<double> kernel,int filterSize){

	//Creates a new image with same size as the input initialized to all 0s (black)
	CImg<double> output(input, "xyzc", 0); 
int fsize,sum=0;
fsize = (filterSize - 1)/2;

CImg<double> inter = reflect_edge(input,filterSize);

	for(int i=fsize; i+fsize < inter.width();i++)
	{
                for(int j=fsize; j+fsize < inter.height(); j++)
		{
			sum=0;
                      for(int u = 0; u < filterSize; u++)
			for(int v = 0; v < filterSize; v++)
				sum += kernel(u,v)*input((i+u-fsize),(j+v-fsize),0);

			output(i-fsize,j-fsize,0) = sum;
		}

	}

	//FIXME Convolve with filter
	return output;

}

CImg<double> meanFilter(CImg<double> input, int filterSize){

	//Creates a new grayscale image (just a matrix) to be our filter
	CImg<double> H(filterSize, filterSize, 1, 1); 

	//FIXME Fill filter values
double value=0,norm=0;
value = (filterSize * filterSize);
norm = 1/value;

	for(int i=0;i<filterSize;i++)
		for(int j=0; j< filterSize; j++)
			H(i,j) = norm;             // in a box mean filter, all the values are equal to 1/(dim*dim)

	//Convole with filter and return
	return filter(input, H, filterSize);

}

CImg<double> gaussianFilter(CImg<double> input, double sigma){

	int filterSize;
	filterSize	= 3*sigma - 1;
	double rad = 0,sum = 0;
	double iter = (filterSize - 1)/2;
	double sig = 2*sigma*sigma;

 CImg<double> H(filterSize, filterSize, 1, 1);

for(int i = -iter; i <= iter; i++)
{
	for(int j = -iter; j <= iter; j++) 
	{
		rad = (i*i + j*j);
		H(i+iter,j+iter) = (exp(-rad/sig))/(PI*sig);
		sum += H(i+iter,j+iter);		
	}
} 
	
for(int x = 0; x < filterSize; x++)
	for(int y = 0; y < filterSize; y++)
		H(x,y) /= sum;
	
	return filter(input, H, filterSize);
}

CImg<double> medianFilter(CImg<double> input, int size){

	//Creates a new image with same size as the input initialized to all 0s (black)
	CImg<double> output(input, "xyzc", 0); 
CImg<double> inter = reflect_edge(input,size);

int fsize=0,array_size=0;
int *array,temp,pos;

fsize = (size-1)/2;
array_size = size*size;

	array = new int[array_size];
	
for(int i=fsize; i+fsize < inter.width();i++)
{ 
	for(int j=fsize; j+fsize < inter.height();j++)
	{

	temp = 0;

		for(int u=0; u<size; u++)
		{
			for(int v=0; v<size; v++)
			{
				array[u*size+v] = input((i+u-fsize),(j+v-fsize),0); 
			}
		}    // copy the contents of 3*3 segment of the image to a 1-d array

		for(int a=0; a<array_size-1; a++)                      //selection sort - sort the array 
		{
		pos = a;			
			for(int b=a+1; b<array_size; b++)
			{
				if(array[pos] > array[b])
					pos = b;
			}
			if(pos!= a)
			{
				temp = array[a];
				array[a] = array[pos];
				array[pos] = temp;
			}
		}
	
		output(i-fsize,j-fsize,0) = array[(array_size - 1)/2];//median value is the 5th element of the sorted array-replace the center value of the output 3*3 segment with the median value
	}	
}
	return output;
}






