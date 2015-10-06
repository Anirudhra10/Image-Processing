/*
	B490/B659 Project 2 Skeleton Code    (2/2015)
	
	Be sure to read over the project document and this code (should you choose to use it) before 
	starting to program. 
	

	Compiling:
		A simple console command 'make' will excute the compilation instructions
		found in the Makefile bundled with this code. It will result in an executable
		named p2.

	Running:
		The executable p2 should take commands of the form:
			./p2 problem_ID input_File ouput_File additional_Arguments
	
*/


//Link to the header file
#include "CImg.h"
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <fft.h>
#include <cmath>
#define PI 3.14159
//Use the cimg namespace to access the functions easily
using namespace cimg_library;
using namespace std;


// This code requires that input be a *square* image, and that each dimension
//  is a power of 2; i.e. that input.width() == input.height() == 2^k, where k
//  is an integer. You'll need to pad out your image (with 0's) first if it's
//  not a square image to begin with. (Padding with 0's has no effect on the FT!)
//
// Forward FFT transform: take input image, and return real and imaginary parts.
//
void fft(const CImg<double> &input, CImg<double> &fft_real, CImg<double> &fft_imag)
{
  fft_real = input;
  fft_imag = input;
  fft_imag = 0.0;

  FFT_2D(1, fft_real, fft_imag);
}

// Inverse FFT transform: take real and imaginary parts of fourier transform, and return
//  real-valued image.
//
void ifft(const CImg<double> &input_real, const CImg<double> &input_imag, CImg<double> &output_real)
{
  output_real = input_real;
  CImg<double> output_imag = input_imag;

  FFT_2D(0, output_real, output_imag);
}

// Write this in Part 3.1
CImg<double> fft_magnitude(const CImg<double> &fft_real, const CImg<double> &fft_imag);

// Write this in Part 3.2
CImg<double> remove_interference(const CImg<double> &input);

// Write this in Part 3.3
CImg<double> fft_filter(const CImg<double> &input, const CImg<double> &filter);

// Write this in Part 3.4
CImg<double> fft_defilter(const CImg<double> &input, const CImg<double> &filter);

// Write this in Part 4 -- add watermark N to image
CImg<double> mark_image(const CImg<double> &input, int N);

// Write this in Part 4 -- check if watermark N is in image
CImg<double> check_image(const CImg<double> &input, int N);

//Gaussian filter function
CImg<double> gaussianFilter(double sigma){

	int filterSize;
	filterSize = 3*sigma - 1;
	double rad = 0,sum = 0;
	double iter = (filterSize - 1)/2;
	double sig = 2*sigma*sigma;

 	CImg<double> H(filterSize, filterSize, 1, 1);
	
	for(int i = -iter; i <= iter; i++){
		for(int j = -iter; j <= iter; j++){
			rad = (i*i + j*j);
			H(i+iter,j+iter) = (exp(-rad/sig))/(PI*sig);
			sum += H(i+iter,j+iter);		
		}
	} 
	
	for(int x = 0; x < filterSize; x++)
		for(int y = 0; y < filterSize; y++)
			H(x,y) /= sum;
	
	return H;
}
//End of function


//Finds the maximum dimension and returns the value.
double findMaxDimension(CImg<double>input)
{
		double max_val = max(input.width(),input.height());
		double max_log = log(max_val)/log(2);

		if(max_log != (int)max_log)
			return pow(2,floor(max_log)+1);
		else
			return pow(2,max_log);
}
//End of function


//Main function
int main(int argc, char **argv)
{
  try {

    if(argc < 4)
      {
	cout << "Insufficent number of arguments; correct usage:" << endl;
	cout << "    p2 problemID inputfile outputfile" << endl;
	return -1;
      }
    
    string part = argv[1];
    char* inputFile = argv[2];
    char* outputFile = argv[3];
    cout << "In: " << inputFile <<"  Out: " << outputFile << endl;
    
    CImg<double> input(inputFile);
    CImg<double> output;

    double init_w = input.width();
    double init_h = input.height();

    
	if(part == "2.1")
	{
 			
	}
	
	
    	if(part == "3.1")
     	{
		double max_dim = findMaxDimension(input);
                input.resize(max_dim,max_dim,1,1,0);		

		CImg<double> fft_real(input,"abc",0);
		CImg<double> fft_imag(input,"def",0);
	 
		 fft(input,fft_real,fft_imag);

		output = fft_magnitude(fft_real,fft_imag);
		output.save(outputFile);	
     	}


	if(part == "3.2")
	{
		output = remove_interference(input);
		output.save(outputFile);	
	}


        if(part == "3.3")
	{
 		if(argc != 5){cout<<"Error: Provide sigma as additional argument"<<endl;return -1;}
		 CImg<double> H = gaussianFilter(atof(argv[4]));
		
		double max_dim = findMaxDimension(input);

		input.resize(max_dim,max_dim,1,1,0);
		H.resize(max_dim,max_dim,1,1,0);		
		
		CImg<double> output(H.width(),H.height(),1,1,0);
		
		output =  fft_filter(input,H);
		output.resize(init_w,init_h,1,1,0);
		output.save(outputFile);
	}


	if(part == "3.4")
	{
 		if(argc != 5){cout<<"Error: Provide sigma as additional argument"<<endl;return -1;}
		 CImg<double> H = gaussianFilter(atof(argv[4]));
		
		double max_dim = findMaxDimension(input);

		input.resize(max_dim,max_dim,1,1,0);
		H.resize(max_dim,max_dim,1,1,0);		
		
		CImg<double> output(H.width(),H.height(),1,1,0);
        	CImg<double> op(H.width(),H.height(),1,1,0);
		
		op =  fft_filter(input,H);
		output = fft_defilter(op,H);	
		
		output.resize(init_w,init_h,1,1,0);
		output.save(outputFile);
	}

  } 
  catch(const string &err) {
    cerr << "Error: " << err << endl;
  }
}
//End of main function


//Part 3.1
CImg<double> fft_magnitude(const CImg<double>&fft_real,const CImg<double>&fft_imag)
{

	CImg<double> magnitude(fft_real,"xyz",0);

		for(int i=0;i<magnitude.width();i++)
			for(int j=0;j<magnitude.height();j++)
				magnitude(i,j,0) = log(sqrt(fft_real(i,j,0)*fft_real(i,j,0) + fft_imag(i,j,0)*fft_imag(i,j,0)));	

	return magnitude;
}
//end of Part 3.1


//Part 3.2
CImg<double> remove_interference(const CImg<double>& input)
{

	CImg<double> ip_real(input.width(),input.height(),1,1,0);
	CImg<double> ip_imag(input.width(),input.height(),1,1,0);

	CImg<double> op(input.width(),input.height(),1,1,0);
	CImg<double> out(input.width(),input.height(),1,1,0);
	fft(input,ip_real,ip_imag);

	op = fft_magnitude(ip_real,ip_imag);

	for(int i=0;i<op.width();i++){
		for(int j=0;j<op.height();j++){
			if(((i>150 && i< 165) && (j>150 && j<165))|| ((i>349 && i< 360) && (j>349 && j<360)) && op(i,j,0)>0 || op(i,j,0)==0){
			ip_real(i,j,0) = -0.0005;
		}
		}
	}
	
	ifft(ip_real,ip_imag,out);
	return out;

}
//end of Part 3.2


//Part 3.3
CImg<double> fft_filter(const CImg<double>&input,const CImg<double>&filter)
{

	CImg<double>ip_real(input.width(),input.height(),1,1,0);
	CImg<double>ip_imag(input.width(),input.height(),1,1,0);

	CImg<double>filter_real(input.width(),input.height(),1,1,0);
	CImg<double>filter_imag(input.width(),input.height(),1,1,0);

	CImg<double>op_real(input.width(),input.height(),1,1,0);
	CImg<double>op_imag(input.width(),input.height(),1,1,0);

	CImg<double>output_real(input.width(),input.height(),1,1,0);

	fft(input,ip_real,ip_imag);
	fft(filter,filter_real,filter_imag);

	for(int i=0;i<input.width();i++){
		for(int j=0;j<input.height();j++){	
			op_real(i,j,0,0) = ip_real(i,j,0,0)*filter_real(i,j,0,0) - ip_imag(i,j,0,0)*filter_imag(i,j,0,0);
			op_imag(i,j,0,0) = (ip_real(i,j,0,0)*filter_imag(i,j,0,0) + ip_imag(i,j,0,0)*filter_real(i,j,0,0));
		}
	}		

	ifft(op_real,op_imag,output_real);

	output_real.normalize(0,255);
	return output_real;
}
//End of Part 3.3


//Part 3.4
CImg<double> fft_defilter(const CImg<double>&input,const CImg<double>&filter)
{

        CImg<double>ip_real(input.width(),input.height(),1,1,0);
        CImg<double>ip_imag(input.width(),input.height(),1,1,0);

        CImg<double>filter_real(input.width(),input.height(),1,1,0);
        CImg<double>filter_imag(input.width(),input.height(),1,1,0);

        CImg<double>op_real(input.width(),input.height(),1,1,0);
        CImg<double>op_imag(input.width(),input.height(),1,1,0);
	
	CImg<double>denom(input.width(),input.height(),1,1,0);
        CImg<double>output_real(input.width(),input.height(),1,1,0);

        fft(input,ip_real,ip_imag);
        fft(filter,filter_real,filter_imag);

        for(int i=0;i<input.width();i++){
                for(int j=0;j<input.height();j++){
			denom(i,j,0,0) = filter_real(i,j,0,0)*filter_real(i,j,0,0) + filter_imag(i,j,0,0)*filter_imag(i,j,0,0);
			op_real(i,j,0,0) = (ip_real(i,j,0,0)*filter_real(i,j,0,0) + ip_imag(i,j,0,0)*filter_imag(i,j,0,0))/denom(i,j,0,0);
			op_imag(i,j,0,0) = (ip_imag(i,j,0,0)*filter_real(i,j,0,0) - ip_real(i,j,0,0)*filter_imag(i,j,0,0))/denom(i,j,0,0);
                }
        }

        ifft(op_real,op_imag,output_real);

        output_real.normalize(0,255);	
        return output_real;
}
//end of part 3.4
