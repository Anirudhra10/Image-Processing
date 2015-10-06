// B490/B659 Project 4 skeleton code
//
// Compile with: "make"
//
// See assignment handout for command line and project specifications.


//Link to the header file
#include "CImg.h"
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Sift.h>
#include<ctime>
#include<math.h>
#include<dirent.h>
#include<fstream>
#include<string.h>

//Use the cimg namespace to access the functions easily
using namespace cimg_library;
using namespace std;


//Function declarations
CImg<double> Image_Warp(CImg<double> input);
CImg<double> Image_Match(CImg<double> ,CImg<double>, string);
CImg<double> closest_Image_Match(CImg<double> ,CImg<double>, string);
void calculate_mean(vector<double>, vector<double>, double*,double*);  
CImg<double> Image_Stitch(CImg<double>, CImg<double>, CImg<double>, double, double);

int desc_count = 0;

//Used for sorting images in 3.2
class image{
	public: string name;
	        double desc;
};

//Used to compare images in 3.2
struct compare {
    bool operator ()(image const *a, image const *b) const {
        
        return (a->desc > b->desc);
    }
};


int main(int argc, char **argv)
{
  try {
	srand(time(NULL));
    if(argc < 2)
      {
	cout << "Insufficent number of arguments; correct usage:" << endl;
	cout << "    p4 part_id ..." << endl;
	return -1;
      }
    
    string part = argv[1];
    
    if(part == "part2")
      {
      	char* inputFile = argv[2];
		char* outputFile = argv[3];

	    CImg<double> input(inputFile);
      	CImg<double> output;
      	
		output = Image_Warp(input);	
		output.save(outputFile);
      }
    
    else if(part == "part3.1")
      {
      	cout<<"Usage for 3.1: ./p4 part3.1 input_1.png input_2.png "<<endl;
      	if(argc != 5) { cout<<"Error: Invalid number of arguments"<<endl; return false;}
    	
    	char* inputFile = argv[2];
		char* inputFile2 = argv[3];
		char* outputFile = argv[4];
		
		CImg<double> input(inputFile);
		CImg<double> input2(inputFile2);
			
		CImg<double> output = Image_Match(input, input2,part);	
		output.save(outputFile);
		
     }   
     else if(part == "part3.2")
     {
     	cout<<"Usage for 3.2: ./p4 part3.2 query.png img1.png img2.png img3.png.... "<<endl;
     	if(argc < 4) { cout<<"Error: Invalid number of arguments"<<endl; return false;}
     	
     	char* inputFile = argv[2];
     	CImg<double> input(inputFile);
     	
     	vector<image*> vec;
     	image* im;
     	
     	for(int i = 3; i < argc; i++){
     		 
     		 im = new image();
     		 im->name = argv[i];
     		 
     		char* inputFile2 = argv[i];	 
     		CImg<double> input2(inputFile2);
     		
     		Image_Match(input,input2,part);	
     		
     		im->desc = desc_count;
     		vec.push_back(im);
     		
     	}
     	vector<image*> :: iterator iter;
     	
     	sort(vec.begin(),vec.end(),compare());
     	
     	for(iter = vec.begin(); iter != vec.end(); iter++){
     		cout<<(*iter)->name<<"     ";
     		cout<<(*iter)->desc<<endl;
     		
     	}
     	
     	
     }
     
     else if(part == "part3.3")
     {
     	cout<<"Usage for 3.3: ./p4 part3.3 ./attraction_images/image-name.png "<<endl;
     	if(argc < 3) { cout<<"Error: Invalid number of arguments"<<endl; return false;}
     	
     	char* inputFile = argv[2];
     	CImg<double> input(inputFile);
     	
     	vector<image*> vec;
     	image* im;
     	
     	DIR *pDIR;
        struct dirent *entry;
        int ct = 0;
        if( pDIR=opendir("/u/anramesh/anramesh-gpoornac-p4/attraction_images") ){
        	
        	std:: clock_t begin = clock();
        	
                while(entry = readdir(pDIR)){
                        if( strcmp(entry->d_name, ".") != 0 && strcmp(entry->d_name, "..") != 0 ){
                                
                                char* name = entry->d_name;
   									
     							string path = "./attraction_images/";
                                string full = path + name;
                                
                                string inputFile2 = full;
                                
                                im = new image();
     		 					im->name = name;
     		 					
     		 					cout << name<<"\t";
     		 					
                                CImg<double> input2(inputFile2.c_str());
                                Image_Match(input,input2,part);
                                ct++ ;
                                cout<<ct<<endl;
     							im->desc = desc_count;
     							vec.push_back(im);
   								
                        }
                }
                cout<<"Total time elapsed:"<<float(std::clock() - begin )/ CLOCKS_PER_SEC<<endl;
                closedir(pDIR);
        }
        vector<image*> :: iterator iter;
     	
     	sort(vec.begin(),vec.end(),compare());
     	
     	for(iter = vec.begin(); iter != vec.end(); iter++){
     		cout<<(*iter)->name<<"     ";
     		cout<<(*iter)->desc<<endl;
     		
     	}	
        	
     }  
      
    else if(part == "part4")
      {
		cout<<"Usage for 4: ./p4 part4 output.png input_1.png input_2.png"<<endl;
		if(argc != 5) { cout<<"Error: Invalid number of arguments"<<endl; return false;}
		
		char* inputFile = argv[3];
		char* inputFile2 = argv[4];
		char* outputFile = argv[2];
		
		CImg<double> input(inputFile);
		CImg<double> input2(inputFile2);
	
		CImg<double> output = closest_Image_Match(input, input2, part);	
		output.save(outputFile);
	
      }

  } 
  catch(const string &err) {
    cerr << "Error: " << err << endl;
  }
}

//Start of function bodies

//Function for Q2
CImg<double>Image_Warp(CImg<double>input){
	
	double width = input.width();
	double height = input.height();
	int x,y;
	double w;
	CImg<double> output(input,"xyzc",255);

	double m[3][3];
	double matrix[3][3];
	
	cout << "Enter the values of matrix row wise"<<endl;
	
	for(int i = 0; i < 3 ;i++)
		for(int j = 0; j < 3; j++)	
			cin >> m[i][j];
			
	
	/*double m[3][3] = {{0.907,0.258,-182},
				      {-0.153,1.44,58},
					  {-0.000306,0.000731,1}};*/
					  
//Method 2 - loop over destination pixels and copy from source			  
double determinant = m[0][0]*(m[1][1]*m[2][2] - m[2][1]*m[1][2]) -
					 m[0][1]*(m[1][0]*m[2][2] - m[2][0]*m[1][2]) +
					 m[0][2]*(m[1][0]*m[2][1] - m[2][0]*m[1][1]);

double inv_det = 1/determinant;

//Inverse of the matrix
matrix[0][0] = (m[1][1] * m[2][2] - m[2][1] * m[1][2]) * inv_det;
matrix[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_det;
matrix[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det;
matrix[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_det;
matrix[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det;
matrix[1][2] = (m[1][0] * m[0][2] - m[0][0] * m[1][2]) * inv_det;
matrix[2][0] = (m[1][0] * m[2][1] - m[2][0] * m[1][1]) * inv_det;
matrix[2][1] = (m[2][0] * m[0][1] - m[0][0] * m[2][1]) * inv_det;
matrix[2][2] = (m[0][0] * m[1][1] - m[1][0] * m[0][1]) * inv_det;


	for(int i = 0; i < width; i++){
		for(int j = 0; j < height; j++){
			
			w = (i * matrix[2][0] + j * matrix[2][1] + matrix[2][2]);
			x = (i * matrix[0][0] + j * matrix[0][1] + matrix[0][2])/w;			
			y = (i * matrix[1][0] + j * matrix[1][1] + matrix[1][2])/w;	
			
				
			if( x >= 0 && x < width && y >= 0 && y < height){
				output(i,j,0,0) = input(x,y,0,0);
				output(i,j,0,1) = input(x,y,0,1);
				output(i,j,0,2) = input(x,y,0,2);
			}
	
		}
	}

	return output;

}

//Function for Q3.1 and Q3.2
CImg<double> Image_Match(CImg<double> input, CImg<double> input2,string part){
		
	double color[] = { 255,0,0 }, num_of_desc = 0;			
	double sum = 0, val = 0, x1, y1, x2, y2;
	
	CImg<double> output(input,"xyzc",0);
	double max_height = max(input.height(),input2.height());
		
	output.resize(input.width()+input2.width(),max_height);

	for(int i = 0; i < input.width() + input2.width(); i++){
		for(int j = 0;j < max_height; j++){
				
			if(i < input.width() && j < input.height()){
				
				if(input.spectrum() == 3){
					output(i,j,0,0) = input(i,j,0,0);
					output(i,j,0,1) = input(i,j,0,1);
					output(i,j,0,2) = input(i,j,0,2);
				}
				else if(input.spectrum() == 1){
					output(i,j,0) = input(i,j,0);
				}
					
			}
			else if(i < input.width() && j >= input.height()){
				output(i,j,0,0) = 0;
				output(i,j,0,1) = 0;
				output(i,j,0,2) = 0;			
			
			}
			else if(i >= input.width() && j < input2.height()){
				
				if(input2.spectrum() == 3){
					output(i,j,0,0) = input2(i - input.width(),j,0,0);
					output(i,j,0,1) = input2(i - input.width(),j,0,1);
					output(i,j,0,2) = input2(i - input.width(),j,0,2);
				}
				else if(input2.spectrum() == 1){
					output(i,j,0) = input2(i - input.width(),j,0);
				}	
				
			}
			else{
				output(i,j,0,0) = 0;
				output(i,j,0,1) = 0;
				output(i,j,0,2) = 0;
			}
		}
	}
	
	vector<SiftDescriptor> descriptors, descriptors2;
	
		if(input.spectrum() != 1){
			CImg<double> gray = input.get_RGBtoHSI().get_channel(2);
			descriptors = Sift::compute_sift(gray);
		}
		
		if(input2.spectrum() != 1){
			CImg<double> gray2 = input2.get_RGBtoHSI().get_channel(2);
			descriptors2 = Sift::compute_sift(gray2);
		}
		
		for(int i=0; i<descriptors.size(); i++){
		   	for(int j=0; j<descriptors2.size(); j++){
					
				sum = 0;
				
				for(int l=0; l<128; l++){
					sum += pow(descriptors[i].descriptor[l] - descriptors2[j].descriptor[l],2); 	
				}
						
				val = sqrt(sum);	
				
				if(val < 100){	
					num_of_desc++;
							
					x1 = descriptors[i].col;     x2 = input.width() + descriptors2[j].col;
					y1 = descriptors[i].row;     y2 = descriptors2[j].row;
					        
					output.draw_line(x1,y1,x2,y2,color);
				}
			}	
        }
        
        if(part == "part3.1"){
        	cout<<"Number of matching descriptors are: "<<num_of_desc<<endl;
        	output.get_normalize(0,255);
			return output;
		}
		else if(part == "part3.2" || part == "part3.3"){
			
			desc_count = num_of_desc;
			
		}
		
}


//Function for Q4

CImg<double> closest_Image_Match(CImg<double>input, CImg<double>input2, string part){

	double color[] = { 255,50,60 };			
	double sum = 0,final_j = 0, val = 0, x1, y1, x2, y2;

	//For 4.2, 4.3, 4.4
	vector<double>x_translation;
	vector<double>y_translation; 
	
	CImg<double> output(input,"xyzc",255);
	double max_height = max(input.height(),input2.height());
		
	CImg<double> gray = input.get_RGBtoHSI().get_channel(2);
		vector<SiftDescriptor> descriptors = Sift::compute_sift(gray);
		
	CImg<double> gray2 = input2.get_RGBtoHSI().get_channel(2);
		vector<SiftDescriptor> descriptors2 = Sift::compute_sift(gray2);

for(int i = 0; i < descriptors.size(); i++){
	
	double min1 = 10000, min2 = 10500; 
	final_j = 0;
	
	for(int j = 0; j < descriptors2.size(); j++){
		sum = 0; 
		
		for(int l = 0; l < 128; l++){
			sum += pow(descriptors[i].descriptor[l] - descriptors2[j].descriptor[l],2); 	
		}
		val = sqrt(sum);
		
		if( val < min1){
			min2 = min1;
			min1 = val;
			final_j = j;
		}
		else if( val >= min1 && val < min2){
			min2 = val;
		}
	}
		double quotient = min1/min2;
		if( quotient < 0.5 ){
			
			x1 = descriptors[i].col;     x2 = input.width() + descriptors2[final_j].col;
			y1 = descriptors[i].row;     y2 = descriptors2[final_j].row;
			
			x_translation.push_back(descriptors2[final_j].col - x1);
			y_translation.push_back(y2 - y1);		        
			
		}
}

//part 4.2, 4.3 and 4.4
double dx = 0, dy = 0;

calculate_mean(x_translation, y_translation, &dx, &dy); 

//4.4
output.resize(input.width() + fabs(dx), input.height() + floor(fabs(dy)));

cout<<"Part 4.4 - Image stitching "<<endl; 
output = Image_Stitch(input, input2, output, floor(dx), floor(dy));

return output;
}

//function for 4.2 and 4.3
void calculate_mean(vector<double>x_translation, vector<double>y_translation, double *dx, double *dy){
	
	double x_sum = 0, y_sum = 0,count = 0;
	
	for(int i = 0; i < x_translation.size(); i++){
		count++;
		x_sum += x_translation[i];
		y_sum += y_translation[i];
	}
	
	//4.2 - Mean
	cout<<"Mean displacement: "<<x_sum/count <<"    "<< y_sum/count <<endl;
	
	int inlier_ct = 0, max_inlier_ct_pos = 0;
	
	for(int i =0; i < 10; i++){
	
		int random_num = rand() % (int)count;	
		inlier_ct = 0;

		for(int j = 0; j < x_translation.size(); j++){
	
			if( j != random_num && x_translation.at(random_num) == x_translation.at(j) && y_translation.at(random_num) == y_translation.at(j)  )
				inlier_ct++;
			
		}
		
		if(inlier_ct >= max_inlier_ct_pos)
			max_inlier_ct_pos = random_num;
				
	}
	//4.3 - position is obtained
		cout<< "final position: "<<max_inlier_ct_pos <<endl;
		
		*dx = x_translation[max_inlier_ct_pos];
		*dy = y_translation[max_inlier_ct_pos];
				
} 

CImg<double> Image_Stitch(CImg<double> input, CImg<double> input2, CImg<double> output, double dx, double dy){

double width = input.width();
double height = input.height();

if(dx > 0){
	CImg<double> temp = input;
	input = input2;
	input2 = temp;
	
	dy = -1 * dy;
}

dx = fabs(dx);
for(int i = 0; i < width; i++){
	for(int j = 0; j < height; j++){
		
		if(dy >= 0){ 
			output(i, j + dy, 0, 0) = input(i, j, 0, 0);
			output(i, j + dy, 0, 1) = input(i, j, 0, 1);
			output(i, j + dy, 0, 2) = input(i, j, 0, 2);
				
			output(i + dx, j, 0, 0) = input2(i, j, 0, 0);
			output(i + dx, j, 0, 1) = input2(i, j, 0, 1);
			output(i + dx, j, 0, 2) = input2(i, j, 0, 2);	
		}
		else{ 
			output(i, j, 0, 0) = input(i, j, 0, 0);
			output(i, j, 0, 1) = input(i, j, 0, 1);
			output(i, j, 0, 2) = input(i, j, 0, 2);
				
			output(i + dx, j + fabs(dy), 0, 0) = input2(i, j, 0, 0);
			output(i + dx, j + fabs(dy), 0, 1) = input2(i, j, 0, 1);
			output(i + dx, j + fabs(dy), 0, 2) = input2(i, j, 0, 2);
		}
	}
}

if(dx != 0 || dy != 0){

	for(int i = dx; i < width; i++){
		for(int j = fabs(dy); j < height; j++){
		
			if( dy >= 0){
				output(i, j, 0, 0) = (input(i, j - dy, 0, 0) + input2(i - dx, j, 0, 0))/2;
				output(i, j, 0, 1) = (input(i, j - dy, 0, 1) + input2(i - dx, j, 0, 1))/2;
				output(i, j, 0, 2) = (input(i, j - dy, 0, 2) + input2(i - dx, j, 0, 2))/2;
			}
			else{
				output(i, j, 0, 0) = (input(i, j, 0, 0) + input2(i - dx, j - fabs(dy), 0, 0))/2;
				output(i, j, 0, 1) = (input(i, j, 0, 1) + input2(i - dx, j - fabs(dy), 0, 1))/2;
				output(i, j, 0, 2) = (input(i, j, 0, 2) + input2(i - dx, j - fabs(dy), 0, 2))/2;
			}
		}
	
	}
}
return output;
}
