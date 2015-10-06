#include "CImg.h"
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <string>
#include<vector>
#include<algorithm>

using namespace cimg_library;

using namespace std;


//Part 2 - Detecting edges and circles

CImg<double> detectEdge(CImg<double> input);

CImg<double> detectCircle(CImg<double> input);

CImg<double> sobelFilter(CImg<double> input);

CImg<double> filter(CImg<double> input, CImg<double> H);


//Part 3 - Finding circular regions

CImg<double> findCircularRegion(CImg<double> input,char* sub_part);

int main(int argc, char **argv){
    
    
    if(argc < 4){
        
        cout << "Insufficent number of arguments. Please see documentation" << endl;
        
        cout << "p3 problemID inputfile outputfile" << endl;
        
        return -1;
        
    }
    
    char* inputFile = argv[2];
    
    char* outputFile = argv[3];
    
    float sub_part = atof(argv[1]);
    int part = (int)sub_part;
    
    cout << "In: " << inputFile <<"  Out: " << outputFile << endl;
    
    
    
    CImg<double> input(inputFile);
    
    CImg<double> output;
    
    
    
    if(part == 2){
        
        cout << "# Problem 2: Detecting edges and circles" << endl;
        
        //if(input.spectrum() != 3){	cout << "INPUT ERROR: Input image is not a color image!" << endl;return -1;}
        
        output = detectEdge(input);
        
        output.save(outputFile);
        
    }
    
    
    else if(part == 3){
        
        cout << "# Problem 3 - Finding Circular regions" << endl;
        
        if(input.width() > 1000 || input.height() > 1000)
        {
            input.resize(0.1*input.width(),0.1*input.height());
            cout<<"here";
        }
        output = findCircularRegion(input,argv[1]);
        
        output.save(outputFile);
        
    }
    
    
    else if(part == 4){
        
        cout << "# Problem 4 - Coin counting" << endl;
        
        //output = countCoin(input, atoi(argv[4]));
        
        //output.save(outputFile);
        
    }
    
    
    
    else{
        
        cout << "Unknown input command" << endl;
        
    }
    
    return 0;
    
}



//Part 2 - Basic image operations

CImg<double> detectEdge(CImg<double> input){
    
    //Creates a new image with same size as the input initialized to all 0s (black)
    CImg<double> input_image1(input);
    CImg<double> input_image2(input);
    CImg<double> output(input.width(), input.height(), 1, 1);
    //CImg<double> binaryImage(input.width(), input.height(), 1, 1);
    
    
    CImg<double> filt;
    
    //Apply Sobel filter
    CImg<double> Kx(3, 3, 1, 1);
    CImg<double> Ky(3,3, 1, 1);
    CImg<double> Gx(input.width(),input.height(), 1, 1);
    CImg<double> Gy(input.width(),input.height(), 1, 1);
    CImg<double> G(input.width(),input.height(), 1, 1);
    
    int Mx[3][3]={{-1,0,1},{-2,0,2},{-1,0,1}};
    
    int My[3][3]={{-1,-2,-1},{0,0,0},{1,2,1}};
    
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            Kx(i,j,0,0)=Mx[i][j];
            Ky(i,j,0,0)=My[i][j];
            
            
        }
        
    }
    
    //Note: Normalize if necessary!
    
    
    Gx=input_image1.convolve(Kx);
    Gx.save("xcomp.jpg");
    Gy=input_image2.convolve(Ky);
    Gy.save("ycomp.jpg");
    
    
    //input to calculate gradient
    static int countmore=0,countless=0;
    for (int i = 0; i < input.width(); i++) {
        for (int j = 0; j < input.height(); j++){
            
            G(i,j,0,0)=sqrt(Gx(i,j,0,0)*Gx(i,j,0,0)+Gy(i,j,0,0)*Gy(i,j,0,0));
            
            if(G(i,j,0,0)<127){
                G(i,j,0,0)=0;
                countless++;
            }
            else{
                G(i,j,0,0)=255;
                countmore++;
            }
        }
    }
    cout<<endl<<"Less than 127:"<<countless<<endl<<"More than 127:"<<countmore;
    
    G.normalize(0,255);
    
    return G;
    
}

CImg<double> detectCircle(CImg<double> input){
    
    //Creates a new grayscale image with same size as the input initialized to all 0s (black)
    
    CImg<double> output(input, "xyzc", 0);
    
    
    
    //FIXME Do stuff
    return output;
    
}


// ------------------------------------------------------------------------------------------------------------------------


//Part 3 - Finding circular regions

class vertex_class{
public:
    double red;
    double green;
    double blue;
    int x;
    int y;
    int compid;
    
    //constructor
    vertex_class(double red,double blue,double green,int x,int y,int compid){
        this->red = red;
        this->green = green;
        this->blue = blue;
        this->x = x;
        this->y = y;
        this->compid = compid;
    }
};

class edge_class{
public: double weight;
    vertex_class *vd1;
    vertex_class *vd2;
    
    edge_class(vertex_class* v1,vertex_class* v2) {
        vd1 = v1;
        vd2 = v2;
        weight = sqrt(pow(vd1->red - vd2->red, 2) + pow(vd1->green - vd2->green, 2) + pow(vd1->blue - vd2->blue, 2));
    }
};

class graph_class{
public: vector<vertex_class*> v;
    vector<edge_class*> e;
    
    graph_class(int vertex_size) {
        v.resize(vertex_size);
    }
    
};

struct compare {
    bool operator ()(edge_class const *a, edge_class const *b) const {
        
        return (a->weight < b->weight);
    }
};

bool push_back(vector<vertex_class*> *vect, vertex_class* element);

CImg<double> gaussian(double sigma);

CImg<double> findCircularRegion(CImg<double> input,char* sub_part){
    
    //Creates a new image with same size as the input initialized to all 0s (black)
    int width = input.width();
    int height = input.height();
    
    
    if(sub_part[2] == 2)
        input.convolve(gaussian(1));
    input.save("convolved.jpg");
    
    CImg<double> output(input, "xyzc", 0);
    
    edge_class *e;
    graph_class g(width*height);
    
    // List of components
    vector<vector<vertex_class*> > v_list;
    v_list.resize(width*height);
    
    // List of edge weights corresponding to each component
    double v_list_max_edge_weight[width*height];
    
    for(int i = 0; i < width; i++)
    {
        for(int j = 0; j < height; j++)
        {
            // Creating a vertex for every pixel
            g.v[i+j*width] = new vertex_class(input(i,j,0,0),input(i,j,0,1),input(i,j,0,2),i,j,i+j*width);
            
            // Adding vertex to a component
            vector<vertex_class*> v;
            v.push_back(g.v[i+j*width]);
            
            // Adding component to component list
            v_list[j*width+i] = v;
            
            // Initializing max edge weight to 0 for every component
            v_list_max_edge_weight[i+j*width] = 0;
        }
    }
    
    //------------------------------------------------------------
    
    //Initialize the edge_class members------------------
    for(int i = 0; i < width; i++)
    {
        for(int j = 0;j < height; j++)
        {
            if(i+1 < width) {
                e = new edge_class(g.v[width*j + i], g.v[width*j + i+1]);
                g.e.push_back(e);
            }
            
            if(j+1 < height) {
                e = new edge_class(g.v[width*j + i], g.v[width*(j+1) + i]);
                g.e.push_back(e);
            }
        }
        
    }
    
    //----------------------------------------------------------------
    
    
    //Sort the edges
    sort(g.e.begin(), g.e.end(), compare());
    
    
    int K = 3;
    
    for (int iter_edge = 0;iter_edge < g.e.size();iter_edge++) {
        
        vertex_class *v1, *v2;
        
        int count = 0;
        double maxEdgeWtComp1, maxEdgeWtComp2, kValComp1, kValComp2;
        
        v1 = g.e[iter_edge]->vd1;
        v2 = g.e[iter_edge]->vd2;
        
        if(v1->compid != v2->compid && v_list[v1->compid].size() > 0 && v_list[v2->compid].size() > 0) {
            
            maxEdgeWtComp1 = v_list_max_edge_weight[v1->compid];
            maxEdgeWtComp2 = v_list_max_edge_weight[v2->compid];
            
            int compid1 = v1->compid;
            int compid2 = v2->compid;
            
            if(sub_part[2] == '1'){
                
                double threshold = 30;
                
                if (g.e[iter_edge]->weight <= threshold) {
                    
                    while(!v_list[compid2].empty())	{
                        vertex_class* v = v_list[compid2].back();
                        v->compid = compid1;
                        v_list[compid1].push_back(v);
                        v_list[compid2].pop_back();
                    }
                }
            }
            
            else if(sub_part[2] =='2'){
                
                kValComp1 = maxEdgeWtComp1 + K / v_list[v1->compid].size();
                kValComp2 = maxEdgeWtComp2 + K / v_list[v2->compid].size();
                
                if (g.e[iter_edge]->weight <= min(kValComp1, kValComp2)) {
                    
                    v_list_max_edge_weight[v1->compid] = max(g.e[iter_edge]->weight, max(maxEdgeWtComp1, maxEdgeWtComp2));
                    
                    while(!v_list[compid2].empty())	{
                        
                        vertex_class* v = v_list[compid2].back();
                        v->compid = compid1;
                        v_list[compid1].push_back(v);
                        v_list[compid2].pop_back();
                    }
                }
            }		
        }
    }	
    
    
    int count = 0;
    
    int r_out, g_out, b_out, x, y;
    
    for(int i = 0; i < v_list.size(); i++) {
        if(!v_list[i].empty())
            count++;
    }
    
    cout << "No of components: " << count << endl;
    
    
    //Randomly color each segment and output the image
    for(vector<vector<vertex_class*> >::size_type i = 0; i < v_list.size(); i++) {
        r_out = rand()%255;
        g_out = rand()%255;
        b_out = rand()%255;
        
        for(vector<vertex_class*>::size_type j = 0; j < v_list[i].size(); j++) {
            
            x = v_list[i][j]->x;
            y = v_list[i][j]->y;
            
            output(x,y,0,0) = r_out;
            output(x,y,0,1) = g_out;
            output(x,y,0,2) = b_out;
        }
        
    }
    
    return output;
    
}

bool push_back(vector<vertex_class*> *vect, vertex_class* element) {
    
    (*vect).resize((*vect).size()+1);
    (*vect)[(*vect).size()-1] = element;
    return true;
}


CImg<double> gaussian(double sigma){
    double PI = 3.1415;
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