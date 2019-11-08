#include<iostream>
#include<cstdio>
#include<fstream>
#include<vector>
using namespace std;
// what I need to do
// filter vertically and horzionatlly using a 3x3 and 5x5 using a Soble operator
// vertical and horzionatl images need to be stored separatedly
//
//
vector<unsigned char> applySobel(vector<vector<unsigned char>> img, vector<vector<char>> sobel_op, int max_val){

    int offset = sobel_op.size()/2;

    int img_dimension = img[0].size();
    vector<unsigned char> output_img;
    // I am going to use a linear array, since that is easier to write to file
    output_img.resize(img_dimension*img_dimension);

    for(int i = 0; i < img_dimension; ++i){

        for(int j = 0; j < img_dimension; ++j){
            // for each pixel, get surrounding piels and find the wieghted sum
            int result = 0; 
            for (int k = 0; k < sobel_op.size(); ++k){
                for (int l = 0; l < sobel_op[k].size(); ++l){

                    // the center pixel gets multiplied by sobel_op[1][1], so check if the bounds are fine
                    if ((i + k  - offset> 0 )&& (i + k - offset < img_dimension)){
                        if ((l + j  - offset> 0 )&&( i + k - offset < img_dimension)){

                            result += img[i - offset + k][j - offset + l]*sobel_op[k][l];

                        }

                    }
                }
            }
            // now that we have the wieghted sum, we know that no difference will be 128, so we can 
            // just check the bounds (< and > than 128) and use branching, or we could do something else
            // if we want to multiply by a coefficient, we could multiply by 127/max possible difference

            // the max difference is either 2*255 + 1*255 + 1*255 = 4*255 = 1020
            // min difference = 2*-255 - 1*255 - 1*255 = -1020
            // so divide by 1020 and add 128, and see what that does
            int positive_result = (result*255)/max_val;
            if (positive_result < 0)
                positive_result = -positive_result;

            output_img[i*img_dimension + j] = positive_result; 
            // if the difference result = 0, then the 128 value will be output to a file
            //
            // using the i*number + j indexing makes sure that the number ends up in the right place
        }
    }
    return output_img;
}

void toFile(vector<unsigned char> img, const char * filename ){

    // make a vector with all the rgb components
    int rgb_size = 3*img.size();
    
    vector<unsigned char> rgb_out(rgb_size, 0);
    for(int i = 0; i <  rgb_size; ++i){
        rgb_out[i] = img[i/3];
    }


    ofstream outFile(filename, ios::binary);
    if (outFile.is_open()){
        outFile.write((char*)rgb_out.data(), rgb_out.size());
        outFile.close();
    }
    else
        cerr << "File " << filename << " Failed to open for writing\n";


}

vector<unsigned char> chooseHighIntensity(vector<unsigned char> first, vector<unsigned char> second){
    int size = first.size();
    if (second.size() < size)
        size = second.size();

    vector<unsigned char> result;
    result.resize(size);

    for (int i = 0; i < size; ++i){
        if (first[i] > second[i])
            result[i] = first[i];
        else
            result[i] = second[i];

    }
    return result;
}
// takes any number above a threshhold (in the array) and sets it to 255
vector<unsigned char> raiseAboveThresh(vector<unsigned char> input, unsigned char thresh){

    vector<unsigned char> output;
    int size = input.size();

    output.resize(size);
    for (int i = 0; i < size; ++i){
        if (input[i] > thresh)
            output[i] = 255;
        else 
            output[i] = input[i];
    }
 return output;
}

int main(){

const char * fname = "input_image.raw";
    ifstream input_image (fname, ios::binary);
    int img_dimension = 500;
    if (input_image.is_open()){
        //
        // get file zie
        input_image.seekg(0, input_image.end);
        int img_size= input_image.tellg();
        input_image.seekg(0, input_image.beg);

        // 
        // the image I exported does not have rgb components, it just has the intensity for each pixel
        vector<vector<unsigned char>> img_matrix;
        img_matrix.resize(img_dimension);

        for(int i = 0; i < img_dimension; ++i){
            // set size for all image parts

            img_matrix[i].resize(img_dimension);
            
            // read in rgb parts
            input_image.read((char *)img_matrix[i].data(),img_dimension);
            //
        }

        input_image.close();

        unsigned char img_thresh = 20;

        // initialize Sobel operator
        vector<vector<char>> sobel_op;
        sobel_op.resize(3);

        // 0,0 is the top left of the operator matrix
            sobel_op[0].assign ( { 1, 0, -1  });

            sobel_op[1].assign({ 2, 0, -2  });

            sobel_op[2].assign({ 1, 0, -1  });

        // this loop multiplies the sobel_op with the image
        int max = (1+1+2)*255;
        vector<unsigned char> horiz = applySobel(img_matrix, sobel_op, max);


        // generate the vertical operator
        // 0,0 is the top left of the operator matrix
            sobel_op[0].assign({ 1, 2, 1  });

            sobel_op[1].assign({ 0, 0, 0  });

            sobel_op[2].assign({ -1, -2, -1  });

        vector<unsigned char> vert = applySobel(img_matrix, sobel_op, max);
        vector<unsigned char> edge = chooseHighIntensity(horiz, vert);
        toFile(edge, "3x3gradient.data");
        edge = raiseAboveThresh(edge, img_thresh);
        toFile(edge, "3x3edge.data");


        // generate the5x5 vertical operator
        sobel_op.resize(5);
            sobel_op[0].assign({ 1, 2, 4, 2,  1  });

            sobel_op[1].assign({ 2, 4, 8, 4, 2  });

            sobel_op[2].assign({0, 0, 0, 0 ,0 });

            sobel_op[3].assign({ -2, -4, -8, -4, -2  });

            sobel_op[4].assign({ -1, -2, -4, -2,  -1  });

        max = 0;
        for (int i = 0; i < sobel_op.size(); ++i){
            max += (sobel_op[0][i] + sobel_op[1][i])*255; 
        }
 vert = applySobel(img_matrix, sobel_op, max);
            
        sobel_op.resize(5);
            sobel_op[0].assign({ 1, 2, 0, -2,  -1  });

            sobel_op[1].assign({ 2, 4, 0, -4, -2  });

            sobel_op[2].assign({4, 8, 0, -8 ,-4 });

            sobel_op[3].assign({ 2, 4, 0, -4, -2  });

            sobel_op[4].assign({ 1, 2, 0, -2,  -1  });


        horiz= applySobel(img_matrix, sobel_op, max);
        edge = chooseHighIntensity(horiz, vert);

        toFile(edge, "5x5gradient.data");
        edge = raiseAboveThresh(edge, img_thresh);
        toFile(edge, "5x5edge.data");
    }
    else {
        cerr << "Couldn't open " << fname << " for reading\n";
    }
    
}


