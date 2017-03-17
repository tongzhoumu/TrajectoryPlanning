#include <cstdio>
#include <cstdlib>
#include <vector>
#include "spline.h"

int main(int argc, char** argv) {

   std::vector<double> X(5), Y(5);
   X[0]=0.1; X[1]=0.4; X[2]=1.2; X[3]=1.8; X[4]=1.4;
   Y[0]=0.1; Y[1]=0.7; Y[2]=0.6; Y[3]=1.1; Y[4]=0.9;
/*
   av_planning::spline s;
   s.setPoints(X,Y);    // currently it is required that X is already sorted

   double x=1.5;
   printf("spline at %f is %f\n", x, s(x));
*/
	av_planning::spline curve_x, curve_y;
	double curve_len;
	av_planning::spline::fitCurve(X, Y, curve_x, curve_y, curve_len); 
	//double len = 1.5;
	//printf("spline at %f is (%f, %f)\n", len, curve_x(len), curve_y(len));

	double stride = 0.03;
	FILE* fp = fopen("tmp.txt","w");
	fprintf(fp, "x = [%.3f", curve_x(0));
	for (double i=stride;i<curve_len;i+=stride) {
		fprintf(fp, ",%.3f", curve_x(i));
	}
	fprintf(fp, "];\n");
	fprintf(fp, "y = [%.3f", curve_y(0));
	for (double i=stride;i<curve_len;i+=stride) {
		fprintf(fp, ",%.3f", curve_y(i));
	}
	fprintf(fp, "];\n");
	fclose(fp);


   return EXIT_SUCCESS;
}