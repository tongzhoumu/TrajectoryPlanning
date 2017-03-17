#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

#define AREA_LEN 750
#define AREA_WIDTH 500

#define OB_LEN 10
#define N_SIDE_POINTS 20

double leftSide(double x) {
	return -0.0025*(x-600)*(x-600)+750;
}

double rightSide(double x) {
	x -= 20;
	return -0.0025*(x-600)*(x-600)+750;
}

int main() {
	FILE *fp;
	fp = fopen("obstacle.txt","w");

	int vehLen = 18, vehWidth = 8;
	printf("%d %d\n",vehLen,vehWidth);
	
	bool area[AREA_WIDTH][AREA_LEN];
	memset(area,0,sizeof(area));
	for (int st_i=0;st_i+OB_LEN<AREA_WIDTH;st_i+=90)
		for (int st_j=0;st_j+OB_LEN<AREA_LEN;st_j+=90)
			for (int i=0;i<OB_LEN;i++)
				for (int j=0;j<OB_LEN;j++) {
					int x = st_i+i, y = st_j+j;
					area[x][y] = true;
					//if (x%2==0 && y%2==0)
					if (y<1.8*(x-50)+150 && y>1.8*(x-50)-150)
						fprintf(fp,"plot(%d,%d,\'k.\');\n",x,y);
				}
	fclose(fp);
	for (int i=0;i<AREA_WIDTH;i++) {
		for (int j=0;j<AREA_LEN;j++)
			printf("%d ",(int)area[i][j]);
		printf("\n");
	}

	double x_stride = AREA_WIDTH/(double)N_SIDE_POINTS;
	std::vector<double> xc,yc;
	for (int i=0;i<N_SIDE_POINTS;i++){
		double x = i*x_stride;
		double y = leftSide(x);
		if (y>0 && y<AREA_LEN) {
			xc.push_back(x);
			yc.push_back(y);
		}
	}
	printf("%d\n",xc.size());
	for (int i=0;i<xc.size();i++){
		printf("%.6lf %.6lf\n",xc[i],yc[i]);
	}
	xc.clear();yc.clear();
	for (int i=0;i<N_SIDE_POINTS;i++){
		double x = i*x_stride;
		double y = rightSide(x);
		if (y>0 && y<AREA_LEN) {
			xc.push_back(x);
			yc.push_back(y);
		}
	}
	printf("%d\n",xc.size());
	for (int i=0;i<xc.size();i++){
		printf("%.6lf %.6lf\n",xc[i],yc[i]);
	}


	return 0;
}