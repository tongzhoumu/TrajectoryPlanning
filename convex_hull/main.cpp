#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "spline.h"
using namespace std;
using namespace av_planning;

#define AREA_LEN 750
#define AREA_WIDTH 500

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))

#define EPS 1e-12


struct Vec{ //double
	double x,y;
	Vec(){}
	Vec(double p1,double p2): x(p1),y(p2) {} 
	inline Vec operator - (const Vec &a) const {
		return Vec(x - a.x, y - a.y); 
	}
	inline Vec operator + (const Vec &a) const {
		return Vec(x + a.x, y + a.y); 
	}
	inline double operator * (const Vec &a) const {
		return x * a.y - y * a.x; //xmult
	}
	inline double operator ^ (const Vec &a) const {
		return x * a.x + y * a.y; //dmult
	}
	void resizeTo(double t) {
		double d = sqrt(x*x+y*y);
		double r = t/d;
		x *= r;
		y *= r;
	}
	void rotate(double ang) {  //anti-clock
		double sin_a = sin(ang), cos_a = cos(ang);
		x = x*cos_a - y*sin_a;
		y = x*sin_a + y*cos_a;
	}
	void prt() {
		cerr<<x<<", "<<y<<endl;
	}
};

struct Point { // int
	int x,y;
	Point(int p1, int p2): x(p1),y(p2) {}
	inline Point operator - (const Point &a) const {
		return Point(x - a.x, y - a.y); 
	}
	inline int operator * (const Point &a) const {
		return x * a.y - y * a.x; //xmult
	}
	inline bool operator <(const Point &a) const {
		return x==a.x ? y<a.y : x<a.x;
	}
};
typedef vector<Point> PointSet;

struct RoadSide {
	spline x, y;
	double len;
	double dist(Point a) const {
		double s,d;
		spline::getClosestPointOnCurveWithExtension(x,y,a.x,a.y,s,d);
		return d;
	}
	double closestPointS(Point a) const {
		double s,d;
		spline::getClosestPointOnCurveWithExtension(x,y,a.x,a.y,s,d);
		return s;
	}
} leftSide, rightSide;

struct ObstaclePolygon {
	int side; //-1->left, 0->center, 1->right
	Vec backHL, frontHL; // HL = half-line
	PointSet ch; //convex hull, from back to front
};


vector<PointSet> pointSets, convexHulls;
vector<ObstaclePolygon> obsPolys;
bool *fConnectLeftSide, *fConnectRightSide;

double vehLen,vehWidth;
bool area[AREA_WIDTH][AREA_LEN]; //x is width, y is len
int safeDist,safeDist_2;


void testCurve(char* name, spline &curve_x, spline &curve_y, double curve_len){
	double stride = 3;
	FILE* fp = fopen(name,"w");
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
}

void getCurve(RoadSide &side) {
	int nCurvePoints;
	scanf("%d",&nCurvePoints);
	//std::cerr<<"!!!!!!"<<nCurvePoints<<std::endl;
	std::vector<double> X(nCurvePoints), Y(nCurvePoints);
	for (int i=0;i<nCurvePoints;i++) {
		scanf("%lf%lf", &X[i],&Y[i]);
		//std::cerr<<X[i]<<" "<<Y[i]<<std::endl;
	}
	spline::fitCurve(X, Y, side.x, side.y, side.len); 
}

void getInput() {
	scanf("%lf%lf",&vehLen,&vehWidth);
	int tmp;
	for (int i=0;i<AREA_WIDTH;i++)
		for (int j=0;j<AREA_LEN;j++) {
			scanf("%d",&tmp);
			area[i][j] = tmp>0;
		}
	getCurve(leftSide);
	getCurve(rightSide);
	//testCurve("lside.txt",leftSide.x,leftSide.y,leftSide.len);
	//testCurve("rside.txt",rightSide.x,rightSide.y,rightSide.len);
	
	safeDist = (int)(vehWidth/2+0.5);
	safeDist_2 = safeDist*safeDist;
}

inline double dist_2(int x1,int y1,int x2,int y2) {
	return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
}

inline bool checkSafeDist(int x1,int y1,int x2,int y2) {
	return dist_2(x1,y1,x2,y2)>=safeDist*safeDist;
}

int binSearch(vector<int> &a, int x) { //find the index of minimal number >=x
	int l = 0, r = a.size()-1, m;
	if (r<0 || a[r]<x) return r+1;
	while (l<r) {
		m = (l+r)>>1;
		if (a[m]<x) {
			l = m+1;
		} else {
			r = m;
		}
	}
	return l;
}


void flood(int x,int y,bool (&visited)[AREA_WIDTH][AREA_LEN], vector<int> (&rows)[AREA_WIDTH]) {
	//cerr<<"flood: "<<x<<" "<<y<<endl;
	PointSet q;
	q.push_back(Point(x,y));
	visited[x][y] = 1;

	int t,p;
	for (int h=0;h<q.size();h++) {
		x = q[h].x;
		y = q[h].y;
		//cerr<<"----h = "<<h<<" ------"<<x<<" "<<y<<endl;
		for (int i=x;i<min(AREA_WIDTH-1,x+safeDist);i++) {
			//if (h==67) cerr<<i<<endl;
			t = binSearch(rows[i],max(0,y-safeDist));
			//if (h==67) {cerr<<"i = "<<i<<", t = "<<t<<endl;}
			for (int j=t;j<rows[i].size();j++) {
				if (rows[i][j]>y+safeDist) break;
				p = rows[i][j];
				if (!visited[i][p] && !checkSafeDist(x,y,i,p)) {
					//cerr<<"insert "<<i<<" "<<p<<endl;
					q.push_back(Point(i,p));
					visited[i][p] = 1;
				}
			}
		}
		//if (h==67) {cerr<<"this h end "<<t<<endl;}
	}
	pointSets.push_back(q);
}

bool insideRoad(int x, int y) {
	double s,dist;
	spline::getClosestPointOnCurveWithExtension(leftSide.x,leftSide.y,x,y,s,dist);
	Vec v1(leftSide.x(s)-x, leftSide.y(s)-y);
	spline::getClosestPointOnCurveWithExtension(rightSide.x,rightSide.y,x,y,s,dist);
	Vec v2(rightSide.x(s)-x, rightSide.y(s)-y);    
	/*
	cerr<<"P: "<<x<<" "<<y<<endl;
	cerr<<"L: "<<leftSide.x(s)<<" "<<leftSide.y(s)<<endl;
	cerr<<"R: "<<rightSide.x(s)<<" "<<rightSide.y(s)<<endl;
	v1.prt(); v2.prt();
	*/
    return (v1^v2)<=0;
}

void testPoints(char* name) {
	FILE *fp;
	fp = fopen(name,"w");
	//int tmp[AREA_WIDTH][AREA_LEN];
	//memset(tmp,0,sizeof(tmp));
	fprintf(fp,"hold on;\n");
	/*
	for (int i=0;i<pointSets.size();i++)
		for (int j=0;j<pointSets[i].size();j++){
			Point v = pointSets[i][j];
			fprintf(fp,"plot(%d,%d,\'k.\');\n",v.x,v.y);
		}
	*/
	for (int i=0;i<AREA_WIDTH;i++)
		for (int j=0;j<AREA_LEN;j++) {
			if (area[i][j] && insideRoad(i,j)) {
				fprintf(fp,"plot(%d,%d,\'k.\');\n",i,j);
			}
		}
	fclose(fp);
}

void clusterPoints() {
	vector<int> rows[AREA_WIDTH];
	bool visited[AREA_WIDTH][AREA_LEN];
	memset(visited,0,sizeof(visited));

	for (int i=0;i<AREA_WIDTH;i++)
		for (int j=0;j<AREA_LEN;j++) {
			if (area[i][j] && insideRoad(i,j)) {
				rows[i].push_back(j);
				//cerr<<i<<" "<<j<<endl;
			}
		}

	for (int i=0;i<AREA_WIDTH;i++)
		for (int j=0;j<rows[i].size();j++) {
			if (!visited[i][rows[i][j]]) 
				flood(i,rows[i][j],visited,rows);
		}

	//cerr<<"here"<<endl;
	//testPoints("points.txt");
}

inline bool turnLeft(const Point &a, const Point &b) { // b in the left of a
	return a*b>0; // change to >=0 means keep the points in the same line
}

inline bool turnLeft(const Vec &a, const Vec &b) {
	return a*b>0;
}

inline bool turnRight(const Vec &a, const Vec &b) {
	return a*b<0;
}

void graham(PointSet &p, PointSet &ch) {
	sort(p.begin(),p.end());
	if (p.size()<3) {
		ch = p;
		return;
	}
	for (int i=0,stepSize=1,m=2;i>=0;i+=stepSize) {
		while (ch.size()>=m){
			Point a = ch.back()-ch[ch.size()-2];
			Point b = p[i]-ch.back();
			if (turnLeft(a,b)) break;
			ch.pop_back();
		}
		ch.push_back(p[i]);
		if (i>=p.size()-1) {
			m = ch.size()+1;
			stepSize = -1; // switch to upper hull
		}
	}
	ch.pop_back();
}

void testCH(char *name) {
	FILE* fp = fopen(name,"w");
	for (int i=0;i<convexHulls.size();i++) {
		/*
		for (int j=1;j<convexHulls[i].size();j++) {
			fprintf(fp,"plot([%d %d],[%d %d],\'r\');\n",convexHulls[i][j-1].x,convexHulls[i][j].x,convexHulls[i][j-1].y,convexHulls[i][j].y);
		}
		fprintf(fp,"plot([%d %d],[%d %d],\'r\');\n\n",convexHulls[i].back().x,convexHulls[i][0].x,convexHulls[i].back().y,convexHulls[i][0].y);
		*/
		for (int j=0;j<convexHulls[i].size();j++) {
			fprintf(fp,"plot(%d,%d,\'r+\');\n",convexHulls[i][j].x,convexHulls[i][j].y);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}

void calcConvexHulls() {
	for (int i=0;i<pointSets.size();i++) {
		convexHulls.push_back(PointSet());
		graham(pointSets[i],convexHulls.back());
	}
	testCH("ch.txt");
}

void addStopLine() {
	cout << "Cannot pass!" <<endl;
	//something to do here in the future
}

bool connectRoadSide(RoadSide &s, PointSet &p) {
	for (int i=0;i<p.size();i++) {
		if (s.dist(p[i])<safeDist) return 1;
	}
	return 0;
}

void testStop(char* name) {
	FILE* fp = fopen(name,"w");
	for (int i=0;i<pointSets.size();i++) 
		if (fConnectLeftSide[i] && fConnectRightSide[i]) {
			for (int j=1;j<pointSets[i].size();j++) {
				fprintf(fp,"plot(%d,%d,\'r.\');\n",pointSets[i][j].x,pointSets[i][j].y);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}


bool checkPass() { //if true, fCLS & fCRS were filled
	bool flag = 1;
	fConnectLeftSide = new bool[pointSets.size()];
	fConnectRightSide = new bool[pointSets.size()];
	for (int i=0;i<pointSets.size();i++) {
		fConnectLeftSide[i] = connectRoadSide(leftSide,pointSets[i]);
		fConnectRightSide[i] = connectRoadSide(rightSide,pointSets[i]);
		if (fConnectLeftSide[i] && fConnectRightSide[i]) {
			flag = 0;
			//return 0;
		}
	}
	//return 1;
	//testStop("stop.txt");
	return flag;
}

void addHL(const PointSet &ch, const RoadSide &rs, bool (*turnBack)(const Vec&, const Vec&), ObstaclePolygon &obsPoly) {
	//----------find the back point & front point
	double s, min_s = rs.closestPointS(ch[0]), max_s = min_s;
	int min_idx = 0, max_idx = 0;
	for (int i=1;i<ch.size();i++) {
		s = rs.closestPointS(ch[i]);
		if (s<min_s) {
			min_s = s;
			min_idx = i;
		}
		if (s>max_s) {
			max_s = s;
			max_idx = i;
		}
	}
	//-----------calculate appropriate angles
	Vec vecBack(rs.x(min_s)-ch[min_idx].x,rs.y(min_s)-ch[min_idx].y);
	Vec vecFront(rs.x(max_s)-ch[max_idx].x,rs.y(max_s)-ch[max_idx].y);
	
	if (!turnBack(vecFront, vecBack)) {
		vecBack.resizeTo(1);
		vecFront.resizeTo(1);
		Vec vecMedium = vecBack + vecFront;
		vecFront = vecBack = vecMedium;
		vecFront.rotate(obsPoly.side * (M_PI/9)); //should keep convex
		vecBack.rotate(-1 * obsPoly.side * (M_PI/9));
	}
	
	obsPoly.frontHL = vecFront;
	obsPoly.backHL = vecBack;
	//-----------sort the points inside the road from back to front
	for (int i=min_idx;i!=max_idx;i-=obsPoly.side) {
		if (i<0) i = ch.size()-1;
		if (i>=ch.size()) i = 0;
		obsPoly.ch.push_back(ch[i]);
	}
	obsPoly.ch.push_back(ch[max_idx]);
}

void testPoly() {
	FILE* fp = fopen("poly.txt","w");
	for (int i=0;i<obsPolys.size();i++) {
		for (int j=1;j<obsPolys[i].ch.size();j++) {
			fprintf(fp,"plot([%d %d],[%d %d],\'r\');\n",obsPolys[i].ch[j-1].x,obsPolys[i].ch[j].x,obsPolys[i].ch[j-1].y,obsPolys[i].ch[j].y);
		}
		obsPolys[i].backHL.resizeTo(50);
		obsPolys[i].frontHL.resizeTo(50);
		fprintf(fp,"plot([%d %.2lf],[%d %.2lf],\'r\');\n",obsPolys[i].ch[0].x,obsPolys[i].ch[0].x+obsPolys[i].backHL.x,obsPolys[i].ch[0].y,obsPolys[i].ch[0].y+obsPolys[i].backHL.y);
		fprintf(fp,"plot([%d %.2lf],[%d %.2lf],\'r\');\n\n",obsPolys[i].ch.back().x,obsPolys[i].ch.back().x+obsPolys[i].frontHL.x,obsPolys[i].ch.back().y,obsPolys[i].ch.back().y+obsPolys[i].frontHL.y);
	}
	fclose(fp);
}

void buildPolygons() {
	for (int i=0;i<convexHulls.size();i++) {
		obsPolys.push_back(ObstaclePolygon());
		obsPolys[i].side = 0;
		if (fConnectLeftSide[i]) { //assert no stop
			obsPolys[i].side = -1;
			addHL(convexHulls[i], leftSide, turnLeft, obsPolys[i]);
		}
		else if (fConnectRightSide[i]) {
			obsPolys[i].side = 1;
			addHL(convexHulls[i], rightSide, turnRight, obsPolys[i]);
		}
		else {
			obsPolys[i].side = 0;
			obsPolys[i].ch = convexHulls[i];
		}
	}
	testPoly();
}


void test() {
	
	int x = 97, y = 97;
	double s,dist;
	spline::getClosestPointOnCurveWithExtension(leftSide.x,leftSide.y,x,y,s,dist);
	Vec v1(leftSide.x(s), leftSide.y(s));
	cerr<<s<<endl;
	v1.prt();
	spline::getClosestPointOnCurveWithExtension(rightSide.x,rightSide.y,x,y,s,dist);
	Vec v2(rightSide.x(s), rightSide.y(s));
	cerr<<s<<endl;
	v2.prt();
	
	/*
	int x = 100, y = 100;
	double s,dist;
	spline::getClosestPointOnCurveWithExtension(leftSide.x,leftSide.y,x,y,s,dist);
	Vec v1(leftSide.x(s), leftSide.y(s));
	v1.prt();
	spline::getClosestPointOnCurveWithExtension(rightSide.x,rightSide.y,x,y,s,dist);
	Vec v2(rightSide.x(s), rightSide.y(s));
	v2.prt();
	*/
}

int main() {

	getInput();
	clusterPoints();
	if (!checkPass()) {
		addStopLine();
		//return 0;
	}
	calcConvexHulls();
	buildPolygons();
	/*
	generateConstraints();
	*/
	test();

	return 0;
}