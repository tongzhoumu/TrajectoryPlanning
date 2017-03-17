#include <cmath>
#include <algorithm>
using namespace std;
const int MAXN = 1000;
const double eps = 1e-8, PI = atan2(0, -1);
inline double sqr(double x){ return x * x; }
inline bool zero(double x){ return (x > 0 ? x : -x) < eps; }
inline int sgn(double x){ return (x > eps ? 1 : (x + eps < 0 ? -1 : 0)); }
struct point{
    double x, y;
    point(double x, double y):x(x), y(y) {}
    point() {}
    bool operator == (const point & a) const{ return sgn(x - a.x) == 0 && sgn(y - a.y) == 0; }
    bool operator != (const point & a) const{ return sgn(x - a.x) != 0 || sgn(y - a.y) != 0; }
    bool operator < (const point & a) const{ return sgn(x - a.x) < 0 || sgn(x - a.x) == 0 && sgn(y - a.y) < 0; }
    point operator + (const point & a) const{ return point(x + a.x, y + a.y); }
    point operator - (const point & a) const{ return point(x - a.x, y - a.y); }
    point operator * (const double & a) const{ return point(x * a, y * a); }
    point operator / (const double & a) const{ return point(x / a, y / a); }
    double operator * (const point & a) const{ return x * a.y - y * a.x; }//xmult
    double operator ^ (const point & a) const{ return x * a.x + y * a.y; }//dmult
    double length() const{ return sqrt(sqr(x) + sqr(y)); }
    point trunc(double a) const{ return (*this) * (a / length()); }
    point rotate(double ang) const{ point p(sin(ang), cos(ang)); return point((*this) * p, (*this) ^ p); }
    point rotate(const point & a) const{ point p(-a.y, a.x); p = p.trunc(1.0); return point((*this) * p, (*this) ^ p); }
};
bool isConvex(int n, const point * p){
    int i, s[3] = {1, 1, 1};
    for(i = 0; i < n && /*s[1] && */s[0] | s[2]; i++)
        s[sgn((p[(i + 1) % n] - p[i]) * (p[(i + 2) % n] - p[i])) + 1] = 0;
    return /*s[1] && */s[0] | s[2];
}//去掉注释即不允许相邻边共线
bool insideConvex(const point & q, int n, const point * p){
    int i, s[3] = {1, 1, 1};
    for(i = 0; i < n && /*s[1] && */s[0] | s[2]; i++)
        s[sgn((p[(i + 1) % n] - p[i]) * (q - p[i])) + 1] = 0;
    return /*s[1] && */s[0] | s[2];
}//去掉注释即严格在形内
inline bool dotsInline(const point & p1, const point & p2, const point & p3){
    return zero((p1 - p3) * (p2 - p3));
}//三点共线
inline int decideSide(const point & p1, const point & p2, const point & l1, const point & l2){
    return sgn((l1 - l2) * (p1 - l2)) * sgn((l1 - l2) * (p2 - l2));
}//点p1和p2,直线l1-l2,-1表示在异侧,0表示在线上,1表示同侧
inline bool dotOnlineIn(const point & p, const point & l1, const point & l2){
    return zero((p - l2) * (l1 - l2)) && (l1.x - p.x) * (l2.x - p.x) < eps && (l1.y - p.y) * (l2.y - p.y) < eps;
}//判点是否在线段及其端点上
inline bool parallel(const point & u1, const point & u2, const point & v1, const point & v2){
    return zero((u1 - u2) * (v1 - v2));
}//判直线平行
inline bool perpendicular(const point & u1, const point & u2, const point & v1, const point & v2){
    return zero((u1 - u2) ^ (v1 - v2));
}//判直线垂直
inline bool intersectIn(const point & u1, const point & u2, const point & v1, const point & v2){
    if(!dotsInline(u1, u2, v1) || !dotsInline(u1, u2, v2))
        return decideSide(u1, u2, v1, v2) != 1 && decideSide(v1, v2, u1, u2) != 1;
    else
        return dotOnlineIn(u1, v1, v2) || dotOnlineIn(u2, v1, v2) || dotOnlineIn(v1, u1, u2) || dotOnlineIn(v2, u1, u2);
}//判两线段相交,包括端点和部分重合
inline bool intersectEx(const point & u1, const point & u2, const point & v1, const point & v2){
    return decideSide(u1, u2, v1, v2) < 0 && decideSide(v1, v2, u1, u2) < 0;
}//判两线段相交,不包括端点和部分重合
inline bool insidePolygon(const point & q, int n, const point * p, bool onEdge = true){
    if(dotOnlineIn(q, p[n - 1], p[0])) return onEdge; for(int i = 0;i + 1 < n;i++) if(dotOnlineIn(q, p[i], p[i + 1])) return onEdge;
    #define getq(i) Q[(sgn(p[i].x-q.x)>0)<<1|sgn(p[i].y-q.y)>0]
    #define difq(a,b,i,j) (a==b?0:(a==((b+1)&3)?1:(a==((b+3)&3)?-1:(sgn((p[i]-q)*(p[j]-q))<<1))))
    int Q[4] = {2, 1, 3, 0}, oq = getq(n-1), nq = getq(0), qua = difq(nq, oq, n - 1, 0); oq = nq;
    for(int i = 1;i < n;i++){ nq = getq(i); qua += difq(nq, oq, i - 1, i); oq = nq; }
    return qua != 0;//象限环顾法,较好
    /*point q1; int i = 0, cnt = 0; const double OFFSET = 1e6;//坐标上限
    for(q1 = point(rand() + OFFSET, rand() + OFFSET);i < n;) for(i = cnt = 0;i < n && !dotsInline(q, q1, p[i]);i++) cnt += intersectEx(q, q1, p[i], p[(i + 1) % n]);
    return cnt & 1;*///考验rp的射线法
}//判点在任意多边形内
inline point intersection(const point & u1, const point & u2, const point & v1, const point & v2){
    return u1 + (u2 - u1) * (((u1 - v1) * (v1 - v2)) / ((u1 - u2) * (v1 - v2)));
}//求两直线交点,须预判是否平行
inline point ptoline(const point & p, const point & l1, const point & l2){
    point t = p; t.x += l1.y - l2.y; t.y += l2.x - l1.x;
    return intersection(p, t, l1, l2);
}//点到直线的最近点,注意l1不能等于l2
inline double disptoline(const point & p, const point & l1, const point & l2){
    return fabs((p - l2) * (l1 - l2)) / (l1 - l2).length();
}//点到直线距离,注意l1不能等于l2
inline point ptoseg(const point & p, const point & l1, const point & l2){
    point t = p; t.x += l1.y - l2.y; t.y += l2.x - l1.x;
    if(sgn((l1 - p) * (t - p)) * sgn((l2 - p) * (t - p)) > 0)
        return (p - l1).length() < (p - l2).length() ? l1 : l2;
    else
        return intersection(p, t, l1, l2);
}//点到线段的最近点,注意l1不能等于l2
inline double disptoseg(const point & p, const point & l1, const point & l2){
    point t = point(l1.y - l2.y, l2.x - l1.x);
    if(sgn((l1 - p) * t) * sgn((l2 - p) * t) > 0)
        return min((p - l1).length(), (p - l2).length());
    else
        return disptoline(p, l1, l2);
}//点到线段距离,注意l1不能等于l2
double fermentpoint(int m, point p[]){
    point u(0, 0), v;
    double step = 0, nowbest = 0, now, maxx = 0, maxy = 0;
    for(int i = 0; i < m; ++i) {
        u = u + p[i];
        maxx = max(maxx, fabs(p[i].x));
        maxy = max(maxy, fabs(p[i].y));
    }
    u = u / m;
    for(int i = 0; i < m; ++i) nowbest += (u - p[i]).length();
    for(step = maxx + maxy;step > 1e-10;step *= 0.97)//对结果有影响,注意调整
        for(int i = -1; i <= 1; i++)
            for(int j = -1; j <= 1; j++){
                v = u + point(i, j) * step; now = 0;
                for(int i = 0; i < m; ++i) now += (v - p[i]).length();
                if(now < nowbest){ nowbest = now; u = v; }
            }
    return nowbest;
}//模拟退火求费马点
void polygonCut(int & n, point * p, const point & l1, const point & l2, const point & side){
    int m = 0, i; point pp[MAXN];//尽量定义成全局变量
    for(i = 0; i < n; i++){
        if(decideSide(p[i], side, l1, l2) == 1) pp[m++] = p[i];
        if(decideSide(p[i], p[(i + 1) % n], l1, l2) < 1 && !(zero((p[i] - l2) * (l1 - l2)) && zero((p[(i + 1) % n] - l2) * (l1 - l2)))) 
            pp[m++] = intersection(p[i], p[(i + 1) % n], l1, l2);
    }
    for(n = i = 0; i < m; i++)
        if(!i || !zero(pp[i].x - pp[i - 1].x) || !zero(pp[i].y - pp[i-1].y)) p[n++] = pp[i];
    if(zero(p[n - 1].x - p[0].x) && zero(p[n - 1].y - p[0].y)) n--;
    if(n < 3) n = 0;
}//将多边形沿l1,l2确定的直线在side侧切割,保证l1,l2,side不共线
inline double Seg_area(const point & p1, const point & p2, const point & p0, double R){
    point tmp = (p0 - p1).rotate(p2 - p1);
    double d = -tmp.y, h1 = -tmp.x, h2 = h1 + (p2 - p1).length();
    if(d >= R || d <= -R) return R * R * (atan2(d, h1) - atan2(d, h2));
    double dh = sqrt(R * R - d * d);
    if(h2 < -dh || dh < h1) return R * R * (atan2(d, h1) - atan2(d, h2));
    double ret = 0;
    if(h1 < -dh) ret += atan2(d, h1) - atan2(d, -dh);
    if(h2 > dh) ret += atan2(d, dh) - atan2(d, h2);
    return ret * R * R + d * (min(h2, dh) - max(h1, -dh));
}//圆与线段交的有向面积
int graham(int n, point * p, point * ch, bool comEdge = false){
    if(n < 3){ for(int i = 0; i < n; i++) ch[i] = p[i]; return n; }
    const double e1 = comEdge ? eps : -eps; int i, j, k;
    sort(p, p + n); ch[0] = p[0]; ch[1] = p[1];
    for(i = j = 2; i < n; ch[j++] = p[i++]) while(j > 1 && (ch[j - 2] - ch[j - 1]) * (p[i] - ch[j - 1]) > e1) j--;
    ch[k = j++] = p[n - 2];
    for(i = n - 3; i > 0; ch[j++] = p[i--]) while(j > k && (ch[j - 2] - ch[j - 1]) * (p[i] - ch[j - 1]) > e1) j--;
    while (j > k && (ch[j - 2] - ch[j - 1]) * (ch[0] - ch[j - 1]) > e1) j--;
    return j;
}//求凸包,p会被打乱顺序,ch为逆时针,comEdge为true时保留共线点,重点会导致不稳定
 