// 计算几何 By 猛犸也钻地 @ 2012.08.21

/* 命名约定 //
   圆：圆心在u，一般情况下半径r大于等于0
   直线：经过点u和v的直线，u不重合于v
   射线：起点在u，途经点v，u不重合于v
   线段：起点在u，终点在v，u不重合于v
   散点集：点的可空集合
   多边形：至少有三个点，沿多边形的边依次排列，边不重合，图形不自交
   凸多边形：各内角均小于180度的多边形
   平面：由不共线的三点uvw所表示
// 所有函数都会默认传入的参数已满足上面的命名约定 */

#include <vector>
#include <cmath>
#include <utility>
#include <algorithm>
using namespace std;
using namespace rel_ops;

// typedef long long NUM;
typedef double NUM;
const NUM EPS = 1e-12, MAGIC = 2.71828e18;
// 因为有相对误差判断，所以EPS不要设得太宽

inline NUM sqr(NUM a){return a*a;}
inline NUM cmp(NUM a, NUM b){
    // return a-b; // 坐标为浮点数时，使用下面这行
    return fabs(a-b)>=EPS+fabs(a)*EPS?a-b:0;
}

//--------------------------------------------------------------------//

struct VEC {NUM x,y;} NOVEC = {MAGIC,MAGIC};
struct RAY {VEC u,v;} NORAY = {NOVEC,NOVEC};
struct CIR {VEC u; NUM r;} NOCIR = {NOVEC,MAGIC};

inline NUM sqr(const VEC& a){return sqr(a.x)+sqr(a.y);}
inline double abs(const VEC& a){return sqrt(sqr(a));}
inline NUM cmp(const VEC& a, const VEC& b){
    NUM at=cmp(a.x,b.x);
    return !at?cmp(a.y,b.y):at;
}

inline VEC operator +(const VEC& a, const VEC& b)
    {return (VEC){a.x+b.x,a.y+b.y};}
inline VEC operator -(const VEC& a, const VEC& b)
    {return (VEC){a.x-b.x,a.y-b.y};}
inline NUM operator *(const VEC& a, const VEC& b)
    {return a.x*b.y-a.y*b.x;}
inline NUM operator %(const VEC& a, const VEC& b)
    {return a.x*b.x+a.y*b.y;}
inline VEC operator -(const VEC& a){return (VEC){-a.x,-a.y};}
inline VEC operator ~(const VEC& a){return (VEC){-a.y,+a.x};}
inline VEC operator *(NUM u, const VEC& a){return (VEC){u*a.x,u*a.y};}
inline VEC operator *(const VEC& a, NUM u){return (VEC){a.x*u,a.y*u};}
inline VEC operator /(const VEC& a, NUM u){return (VEC){a.x/u,a.y/u};}
inline VEC operator /(const VEC& a, const VEC& b){return a%b/sqr(b)*b;}
inline bool operator ==(const VEC& a, const VEC& b){return !cmp(a,b);}
inline bool operator <(const VEC& a, const VEC& b){return cmp(a,b)<0;}

//  返回值        cmp_side             cmp_axis
//  ==  0       a和b相互平行    /    a和b相互垂直
//  <= -EPS     a在b的左手侧    /    a和b朝向相反（内角大于90度）
//  >= +EPS     a在b的右手侧    /    a和b朝向相同（内角小于90度）
NUM cmp_side(const VEC& a, const VEC& b){return cmp(a.x*b.y,+a.y*b.x);}
NUM cmp_axis(const VEC& a, const VEC& b){return cmp(a.x*b.x,-a.y*b.y);}

//--------------------------------------------------------------------//

// 求向量a长度缩放至u单位后的新向量，a不能是零向量
// 求向量a绕坐标原点o，逆时针转u度后的新向量
VEC resize(const VEC& a, NUM u){return u/abs(a)*a;}
VEC rotate(const VEC& a, NUM u)
    {return (VEC){cos(u)*a.x-sin(u)*a.y,sin(u)*a.x+cos(u)*a.y};}

// 点在直线上的投影(到直线的最近点)
// 点在圆周上的投影(到圆周的最近点)
VEC project(const VEC& p, const RAY& l){
    return (p-l.u)/(l.v-l.u)+l.u;
}
VEC project(const VEC& p, const CIR& c){
    if(!cmp(p,c.u)) return NOVEC;
    return resize(p-c.u,c.r)+c.u;
}

// 求两直线的交点
// 求直线与圆的交点，交线段的方向与原先直线相同
// 求两圆相交的交点，交线段的方向为圆心a到b连线方向逆指针转90度
// 求直线与凸多边形的交点，交线段的方向与原先直线相同，复杂度O(logn)
VEC intersect(const RAY& a, const RAY& b){
    VEC s=a.u-a.v,t=b.u-b.v;
    NUM at=cmp_side(s,t);
    if(!at) return NOVEC;
    return a.u+(b.u-a.u)*t/at*s;
}
RAY intersect(const RAY& l, const CIR& c){
    VEC s=l.u+(c.u-l.u)/(l.v-l.u);
    NUM at=cmp(c.r*c.r,sqr(s-c.u));
    if(at<0) return NORAY;
    VEC t=resize(l.v-l.u,sqrt(at));
    return (RAY){s-t,s+t};
}
RAY intersect(const CIR& a, const CIR& b){
    NUM l=sqr(b.u-a.u);
    NUM w=(1+(a.r*a.r-b.r*b.r)/l)*0.5;
    NUM e=cmp(a.r*a.r/l,w*w);
    if(e<0) return NORAY;
    VEC t=sqrt(e)*~(b.u-a.u);
    VEC s=a.u+w*(b.u-a.u);
    return (RAY){s-t,s+t};
}

// 判断三点是否共线
// 判断点在直线上的投影点，是否在线段上
// 判断点和线的位置关系，在外侧为0，在直线上为1或2(在线段上时为2)
// 判断点和圆的位置关系，在外侧为0，内部为1，边上为2
// 判断点和任意简单多边形的位置关系，在外侧为0，内部为1，边上为2
// 快速地判断点和凸多边形的位置关系，在外侧为0，内部为1，边上为2
// 判断两条线的位置关系，斜相交为0，垂直为1，平行为2，重合为3
bool collinear(const VEC& a, const VEC& b, const VEC& c){
    return !cmp_side(a-b,b-c);
}
bool seg_range(const VEC& p, const RAY& l){
    return cmp_axis(p-l.u,p-l.v)<=0;
}
int relation(const VEC& p, const RAY& l){
    if(cmp_side(p-l.u,p-l.v)) return 0;
    return cmp_axis(p-l.u,p-l.v)>0?1:2;
}
int relation(const VEC& p, const CIR& c){
    NUM at=cmp(sqr(c.r),sqr(c.u-p));
    return at?at<0?0:1:2;
}
int relation(const VEC& p, const vector<VEC>& u){
    int n=u.size(),ret=0;
    for(int i=0;i<n;i++){
        VEC s=u[i]-p,t=u[(i+1)%n]-p;
        if(t<s) swap(s,t);
        if(!cmp_side(s,t) && cmp_axis(s,t)<=0) return 2;
        if(cmp(s.x+p.x,p.x)<=0 && cmp(t.x+p.x,p.x)>0
        && cmp_side(s,t)>0) ret^=1;
    }
    return ret;
}
int relation_convex(const VEC& p, const vector<VEC>& u){
    int n=u.size(),l=0,r=n-1,o=cmp_side(u[1]-u[0],u[r]-u[0])<0?-1:1;
    if(relation(p,(RAY){u[0],u[1]})==2
    || relation(p,(RAY){u[0],u[r]})==2) return 2;
    while(l<r){
        int m=(l+r+1)/2;
        if(cmp_side(p-u[0],u[m]-u[0])*o<=0) l=m; else r=m-1;
    }
    if(!r || r==n-1) return 0;
    NUM at=cmp_side(p-u[r],u[r+1]-u[r])*o;
    return at?at<0:2;
}
int relation(const RAY& a, const RAY& b){
    NUM at=cmp_side(a.u-a.v,b.u-b.v);
    return at?!cmp_axis(a.u-a.v,b.u-b.v):!cmp_side(a.u-b.u,a.u-b.v)+2;
}

// 由ax+by+c=0构造直线
// 由直径上的两点构造一个圆
// 由三角形的顶点构造外接圆
RAY make_line(NUM a, NUM b, NUM c){
    if(!cmp(a,0) && !cmp(b,0)) return NORAY;
    else if(!cmp(a,0)) return (RAY){{0,-c/b},{1,-c/b}};
    else if(!cmp(b,0)) return (RAY){{-c/a,0},{-c/a,1}};
    return (RAY){{0,-c/b},{-c/a,0}};
}
CIR make_circle(const VEC& a, const VEC& b){
    return (CIR){(a+b)/2,abs(a-b)/2};
}
CIR make_circle(const VEC& a, const VEC& b, const VEC& c){
    if(!cmp_side(a-b,a-c)) return NOCIR;
    NUM x=(c-b)%(a-c),y=(c-b)*(a-b);
    VEC m=(x/y*~(a-b)+a+b)/2;
    return (CIR){m,abs(a-m)};
}

// 求三点的内切圆
// 求点到圆的两个切点，返回的切点分别在点到圆心连线方向的左侧和右侧
// 求两圆的两条公切线，切线段的方向与圆心a到b连线方向相同
//     默认是外公切线，若将其中的一个圆半径设为负数，则求出的是内公切线
CIR tangent_circle(const VEC& a, const VEC& b, const VEC& c){
    if(!cmp_side(a-b,a-c)) return NOCIR;
    NUM x=abs(b-c),y=abs(c-a),z=abs(a-b);
    VEC m=(a*x+b*y+c*z)/(x+y+z);
    return (CIR){m,fabs((m-a)*(a-b)*1.0/z)};
}
RAY tangent(const VEC& p, const CIR& c){
    NUM l=sqr(p-c.u),e=cmp(l,c.r*c.r);
    if(e<0) return NORAY;
    NUM x=c.r/sqrt(l),y=sqrt(e/l);
    VEC s=resize(p-c.u,1),t=~s;
    RAY lr={c.u+c.r*x*s-c.r*y*t,
            c.u+c.r*x*s+c.r*y*t};
    return lr;
}
pair<RAY,RAY> tangent(const CIR& a, const CIR& b){
    NUM o=a.r-b.r,l=sqr(b.u-a.u),e=cmp(l,o*o);
    if(e<0) return make_pair(NORAY,NORAY);
    NUM x=o/sqrt(l),y=sqrt(e/l);
    VEC s=resize(b.u-a.u,1),t=~s;
    RAY ll={a.u+a.r*x*s+a.r*y*t,
            b.u+b.r*x*s+b.r*y*t};
    RAY rr={a.u+a.r*x*s-a.r*y*t,
            b.u+b.r*x*s-b.r*y*t};
    return make_pair(ll,rr);
}

// 由散点集构造一个最小覆盖圆，期望复杂度O(n)
CIR min_covering_circle(vector<VEC> u){
    random_shuffle(u.begin(),u.end());
    int n=u.size(),i,j,k,z=1%n;
    CIR ret;
    for(ret=make_circle(u[0],u[z]),i=2;i<n;i++) if(!relation(u[i],ret))
    for(ret=make_circle(u[0],u[i]),j=1;j<i;j++) if(!relation(u[j],ret))
    for(ret=make_circle(u[i],u[j]),k=0;k<j;k++) if(!relation(u[k],ret))
        ret=make_circle(u[i],u[j],u[k]);
    return ret;
}

// 求散点集的二维凸包，并按逆时针顺序排列
// 若传入的点集不足以构成凸多边形，则返回的点集是退化后的点或线段
vector<VEC> convex_hull(vector<VEC> u){
    sort(u.begin(),u.end()); // 这两行是排序+去重，如果数据已经有保证
    u.erase(unique(u.begin(),u.end()),u.end()); // 则可省略相应的操作
    if(u.size()<3) return u;
    vector<VEC> c;
    for(size_t i=0,o=1,m=1;~i;i+=o){
        while(c.size()>m){
            VEC a=c.back()-c[c.size()-2];
            VEC b=c.back()-u[i];
            if(cmp_side(a,b)<0) break; // 改成<=0则保留共线点
            c.pop_back();
        }
        c.push_back(u[i]);
        if(i+1==u.size()) m=c.size(),o=-1; // 条件成立时切换至上凸壳
    }
    c.pop_back();
    return c;
}

/* 警告：下面这两个函数没有被正确地实现，等待修正中
   比赛时请不要使用这两个函数，有较高几率出错

// 求凸多边形上，朝某个方向看过去的最远点的编号，复杂度O(logn)
// 如果有多解，则返回相对于观测向量，在凸多边形上相对顺序更靠前的点
int apoapsis(const VEC& v, const vector<VEC>& u){
    if(!cmp((VEC){0,0},v)) return -1;
    int l=0,r=u.size()-1;
    NUM s=cmp_axis(u[r]-u[0],v);
    NUM t=cmp_axis(u[1]-u[0],v);
    if(s<=0 && t<=0) return !s?r:0;
    while(l<r){
        int m=(l+r)/2,e=cmp_axis(u[m]-u[0],v);
        if((e>=0 && e<cmp_axis(u[m+1]-u[0],v))
        || (e< 0 && t<0)) l=m+1; else r=m;
    }
    return r;
}

// 求直线与凸多边形的交点，交线段的方向与原先直线相同，复杂度O(logn)
RAY intersect(RAY l, const vector<VEC>& u){
    int n=u.size(),p,q,lo,hi;
    VEC o=l.v-l.u;
    if(cmp_side(u[1]-u[0],u[2]-u[0])<0) o=-o;
    NUM pt=cmp_side(o,u[p=apoapsis(~ o,u)]-l.u);
    NUM qt=cmp_side(o,u[q=apoapsis(~-o,u)]-l.u);
    if(pt*qt>0) return NORAY;
    for(;p<n+n;o=-o){ // 只执行两次，分别计算(p,q]和(q,p]段和直线的交点
        lo=p,hi=q+=n;
        swap(p+=n,q);
        while(lo<hi){
            int at=(lo+hi+1)/2;
            if(cmp_side(o,u[at%n]-l.u)>=0) lo=at; else hi=at-1;
        }
        if(!cmp_side(o,u[lo%n]-l.u)) l.u=u[lo%n];
        else l.u=intersect((RAY){u[lo%n],u[(lo+1)%n]},l);
        swap(l.u,l.v);
    }
    return l;
}

/*--------------------------------------------------------------------*/

struct TOR{NUM x,y,z;} NOTOR = {MAGIC,MAGIC,MAGIC};
struct SIG{TOR u,v;} NOSIG = {NOTOR,NOTOR};
struct PLN{TOR u,v,w;} NOPLN = {NOTOR,NOTOR,NOTOR};

inline NUM sqr(const TOR& a){return sqr(a.x)+sqr(a.y)+sqr(a.z);}
inline double abs(const TOR& a){return sqrt(sqr(a));}
inline NUM cmp(const TOR& a, const TOR& b){
    NUM at=cmp(a.x,b.x);
    if(!at) at=cmp(a.y,b.y);
    return !at?cmp(a.z,b.z):at;
}

inline TOR operator +(const TOR& a, const TOR& b)
    {return (TOR){a.x+b.x,a.y+b.y,a.z+b.z};}
inline TOR operator -(const TOR& a, const TOR& b)
    {return (TOR){a.x-b.x,a.y-b.y,a.z-b.z};}
inline TOR operator *(const TOR& a, const TOR& b)
    {return (TOR){a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
inline NUM operator %(const TOR& a, const TOR& b)
    {return a.x*b.x+a.y*b.y+a.z*b.z;}
inline TOR operator -(const TOR& a) {return (TOR){-a.x,-a.y,-a.z};}
inline TOR operator *(NUM u, const TOR& a)
    {return (TOR){u*a.x,u*a.y,u*a.z};}
inline TOR operator *(const TOR& a, NUM u)
    {return (TOR){a.x*u,a.y*u,a.z*u};}
inline TOR operator /(const TOR& a, NUM u)
    {return (TOR){a.x/u,a.y/u,a.z/u};}
inline TOR operator /(const TOR& a, const TOR& b){return a%b/sqr(b)*b;}
inline bool operator ==(const TOR& a, const TOR& b){return !cmp(a,b);}
inline bool operator <(const TOR& a, const TOR& b){return cmp(a,b)<0;}

// 下面两个函数类似于它们的二维版本，但cmp_side只能用于判定向量是否平行
int cmp_side(const TOR& a, const TOR& b){
    return cmp(a.y*b.z,a.z*b.y)
        || cmp(a.z*b.x,a.x*b.z)
        || cmp(a.x*b.y,a.y*b.x);
}
NUM cmp_axis(const TOR& a, const TOR& b){
    NUM x=a.x*b.x,y=a.y*b.y,z=a.z*b.z;
    if((x<0)==(y<0)) return cmp(x+y,-z);
    if((x<0)==(z<0)) return cmp(x+z,-y);
    // 注释掉上面两行可以提升程序速度，但有极小概率出现精度问题
    return cmp(y+z,-x);
}

//--------------------------------------------------------------------//

// 求平面c的法向量
// 求向量a长度缩放至u单位后的新向量，a不能是零向量
// 求向量a绕转轴向量o，逆时针转u度后的新向量
inline TOR normal(const PLN& c){return (c.v-c.u)*(c.w-c.u);}
TOR resize(const TOR& a, NUM u){return u/abs(a)*a;}
TOR rotate(const TOR& a, NUM u, const TOR& o)
    {return a*cos(u)+resize(o,1)*a*sin(u);}

// 点在直线上的投影(到直线的最近点)
// 点在平面上的投影(到平面的最近点)
TOR project(const TOR& p, const SIG& l){
    return (p-l.u)/(l.v-l.u)+l.u;
}
TOR project(const TOR& p, const PLN& c){
    return (c.u-p)/normal(c)+p;
}

// 求两直线的交点
// 求直线与平面的交点
// 求两平面的交线
TOR intersect(const SIG& a, const SIG& b){
    TOR s=b.u-b.v,p=s*(b.u-a.u);
    TOR t=a.u-a.v,q=s*t;
    if(cmp_axis(p,t) || !cmp_side(s,t)) return NOTOR;
    NUM at=cmp_axis(p,q);
    return a.u+(at?at<0?-1:1:0)*sqrt(sqr(p)/sqr(q))*t;
}
TOR intersect(const SIG& l, const PLN& c){
    TOR at=l.v-l.u,o=normal(c);
    if(!cmp_axis(o,at)) return NOTOR;
    return l.u+(c.u-l.u)%o/(at%o)*at;
}
SIG intersect(const PLN& a, const PLN& b){
    TOR o=normal(a);
    SIG s={b.u,b.v},t={b.u,b.w},r={b.v,b.w};
    s.u=intersect(cmp_axis(s.u-s.v,o)?s:r,a);
    t.u=intersect(cmp_axis(t.u-t.v,o)?t:r,a);
    return (SIG){s.u,t.u};
}

// 判断四点是否共面
// 判断三点是否共线
// 判断点在直线上的投影点，是否在线段上
// 判断点和线的位置关系，在外侧为0，在直线上为1或2(在线段上时为2)
// 判断点和面的位置关系，在面内为0，在正方向为1，在负方向为-1
// 判断两条线的位置关系，其他情况下为0，垂直为1，平行为2，重合为3
// 判断线和面的位置关系，斜相交为0，垂直为1，平行为2，线在面内为3
// 判断两平面的位置关系，斜相交为0，垂直为1，平行为2，重合为3
bool coplanar(const TOR& a, const TOR& b, const TOR& c, const TOR& d){
    return !cmp_axis(a-b,(a-c)*(a-d));
}
bool collinear(const TOR& a, const TOR& b, const TOR& c){
    return !cmp_side(a-b,b-c);
}
bool seg_range(const TOR& p, const SIG& l){
    return cmp_axis(p-l.u,p-l.v)<=0;
}
int relation(const TOR& p, const SIG& l){
    if(cmp_side(p-l.u,p-l.v)) return 0;
    return cmp_axis(p-l.u,p-l.v)>0?1:2;
}
int relation(const TOR& p, const PLN& c){
    NUM at=cmp_axis(p-c.u,normal(c));
    return at?at<0?-1:1:0;
}
int relation(const SIG& a, const SIG& b){ // 注意，异面垂直也算垂直
    NUM at=cmp_side(a.u-a.v,b.u-b.v);
    return at?!cmp_axis(a.u-a.v,b.u-b.v):!cmp_side(a.u-b.u,a.u-b.v)+2;
}
int relation(const SIG& l, const PLN& c){
    TOR o=normal(c),e=l.v-l.u;
    return cmp_axis(e,o)?!cmp_side(e,o):!cmp_axis(c.u-l.u,o)+2;
}
int relation(const PLN& a, const PLN& b){
    TOR p=normal(a),q=normal(b);
    return cmp_side(p,q)?!cmp_axis(p,q):!cmp_axis(a.u-b.u,p)+2;
}

// 由ax+by+cz+d=0构造平面
PLN make_plane(NUM a, NUM b, NUM c, NUM d){
    if(cmp(a,0)) return (PLN){{-d/a,0,0},{(-b-d)/a,1,0},{(-c-d)/a,0,1}};
    if(cmp(b,0)) return (PLN){{0,-d/b,0},{1,(-a-d)/b,0},{0,(-c-d)/b,1}};
    if(cmp(c,0)) return (PLN){{0,0,-d/c},{1,0,(-a-d)/c},{0,1,(-b-d)/c}};
    return NOPLN;
}

// 求散点集的三维凸包，返回每个三角面的顶点编号，复杂度O(nlogn)
// 从凸包外看，每个面的顶点都按逆时针的顺序排列，edge存储了邻面编号
// 比如edge::u[0]表示的是：face::u[0]至face::u[1]这条边所对应的邻面
// 传入的点集不能含有重点，若返回值为空集，则说明所有的点共面
// NUM的类型为long long时，坐标的范围不要超过10^6，建议使用浮点类型
struct TPL{int u[3];};
vector<TPL> convex_hull(const vector<TOR>& p){
    vector<TPL> face,edge;
    static vector<int> F[100005],G[100005*7]; // 注意设置最大结点数
    int n=p.size(),i,j,k;
    if(n<=3) return face;
    vector<int> u(n),v(4),at(n),go(n),by(n);
    for(i=0;i<n;i++) u[i]=i;
    random_shuffle(u.begin(),u.end());
    TOR a=p[u[0]]-p[u[1]],b;
    for(i=2;i<n;i++) if(cmp_side(a,b=p[u[0]]-p[u[i]])) break;
    for(j=i;j<n;j++) if(cmp_axis(a*b,p[u[0]]-p[u[j]])) break;
    if(i>=n || j>=n) return face;
    swap(u[i],u[2]),swap(u[j],u[3]);
    b=p[u[0]]+p[u[1]]+p[u[2]]+p[u[3]];
    for(i=0;i<4;i++){
        a=(p[u[i]]-p[u[j=(i+1)%4]])*(p[u[i]]-p[u[k=(i+2)%4]]);
        if(cmp_axis(p[u[i]]*4-b,a)<0) swap(j,k),a=-a;
        face.push_back((TPL){{u[i],u[j],u[k]}});
        edge.push_back((TPL){{(k+1)%4,(i+1)%4,(j+1)%4}});
        for(j=4;j<n;j++) if(cmp_axis(p[u[j]]-p[u[i]],a)>0)
            F[j].push_back(i),G[i].push_back(j);
    }
    for(i=4;i<n;F[i++].clear()){
        int x=n,m=F[i].size(),c=v.size();
        for(j=0;j<m;j++) v[F[i][j]]++;
        for(j=0;j<m;j++) if(v[F[i][j]]>0){
            v[F[i][j]]=-1234567890;
            for(k=0;k<3;k++){
                if(v[edge[F[i][j]].u[k]]) continue;
                at[x=face[F[i][j]].u[k]]=F[i][j];
                go[x]=k;
            }
        }
        if(x==n) continue;
        for(j=x,k=-1;k!=x;j=k){
            k=face[at[j]].u[(go[j]+1)%3];
            a=(p[j]-p[k])*(p[j]-p[u[i]]);
            int t=v.size(),w=edge[at[j]].u[go[j]];
            v.push_back(0);
            face.push_back((TPL){{j,k,u[i]}});
            edge.push_back((TPL){{w,t+1,t-1}});
            *find(edge[w].u,edge[w].u+3,at[j])=t;
            vector<int>::const_iterator o,z;
            z=set_union(G[at[j]].begin(),G[at[j]].end(),
                        G[w].begin(),G[w].end(),by.begin());
            for(o=by.begin();o!=z;++o)
                if(*o>i && cmp_axis(p[u[*o]]-p[u[i]],a)>0)
                    F[*o].push_back(t),G[t].push_back(*o);
        }
        edge[edge.back().u[1]=c].u[2]=edge.size()-1;
    }
    int m=v.size();
    for(i=j=0;i<m;G[i++].clear())
        if(!v[i]) face[j]=face[i],edge[j++]=edge[i];
    face.erase(face.begin()+j,face.end());
    edge.erase(edge.begin()+j,edge.end());
    return face;
}