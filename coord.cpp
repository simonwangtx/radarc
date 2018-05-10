#include "coord.h"


void RadarConv::init(double B0,double L0,double H0)
{
    B=B0;L=L0;H=H0;
    alpha = 1/298.257222101;
    a = 6378137.0000;
    b = a*(1-alpha);
    ee = (a*a-b*b)/(a*a);
    ee2 = (a*a-b*b)/(b*b);
    theR = a- 2*a*alpha; //WGS84的基本参数

    siteR = a*sqrt(1-ee)/(1-ee*sin(B0)*sin(B0));
    N=a/(1-ee*sin(B0)*sin(B0));  //各纬度的地球半径

    //--地心坐标至雷达站坐标转换矩阵
    G00=(-sin(B0)*cos(L0));
    G01=(-sin(B0)*sin(L0));
    G02=(cos(B0));

    G10=(-sin(L0));
    G11=(cos(L0));
    G12=0;

    G20=(cos(B0)*cos(L0));
    G21=(cos(B0)*sin(L0));
    G22=(sin(B0));

    //--雷达站坐标至地心坐标转换矩阵
    R00=(-sin(B0)*cos(L0));
    R01=(-sin(L0));
    R02=(cos(B0)*cos(L0));

    R10=(-sin(B0)*sin(L0));
    R11=(cos(L0));
    R12=(cos(B0)*sin(L0));

    R20=(cos(B0));
    R21=0;
    R22=(sin(B0));

    cosB0=cos(B0);
    sinB0=sin(B0); //预计算值
    Bc=B0-(ee/2+5*ee*ee/24+3*ee*ee*ee/32)*sin(2.0*B0)+(5*ee*ee/48+7*ee*ee*ee/80)*sin(4.0*B0);
    Re=a*cos(B0)/((1-ee*sin(B0)*sin(B0))*cos(Bc));

    //-------------------
    RadarXYZ = Geo2XYZ((GEO)(*this));
};

//大地坐标转地心直角坐标
XYZ RadarConv::Geo2XYZ(GEO in) //X指向本初子午线
{
    XYZ GeoXYZ;
    double Nr = a/sqrt(1-ee*sin(in.B)*sin(in.B));
    GeoXYZ.x = (Nr+in.H)*cos(in.B)*cos(in.L);
    GeoXYZ.y = (Nr+in.H)*cos(in.B)*sin(in.L);
    GeoXYZ.z = (Nr*(1-ee)+in.H)*sin(in.B);
    return GeoXYZ;
};

//大地坐标转站心坐标
XYZ RadarConv::Geo2Site(GEO in)//X指向北极
{
    XYZ inXYZ = Geo2XYZ(in);
    XYZ outXYZ;

    double dx=(inXYZ.x-RadarXYZ.x);
    double dy=(inXYZ.y-RadarXYZ.y);
    double dz=(inXYZ.z-RadarXYZ.z);

    outXYZ.x=G00*dx+G01*dy+G02*dz;
    outXYZ.y=G10*dx+G11*dy+G12*dz;
    outXYZ.z=G20*dx+G21*dy+G22*dz;

    return outXYZ;
}

//大地坐标转站心坐标
XYZ RadarConv::Geo2Map(GEO in)//Y指向北极，X指向正东，适合中国地图投影
{
    XYZ inXYZ = Geo2XYZ(in);
    XYZ outXYZ;

    double dx=(inXYZ.x-RadarXYZ.x);
    double dy=(inXYZ.y-RadarXYZ.y);
    double dz=(inXYZ.z-RadarXYZ.z);

    outXYZ.x=G00*dx+G01*dy+G02*dz;
    outXYZ.y=G10*dx+G11*dy+G12*dz;
    outXYZ.z=G20*dx+G21*dy+G22*dz;

    inXYZ.x= outXYZ.y;  //按照中国地图坐标系,变换坐标参数
    inXYZ.y= outXYZ.x;
    inXYZ.z= outXYZ.z;
    return inXYZ;
}


//地心直角坐标转大地坐标
GEO RadarConv::XYZ2Geo(XYZ in) //要求输入的位置x,y不为零
{
    GEO geo;
    double length=sqrt(in.x*in.x+in.y*in.y);
    //if(length==0.0) throw "XYZ2GEO math error!! input length=0";
    double theta = atan(a*in.z/(b*length));
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double tmpvalue=length-a*ee*cos_theta*cos_theta*cos_theta;
    if(tmpvalue==0.0) geo.B=PI/2.0;
    else geo.B = atan((in.z+b*ee2*sin_theta*sin_theta*sin_theta)/tmpvalue);
    if(in.x==0.0) geo.L=(in.y>0.0)?PI/2:3.0*PI/2;
    else geo.L = atan(in.y/in.x);
    double Nr = a/sqrt(1-ee*sin(geo.B)*sin(geo.B)); //分母不会为零
    geo.H = sqrt(in.x*in.x+in.y*in.y)/cos(geo.B)-Nr;
    return geo;
    /*
    GEO geo2;  //按1米精度迭代法计算值与上述解析法差不多，解析法更好
    double preB=PI/4.0;
    //double Nr;
    double delta;
    double r=sqrt(in.x*in.x+in.y*in.y);
    geo2.H=0.0;
    do
    {
    Nr=a/sqrt(1-ee*sin(preB)*sin(preB));
    geo2.B=atan(in.z*(Nr+geo2.H)/(r*(Nr*(1.0-ee)+geo2.H)));
    geo2.H=r/cos(geo2.B)-Nr;
    delta=(geo2.B-preB)*Nr;
    preB=geo2.B;
    }while(delta>=1.0 || delta<=-1.0);
    if(in.x==0.0) geo2.L=(in.y>0.0)?PI/2:3.0*PI/2;
    else geo2.L = atan(in.y/in.x);
    return geo2;
    */
};

//大地坐标转雷达坐标
Radar RadarConv::Geo2Radar(GEO in)
{
    XYZ inXYZ = Geo2XYZ(in);
    XYZ outXYZ;

    double dx=(inXYZ.x-RadarXYZ.x);
    double dy=(inXYZ.y-RadarXYZ.y);
    double dz=(inXYZ.z-RadarXYZ.z);

    outXYZ.x=G00*dx+G01*dy+G02*dz;
    outXYZ.y=G10*dx+G11*dy+G12*dz;
    outXYZ.z=G20*dx+G21*dy+G22*dz;

    //转换为雷达位置的雷达直角坐标,注意，直角坐标系的X轴指向正北，Y轴指向正东
    Radar rt;
    rt.x= outXYZ.x;
    rt.y= outXYZ.y;
    rt.z= outXYZ.z;

    rt.range = sqrt(outXYZ.x*outXYZ.x+outXYZ.y*outXYZ.y+outXYZ.z*outXYZ.z);
    rt.rangeXY = sqrt(outXYZ.x*outXYZ.x+outXYZ.y*outXYZ.y);
    if(outXYZ.x==0.0)
        rt.theta=PI/2;
    else
    {
        rt.theta = atan(outXYZ.y/outXYZ.x);  //正北顺时针角度
        if(outXYZ.y>0.0 && outXYZ.x<0.0) rt.theta+=PI;//第二相限
        if(outXYZ.y<0.0 && outXYZ.x<0.0) rt.theta+=PI;//第三相限
        if(outXYZ.y<0.0 && outXYZ.x>0.0) rt.theta+=2*PI;//第四相限
    }
    rt.level=(int)(in.H/0.3048+0.5); //注意，level的单位是英尺
    return rt;
}


//将雷达斜距,方位角,ModeC转换为雷达站心切平面,x,y,z值
void RadarConv::RadarRange(Radar &in)
{
    /*
    if(in.range<10000) //10Km以内，地球变形较小，使用直接计算避免三角计算误差
    {
    in.z=in.level*0.3048;
    if(in.range>in.z)
    in.rangeXY=sqrt(in.range*in.range-in.z*in.z);
    else
    in.rangeXY=in.range;
    in.x=in.rangeXY*cos(in.theta+RadarBias);
    in.y=in.rangeXY*sin(in.theta+RadarBias);
    return;
    }

    double Rs=N+H;       //雷达站高度
    double Rr=N+in.level*0.3048; //雷达ModeC高度
    double costheta=(Rs*Rs+Rr*Rr-(in.range+RangeDelta)*(in.range+RangeDelta))/(2.0*Rs*Rr);
    double pulldown=(Rs/costheta)-N;
    in.z=(in.level*0.3048-pulldown)*costheta;
    if(in.range>in.z)
    in.rangeXY=sqrt(in.range*in.range-in.z*in.z);
    else
    in.rangeXY=in.range;
    in.x=in.rangeXY*cos(in.theta+RadarBias);
    in.y=in.rangeXY*sin(in.theta+RadarBias);
    */

    //转换为雷达位置的雷达直角坐标,注意，直角坐标系的Y轴指向正北，X轴指向正东
    //雷达送入的theta是与正北，也就是Y轴的夹角

    double F=in.level*0.3048;
    double r=in.range;
    double beta=acos((2.0*N*(H-F)+(H-F)*(H-F)+r*r)/(2.0*(N+H)*r));
    in.rangeXY=r*cos(beta-PI/2);
    in.z=r*sin(beta-PI/2);
    in.x=in.rangeXY*cos(in.theta+RadarBias);
    in.y=in.rangeXY*sin(in.theta+RadarBias);
}

//雷达坐标转大地坐标，转换均是基于雷达站心切平面坐标系的x,y值计算
//对于雷达站心切平面来说，x,y值计算应该考虑pull-down影响
GEO RadarConv::Radar2Geo(Radar &in)
{

    /*
    //极射投影坐标变换精度比地心坐标系低不少(200Km大约能到200-300m)，因此不再采用
    GEO rt;
    double Rs=N+H;       //雷达站高度
    double Rr=N+in.level*0.3048; //雷达ModeC高度
    double theta=(Rs*Rs+Rr*Rr-(in.range+RangeDelta)*(in.range+RangeDelta))/(2.0*Rs*Rr);
    if(theta>1.0) return rt; //代入值错误，无法计算
    theta=acos(theta);
    double pr=2.0*N*tan(theta/2.0);
    double costheta=cos(theta);
    double sintheta=sin(theta);
    double cosintheta=cos(in.theta+RadarBias);
    double sinintheta=sin(in.theta+RadarBias);  //使用极射投影算法，首先得到极射平面参数，再转回大地坐标

    rt.B=costheta*sinB0+cosintheta*sintheta*cosB0;
    if(rt.B<0.0 || rt.B>1.0) return rt; //代入值错误，无法计算
    rt.B=asin(rt.B);
    rt.L=sinintheta*sintheta/(costheta*cosB0-cosintheta*sintheta*sinB0);
    rt.L=L+atan(rt.L);
    rt.H=in.level*0.3048;//雷达ModeC高度
    */
    XYZ outXYZ;
    RadarRange(in); //转换到雷达站心切平面坐标系

    outXYZ.x=R00*in.x+R01*in.y+R02*in.z;
    outXYZ.y=R10*in.x+R11*in.y+R12*in.z;
    outXYZ.z=R20*in.x+R21*in.y+R22*in.z; //矩阵旋转

    outXYZ.x+=RadarXYZ.x; //平移
    outXYZ.y+=RadarXYZ.y;
    outXYZ.z+=RadarXYZ.z;//得到地心直角坐标系


    GEO rt2=XYZ2Geo(outXYZ); //转换到大地坐标系
    if(rt2.L<0) rt2.L+=PI; //解决负经度问题，忘记了为什么加180度，但这样才是对的，也许在西半球会有问题

    return rt2;


}

GEO RadarConv::RadarXY2Geo(Radar &in)
{
    XYZ outXYZ;
    outXYZ.x=R00*in.x+R01*in.y+R02*in.z;
    outXYZ.y=R10*in.x+R11*in.y+R12*in.z;
    outXYZ.z=R20*in.x+R21*in.y+R22*in.z; //矩阵旋转

    outXYZ.x+=RadarXYZ.x; //平移
    outXYZ.y+=RadarXYZ.y;
    outXYZ.z+=RadarXYZ.z;//得到地心直角坐标系


    GEO rt2=XYZ2Geo(outXYZ); //转换到大地坐标系
    if(rt2.L<0) rt2.L+=PI; //解决负经度问题，忘记了为什么加180度，但这样才是对的，也许在西半球会有问题

    return rt2;


}

GEO RadarConv::Site2Geo(XYZ in)
{
    XYZ outXYZ;

    outXYZ.x=R00*in.x+R01*in.y+R02*in.z;
    outXYZ.y=R10*in.x+R11*in.y+R12*in.z;
    outXYZ.z=R20*in.x+R21*in.y+R22*in.z; //矩阵旋转

    outXYZ.x+=RadarXYZ.x; //平移
    outXYZ.y+=RadarXYZ.y;
    outXYZ.z+=RadarXYZ.z;
    //得到地心直角坐标系

    GEO rt=XYZ2Geo(outXYZ); //转换到大地坐标系
    if(rt.L<0) rt.L+=PI; //解决负经度问题
    if(rt.B<0) rt.B+=PI;
    return rt;
}

GEO RadarConv::Map2Geo(XYZ in)
{
    XYZ outXYZ;

    double xx=in.x;
    double yy=in.y;
    in.x=yy;in.y=xx; //将地图坐标转换为内部坐标（X指向北极）

    outXYZ.x=R00*in.x+R01*in.y+R02*in.z;
    outXYZ.y=R10*in.x+R11*in.y+R12*in.z;
    outXYZ.z=R20*in.x+R21*in.y+R22*in.z; //矩阵旋转

    outXYZ.x+=RadarXYZ.x; //平移
    outXYZ.y+=RadarXYZ.y;
    outXYZ.z+=RadarXYZ.z;
    //得到地心直角坐标系

    GEO rt=XYZ2Geo(outXYZ); //转换到大地坐标系
    if(rt.L<0) rt.L+=PI; //解决负经度问题，忘记了为什么加180度，但这样才是对的，也许在西半球会有问题
    //if(rt.B<0) rt.B+=PI;
    return rt;
}

//计算参考高度，由于采用球体模型代替椭球，计算可能带来数十米的高度误差
//使用Radar2Geo函数反算，则消除了这个计算误差
double RadarConv::Altutide(XYZ in)
{
    double corR=siteR;
    double range = sqrt(in.x*in.x+in.y*in.y);
    double theta = atan(range/(siteR+H));
    double nth = (siteR+H+in.z)/cos(theta)-corR;
    return nth;
}

static double f=1/298.257222101; //扁率
static double a=6378137.0000; //地球长半轴半径
static double b=a*(1-f); //地球短半轴半径
static double k0=1.00000254;   //中心经线缩放因子
static double ee=f*(2-f);    //e^2
static double n=f/(2-f);
static double ap=a*(1+n*n/4+n*n*n*n/64)/(1+n);
static double beta1=n/2-2*n*n/3+5*n*n*n/16+41*n*n*n*n/180;
static double beta2=13*n*n/48-3*n*n*n/5+557*n*n*n*n/1440;
static double beta3=61*n*n*n/240-103*n*n*n*n/140;
static double beta4=49561*n*n*n*n/161280;
static double agA=ee;
static double agB=(5*ee*ee-ee*ee*ee)/6;
static double agC=(104*ee*ee*ee-45*ee*ee*ee*ee)/120;
static double agD=1237*ee*ee*ee*ee/1260;

static double n2=n*n;
static double n4=n2*n2;
static double dt1=n/2-2*n2/3+37*n2*n/96-n4/360;
static double dt2=n2/48+n2*n/15-437*n4/1440;
static double dt3=17*n2*n/480-37*n4/840;
static double dt4=4397*n4/161280;
static double Aa=ee+ee*ee+ee*ee*ee+ee*ee*ee*ee;
static double Ab=-(7*ee*ee+17*ee*ee*ee+30*ee*ee*ee*ee)/17;
static double Ac=(224*ee*ee*ee+889*ee*ee*ee*ee)/120;
static double Ad=-4279*ee*ee*ee*ee/1260;


double GEO::length(GEO &bhe)
{
    GEO* pGeo=&bhe;  //计算到指定点大圆线的距离,使用大圆计算，会有一些误差
    double Bb=(pGeo->B+B)/2;
    double N=a/(1-ee*sin(Bb)*sin(Bb));
    double N1=a/(1-ee*sin(B)*sin(B));
    double N2=a/(1-ee*sin(pGeo->B)*sin(pGeo->B));
    double R=a*sqrt(1-ee)/(1-ee*sin(Bb)*sin(Bb));

    double x1=N1*cos(B)*cos(L);
    double x2=N2*cos(pGeo->L)*cos(pGeo->B);
    double y1=N1*cos(B)*sin(L);
    double y2=N2*sin(pGeo->L)*cos(pGeo->B);
    double z1=N1*sin(B)*(1-ee);
    double z2=N2*sin(pGeo->B)*(1-ee);
    double v=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
    double ang=acos((N1*N1+N2*N2-v)/(2*N1*N2));
    return N*ang;//*0.999;
    /*
    //计算到指定点大圆线的距离,使用大椭圆计算，纬度相近的时候会有较大误差
    double detL=bhe.L-L;
    double detB=bhe.B-B;
    double t1=atan((1-ee)*tan(B));
    double t2=atan((1-ee)*tan(bhe.B));


    double A1=atan(sin(detL)/(cos(t1)*(tan(t1)-tan(t2)*cos(detL))));
    double A2=atan(sin(detL)/(cos(t2)*(tan(t2)-tan(t1)*cos(detL))));

    double th1=atan(cos(A1)/tan(t1));
    double th2=atan(cos(A2)/tan(t2));

    double tpe=cos(t1)*cos(t1)*sin(A1)*sin(A1);
    double bb=b/(sqrt(1-ee*tpe));
    double bee=(1-tpe)*ee/(1-ee*tpe);
    double leng=(1+bee/4+13*bee*bee/64+25*bee*bee*bee/256)*(th2-th1);
    leng-=(bee/8+3*bee*bee/32+95*bee*bee*bee/1024)*(sin(2*th2)-sin(2*th1));
    leng=bb*leng;
    if(leng<0) leng=-leng;
    return leng;
    */
}

double GEO::distance(GEO *pGeo) //计算到指定点的空间直线距离
{

    double N1=a/(1-ee*sin(B)*sin(B))+H;
    double N2=a/(1-ee*sin(pGeo->B)*sin(pGeo->B))+pGeo->H;
    double x1=N1*cos(B)*cos(L);
    double x2=N2*cos(pGeo->L)*cos(pGeo->B);
    double y1=N1*cos(B)*sin(L);
    double y2=N2*sin(pGeo->L)*cos(pGeo->B);
    double z1=(N1-H)*sin(B)*(1-ee);
    double z2=(N2-pGeo->H)*sin(pGeo->B)*(1-ee);
    double v=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
    return v;

};

double GEO::arc(GEO *pGeo) //计算到指定点地心连线夹角
{

    double N1=a/(1-ee*sin(B)*sin(B))+H;
    double N2=a/(1-ee*sin(pGeo->B)*sin(pGeo->B))+pGeo->H;

    double x1=N1*cos(B)*cos(L);
    double x2=N2*cos(pGeo->L)*cos(pGeo->B);
    double y1=N1*cos(B)*sin(L);
    double y2=N2*sin(pGeo->L)*cos(pGeo->B);
    double z1=(N1-H)*sin(B)*(1-ee);
    double z2=(N2-pGeo->H)*sin(pGeo->B)*(1-ee);
    double v=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
    double ang=acos((N1*N1+N2*N2-v)/(2*N1*N2));
    return ang;
}

double GEO::length(GEO *pGeo) //计算到指定点大圆线的距离,使用大圆计算，会有一些误差
{
    double Bb=(pGeo->B+B)/2;
    double N=a/(1-ee*sin(Bb)*sin(Bb));
    double N1=a/(1-ee*sin(B)*sin(B));
    double N2=a/(1-ee*sin(pGeo->B)*sin(pGeo->B));
    double R=a*sqrt(1-ee)/(1-ee*sin(Bb)*sin(Bb));

    double x1=N1*cos(B)*cos(L);
    double x2=N2*cos(pGeo->L)*cos(pGeo->B);
    double y1=N1*cos(B)*sin(L);
    double y2=N2*sin(pGeo->L)*cos(pGeo->B);
    double z1=N1*sin(B)*(1-ee);
    double z2=N2*sin(pGeo->B)*(1-ee);
    double v=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
    double ang=acos((N1*N1+N2*N2-v)/(2*N1*N2));
    return N*ang;//*0.999;

    /* //椭球大圆线距离，还有些问题，不使用
    static double eps=4.0E-4;
    double detL=pGeo->L-L;
    double detB=pGeo->B-B;

    double tmp1=(1-ee)*tan(B);
    double tmp2=(1-ee)*tan(pGeo->B);
    double t1=atan(tmp1);
    double t2=atan(tmp2);

    if(fabs(detB)<eps)
    {
    double N=a/(1-ee*sin(B)*sin(B));
    double x1=cos(B)*cos(L);
    double x2=cos(pGeo->L)*cos(pGeo->B);
    double y1=cos(B)*sin(L);
    double y2=sin(pGeo->L)*cos(pGeo->B);
    double z1=sin(B)*(1-ee);
    double z2=sin(pGeo->B)*(1-ee);
    double v=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
    double vr=v*N;
    return vr;
    }

    double ta1=cos(t1)*(tan(t2)-tan(t1)*cos(detL));
    double ta2=cos(t2)*(tan(t1)-tan(t2)*cos(detL));

    double A1,A2;
    A1=(ta1!=0.0)?atan(sin(detL)/ta1):(PI/2);
    A2=(ta2!=0.0)?atan(sin(detL)/ta2):(PI/2);

    double th1,th2;
    th1=(tmp1!=0.0)?atan(cos(A1)/tmp1):(PI/2);
    th2=(tmp2!=0.0)?atan(cos(A2)/tmp2):(PI/2);

    double tpe=cos(t1)*cos(t1)*sin(A1)*sin(A1);
    double bb=b/(sqrt(1-ee*tpe));
    double bee=(1-tpe)*ee/(1-ee*tpe);
    double leng=(1+bee/4+13*bee*bee/64+25*bee*bee*bee/256)*(th2-th1);
    leng-=(bee/8+3*bee*bee/32+95*bee*bee*bee/1024)*(sin(2*th2)-sin(2*th1));
    leng=bb*leng;
    if(leng<0) leng=-leng;

    return leng;
    */
}

double GEO::stereoX(GEO *pGeo)  //stereographic projection;
{
    if(pGeo==NULL) return 0.0;
    double Bc=B-(ee/2+5*ee*ee/24+3*ee*ee*ee/32)*sin(2.0*B)+(5*ee*ee/48+7*ee*ee*ee/80)*sin(4.0*B);
    double R=a*cos(B)/((1-ee*sin(B)*sin(B))*cos(Bc));
    double k=2*R/(1+sin(B)*sin(pGeo->B)+cos(B)*cos(pGeo->B)*cos(pGeo->L-L));
    double X=k*cos(pGeo->B)*sin(pGeo->L-L);
    return X;
}
double GEO::stereoY(GEO *pGeo)  //stereographic projection;
{
    if(pGeo==NULL) return 0.0;
    double Bc=B-(ee/2+5*ee*ee/24+3*ee*ee*ee/32)*sin(2.0*B)+(5*ee*ee/48+7*ee*ee*ee/80)*sin(4.0*B);
    double R=a*cos(B)/((1-ee*sin(B)*sin(B))*cos(Bc));
    double k=2*R/(1+sin(B)*sin(pGeo->B)+cos(B)*cos(pGeo->B)*cos(pGeo->L-L));
    double Y=k*(cos(B)*sin(pGeo->B)-sin(B)*cos(pGeo->B)*cos(pGeo->L-L));
    return Y;
}

void GEO::stereoXY(GEO *pGeo,double &x,double &y)  //stereographic projection;
{
    if(pGeo==NULL) {x=0.0;y=0.0;return;}
    double Bc=B-(ee/2+5*ee*ee/24+3*ee*ee*ee/32)*sin(2.0*B)+(5*ee*ee/48+7*ee*ee*ee/80)*sin(4.0*B);
    double R=a*cos(B)/((1-ee*sin(B)*sin(B))*cos(Bc));
    double k=2*R/(1+sin(B)*sin(pGeo->B)+cos(B)*cos(pGeo->B)*cos(pGeo->L-L));
    x=k*cos(pGeo->B)*sin(pGeo->L-L);
    y=k*(cos(B)*sin(pGeo->B)-sin(B)*cos(pGeo->B)*cos(pGeo->L-L));
}


void GEO::unstereo(double x,double y,GEO *pGeo)
{
    if(pGeo==NULL) return;
    double Bc=B-(ee/2+5*ee*ee/24+3*ee*ee*ee/32)*sin(2.0*B)+(5*ee*ee/48+7*ee*ee*ee/80)*sin(4.0*B);
    double R=a*cos(B)/((1-ee*sin(B)*sin(B))*cos(Bc));

    double p=sqrt(x*x+y*y);
    double c=2*atan(p/(2*R));
    pGeo->B=asin(cos(c)*sin(B)+(y*sin(c)*cos(B)/p));
    pGeo->L=L+atan(x*sin(c)/(p*cos(B)*cos(c)-y*sin(B)*sin(c)));
}

double GEO::LambertX(GEO *pGeo)
{
    if(pGeo==NULL) return 0.0;
    double k=sqrt(2/(1+sin(B)*sin(pGeo->B)+cos(B)*cos(pGeo->B)*cos(pGeo->L-L)));
    double X=k*cos(pGeo->B)*sin(pGeo->L-L);
    return X;
}
double GEO::LambertY(GEO *pGeo)
{
    if(pGeo==NULL) return 0.0;
    double k=sqrt(2/(1+sin(B)*sin(pGeo->B)+cos(B)*cos(pGeo->B)*cos(pGeo->L-L)));
    double Y=k*(cos(B)*sin(pGeo->B)-sin(B)*cos(pGeo->B)*cos(pGeo->L-L));
    return Y;
}
void GEO::LambertXY(GEO *pGeo,double &x,double &y)
{
    if(pGeo==NULL) return;
    double k=sqrt(2/(1+sin(B)*sin(pGeo->B)+cos(B)*cos(pGeo->B)*cos(pGeo->L-L)));
    y=k*(cos(B)*sin(pGeo->B)-sin(B)*cos(pGeo->B)*cos(pGeo->L-L));
    x=k*cos(pGeo->B)*sin(pGeo->L-L);
}

void GEO::unLambert(double x,double y,GEO *pGeo)
{
    if(pGeo==NULL) return;

    double p=sqrt(x*x+y*y);
    double c=2*asin(p/2);
    pGeo->B=asin(cos(c)*sin(B)+(y*sin(c)*cos(B)/p));
    pGeo->L=L+atan(x*sin(c)/(p*cos(B)*cos(c)-y*sin(B)*sin(c)));
}


void GEO::projectXY(GEO *pGeo,double &x,double &y)
{
    double detL=pGeo->L- L;
    double cosdetL=cos(detL);
    double sinB=sin(pGeo->B);
    double cosB=cos(pGeo->B);
    double sinB2=sinB*sinB;
    double sinB4=sinB2*sinB2;
    double B1=pGeo->B-sinB*cosB*(agA+agB*sinB2+agC*sinB4+agD*sinB4*sinB2);
    double eps;

    double ysinB2=sin(B)*sin(B);
    double ysinB4=sinB2*sinB2;
    eps=B-sin(B)*cos(B)*(agA+agB*ysinB2+agC*ysinB4+agD*ysinB4*ysinB2);
    double y0=k0*ap*(eps+beta1*sin(2*eps)+beta2*sin(4*eps)+beta3*sin(6*eps)+beta4*sin(8*eps));

    if(cosdetL==0.0) eps=PI/2.0;
    else eps=atan(tan(B1)/cosdetL);
    double nuu=atanh(cos(B1)*sin(detL));

    x=k0*ap*(nuu+beta1*cos(2*eps)*sinh(2*nuu)+
             beta2*cos(4*eps)*sinh(4*nuu)+
             beta3*cos(6*eps)*sinh(6*nuu)+
             beta4*cos(8*eps)*sinh(8*nuu));

    y=k0*ap*(eps+beta1*sin(2*eps)*cosh(2*nuu)+
             beta2*sin(4*eps)*cosh(4*nuu)+
             beta3*sin(6*eps)*cosh(6*nuu)+
             beta4*sin(8*eps)*cosh(8*nuu))-y0;

}

double GEO::projectX(GEO *pGeo)
{
    double detL=pGeo->L- L;
    double cosdetL=cos(detL);
    double sinB=sin(pGeo->B);
    double cosB=cos(pGeo->B);
    double sinB2=sinB*sinB;
    double sinB4=sinB2*sinB2;
    double B1=pGeo->B-sinB*cosB*(agA+agB*sinB2+agC*sinB4+agD*sinB4*sinB2);
    double eps;
    if(cosdetL==0.0) eps=PI/2.0;
    else eps=atan(tan(B1)/cosdetL);
    double nuu=atanh(cos(B1)*sin(detL));

    double x=k0*ap*(nuu+beta1*cos(2*eps)*sinh(2*nuu)+
                    beta2*cos(4*eps)*sinh(4*nuu)+
                    beta3*cos(6*eps)*sinh(6*nuu)+
                    beta4*cos(8*eps)*sinh(8*nuu));

    return x;
}

double GEO::projectY(GEO *pGeo)
{
    double sinB=sin(B);
    double cosB=cos(B);
    double sinB2=sinB*sinB;
    double sinB4=sinB2*sinB2;
    double B1=B-sinB*cosB*(agA+agB*sinB2+agC*sinB4+agD*sinB4*sinB2);
    double eps=B1;
    double y0=k0*ap*(eps+beta1*sin(2*eps)+beta2*sin(4*eps)+beta3*sin(6*eps)+beta4*sin(8*eps));

    double detL=pGeo->L- L;
    double cosdetL=cos(detL);
    sinB=sin(pGeo->B);
    cosB=cos(pGeo->B);
    sinB2=sinB*sinB;
    sinB4=sinB2*sinB2;
    B1=pGeo->B-sinB*cosB*(agA+agB*sinB2+agC*sinB4+agD*sinB4*sinB2);
    if(cosdetL==0.0) eps=PI/2.0;
    else eps=atan(tan(B1)/cosdetL);
    double nuu=atanh(cos(B1)*sin(detL));
    double y=k0*ap*(eps+beta1*sin(2*eps)*cosh(2*nuu)+
                    beta2*sin(4*eps)*cosh(4*nuu)+
                    beta3*sin(6*eps)*cosh(6*nuu)+
                    beta4*sin(8*eps)*cosh(8*nuu));
    return y-y0;
}

void GEO::unproject(double x,double y,GEO *pGeo)
{
    double sinB=sin(B);
    double cosB=cos(B);
    double sinB2=sinB*sinB;
    double sinB4=sinB2*sinB2;
    double B1=B-sinB*cosB*(agA+agB*sinB2+agC*sinB4+agD*sinB4*sinB2);
    double eps=B1;
    double y0=k0*ap*(eps+beta1*sin(2*eps)+beta2*sin(4*eps)+beta3*sin(6*eps)+beta4*sin(8*eps));

    double ta=(y+y0)/(k0*ap);
    double tb=x/(k0*ap); //k0*ap是常数，不会为0

    double taa=ta-dt1*sin(2*ta)*cosh(2*tb)
               -dt2*sin(4*ta)*cosh(4*tb)
               -dt3*sin(6*ta)*cosh(6*tb)
               -dt4*sin(8*ta)*cosh(8*tb);
    double tbb=tb-dt1*cos(2*ta)*sinh(2*tb)
               -dt2*cos(4*ta)*sinh(4*tb)
               -dt3*cos(6*ta)*sinh(6*tb)
               -dt4*cos(8*ta)*sinh(8*tb);
    double coshtbb=cosh(tbb);
    double bb=asin(sin(taa)/cosh(tbb)); //cosh函数恒不为零，无需处理被0除问题
    double costaa=cos(taa);
    double ll;
    if(costaa==0.0) ll=PI/2.0;
    else ll=atan(sinh(tbb)/costaa);
    double sinb=sin(bb);

    pGeo->L=L+ll;
    pGeo->B=bb+sinb*cos(bb)*(Aa+Ab*sinb*sinb+Ac*sinb*sinb*sinb*sinb+Ad*sinb*sinb*sinb*sinb*sinb*sinb);

    //cout<<"c++ B="<<pGeo->B/PI*180.0;
    //cout<<",c++ L="<<pGeo->L/PI*180.0;
}

double GEO::avgR() //地球在特定纬度的曲率半径
{
    return (a/sqrt(1-ee*sin(B)));
}

