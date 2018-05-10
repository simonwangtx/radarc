//
// Created by WangSiyuan on 5/6/18.
//

#ifndef ZYS_COORD_H
#define ZYS_COORD_H

#include <stdio.h>
#include <cmath>
#include <string>
#include <list>
#include <vector>
#include <map>
#include <deque>
#include <time.h>
#include <float.h>


using namespace std;
#define PI 3.14159265358979323846
struct GEO //经纬度以弧度表示，高度以米表示
{
    double B;  //纬度
    double L;  //经度
    double H;  //高程

    GEO(const GEO &v)
    {
        B = v.B;
        L = v.L;
        H = v.H;
    }

    GEO(const char* pS,double iH)
    {
        init(pS,iH);
    }

    void init(const char* pS,double iH=0)
    {
        int wpx=0;
        int wpy=0;
        int flag_NS=-1;
        int flag_EW=-1;
        H=iH;

        if(*pS=='N' || *pS=='n') {flag_NS=1;++pS;}
        if(*pS=='S' || *pS=='s') {flag_NS=0;++pS;}

        while(*pS!=0 && *pS>='0' && *pS<='9')
        {
            wpx=wpx*10+(*pS-48);
            ++pS;
        }
        int bf=0;
        float bl=1.0;
        if(*pS=='.')
        {
            ++pS;
            while(*pS>='0' && *pS<='9')
            {
                bf=bf*10+(*pS-48);
                bl=bl*10.0;
                ++pS;
            }
        }
        if(flag_NS==-1)
        {
            if(*pS=='N' || *pS=='n') {flag_NS=1;++pS;}
            if(*pS=='S' || *pS=='s') {flag_NS=0;++pS;}
        }
        else
        {
            if(*pS=='E' || *pS=='e') {flag_EW=1;++pS;}
            if(*pS=='W' || *pS=='w') {flag_EW=0;++pS;}
        }
        if(flag_NS==1) set_DMS_B(wpx);
        if(flag_NS==0) {set_DMS_B(wpx);B=-B;}

        if(bf>0) B=B+4.8481368E-6*bf/bl;
        while(*pS!=0 && *pS>='0' && *pS<='9')
        {
            wpy=wpy*10+(*pS-48);
            ++pS;
        }
        int lf=0;
        float ll=1.0;
        if(*pS=='.')
        {
            ++pS;
            while(*pS>='0' && *pS<='9')
            {
                lf=lf*10+(*pS-48);
                ll=ll*10.0;
                ++pS;
            }
        }
        if(flag_EW==-1)
        {
            if(*pS=='E' || *pS=='e') {set_DMS_L(wpy);++pS;}
            if(*pS=='W' || *pS=='w') {set_DMS_L(wpy);L=-L;++pS;}
        }
        else
        {
            if(flag_EW==1) set_DMS_L(wpy);
            if(flag_EW==0) {set_DMS_L(wpy);L=-L;}
        }
        if(lf>0) L=L+4.8481368E-6*lf/ll;
    }

    void init(const char* pS,const char* pE,double iH=0)
    {
        int wpx=0;
        int wpy=0;
        int flag_NS=-1;
        int flag_EW=-1;
        H=iH;

        if(*pS=='N' || *pS=='n') {flag_NS=1;++pS;}
        if(*pS=='S' || *pS=='s') {flag_NS=0;++pS;}

        while(pS!=pE && *pS>='0' && *pS<='9')
        {
            wpx=wpx*10+(*pS-48);
            ++pS;
        }
        int bf=0;
        float bl=1.0;
        if(*pS=='.')
        {
            ++pS;
            while(*pS>='0' && *pS<='9')
            {
                bf=bf*10+(*pS-48);
                bl=bl*10.0;
                ++pS;
            }
        }
        if(flag_NS==-1)
        {
            if(*pS=='N' || *pS=='n') {flag_NS=1;++pS;}
            if(*pS=='S' || *pS=='s') {flag_NS=0;++pS;}
        }
        else
        {
            if(*pS=='E' || *pS=='e') {flag_EW=1;++pS;}
            if(*pS=='W' || *pS=='w') {flag_EW=0;++pS;}
        }
        if(flag_NS==1) set_DMS_B(wpx);
        if(flag_NS==0) {set_DMS_B(wpx);B=-B;}

        if(bf>0) B=B+4.8481368E-6*bf/bl;
        while(*pS!=0 && *pS>='0' && *pS<='9')
        {
            wpy=wpy*10+(*pS-48);
            ++pS;
        }
        int lf=0;
        float ll=1.0;
        if(*pS=='.')
        {
            ++pS;
            while(*pS>='0' && *pS<='9')
            {
                lf=lf*10+(*pS-48);
                ll=ll*10.0;
                ++pS;
            }
        }
        if(flag_EW==-1)
        {
            if(*pS=='E' || *pS=='e') {set_DMS_L(wpy);++pS;}
            if(*pS=='W' || *pS=='w') {set_DMS_L(wpy);L=-L;++pS;}
        }
        else
        {
            if(flag_EW==1) set_DMS_L(wpy);
            if(flag_EW==0) {set_DMS_L(wpy);L=-L;}
        }
        if(lf>0) L=L+4.8481368E-6*lf/ll;
    }

    void set_value(float fB,float fL)
    {
        H=0.0;
        B=(double)fB;
        L=(double)fL;
    };

    GEO(int iB,int iL)
    {
        H=0.0;
        set_DMS_L(iL);
        set_DMS_B(iB);
    };

    GEO(double dB,double dL,double dH)
    {
        H=dH;
        B=dB*PI/180.0;
        L=dL*PI/180.0;
    };

    GEO(int iB,int iL,double iH)
    {
        H=iH;
        set_DMS_L(iL);
        set_DMS_B(iB);
    };

    GEO(){L=0;B=0;H=0;};  //构造函数，清零
    ~GEO(){};

    double offset(GEO* pos) //在小距离的情况下，直接按照经纬度计算的距离
    {
        return 6378137.0*sqrt((B-pos->B)*(B-pos->B)+(L-pos->L)*(L-pos->L)*cos(B)*cos(B));
    }
    void offsetX(double ofx){L+=ofx/(6378137.0000*cos(B));} //在小距离的情况下，直接按照距离偏移经纬度
    void offsetY(double ofy){B+=ofy/(6378137.0000*277.0/278.0);}  //在小距离的情况下，直接按照距离偏移经纬度
    void offsetXY(double ofx,double ofy){L+=ofx/(6378137.0000*cos(B));B+=ofy/(6378137.0000*277.0/278.0);} //在小距离的情况下，直接按照距离偏移经纬度
    //注意，直接偏移存在误差，建议直接偏移值控制在100KM范围内，否则会引起较大的误差
    double L_deg(){return L*180/PI;};  //返回经度的角度表示
    double B_deg(){return B*180/PI;};  //返回纬度的角度表示
    void set_DEG_L(double l){L=l*PI/180;}; //按度设置经度参数
    void set_DEG_B(double b){B=b*PI/180;}; //按度设置纬度度参数
    void set_DEG_L(float l){L=l*PI/180;}; //按度设置经度参数
    void set_DEG_B(float b){B=b*PI/180;}; //按度设置纬度度参数
    void set_DEG(double b,double l){L=l*PI/180;B=b*PI/180;}
    void set_DEG2(double l,double b){L=l*PI/180;B=b*PI/180;}
    void setH(float h)
    {
        H=h;
    };

    void set_DEG_L(double Ld,double  Lm,double  Ls) //按度分秒设置经度参数
    {
        L=Ld+Lm/60+Ls/3600;
        L=L*PI/180;
    };

    void set_DEG_B(double  Bd,double  Bm,double  Bs) //按度分秒设置纬度度参数
    {
        B=Bd+Bm/60+Bs/3600;
        B=B*PI/180;
    };

    void set_DMS_L(int l)
    {
        int s=l % 100;
        int m=l % 10000;
        m=(m-s)/100;
        int d=l/10000;
        L=d+(m/60.00)+s/3600.00;
        L=L*PI/180;
    };

    void set_DMS_B(int b)
    {
        int s=b % 100;
        int m=b % 10000;
        m=(m-s)/100;
        int d=b/10000;
        B=d+m/60.00+s/3600.00;
        B=B*PI/180;
    };

    //-----------------------------------------------------
    //设置坐标无效
    void set_null()
    {
        L=0;
        B=0;
        H=0;
    };
    //判明该坐标是否无效
    bool is_null()
    {
        if(L==0 && B==0 && H==0)
            return true;
        return false;
    };
    //计算本点到指定点的角度，此角度并非球面角度，而是直接投影后的平面角度
    //这里使用sterographics投影计算航向，误差状况未分析
    //航向返回值为弧度,为与正北，顺时针的方向与经线的夹角
    double heading(GEO &v)
    {
        double x=cos(v.B)*sin(v.L-L);
        double y=cos(B)*sin(v.B)-sin(B)*cos(v.B)*cos(v.L-L);
        double arg;
        if(x==0.0 && y==0.0) return 0.0;
        if(y==0.0) arg=PI/2;
        else arg=fabs(atan(x/y));

        if(x<0 && y>0) arg=2*PI-arg;
        if(x<0 && y<0) arg=PI+arg;
        if(x>0 && y<0) arg=PI-arg;
        return arg;
    };

    double heading(GEO *pGeo)
    {
        double x=cos(pGeo->B)*sin(pGeo->L-L);
        double y=cos(B)*sin(pGeo->B)-sin(B)*cos(pGeo->B)*cos(pGeo->L-L);
        double arg;
        if(x==0.0 && y==0.0) return 0.0;
        if(y==0.0) arg=PI/2;
        else arg=fabs(atan(x/y));

        if(x<0 && y>0) arg=2*PI-arg;
        if(x<0 && y<0) arg=PI+arg;
        if(x>0 && y<0) arg=PI-arg;
        return arg;
    };

    void copyvalue(GEO &obj)
    {
        memcpy(this,&obj,sizeof(GEO));
    };

    void copy(GEO *pobj)
    {
        memcpy(this,pobj,sizeof(GEO));
    };

    //重载的运算符，判明是否相等
    int operator ==(GEO &obj)
    {
        if(obj.B==B && obj.L==L && obj.H==H) return true;
        return false;
    };

    string L_DMS()
    {
        char t[16];
        int d,m,s;
        double tmp=L*180/PI;
        d=int(tmp);
        tmp=tmp-d;
        tmp=tmp*60;
        m=int(tmp);
        tmp=(tmp-m)*60+0.5;
        s=int(tmp);
        if(s>=60) {s=0;m+=1;if(m>=60){m=0;d+=1;}}
        sprintf(t,"%d'%d'%d",d,m,s);
        string rt(t);
        return rt;
    };

    int LN_DMS()
    {
        int t;
        int d,m,s;
        double tmp=L*180/PI;
        d=int(tmp);
        tmp=tmp-d;
        tmp=tmp*60;
        m=int(tmp);
        tmp=(tmp-m)*60+0.5;
        s=int(tmp);
        if(s>=60) {s=0;m+=1;if(m>=60){m=0;d+=1;}}
        t=d*10000+m*100+s;
        return t;
    };
    double LN_DMSS()
    {
        double rt;
        int t;
        int d,m,s;
        double tmp=L*180/PI;
        d=int(tmp);
        tmp=tmp-d;
        //------------
        tmp=tmp*60;
        m=int(tmp);
        tmp=(tmp-m)*60;
        s=int(tmp);
        //-----------------
        rt=tmp-s;
        if(s>=60) {s=0;m+=1;if(m>=60){m=0;d+=1;}}
        t=d*10000+m*100+s;
        rt+=t;
        return rt;
    };
    //---------------------------------------------------------
    //返回纬度,度，分，秒表达方式
    string B_DMS()
    {
        char t[16];
        int d,m,s;
        double tmp=B*180/PI;
        d=int(tmp);
        tmp=tmp-d;
        tmp=tmp*60;
        m=int(tmp);
        tmp=(tmp-m)*60+0.5;
        s=int(tmp);
        if(s>=60) {s=0;m+=1;if(m>=60){m=0;d+=1;}}
        sprintf(t,"%d'%d'%d",d,m,s);
        string rt(t);
        return rt;
    };

    int BN_DMS()
    {
        int t;
        int d,m,s;
        double tmp=B*180/PI;
        d=int(tmp);
        tmp=tmp-d;
        tmp=tmp*60;
        m=int(tmp);
        tmp=(tmp-m)*60+0.5;
        s=int(tmp);
        if(s>=60) {s=0;m+=1;if(m>=60){m=0;d+=1;}}
        t=d*10000+m*100+s;
        return t;
    };

    double BN_DMSS()
    {
        double rt;
        int t;
        int d,m,s;
        double tmp=B*180/PI;

        d=int(tmp);
        tmp=tmp-d;
        //------------
        tmp=tmp*60;
        m=int(tmp);
        tmp=(tmp-m)*60;
        s=int(tmp);
        //-----------------
        rt=tmp-s;
        if(s>=60) {s=0;m+=1;if(m>=60){m=0;d+=1;}}
        t=d*10000+m*100+s;
        rt+=t;
        return rt;
    };

    //两种LBH的格式的输出
    bool write1(char *buf)
    {
        sprintf(buf,"LBH(L=%d,B=%d,H=%f)",LN_DMS(),BN_DMS(),H);
        return true;
    };

    //两种LBH的格式的输出
    bool write2(char *buf)
    {
        sprintf(buf,"LBH(LB=E%dN%d,H=%f)",LN_DMS(),BN_DMS(),H);
        return true;
    };

    //两种LBH的格式的输出
    bool write3(char *buf)
    {
        sprintf(buf,"N%dE%d",BN_DMS(),LN_DMS());
        return true;
    };

    bool write4(char *buf)
    {
        sprintf(buf,"%dN%dE",BN_DMS(),LN_DMS());
        return true;
    };

    bool write5(char *buf)
    {
        sprintf(buf,"%6.2fN%7.2fE",BN_DMSS(),LN_DMSS());
        return true;
    };

    string write()
    {
        char buf[32];
        sprintf(buf,"%dN%dE",BN_DMS(),LN_DMS());
        return string(buf);
    };

    string geo_str()
    {
        char buf[32];
        sprintf(buf,"%dN%dE",BN_DMS(),LN_DMS());
        return string(buf);
    };

    bool in_range(const char* pos,double range) //给定的点是否在本点的距离范围内
    {
        GEO npt(pos,0.0);
        return (length(npt)<range);
    }
    void directXY(GEO *pGeo,double &x,double &y)  //两个点之间的直接经纬度偏差;
    {
        y=(pGeo->B-B)*6371393.0;
        x=(pGeo->L-L)*6371393.0*cos(B);
    }
    double length(GEO &bhe); //计算到指定点大圆线的距离,使用大椭圆计算，纬度相近的时候会有较大误差

    double distance(GEO *pGeo); //计算到指定点的空间直线距离

    double arc(GEO *pGeo); //计算到指定点地心连线夹角

    double length(GEO *pGeo); //计算到指定点大圆线的距离,使用大圆计算，会有一些误差

    double stereoX(GEO *pGeo);  //极射 projection;

    double stereoY(GEO *pGeo);  //极射 projection;

    void stereoXY(GEO *pGeo,double &x,double &y);  //stereographic projection;

    void unstereo(double x,double y,GEO *pGeo);

    double LambertX(GEO *pGeo);  //兰伯特 projection;

    double LambertY(GEO *pGeo);  //兰伯特 projection;

    void LambertXY(GEO *pGeo,double &x,double &y);  //兰伯特 projection;

    void unLambert(double x,double y,GEO *pGeo);

#ifdef WIN32
    //window没有atanh这个破函数，恶心!只能自己写了
	inline double atanh(const double &x)
	{
		return log((1.0 + x) / (1.0 - x)) / 2.0;
	};
#endif

    void projectXY(GEO *pGeo,double &x,double &y); //高斯投影

    double projectX(GEO *pGeo);

    double projectY(GEO *pGeo);

    void unproject(double x,double y,GEO *pGeo);

    double avgR(); //地球在特定纬度的曲率半径

    inline bool is_in_grid(GEO *pGeo,double grid_size)
    {
        return(pGeo->B>B-grid_size && pGeo->B<B+grid_size
               && pGeo->L>L-grid_size && pGeo->L<L+grid_size);
    }

};

struct XYZ //地心直角坐标或者雷达站平面坐标
{
    double x;
    double y;
    double z;
};

struct Radar
{
    double x; //x,y,z是雷达站心切平面直角坐标系，与雷达坐标相关，不应该用雷达解码器填入，ADS-B的x,y,z无意义，通常填入经纬度值
    double y;
    double z; //注意，z不是雷达报告的ModeC,而是一个需要有pull-down计算的站心切平面的垂直高度，这个值的计算与雷达的地理位置有关
    //-------------------------------------------------
    double theta;		//到正北的偏转角，弧度单位
    double range;		//到站心的距离，这不是水平投影距离，而是直接距离
    double rangeXY;		//到站心的距离，这是水平投影距离，不是直接距离
    int    level;		//这是雷达报告的ModeC值高度,注意，这里使用英尺单位,使用英尺单位是保持与雷达头一致，便于高度跟踪
};

class RadarConv:public GEO
{
private:
    double G00,G01,G02;
    double G10,G11,G12;
    double G20,G21,G22;

    double R00,R01,R02;
    double R10,R11,R12;
    double R20,R21,R22;

    double siteR; // 雷达站位置法平面地球半径
    double alpha; //偏心率
    double a;  //短半轴
    double b;  //长半轴
    double ee;//
    double ee2;
    double theR; //地球平均半径
    double RadarBias; //弧度单位,角度修正量
    double RangeDelta;  //距离修正量，米为单位
    //double Range;  //距离修正量，米为单位

    double cosB0;  //站点纬度的正余弦预计算值
    double sinB0;
    double Bc;     //保角纬度，用于极射投影计算
    double Re;     //极射投影用的地球半径
    double N;      //站点附近的地球半径
    int transfer_delay;  //雷达传输延迟，毫秒单位
public:
    XYZ RadarXYZ; //这是计算以后的雷达站地心直角坐标的值，在类初始化的时候根据地理位置自动计算

    RadarConv(){RadarBias=0.0;transfer_delay=0;RangeDelta=0.0;}

    RadarConv(GEO site,double bias=0.0):GEO(site)
    {
        init(GEO::B,GEO::L,GEO::H);
        RadarBias=bias/180.0*PI;
        RangeDelta=0.0;
        transfer_delay=0;
    }

    void set_bias(double bias) //度为单位
    {
        RadarBias=bias/180.0*PI;
    }

    void set_range_offset(double range) //米为单位
    {
        RangeDelta=range;
    }

    void set_transfer_delay(int d) //毫秒单位的延迟量
    {
        transfer_delay=d;
    }

    GEO RadarXY2Geo(Radar &in);

    void init(double B0,double L0,double H0=0.0);

    //大地坐标转地心直角坐标
    XYZ Geo2XYZ(GEO in);

    //大地坐标转站心坐标,注意，X指向正北，Y指向正东，Z指向法线
    XYZ Geo2Site(GEO in);

    //大地坐标转站心坐标,注意，X指向正东，Y指向正北，Z指向法线
    XYZ Geo2Map(GEO in);

    //地心直角坐标转大地坐标
    GEO XYZ2Geo(XYZ in); //要求输入的位置x,y不为零

    //站心直角坐标转大地坐标
    GEO Site2Geo(XYZ in);

    //站心直角坐标转大地坐标
    GEO Map2Geo(XYZ in);

    //大地坐标转雷达坐标
    Radar Geo2Radar(GEO in);


    //将雷达斜距,方位角,ModeC转换为雷达站心切平面,x,y,z值
    void RadarRange(Radar &in);

    //雷达坐标转大地坐标，转换均是基于雷达站心切平面坐标系的x,y值计算
    //对于雷达站心切平面来说，x,y值计算应该考虑pull-down影响
    GEO Radar2Geo(Radar &in);

    //计算参考高度，由于采用球体模型代替椭球，计算可能带来数十米的高度误差
    //使用Radar2Geo函数反算，则消除了这个计算误差
    double Altutide(XYZ in);

};

#endif //ZYS_COORD_H
