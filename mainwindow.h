#pragma execution_character_set("utf-8")
#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include <QWebEngineView>
#include <QWebEngineProfile>
#include <QWebEngineSettings>
#include <QWebEnginePage>
#include <QWebEngineHttpRequest>

#include "rtklib.h"
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <unordered_set>
#include <unordered_map>
#include <QMessageBox>

#include <QMainWindow>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_openN_clicked();

    void on_openO_clicked();

    void on_startBtn_clicked();


    void on_saveBtn_clicked();

    void on_downloadBtn_clicked();

private:
    Ui::MainWindow *ui;
    QWebEngineView* view;
    obs_t obs = { 0 }; // 观测数据结构
    nav_t nav = { 0 }; // 导航数据结构
    double APPROX_POSITION_XYZ[3];//接收机坐标
    //定义一个结构体，名字叫做useableSat
    struct useableSat
    {
        //时刻
        gtime_t time;
        //伪距
        double C1;
        //卫星钟差
        double tgd;//卫星钟差
        //卫星坐标
        double X;
        double Y;
        double Z;
        //卫星名称：如G01
        char sat[4];


    };
    //定义一个结构体，名字叫做resultReciver
    struct resultReciver{
        //时刻
        gtime_t time;
        //接收机坐标
        double X;
        double Y;
        double Z;
        double B;
        double L;
        double H;
        int ns;//可用卫星数
        double PDOP;
        double GDOP;
        double TDOP;
        double sdx;//接收机坐标标准差
        double sdy;
        double sdz;

    };
    std::vector<resultReciver> resultReciverList;
    //Saastamoinen 模型
    //azel[0], azel[1]分别是方位角和仰角，pos[0], pos[1]分别是接收机纬度和经度，它们的单位都是弧度，pos[2]为接收机高度，单位为米。
    double T_Saastamoinen(const double *pos, const double *azel)
    {
        const double humi = 0.5;
        const double temp0=15.0; /* temparature at sea level */
        double hgt,pres,temp,e,z,trph,trpw;

        if (pos[2]<-100.0||1E4<pos[2]||azel[1]<=0) return 0.0;

        /* standard atmosphere */
        hgt=pos[2]<0.0?0.0:pos[2];

        pres=1013.25*pow(1.0-2.2557E-5*hgt,5.2568);
        temp=temp0-6.5E-3*hgt+273.16;
        e=6.108*humi*exp((17.15*temp-4684.0)/(temp-38.45));

        /* saastamoninen model */
        z=PI/2.0-azel[1];
        trph=0.0022768*pres/(1.0-0.00266*cos(2.0*pos[0])-0.00028*hgt/1E3)/cos(z);
        trpw=0.002277*(1255.0/temp+0.05)*e/cos(z);
        return trph+trpw;
    }
    //klobuchar模型
    //t为GPS时刻，ion为电离层参数，pos为接收机坐标，azel为卫星方位角和仰角，单位都是弧度。
    double I_klobuchar(gtime_t t, const double *ion, const double *pos,
                       const double *azel)
    {
        //        const double ion_default[]={ /* 2004/1/1 */
        //            0.1118E-07,-0.7451E-08,-0.5961E-07, 0.1192E-06,
        //            0.1167E+06,-0.2294E+06,-0.1311E+06, 0.1049E+07
        //        };
        double tt,f,psi,phi,lam,amp,per,x;
        int week;

        if (pos[2]<-1E3||azel[1]<=0) return 0.0;
        if (norm(ion,8)<=0.0) return 0.0;//ion=ion_default;  //若没有电离层参数，用默认参数


        /* earth centered angle (semi-circle) */    //地球中心角
        psi=0.0137/(azel[1]/PI+0.11)-0.022;         //计算地心角(E.5.6)

        /* subionospheric latitude/longitude (semi-circle) */
        phi=pos[0]/PI+psi*cos(azel[0]);                 //计算穿刺点地理纬度(E.5.7)
        if      (phi> 0.416) phi= 0.416;        //phi不超出(-0.416,0.416)范围
        else if (phi<-0.416) phi=-0.416;
        lam=pos[1]/PI+psi*sin(azel[0])/cos(phi*PI);     //计算穿刺点地理经度(E.5.8)

        /* geomagnetic latitude (semi-circle) */
        phi+=0.064*cos((lam-1.617)*PI);                 //计算穿刺点地磁纬度(E.5.9)

        /* local time (s) */
        tt=43200.0*lam+time2gpst(t,&week);              //计算穿刺点地方时(E.5.10)
        tt-=floor(tt/86400.0)*86400.0; /* 0<=tt<86400 */

        /* slant factor */
        f=1.0+16.0*pow(0.53-azel[1]/PI,3.0);            //计算投影系数(E.5.11)

        /* ionospheric delay */
        amp=ion[0]+phi*(ion[1]+phi*(ion[2]+phi*ion[3]));
        per=ion[4]+phi*(ion[5]+phi*(ion[6]+phi*ion[7]));
        amp=amp<    0.0?    0.0:amp;
        per=per<72000.0?72000.0:per;
        x=2.0*PI*(tt-50400.0)/per;                      //(E.5.12)

        return CLIGHT*f*(fabs(x)<1.57?5E-9+amp*(1.0+x*x*(-0.5+x*x/24.0)):5E-9);     //(E.5.13)
    }
};
#endif // MAINWINDOW_H
