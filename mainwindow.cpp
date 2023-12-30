#include "mainwindow.h"
#include "ui_mainwindow.h"


#include <QFileDialog>
#include <QDebug>
#include <QMessageBox>
#include <iostream>
#include <QDir>
#include <QWebChannel>
//QDesktopServices
#include <QDesktopServices>
MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    //设置窗口标题
    this->setWindowTitle("GNSS伪距单点定位软件");
    //禁用全屏
    this->setWindowFlags(Qt::WindowCloseButtonHint);
    //禁止放大缩小
    this->setFixedSize(this->width(), this->height());
     //表格设置
    ui->resultTable->setColumnCount(14);
    QStringList header;
    header<<"时间(GPST)"<<"卫星数目"<<"X"<<"Y"<<"Z"<<"B"<<"L"<<"H"<< "sdx"<<"sdy"<<"sdz"<<"GDOP"<< "PDOP"<<"TDOP";
    ui->resultTable->setHorizontalHeaderLabels(header);
    //不可编辑
    ui->resultTable->setEditTriggers(QAbstractItemView::NoEditTriggers);

    QWidget* web = ui->webview;
    //QWebEngineView
    view = new QWebEngineView(web);
    QString currentPath = QCoreApplication::applicationDirPath();
    QString filePath = QDir::toNativeSeparators(currentPath + "/map/map.html");




    //设置窗口大小
    view->resize(web->width(),web->height());

    view->page()->setUrl(QUrl::fromLocalFile(filePath));



    ui->progressBar->setValue(0);
    //不可编辑
    ui->nPath->setReadOnly(true);
    ui->oPath->setReadOnly(true);






}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_openN_clicked()
{

    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"), "", tr("RINEX文件(*.rnx brdc*.*n brdc*.*g brdc*.*b *.*p *.*n)"));
    if(fileName.isEmpty())
    {
        return;
    }
    //读取文件的第一行
    QFile file(fileName);
    if(!file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        QMessageBox::warning(this, "提示", "读取导航数据失败", QMessageBox::Ok);
        return;
    }
    //读取第一行
    QByteArray line = file.readLine();
    //读取rnx版本:     4.00           NAVIGATION DATA     M                   RINEX VERSION / TYPE
    QString rnxVersion = QString(line.mid(5, 1));
    //判断是否为2或3版本
    if(rnxVersion != "2" && rnxVersion != "3")
    {
        QMessageBox::warning(this, "提示", "只支持v2/v3版本RINEX文件", QMessageBox::Ok);
        return;
    }
    ui->nPath->setText(fileName);
    //nav清空
    nav = {0};
    int rcv = 1; // 接收机编号
    const char* navOpt = ""; // 导航文件选项
    sta_t sta = { 0 }; // 基准站信息结构
    char* navFile;
    QByteArray ba = fileName.toLatin1();
    navFile = ba.data();
    for(int i = 0; i < fileName.length(); i++)
    {
        if(navFile[i] == '/')
        {
            navFile[i] = '\\';
        }
    }
    ui->progressBar->setValue(50);

    // 读取导航数据
    int ret = readrnx(navFile, rcv, navOpt, &obs, &nav, &sta);


    if(nav.n != 0 || nav.ng != 0){
        //qdebug：共读取多少条导航数据，正在导入
        qDebug()<<"共读取到"<<nav.n+nav.ng<<"条导航数据，正在导入";
        ui->progressBar->setValue(100);
        QMessageBox::information(this, "提示", "共读取到"+QString::number(nav.n+nav.ng)+"条导航数据", QMessageBox::Ok);
        ui->progressBar->setValue(0);
    }else{
        //状态栏显示
        ui->statusbar->showMessage("读取导航数据失败");
        QMessageBox::warning(this, "提示", "读取导航数据失败", QMessageBox::Ok);
        return;
    }
}


void MainWindow::on_openO_clicked()
{
    //判断nav是否为空
    if(nav.n == 0 && nav.ng == 0)
    {
        QMessageBox::warning(this, "提示", "请先导入导航数据", QMessageBox::Ok);
        return;
    }
    //打开文件后缀为o或者rnx的文件
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"), "", tr("RINEX文件(*.rnx *.*o)"));
    if(fileName.isEmpty())
    {
        return;
    }
    //读取文件的第一行
    QFile file(fileName);
    if(!file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        QMessageBox::warning(this, "提示", "读取观测数据失败", QMessageBox::Ok);
        return;
    }
    //读取第一行
    QByteArray line = file.readLine();
    //读取rnx版本:     4.00           NAVIGATION DATA     M                   RINEX VERSION / TYPE
    QString rnxVersion = QString(line.mid(5, 1));
    //判断是否为2或3版本
    if(rnxVersion != "2" && rnxVersion != "3")
    {
        QMessageBox::warning(this, "提示", "只支持v2/v3版本RINEX文件", QMessageBox::Ok);
        return;
    }
    ui->oPath->setText(fileName);
    //obs清空
    obs = {0};
    //接着读取下一行，直至包含“END OF HEADER”
    while(!line.contains("END OF HEADER"))
    {
        line = file.readLine();
        // -3120042.4825  4084614.6856  3764025.7815                  APPROX POSITION XYZ
        if(line.contains("APPROX POSITION XYZ"))
        {
            //获取接收机坐标
            double x = line.mid(0, 14).toDouble();
            double y = line.mid(14, 14).toDouble();
            double z = line.mid(28, 14).toDouble();
            APPROX_POSITION_XYZ[0] = x;
            APPROX_POSITION_XYZ[1] = y;
            APPROX_POSITION_XYZ[2] = z;
            //打印接收机坐标
            qDebug()<<"接收机坐标："<<APPROX_POSITION_XYZ[0]<<APPROX_POSITION_XYZ[1]<<APPROX_POSITION_XYZ[2];
                break;
        }
    }
    char* obsFile;
    QByteArray ba = fileName.toLatin1();
    obsFile = ba.data();
    for(int i = 0; i < fileName.length(); i++)
    {
        if(obsFile[i] == '/')
        {
                obsFile[i] = '\\';
        }
    }

    ui->progressBar->setValue(25);



    int rcv = 1; // 接收机编号
    const char* obsOpt = ""; // 观测文件选项
    sta_t sta = { 0 }; // 基准站信息结构
    // 读取观测数据
    int ret = readrnx(obsFile, rcv, obsOpt, &obs, &nav, &sta);
    if(obs.n != 0){
        //qdebug：共读取多少条观测数据，正在导入
        qDebug()<<"共读取到"<<obs.n<<"条观测数据，正在导入";
        ui->progressBar->setValue(50);
   }
    else{
        //状态栏显示
        ui->statusbar->showMessage("读取观测数据失败");
        QMessageBox::warning(this, "提示", "读取观测数据失败", QMessageBox::Ok);
        return;
    }

   // 提取nav中各个卫星的名称
   std::unordered_set<int> satSet;

   // 遍历nav
   for (int i = 0; i < nav.n; i++) {
        // 获取导航数据的卫星编号
        int sat = nav.eph[i].sat;
        // 添加到哈希集合中
        satSet.insert(sat);
   }

   // 遍历glonass
   for (int i = 0; i < nav.ng; i++) {
        // 获取导航数据的卫星编号
        int sat = nav.geph[i].sat;
        // 添加到哈希集合中
        satSet.insert(sat);
   }
   ui->progressBar->setValue(75);
   // 遍历观测数据，删除没有导航数据的卫星
   // 创建一个临时的观测数据数组
   obsd_t* tempData = new obsd_t[obs.n];
   int tempIndex = 0;

   // 遍历观测数据，将存在导航数据的卫星复制到临时数组中
   for (int i = 0; i < obs.n; i++) {
        // 获取观测数据的卫星编号
        int sat = obs.data[i].sat;
        // 判断是否存在于哈希集合中
        if (satSet.find(sat) != satSet.end()) {
                // 将观测数据复制到临时数组中
                tempData[tempIndex++] = obs.data[i];
        }
   }

   // 删除原始的观测数据数组
   delete[] obs.data;

   // 更新观测数据数组为临时数组
   obs.data = tempData;
   obs.n = tempIndex;
      ui->progressBar->setValue(100);
    //提示信息
            //qdebug：共读取多少条观测数据，正在导入
            qDebug()<<"共读取到"<<obs.n<<"条观测数据，正在导入";
            QMessageBox::information(this, "提示", "共读取到"+QString::number(obs.n)+"条可用观测数据", QMessageBox::Ok);
         ui->progressBar->setValue(0);

}


void MainWindow::on_startBtn_clicked()
{

    //伪距单点定位
    //1.获取观测数据
    //2.获取导航数据
    //3.计算卫星位置
    //4.计算接收机位置
    //5.输出结果

    //1.获取观测数据
    //检查观测数据是否为空
    if(obs.n == 0)
    {
        QMessageBox::warning(this, "提示", "观测数据为空", QMessageBox::Ok);
        return;
    }
    //检查导航数据是否为空
    if(nav.n == 0 && nav.ng == 0)
    {
        QMessageBox::warning(this, "提示", "导航数据为空", QMessageBox::Ok);
        return;
    }

    //处理每一条观测数据
    std::vector<useableSat> useableSatList;
    for(int i = 0; i < obs.n; i++)
    {

        ui->progressBar->setValue( (double(i) / obs.n)*25);

        // 获取观测数据的时间
        gtime_t time = obs.data[i].time;
        // 获取观测数据的卫星编号
        int sat = obs.data[i].sat;
        // 获取观测数据的卫星名称
        char satName[16] = { 0 };
        satno2id(sat, satName);
        double rs[3] = { 0 };
        double dts = 0;//卫星钟差
        double var = 0;//卫星钟差，单位s
        int svh;
        satposs(time, &obs.data[i], 1, &nav, 0, rs, &dts, &var, &svh);
        //获取卫星钟差
        double tgd = dts;//卫星钟差


        //定义结构体
        useableSat useableSat;
        useableSat.time = time;
        useableSat.C1 = obs.data[i].P[0];
        useableSat.tgd = tgd;
        useableSat.X = rs[0];
        useableSat.Y = rs[1];
        useableSat.Z = rs[2];
        strcpy(useableSat.sat, satName);

        //添加到数组
        useableSatList.push_back(useableSat);
    }
    ui->progressBar->setValue(25);
    //qdebug：共有多少个可用卫星
    qDebug()<<"共有"<<useableSatList.size()<<"个可用数据";
    //按useableSatList的时间分组成多个时间相同的数组
    //定义一个数组存储时间
    std::vector<gtime_t> timeList;
    //遍历useableSatList
    for(int i = 0; i < useableSatList.size(); i++)
    {
        ui->progressBar->setValue(25+(double(i) / useableSatList.size())*25);

        //获取时间
        gtime_t time = useableSatList[i].time;
        //判断是否已经存在
        bool isExist = false;
        for(int j = 0; j < timeList.size(); j++)
        {
                if(timediff(time, timeList[j]) == 0)
                {
                    isExist = true;
                    break;
                }
        }
        if(!isExist)
        {
                timeList.push_back(time);
        }
    }
    ui->progressBar->setValue(50);

    //按时间分组
    std::vector<std::vector<useableSat>> useableSatListList;
    //遍历timeList
    for(int i = 0; i < timeList.size(); i++)
    {
        ui->progressBar->setValue(50+(double(i) / timeList.size())*25);
        //定义数组
        std::vector<useableSat> useableSatList1;
        //遍历useableSatList
        for(int j = 0; j < useableSatList.size(); j++)
        {
                //获取时间
                gtime_t time = useableSatList[j].time;
                //判断是否相同
                if(timediff(time, timeList[i]) == 0)
                {
                    //添加到数组
                    useableSatList1.push_back(useableSatList[j]);
                }
        }
        //添加到数组
        useableSatListList.push_back(useableSatList1);

    }
    qDebug()<<"共有"<<useableSatListList.size()<<"个时间段";
    ui->progressBar->setValue(75);

    //定义resultReciverList数组
    //清空
    resultReciverList.clear();
    //处理每个时间段
    for(int i=0;i<  useableSatListList.size();
         i++){
        ui->progressBar->setValue(75+(double(i) / useableSatListList.size())*15);
        std::vector<useableSat> now_useableSatList = useableSatListList[i];
        Eigen::MatrixXd x1 = Eigen::MatrixXd::Zero(4, 1);
        x1.setOnes();
        int m = 0;
        //初始xyz
        double x = APPROX_POSITION_XYZ[0];
        double y = APPROX_POSITION_XYZ[1];
        double z = APPROX_POSITION_XYZ[2];
        Eigen::MatrixXd _B(now_useableSatList.size(), 4);
        _B.setZero();
        Eigen::MatrixXd _L(now_useableSatList.size(), 1);
        _L.setZero();

        //迭代下面的代码，直至x1(0,0)小于0.0001
        while(abs(x1(0,0)) > 10e-4 || abs(x1(1,0)) > 10e-4  || abs(x1(2,0)) > 10e-4)
        {
                m++;
                if(m>999){
                    break;
                }
                //初始B矩阵
                Eigen::MatrixXd B(now_useableSatList.size(), 4);
                B.setZero();
                //初始L矩阵
                Eigen::MatrixXd L(now_useableSatList.size(), 1);
                L.setZero();
                //初始P矩阵
                Eigen::MatrixXd P(now_useableSatList.size(), now_useableSatList.size());
                P.setZero();
                //对BLP进行赋值
                for(int k = 0; k < now_useableSatList.size(); k++)
                {
                    //计算距离
                    double r = sqrt(pow(now_useableSatList[k].X-x, 2) + pow(now_useableSatList[k].Y-y, 2) + pow(now_useableSatList[k].Z-z, 2));//卫星到接收机的距离

                    //初始l
                    double l = (now_useableSatList[k].X-x)/r;
                    //初始m
                    double m = (now_useableSatList[k].Y-y)/r;
                    //初始n
                    double n = (now_useableSatList[k].Z-z)/r;

                    //将l,m,n放入B矩阵
                    B(k, 0) = -l;
                    B(k, 1) = -m;
                    B(k, 2) = -n;
                    B(k, 3) = 1;

                    //获取卫星高度截止角
                    //ecef2pos
                    double blh[3] = { 0 };
                    ecef2pos(Eigen::Vector3d(x,y,z).data(), blh);
                    //ecef2enu
                    double enu[3] = { 0 };
                    double rr[3];
                    rr[0] = now_useableSatList[k].X-x;
                    rr[1] = now_useableSatList[k].Y-y;
                    rr[2] = now_useableSatList[k].Z-z;
                    ecef2enu(blh, rr, enu);
                    //获取卫星高度截止角
                    double el = atan2(enu[2], sqrt(pow(enu[0], 2) + pow(enu[1], 2)));
                    double az = atan2(enu[0], enu[1]);
                    double azel[2];
                    azel[0] = az;
                    azel[1] = el;
                    //如果高度截止角小于15度，跳过
                    if(el < 15*PI/180)
                    {
                        continue;
                    }
                    //初始化P
                    double p = sin(el)*sin(el);
                    P(k, k) = p;
                    //初始T
                    double T = T_Saastamoinen(blh,azel);
                    //初始I
                    //t为GPS时刻，ion为电离层参数，pos为接收机坐标，azel为卫星方位角和仰角，单位都是弧度。
                    //                double I_klobuchar(gtime_t t, const double *ion, const double *pos,
                    //                                   const double *azel)
                    //判断卫星类型
                    double I =0.0;
                    if(now_useableSatList[k].sat[0] == 'G')
                    {
                        I = I_klobuchar(now_useableSatList[k].time, nav.ion_gps, blh, azel);
                    }else if(now_useableSatList[k].sat[0] == 'C')
                    {
                        I = I_klobuchar(now_useableSatList[k].time, nav.ion_cmp, blh, azel);
                    }else if(now_useableSatList[k].sat[0] == 'E')
                    {
                        I = I_klobuchar(now_useableSatList[k].time, nav.ion_gal, blh, azel);
                    }else if(now_useableSatList[k].sat[0] == 'J')
                    {
                        I = I_klobuchar(now_useableSatList[k].time, nav.ion_qzs, blh, azel);
                    }
                    //初始L
                    double L1 = now_useableSatList[k].C1 - r + now_useableSatList[k].tgd*CLIGHT-T-I;
                    //将L放入L矩阵
                    L(k, 0) = L1;
                }





                //x=(BTPB)-1BTPL
                //计算BTPB
                Eigen::MatrixXd BTB = B.transpose()*B;
                //计算BTPL
                Eigen::MatrixXd BTL = B.transpose()*L;
                //计算(BTPB)-1
                Eigen::MatrixXd BTB_1 = BTB.inverse();
                //计算x
                x1 = BTB_1*BTL;
                x += x1(0, 0);
                y += x1(1, 0);
                z += x1(2, 0);
                //打印x1
                qDebug()<<"第"<<i+1<<"个时间段的x11："<<x1(0,0);
                _B = B;
                _L = L;

        }


        //        if(m>100){
        //            continue;
        //        }
        //        //打印结果
        //        qDebug()<<"第"<<i+1<<"个时间段的结果：";
        //        qDebug()<<"x="<<x;
        //        qDebug()<<"y="<<y;
        //        qDebug()<<"z="<<z;
        //        //转BLH
        double blh[3] = { 0 };
        ecef2pos(Eigen::Vector3d(x, y, z).data(), blh);
        //        qDebug()<<"B="<<blh[0]*R2D;
        //        qDebug()<<"L="<<blh[1]*R2D;
        //        qDebug()<<"H="<<blh[2];
        //初始化Q
        Eigen::MatrixXd Q(4,4);
        Q = (_B.transpose()*_B).inverse();//(BTB)-1
        //计算GDOP
        double GDOP = sqrt(Q(0, 0)+Q(1, 1)+Q(2, 2)+Q(3, 3));
        //计算PDOP
        double PDOP = sqrt(Q(0, 0)+Q(1, 1)+Q(2, 2));
        //计算TDOP
        double TDOP = sqrt(Q(3, 3));
        //初始V
        Eigen::MatrixXd V(now_useableSatList.size(), 1);
        V.setZero();
        //计算V
        V = _B*x1+_L;
        //计算sigma0
        Eigen::MatrixXd VTV = V.transpose()*V;
        double sigma0 = sqrt(VTV(0, 0)/(now_useableSatList.size()-4));
        //计算sdx
        double sdx = sqrt(sigma0*Q(0, 0));
        //计算sdy
        double sdy = sqrt(sigma0*Q(1, 1));
        //计算sdz
        double sdz = sqrt(sigma0*Q(2, 2));





        //定义结构体
        resultReciver resultReciver1;
        resultReciver1.time = now_useableSatList[0].time;
        resultReciver1.X = x;
        resultReciver1.Y = y;
        resultReciver1.Z = z;
        resultReciver1.B = blh[0]*R2D;
        resultReciver1.L = blh[1]*R2D;
        resultReciver1.H = blh[2];
        resultReciver1.ns = now_useableSatList.size();
        resultReciver1.GDOP = GDOP;
        resultReciver1.PDOP = PDOP;
        resultReciver1.TDOP = TDOP;
        resultReciver1.sdx = sdx;
        resultReciver1.sdy = sdy;
        resultReciver1.sdz = sdz;
        //判断是否有效
        if(resultReciver1.H>1000|| resultReciver1.H<0 || resultReciver1.ns<4 || isnan(resultReciver1.H) == 1 || isnan(resultReciver1.GDOP) == 1
            || isnan(resultReciver1.PDOP) == 1
            || isnan(resultReciver1.TDOP) == 1
            || sdx == 0 || sdy == 0 || sdz == 0){
                continue;
        }


        //添加到数组
        resultReciverList.push_back(resultReciver1);
    }

   ui->progressBar->setValue(90);
     //resultTable渲染数据
    //设置行数
    ui->resultTable->setRowCount(resultReciverList.size());
    //resultReciverList
    for(int i = 0; i < resultReciverList.size();
         i++)
    {
        ui->progressBar->setValue(90+(double(i) / resultReciverList.size())*10);
        // 获取接收机坐标
        double X = resultReciverList[i].X;
        double Y = resultReciverList[i].Y;
        double Z = resultReciverList[i].Z;
        // 获取接收机BLH
        double B = resultReciverList[i].B;
        double L = resultReciverList[i].L;
        double H = resultReciverList[i].H;

        // 获取接收机时间
        gtime_t time = resultReciverList[i].time;
        //gpstime转年月日
        //    extern void time2str(gtime_t t, char *s, int n)
        char str_time[32] = { 0 };
        time2str(time, str_time, 0);

        //渲染表格
        ui->resultTable->setItem(i,0,new QTableWidgetItem(QString(str_time)));
        ui->resultTable->setItem(i,1,new QTableWidgetItem(QString::number(resultReciverList[i].ns)));
       ui->resultTable->setItem(i,2,new QTableWidgetItem(QString::number(X,'f',6)));
       ui->resultTable->setItem(i,3,new QTableWidgetItem(QString::number(Y,'f',6)));
       ui->resultTable->setItem(i,4,new QTableWidgetItem(QString::number(Z,'f',6)));
       ui->resultTable->setItem(i,5,new QTableWidgetItem(QString::number(B,'f',6)));
       ui->resultTable->setItem(i,6,new QTableWidgetItem(QString::number(L,'f',6)));
       ui->resultTable->setItem(i,7,new QTableWidgetItem(QString::number(H,'f',6)));
       ui->resultTable->setItem(i,8,new QTableWidgetItem(QString::number(resultReciverList[i].sdx,'f',6)));
       ui->resultTable->setItem(i,9,new QTableWidgetItem(QString::number(resultReciverList[i].sdy,'f',6)));
       ui->resultTable->setItem(i,10,new QTableWidgetItem(QString::number(resultReciverList[i].sdz,'f',6)));
       ui->resultTable->setItem(i,11,new QTableWidgetItem(QString::number(resultReciverList[i].GDOP,'f',6)));
       ui->resultTable->setItem(i,12,new QTableWidgetItem(QString::number(resultReciverList[i].PDOP,'f',6)));
       ui->resultTable->setItem(i,13,new QTableWidgetItem(QString::number(resultReciverList[i].TDOP,'f',6)));
    }
    //自适应列宽
    ui->resultTable->resizeColumnsToContents();
    //将bl存储到bl
    std::string bl;
    for(int i = 0; i < resultReciverList.size(); i++)
    {
        bl += std::to_string(resultReciverList[i].B);
        bl += ",";
        bl += std::to_string(resultReciverList[i].L);
        bl += ";";
    }
    QString jsCode = QString("showOnMap('%1')").arg(bl.c_str());
    view->page()->runJavaScript(jsCode, [](const QVariant &v) { qDebug() << v.toString(); });




    ui->progressBar->setValue(100);

       //提示信息
       //qdebug：共有多少个结果
       qDebug()<<"共有"<<resultReciverList.size()<<"个结果";
       QMessageBox::information(this, "提示", "共有"+QString::number(resultReciverList.size())+"个结果", QMessageBox::Ok);
   ui->progressBar->setValue(0);
}


void MainWindow::on_saveBtn_clicked()
{
   //判断resultReciverList是否为空
   if(resultReciverList.size() == 0)
   {
       QMessageBox::warning(this, "提示", "请先进行伪距单点定位", QMessageBox::Ok);
       return;
   }
   //打开文件
   QString fileName = QFileDialog::getSaveFileName(this, tr("Save File"), "", tr("文本文件(*.txt)"));
   if(fileName.isEmpty())
   {
       return;
   }
   QFile file(fileName);
   if(!file.open(QIODevice::WriteOnly | QIODevice::Text))
   {
       QMessageBox::warning(this, "提示", "保存失败", QMessageBox::Ok);
       return;
   }
   //写入文件
   QTextStream out(&file);
   //遍历resultReciverList
   for(int i = 0; i < resultReciverList.size(); i++)
   {
       // 获取接收机坐标
       double X = resultReciverList[i].X;
       double Y = resultReciverList[i].Y;
       double Z = resultReciverList[i].Z;
       // 获取接收机BLH
       double B = resultReciverList[i].B;
       double L = resultReciverList[i].L;
       double H = resultReciverList[i].H;
       // 获取接收机时间
       gtime_t time = resultReciverList[i].time;
       //gpstime转年月日
       //    extern void time2str(gtime_t t, char *s, int n)
       char str_time[32] = { 0 };
       time2str(time, str_time, 0);
       //写入文件
       out<<str_time<<" "<<resultReciverList[i].ns<<" "<<X<<" "<<Y<<" "<<Z<<" "<<B<<" "<<L<<" "<<H<<" "<<resultReciverList[i].sdx<<" "<<resultReciverList[i].sdy<<" "<<resultReciverList[i].sdz<<" "<<resultReciverList[i].GDOP<<" "<<resultReciverList[i].PDOP<<" "<<resultReciverList[i].TDOP<<"\n";
   }
   //关闭文件
   file.close();
   //提示信息
   QMessageBox::information(this, "提示", "保存成功", QMessageBox::Ok);


}


void MainWindow::on_downloadBtn_clicked()
{
   QDesktopServices::openUrl(QUrl(QLatin1String("https://cddis.nasa.gov/archive/gnss/data/daily/")));//打开网页

}

