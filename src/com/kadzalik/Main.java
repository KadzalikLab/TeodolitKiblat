package com.kadzalik;


import java.util.Locale;
import static java.lang.Math.*;


 //KadzalikLab 2021
public class Main {
    public static void main(String[] args) {


        double tgl=8,bln=3,thn=2012;
        double jam_setempatDeg=11,jam_setempatMin=30,jam_setempatSec=30;

        double derajat_lintang=-7;
        double menit_lintang=11;
        double detik_lintang=53.8;
        double lintang=Konversi.dmsKeDesimal(derajat_lintang,menit_lintang,detik_lintang);

        double derajat_bujur=113;
        double menit_bujur=28;
        double detik_bujur=21.5;
        double bujur=Konversi.dmsKeDesimal(derajat_bujur,menit_bujur,detik_bujur);

        double derajat_lintangKB=21;
        double menit_lintangKB=25;
        double detik_lintangKB=14.7;
        double lintangKB=Konversi.dmsKeDesimal(derajat_lintangKB,menit_lintangKB,detik_lintangKB);

        double derajat_bujurKB=39;
        double menit_bujurKB=49;
        double detik_bujurKB=40.39;
        double bujurKB=Konversi.dmsKeDesimal(derajat_bujurKB,menit_bujurKB,detik_bujurKB);



        double wd = Konversi.dmsKeDesimal(jam_setempatDeg,jam_setempatMin,jam_setempatSec);
        double tz=7;
        double bwd = tz*15;
        //double wd=jam_setempat;
        double jd_ut = Konversi.kmJd(tgl,bln,thn,wd-tz);
        double jd_wd = Konversi.kmJd(tgl,bln,thn,wd);

        double t = (jd_ut-2451545)/36525;
        double s = abs360(( (280.46645 + 36000.76983 * t) / 360%1)*360);
        double n =abs360( ( (125.04  -1934.136 * t) / 360%1)*360);
        double m =abs360(((357.5291+ 35999.0503 * t) / 360%1)*360);

        print("JD ",jd_ut);
        print("T",t);
        print("S",s);
        print("M",m);
        print("N",n);

        double kr1=(17.264/3600)*sin(toRadians(n))+(0.206/3600)*sin(toRadians(2*n));
        double kr2=(-1.264/3600)*sin(toRadians(2*s));
        double kr3 = (9.23/3600)*cos(toRadians(n))-(0.09/3600)*cos(toRadians(2*n));
        double kr4 = (0.548/3600)*cos(toRadians(2*s));

        printDigit("kr1",kr1,9);
        printDigit("kr2",kr2,9);
        printDigit("kr3",kr3,9);
        printDigit("kr4",kr4,9);

        double q1 = 23.43929111+kr3+kr4-(46.815/3600)*t;
        double E = (6898.06/3600)*sin(toRadians(m))+(72.095/3600)*sin(toRadians(2*m))+(0.966/3600)*sin(toRadians(3*m));
        double s1 = s+E+kr1+kr2-(20.47/3600);
        double mail_syamsi = toDegrees(asin(sin(toRadians(s1))*sin(toRadians(q1))));

        print("q1",q1);
        print("E",E);
        print("s1",s1);
        print("Mail S",mail_syamsi);

        double pt_a = toDegrees(atan(tan(toRadians(s1))*cos(toRadians(q1))));
        double pt_b; if (s1>=0&&s1<=90)pt_b=pt_a; else pt_b=0;
        double pt_c; if (s1>=90&&s1<=270)pt_c=pt_a+180; else pt_c=0;
        double pt_d; if (s1>=270&&s1<=360)pt_d=pt_a+360; else pt_d=0;
        double pt = pt_b+pt_c+pt_d;

        print("pt",pt);
        double e = (s-pt)/15;
        double sd = 0.267/(1-0.017*cos(toRadians(m)));
        print("e",e);
        print("sd",sd);
        print("wd",wd);

        System.out.println(Konversi.jd_kmString(jd_wd));
        printDms("Deklinasi Matahari",mail_syamsi);
        printDms("Equation of Time",e);
        printDms("SemiDiameter",sd);

        double wh = wd+e-(bwd-bujur)/15;
        double t1 = (wh-12)*15;
        double h = toDegrees(asin(sin(toRadians(lintang))*sin(toRadians(mail_syamsi))+cos(toRadians(lintang))*cos(toRadians(mail_syamsi))*cos(toRadians(t1))));
        double az1 = toDegrees(atan(tan(toRadians(mail_syamsi))*cos(toRadians(lintang))/sin(toRadians(t1))-sin(toRadians(lintang))/tan(toRadians(t1))));
        double az2; if (wh<12)az2=az1+90; else az2=az1+270;
        double az = 360-az2;

        printDms("wh",wh);
        printDms("t1",t1);
        printDms("h",h);
        printDms("az1",az1);
        printDms("az2",az2);
        printDms("az",az);

        double sudut_b = abs(az2-270);
        double r = 50;
        double jarak_a_c = r/sin(toRadians((180-sudut_b)/2))*sin(toRadians(sudut_b));
        print("sudut b",sudut_b);
        print("r",r);
        print("jarak a-c",jarak_a_c);

        double sisi_a = 90-lintang;
        double sisi_b = 90-lintangKB;
        double sisi_c = bujur-bujurKB;

        double dr = PI/180;
        //     double aq = toDegrees(atan(1/tan(toRadians(sisi_b)))*sin(toRadians(sisi_a))/sin(toRadians(sisi_c))-cos(toRadians(sisi_a))*1/tan(toRadians(sisi_c)));
        double aq =atan(1/tan(sisi_b*dr)*sin(sisi_a*dr)/sin(sisi_c*dr)-cos(sisi_a*dr)*1/tan(sisi_c*dr))*180/PI;
        double az_kb; if (sisi_c<0)az_kb=aq+90; else az_kb=aq+270;

        printDms("lt_m",lintang);
        printDms("bj_m",bujur);
        printDms("lt_kb",lintangKB);
        printDms("bj_kb",bujurKB);

        print("sisi_a",sisi_a);
        print("sisi_b",sisi_b);
        print("sisi_c",sisi_c);
        printDms("aq",aq);
        printDms("az_kb",az_kb);

        double m_jarak_kb = toDegrees(acos(sin(toRadians(lintang))*sin(toRadians(lintangKB))+cos(toRadians(lintang))*cos(toRadians(lintangKB))*cos(toRadians(sisi_c))));
        double km = m_jarak_kb/360*6.283185307*6378.388;
        print("m_jarak_kb",m_jarak_kb);
        print("km",km);


        double sb = toDegrees(atan(1/(1/tan(toRadians(aq))*sin(toRadians(lintang)))));
        double bq = ((toDegrees(acos(1/tan(toRadians(lintang))*tan(toRadians(mail_syamsi))*cos(toRadians(sb))))+sb)/15+(12-e-(bujur-bwd)/15));

        print("sb",sb);
        print("bq",bq);
        printDms("bq",bq);
        double kesalahan_sudut = 5;
        double serong_dari_makkah = km/sin(toRadians((180-kesalahan_sudut)/2))*sin(toRadians(kesalahan_sudut));
        print("Serong dari makkah (km)",serong_dari_makkah);


//        print("",);
//        print("",);
//        print("",);
//        print("",);
//        print("",);
//        print("",);



    }




    static void printDms(String teks, double dms_value){
        String format = "%-25s%5s%n";
        System.out.printf(format,teks,":"+Konversi.dmsMili(dms_value));

    }

    static void printHms(String teks, double hms_value){
        String format = "%-25s%5s%n";
        System.out.printf(format,teks,":"+Konversi.hms(hms_value));

    }

    static void print(String teks, double value){
        String format = "%-25s%5s%n";
        System.out.printf(format,teks,":"+value);

    }

    static void printDigit(String teks, double value,int digit){
        String format = "%-25s%5s%n";
        System.out.printf(format,teks,":"+String.format(Locale.JAPAN,"%."+digit+"f",value));

    }

    static double abs360(double val){
        if (val<0)val+=360;
        return val;
    }


}

