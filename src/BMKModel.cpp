#include "BMKModel.h"
#include <iostream>
using namespace std;

// Global flags and CFF renormalization
bool hasH = true;
bool hasHt = false;
bool hasE = false;
bool hasEt = false;
double renormImag = 1.0;
double renormReal = 1.0;

// Constructor
BMKModel::BMKModel(double rq_beam, double rL_beam, double rL_target,
                   double rEB, double rxB, double rQ2, double rt, double rphi,
                   double rtheta_Tpol, double rphi_Tpol) {
    q_beam = rq_beam;
    L_beam = rL_beam;
    L_target = rL_target;
    EB = rEB;
    xB = rxB;
    Q2 = rQ2;
    t = -TMath::Abs(rt);
    phi = rphi;
    theta_Tpol = PI * rtheta_Tpol / 180.;
    phi_Tpol = PI * rphi_Tpol / 180.;
    VERB = true;
    setSecondaryVars();
}

void BMKModel::setPrimaryVars(double rq_beam, double rL_beam, double rL_target,
                              double rEB, double rxB, double rQ2, double rt, double rphi,
                              double rtheta_Tpol, double rphi_Tpol) {
    q_beam = rq_beam;
    L_beam = rL_beam;
    L_target = rL_target;
    EB = rEB;
    xB = rxB;
    Q2 = rQ2;
    t = -TMath::Abs(rt);
    phi = rphi;
    theta_Tpol = PI * rtheta_Tpol / 180.;
    phi_Tpol = PI * rphi_Tpol / 180.;
    setSecondaryVars();
}

void BMKModel::setSecondaryVars(void) {
    xi = xB * (1 + 0.5 * t / Q2) / (2 - xB + xB * t / Q2);
    phi_BMK = PI * (1 - phi / 180.);
    nu = Q2 / (2. * M * xB);
    y = nu / EB;
    eps = 2 * xB * M / TMath::Sqrt(Q2);
    Jacob = y / Q2;
    eps2 = eps * eps;
    t_min = -Q2 * (2. * (1. - xB) * (1. - TMath::Sqrt(1. + eps2)) + eps2)
            / (4. * xB * (1. - xB) + eps2);

    if (t_min < t) K2 = 0;
    else K2 = -(t - t_min) / Q2 * (1. - xB) * (1. - y - 0.25 * y * y * eps2) *
              (TMath::Sqrt(1. + eps2) + (4. * xB * (1. - xB) + eps2) / (4. * (1. - xB)) * (t - t_min) / Q2);
    K = TMath::Sqrt(K2);

    J = (1. - y - 0.5 * y * eps2) * (1. + t / Q2) - (1. - xB) * (2. - y) * t / Q2;

    if (t_min < t) Ktild2 = 0;
    else Ktild2 = (t_min - t) * ((1 - xB) * (1 + eps2) +
                (t_min - t) * (eps2 + 4 * (1 - xB) * xB) / (4 * Q2));
    Ktilda = TMath::Sqrt(Ktild2);

    F1 = GetF1(t);
    F2 = GetF2(t);
    FF_comb1 = F1 * F1 - t * F2 * F2 / (4 * M * M);
    FF_comb2 = TMath::Power(F1 + F2, 2);
    FF_comb3 = F1 + t * F2 / (4 * M * M);

    ImH = GetImH(xi, t);
    ImHt = GetImHt(xi, t);
    ImE = GetImE(xi, t);
    ImEt = GetImEt(xi, t);
    ReH = GetReH(xi, t);
    ReHt = GetReHt(xi, t);
    ReE = GetReE(xi, t);
    ReEt = GetReEt(xi, t);

    if (VERB) {
        cout << "[DVCS] xB=" << xB << " Q2=" << Q2 << " t=" << t << " phi=" << phi << endl;
    }
}

double BMKModel::CrossSection() {
    double prefactor = 1e6 * hbarc2 * 2 * PI * Jacob * alpha3 * xB * y /
                       (16 * PI * PI * Q2 * TMath::Sqrt(1 + eps2));
    return prefactor * T2();
}

double BMKModel::BSA() {
    q_beam = -1;
    L_target = 0;
    L_beam = 1;
    std::cout<<L_beam<<std::endl;
    double xsec1 = CrossSection();
    L_beam = -1;
    std::cout<<L_beam<<std::endl;
    double xsec2 = CrossSection();
    return (xsec1 - xsec2) / (xsec1 + xsec2);
}

double BMKModel::T2() {
    return BH2() + DVCS2() - q_beam * BHDVCS();
}

double BMKModel::BH2() {
    double denom = xB * xB * y * y * t * BHP1() * BHP2() * TMath::Power(1 + eps2, 2);
    double harmonic = c0_BH() + c1_BH() * TMath::Cos(phi_BMK) + c2_BH() * TMath::Cos(2 * phi_BMK);
    return harmonic / denom;
}

double BMKModel::DVCS2() {
    double denom = y * y * Q2;
    return c0_DVCS() / denom;
}

double BMKModel::BHDVCS() {
    double denom = xB * y * y * y * t * BHP1() * BHP2();
    double harmonic = c0_I() + c1_I() * TMath::Cos(phi_BMK) + s1_I() * TMath::Sin(phi_BMK);
    return harmonic / denom;
}

double BMKModel::c0_BH() {
    double term = (2 + eps2) * (4 * xB * xB * M * M * TMath::Power(1 + t / Q2, 2) / t + 4 * (1 - xB) * (1 + xB * t / Q2));
    term += 4 * xB * xB * (xB + (1 - xB + 0.5 * eps2) * TMath::Power(1 - t / Q2, 2) - xB * (1 - 2 * xB) * t * t / (Q2 * Q2)) * FF_comb2;
    term *= TMath::Power(2 - y, 2);
    term += 8 * K2 * ((2 + 3 * eps2) * Q2 * FF_comb1 / t + 2 * xB * xB * FF_comb2);
    return term;
}

double BMKModel::c1_BH() {
    return 8 * K * (2. - y) * ((4 * xB * xB * M * M / t - 2 * xB - eps2) * FF_comb1
            + 2 * xB * xB * (1 - (1 - 2 * xB) * t / Q2) * FF_comb2);
}

double BMKModel::c2_BH() {
    return 8 * xB * xB * K2 * (4 * M * M * FF_comb1 / t + 2 * FF_comb2);
}

double BMKModel::c0_I() {
    double CFF_unp0 = F1 * ReH - t / (4 * M * M) * F2 * ReE + xB * (F1 + F2) * ReHt / (2 - xB + xB * t / Q2);
    double CFF_unpV = xB * (F1 + F2) * (ReH + ReE) / (2 - xB + xB * t / Q2);
    double CFF_unpA = xB * (F1 + F2) * ReHt / (2 - xB + xB * t / Q2);

    double term0 = Ktild2 / Q2 * TMath::Power(2 - y, 2) / TMath::Sqrt(1 + eps2);
    term0 += t / Q2 * (1 - y - 0.25 * eps2 * y * y) * (2 - xB)
        * (1 + (2 * xB * (2 - xB + 0.5 * (TMath::Sqrt(1 + eps2) - 1) + 0.5 * eps2 / xB) * t / Q2 + eps2)
        / ((2 - xB) * (1 + TMath::Sqrt(1 + eps2))));
    term0 *= -4 * (2 - y) * (1 + TMath::Sqrt(1 + eps2)) / TMath::Power(1 + eps2, 2);

    double termV = (1 - y - 0.25 * y * y * eps2) * 0.5 * (1 + TMath::Sqrt(1 + eps2)) * (1 + t / Q2)
        * (1 + (TMath::Sqrt(1 + eps2) - 1 + 2 * xB) / (1 + TMath::Sqrt(1 + eps2)) * t / Q2);
    termV += TMath::Power(2 - y, 2) * Ktild2 / (TMath::Sqrt(1 + eps2) * Q2);
    termV *= 8 * (2 - y) * xB * t / (TMath::Power(1 + eps2, 2) * Q2);

    double termA = 0.5 * (1 + TMath::Sqrt(1 + eps2)) * (1 + TMath::Sqrt(1 + eps2) - xB +
        (TMath::Sqrt(1 + eps2) - 1 + xB * (3 + TMath::Sqrt(1 + eps2) - 2 * xB)
        / (1 + TMath::Sqrt(1 + eps2))) * t / Q2) - 2 * Ktild2 / Q2;
    termA = 8 * (2 - y) / TMath::Power(1 + eps2, 2) * t / Q2 *
        (TMath::Power(2 - y, 2) * Ktild2 / (TMath::Sqrt(1 + eps2) * Q2) *
         0.5 * (1 + TMath::Sqrt(1 + eps2) - 2 * xB)
         + (1 - y - 0.25 * eps2 * y * y) * termA);

    return term0 * CFF_unp0 + termV * CFF_unpV + termA * CFF_unpA;
}

double BMKModel::c1_I() {
    double CFF_unp0 = F1 * ReH - t / (4 * M * M) * F2 * ReE + xB * (F1 + F2) * ReHt / (2 - xB + xB * t / Q2);
    double CFF_unpV = xB * (F1 + F2) * (ReH + ReE) / (2 - xB + xB * t / Q2);
    double CFF_unpA = xB * (F1 + F2) * ReHt / (2 - xB + xB * t / Q2);

    double term0 = -4 * K * (2 - 2 * y + y * y + 0.5 * eps2 * y * y)
        * (1 + TMath::Sqrt(1 + eps2) - eps2) / TMath::Power(1 + eps2, 2.5)
        * (1 - (1 - 3 * xB) * t / Q2 +
           (1 - TMath::Sqrt(1 + eps2) + 3 * eps2) / (1 + TMath::Sqrt(1 + eps2) - eps2) * xB * t / Q2);

    double term0b = -16 * K * (1 - y - 0.25 * y * y * eps2) / TMath::Power(1 + eps2, 2.5) *
        ((1 + (1 - xB) * 0.5 * (TMath::Sqrt(1 + eps2) - 1) / xB) * xB * t / Q2 - 0.75 * eps2);

    double termV = 16 * K * xB * t / (Q2 * TMath::Power(1 + eps2, 2.5)) *
        (TMath::Power(2 - y, 2) * (1 - (1 - 2 * xB) * t / Q2) +
         (1 - y - 0.25 * eps2 * y * y) * 0.5 * (1 + TMath::Sqrt(1 + eps2) - 2 * xB) * (t - t_min) / Q2);

    double termA = -16 * K * t / (Q2 * TMath::Power(1 + eps2, 2)) *
        ((1 - y - 0.25 * eps2 * y * y) *
         (1 - (1 - 2 * xB) * t / Q2 +
          (4 * xB * (1 - xB) + eps2) * 0.25 / TMath::Sqrt(1 + eps2) * (t - t_min) / Q2)
         - TMath::Power(2 - y, 2) *
         (1 - 0.5 * xB + 0.25 * (1 + TMath::Sqrt(1 + eps2) - 2 * xB) * (1 - t / Q2)
          + (4 * xB * (1 - xB) + eps2) * 0.5 / TMath::Sqrt(1 + eps2) * (t - t_min) / Q2));

    return (term0 + term0b) * CFF_unp0 + termV * CFF_unpV + termA * CFF_unpA;
}

double BMKModel::s1_I() {
    double CFF_unp0 = F1 * ImH - t / (4 * M * M) * F2 * ImE + xB * (F1 + F2) * ImHt / (2 - xB + xB * t / Q2);
    double CFF_unpV = xB * (F1 + F2) * (ImH + ImE) / (2 - xB + xB * t / Q2);
    double CFF_unpA = xB * (F1 + F2) * ImHt / (2 - xB + xB * t / Q2);

    double term0 = 8 * L_beam * K * (2 - y) * y / (1 + eps2)
        * (1 + (1 - xB + 0.5 * (TMath::Sqrt(1 + eps2) - 1)) / (1 + eps2) * (t - t_min) / Q2);

    double termV = -8 * L_beam * K * (2 - y) * y / TMath::Power(1 + eps2, 2) *
        xB * t / Q2 * (TMath::Sqrt(1 + eps2) - 1 + (1 + TMath::Sqrt(1 + eps2) - 2 * xB) * t / Q2);

    double termA = 8 * L_beam * K * (2 - y) * y / (1 + eps2) * t / Q2 *
        (1 - (1 - 2 * xB) * (1 + TMath::Sqrt(1 + eps2) - 2 * xB) * 0.5 / TMath::Sqrt(1 + eps2) * (t - t_min) / Q2);

    return term0 * CFF_unp0 + termV * CFF_unpV + termA * CFF_unpA;
}

double BMKModel::c0_DVCS() {
    double prefactor = 2 * (2 - 2 * y + y * y + 0.5 * eps2 * y * y) / (1 + eps2);
    double C = TMath::Power(2 - xB + xB * t / Q2, -2);

    double term = 4 * (1 - xB) * (1 + xB * t / Q2) * (ReH * ReH + ImH * ImH + ReHt * ReHt + ImHt * ImHt);
    term += (2 + t / Q2) * eps2 * (ReHt * ReHt + ImHt * ImHt);
    term -= 0.25 * t / (M * M) * (ReE * ReE + ImE * ImE);
    term -= xB * xB * (TMath::Power(1 + t / Q2, 2) * (2 * ReH * ImE + 2 * ReE * ImH + ReE * ReE + ImE * ImE)
                       + 2 * ReHt * ImEt + 2 * ReEt * ImHt
                       + 0.25 * t / (M * M) * (ReE * ReE + ImE * ImE));

    return prefactor * C * term;
}

double BMKModel::BHP1() {
    return -(J + 2 * TMath::Sqrt(K2) * TMath::Cos(phi_BMK)) / (y * (1 + eps2));
}

double BMKModel::BHP2() {
    return 1 + t / Q2 - BHP1();
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////

double GetF1(double T){
	T = TMath::Abs(T);
	double tau = T/(4*BMKModel::M*BMKModel::M);
	return (GetGEP(tau)+tau*GetGMP(tau))/(1+tau);
}

double GetF2(double T){
	T = TMath::Abs(T);
	double tau = T/(4*BMKModel::M*BMKModel::M);
	return (GetGMP(tau)-GetGEP(tau))/(1+tau);
}

double GetGMP(double tau){
	double res;
	double Q2 = 4*BMKModel::M*BMKModel::M*tau;
	double tcut = 4*0.13957*0.13957;
	double t0 = -0.7;
	double z = ( TMath::Sqrt(tcut+Q2) - TMath::Sqrt(tcut-t0) ) / ( TMath::Sqrt(tcut+Q2) + TMath::Sqrt(tcut-t0) ) ;
	res = 0;
	res += 0.264142994136;
	res += -1.095306122120*z;
	res += 1.218553781780 * TMath::Power(z,2);
	res += 0.661136493537 * TMath::Power(z,3);
	res += -1.405678925030 * TMath::Power(z,4);
	res += -1.356418438880 * TMath::Power(z,5);
	res += 1.447029155340 * TMath::Power(z,6);
	res += 4.235669735900 * TMath::Power(z,7);
	res += -5.334045653410 * TMath::Power(z,8);
	res += -2.916300520960 * TMath::Power(z,9);
	res += 8.707403067570 * TMath::Power(z,10);
	res += -5.706999943750 * TMath::Power(z,11);
	res += 1.280814375890 * TMath::Power(z,12);
	return BMKModel::muP*res;
}

double GetGEP(double tau){
	double res;
	double Q2 = 4*BMKModel::M*BMKModel::M*tau;
	double tcut = 4*0.13957*0.13957;
	double t0 = -0.7;
	double z = ( TMath::Sqrt(tcut+Q2) - TMath::Sqrt(tcut-t0) ) / ( TMath::Sqrt(tcut+Q2) + TMath::Sqrt(tcut-t0) ) ;
	res = 0;
	res += 0.239163298067;
	res += -1.109858574410 * z;
	res += 1.444380813060 * TMath::Power(z,2);
	res += 0.479569465603 * TMath::Power(z,3);
	res += -2.286894741870 * TMath::Power(z,4);
	res += 1.126632984980 * TMath::Power(z,5);
	res += 1.250619843540 * TMath::Power(z,6);
	res += -3.631020471590 * TMath::Power(z,7);
	res += 4.082217023790 * TMath::Power(z,8);
	res += 0.504097346499 * TMath::Power(z,9);
	res += -5.085120460510 * TMath::Power(z,10);
	res += 3.967742543950 * TMath::Power(z,11);
	res += -0.981529071103 * TMath::Power(z,12);
	return res;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

double GetImH(double xi, double t){
	//return 0;
	double res = 0;
	if(hasH){
		res=(1-xi)/t;
		// VALENCE
		double r = 0.9;
		double alpha = 0.43 + 0.85*t;
		double n = 1.35;
		double b = 0.4;
		double Mm2 = 0.64;
		double P = 1;
		res = TMath::Pi()*5.0/9.0 * n * r /(1+xi) * TMath::Power( 2*xi/(1+xi),-alpha ) * TMath::Power( (1-xi)/(1+xi) , b ) * TMath::Power( 1 - (1-xi)/(1+xi) *t/Mm2 , -P );
		// SEA
		if(false){
			r = 1;
			alpha = 1.13 + 0.15*t;
			n = 1.5;
			b = 0.4;
			Mm2 = 0.5;
			P = 2;
			res += TMath::Pi()*2.0/9.0 * n * r /(1+xi) * TMath::Power( 2*xi/(1+xi),-alpha ) * TMath::Power( (1-xi)/(1+xi) , b ) * TMath::Power( 1 - (1-xi)/(1+xi) *t/Mm2 , -P );
		}
		res *= 2.0;// correction from VGG
	}
	return res * renormImag;
}

double GetImHt(double xi, double t){
	//return 0;
	double res = 0;
	if(hasHt){
		res=0.1*(1-xi)/t;
		double r = 7;
		double alpha = 0.43 + 0.85*t;
		double n = 0.6;
		double b = 2;
		double Mm2 = 0.8;
		double P = 1;
		res = TMath::Pi()*5.0/9.0 * n * r /(1+xi) * TMath::Power( 2*xi/(1+xi),-alpha ) * TMath::Power( (1-xi)/(1+xi) , b ) * TMath::Power( 1 - (1-xi)/(1+xi) *t/Mm2 , -P );
		res *= 0.4;// correction from VGG
	}
	return res * renormImag;
}

double GetImE(double xi, double t){
	//return 0;
	// from HERMES CFF paper 1301.1230 argue ImE = 0.5 x ImH
	double res = 0;
	if(hasE){
		res=1.5*(1-xi)/t;
		// VALENCE
		double r = 0.9;
		double alpha = 0.43 + 0.85*t;
		double n = 1.35;
		double b = 0.4;
		double Mm2 = 0.64;
		double P = 1;
		res = TMath::Pi()*5.0/9.0 * n * r /(1+xi) * TMath::Power( 2*xi/(1+xi),-alpha ) * TMath::Power( (1-xi)/(1+xi) , b ) * TMath::Power( 1 - (1-xi)/(1+xi) *t/Mm2 , -P );
		if(false){
			// SEA
			r = 1;
			alpha = 1.13 + 0.15*t;
			n = 1.5;
			b = 0.4;
			Mm2 = 0.5;
			P = 2;
			res += TMath::Pi()*2.0/9.0 * n * r /(1+xi) * TMath::Power( 2*xi/(1+xi),-alpha ) * TMath::Power( (1-xi)/(1+xi) , b ) * TMath::Power( 1 - (1-xi)/(1+xi) *t/Mm2 , -P );
		}
	}
	return res * renormImag;
}

double GetImEt(double xi, double t){
	return 0;
}

double GetReH(double xi, double t){
	double res = 0;
	if(hasH){
		res  = -12*xi*TMath::Power(1-xi,2) *TMath::Sqrt(TMath::Abs(t))/TMath::Power(1-t/0.7,2);
		res += -3*TMath::Power(1-xi,4) /TMath::Power(1-t/1.1,2);
	}
	return res / ( 1 + TMath::Power( t/0.8 , 4) );
	return res*renormReal / ( 1 + TMath::Power( t , 4) );
}

double GetReHt(double xi, double t){
	double res = 0;
	if(hasHt)res=-12*xi*TMath::Power(1-xi,2) /TMath::Power(1-t/1.5,2);
	return res*renormReal;
}

double GetReE(double xi, double t){
	double res = 0;
	if(hasE){
		res = -7 * xi*TMath::Power(1-xi,2) *TMath::Sqrt(TMath::Abs(t))/TMath::Power(1-t/0.7,2);
		res += -3*TMath::Power(1-xi,2)/TMath::Power(1-t/1.2,2);

	}
	return res*renormReal / ( 1 + TMath::Power( t , 4) );
	return res*renormReal;
}

double GetReEt(double xi, double t){
	//return 0;
	double res = 0;
	if(hasEt)res=10/t * 1/( 1 + TMath::Power( 3*xi , 4) );
	return res*renormReal;
}
