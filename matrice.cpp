//
//  matrice.cpp
//  echantillonnage_preferentiel
//
//  Created by DOUGE Louis on 4/14/16.
//  Copyright © 2016 DOUGE Louis. All rights reserved.
//

#include "matrice.hpp"
#include <cmath>

matrice::matrice(int dl, int dc, double v) {
    m_dim_l=dl;m_dim_c=dc;m_val=NULL;
    int d=dl*dc;
    if(d<=0)return; //aucune allocation matrice vide
    m_val=new double [d];
    for(int k=0;k<d;k++){
        m_val[k]=v;
    } //affectation
}

matrice::matrice(const matrice & V){
    m_dim_l=V.m_dim_l;m_dim_c=V.m_dim_c;m_val=NULL;
    int d=m_dim_c*m_dim_l;
    if (d<=0) {
        return ;
    }
    
    m_val=new double[d];
    for (int k=0; k<d; k++) {
        m_val[k]=V.m_val[k]; //recopie
    }
}

matrice::~matrice(){
    if(m_val!=NULL)delete [] m_val;
}

//operateur d'assignation U=V
matrice& matrice::operator=(const matrice & B) {
    delete [] m_val; //libère la mémoire
    m_dim_l=B.m_dim_l;m_dim_c=B.m_dim_c;m_val=NULL;
    int d=m_dim_l*m_dim_c;
    if (d<=0) {
        return *this;
    }
    m_val=new double[d];
    for (int k=0; k<d; k++) {
        m_val[k]=B.m_val[k]; //recopie
    }
    return *this;
}

//acces a un element
double& matrice::operator()(int i, int j) const {
    return m_val[(j-1)*m_dim_l+i-1];
}
 //le reste a été fait dans le .hpp

//operations algebriques internes
matrice& matrice::operator+=(double x){ //A+x
    double *pA=m_val;
    int d=m_dim_l*m_dim_c;
    for (int i=1; i<=d; i++,pA++) {
        (*pA)+=x;
    }
    return *this;
}

matrice& matrice::operator-=(double x){//A-x
    double *pA=m_val;
    int d=m_dim_l*m_dim_c;
    for (int i=1; i<=d; i++,pA++) {
        (*pA)-=x;
    }
    return *this;
}


matrice& matrice::operator*=(double x){ //A*x
    double *pA=m_val;
    int d=m_dim_l*m_dim_c;
    for (int i=1; i<=d; i++,pA++) {
        (*pA)*=x;
    }
    return *this;
}


matrice& matrice::operator/=(double x){ //A/x
    
    if(x==0){
        cout << "Division par zero dans A/=x";
        exit(-1);
    }
    
    double *pA=m_val;
    int d=m_dim_l*m_dim_c;
    for (int i=1; i<=d; i++,pA++) {
        (*pA)/=x;
    }
    return *this;
}

matrice& matrice::operator+=(const matrice& B){ //A+=B
    double *pA=m_val;
    double *pB=B.m_val;
    int d=m_dim_l*m_dim_c;
    for (int i=1; i<=d; i++,pA++,pB++) {
        (*pA)+=(*pB);
    }
    return *this;
}

matrice& matrice::operator-=(const matrice& B){ //A-=B
    double *pA=m_val;
    double *pB=B.m_val;
    int d=m_dim_l*m_dim_c;
    for (int i=1; i<=d; i++,pA++,pB++) {
        (*pA)-=(*pB);
    }
    return *this;
}

//operateur d'egalite
bool matrice::operator==(const matrice & B){
    if (m_dim_l!=B.m_dim_l || m_dim_c!=B.m_dim_c) {
        return false;
    }
    
    double *pA=m_val;
    double *pB=B.m_val;
    int d=m_dim_l*m_dim_c;
    for (int i=1; i<=d; pA++, pB++) {
        if ((*pA)!=(*pB)) {
            return false;
        }
    }
    return true;
}

bool matrice::operator!=(const matrice & V){
    return !((*this)==V);
}

//impression

void matrice::print() const{
    cout << "matrice de dimension " << m_dim_l << "par" << m_dim_c << endl;
    for (int i=1; i<=m_dim_l; i++) {
        for (int j=1; j<=m_dim_c; j++) {
            cout << (*this)(i,j) << " ";
            cout << endl;
        }
        cout << endl;
    }

}


//---------------------------
// fonctions externes à la classe
//---------------------------

matrice operator-(const matrice & A){ //-A
    matrice R(A.getDim_l(),A.getDim_c(),0);
    return R-=A;
}

matrice operator+(const matrice & A){ //+A
    return A;
}

matrice operator+(const matrice & A, double x) { //A+x
    matrice R(A);
    return R+=x;
}

matrice operator-(const matrice & A, double x) { //A-x
    matrice R(A);
    return R-=x;
}

matrice operator*(const matrice & A, double x) { //A*x
    matrice R(A);
    return R*=x;
}

matrice operator/(const matrice & A, double x) { //A/x
    matrice R(A);
    return R/=x;
}

matrice operator+(double x, const matrice& A){ //x+A
    return A+x;
}

matrice operator-(double x, const matrice & A){ //x-A
    matrice R(A.getDim_l(),A.getDim_c(),x);
    return R-=A;
}

matrice operator*(double x, const matrice &A){ //x*A
    return A*x;
}


matrice operator+(const matrice &A, const matrice & B){ //A+B
    matrice R(A);
    return R=+B;
}

matrice operator-(const matrice &A, const matrice & B){ //A-B
    matrice R(A);
    return R=-B;
}

//operations de flux
ostream& operator <<(ostream &os, const matrice& A) { //os << A
    os<<"matrice de dimension "<<A.getDim_l()<<" par "<<A.getDim_c()<<endl;
    for (int i=1; i<=A.getDim_l(); i++) {
        for (int j=1; j<=A.getDim_c(); j++) {
            os << A(i,j)<<" ";
        }
        os<<endl;
    }
    os<<endl;
    return os;
}

istream & operator>>(istream &is, matrice &A){ //is >>A
    for (int i=1; i<=A.getDim_l(); i++) {
        for (int j=1; j<=A.getDim_c(); j++) {
            is>>A(i,j);
        }
    }
    return is;
}


//produit matrice vecteur AX
vector<double> operator*(const matrice& A, const vector<double> & V){
    vector<double> R(A.getDim_l(),0); //vecteur contenant le resultat
    for (int i=1; i<=A.getDim_l(); i++) {
        for (int j=1; j<=A.getDim_c(); j++) {
            R[i-1]+=A(i,j)*V[j-1];
        }
    }
    return R;
}


//calcul du déterminant d'une matrice carrée
double matrice::det_3d() const{
    
    if (this->m_dim_c!=this->m_dim_l) {
        cout << "Determinant d'une matrice non carree ... Oups !" << endl;
        exit(-1);
    }
    
    if (this->m_dim_c!=3) {
        cout << "Cette fonction ne calcule que les déterminants de matrice d'ordre 3" << endl;
        exit(-1);
    }
    else{
        return (*this)(1,1)*((*this)(2,2)*(*this)(3,3)-(*this)(2,3)*(*this)(3,2))-(*this)(1,2)*((*this)(2,1)*(*this)(3,3)-(*this)(2,3)*(*this)(3,1))+(*this)(1,3)*((*this)(2,1)*(*this)(3,2)-(*this)(2,2)*(*this)(3,1));
    }
}

//inverse d'une matrice 3x3

matrice matrice::inv_3d() const{
    if (this->det_3d()==0) {
        cout << " La matrice suivante : " << endl;
        this->print();
        cout << endl;
        cout << "n'est pas inversible" << endl;
        exit(-1);
    }
    else{
        matrice Inverse (3,3,0); //matrice inverse que l'on veut renvoyer à la fin
        
        Inverse(1,1)=(*this)(2,2)*(*this)(3,3)-(*this)(2,3)*(*this)(3,2);
        Inverse(1,2)=(*this)(1,3)*(*this)(3,2)-(*this)(1,2)*(*this)(3,3);
        Inverse(1,3)=(*this)(1,2)*(*this)(2,3)-(*this)(1,1)*(*this)(2,3);
        
        Inverse(2,1)=(*this)(3,2)*(*this)(3,1)-(*this)(2,1)*(*this)(3,3);
        Inverse(2,2)=(*this)(1,1)*(*this)(3,3)-(*this)(1,3)*(*this)(3,1);
        Inverse(2,3)=(*this)(1,3)*(*this)(2,1)-(*this)(1,1)*(*this)(2,3);
        
        Inverse(3,1)=(*this)(2,1)*(*this)(3,2)-(*this)(2,2)*(*this)(3,1);
        Inverse(3,2)=(*this)(1,2)*(*this)(3,1)-(*this)(1,1)*(*this)(3,2);
        Inverse(3,3)=(*this)(1,1)*(*this)(2,2)-(*this)(1,2)*(*this)(2,1);
        
        return Inverse/(this->det_3d());
       
        
        
    }
}
//decomposition de cholesky 3D, très facilement extensible au cas n quelconque (changer les bornes des boucles for et assurer la concordance des dimensions des matrices en jeu)
matrice matrice::cholesky_3d() const {
    matrice Factorisation(3,3,0);
    //definition de la premiere colonne
    Factorisation(1,1)=sqrt((*this)(1,1));
    for (int j=2; j<=3; j++) { // j est le numero de la ligne
        Factorisation(j,1)=(*this)(1,j)/(sqrt((*this)(1,1)));
    }
    //definition des autres colonnes
    for (int i=2; i<=3; i++) {
        //on calcule d'abord une somme interne pour les termes diagonaux
        double sum=0;
        for (int k=1; k<=i-1; k++) {
            sum+=pow(Factorisation(i,k), 2.0);
        }
        if((*this)(i,i)-sum < 0){
            cout << " Il n'existe pas de décomposition de Cholesky pour cette matrice. " << endl;
            exit(-1);
        }
        Factorisation(i,i)=sqrt((*this)(i,i)-sum);
        //de même pour les termes non diagonaux
        for (int j=i+1; j<=3; j++) {
            //somme préalable
            double sum_2=0;
            for (int k=1; k<=i-1; k++) {
                sum_2+=Factorisation(i,k)*Factorisation(j,k);
            }
            Factorisation(j,i)=((*this)(i,j)-sum_2)/(Factorisation(i,i));
        }
    }
    return Factorisation;
}

//transposée
matrice matrice::transpose_3d() const{
    matrice transpose=*this;
    transpose(1,2)=(*this)(2,1);
    transpose(1,3)=(*this)(3,1);
    transpose(2,1)=(*this)(1,2);
    transpose(2,3)=(*this)(3,2);
    transpose(3,1)=(*this)(1,3);
    transpose(3,2)=(*this)(2,3);
    
    return transpose;
}

