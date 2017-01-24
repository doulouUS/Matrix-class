//
//  matrice.hpp
//  echantillonnage_preferentiel
//
//  Created by DOUGE Louis on 4/14/16.
//  Copyright © 2016 DOUGE Louis. All rights reserved.
//

#ifndef matrice_hpp
#define matrice_hpp

#include <stdio.h>

#include <iostream>
#include <vector>

using namespace std;

//---------------------------------
//      class matrice (définition)
//---------------------------------

class matrice {
private:
    int m_dim_l,m_dim_c;//dimension de la matrice
    double* m_val;
    
public:
    
    //constructeurs/destructeurs
    matrice():m_dim_l(0), m_dim_c(0), m_val(NULL) {} //constructeur par défaut
    matrice(int l, int c,double v=0);//constructeur dimensions et val. initiale
    matrice (const matrice & other); //constructeur par copie
    ~matrice();
    
    //assignation
    matrice& operator=(const matrice & other);
    
    //accès à un élément
    double& operator()(int l, int c) const;//acces au terme (i,j) en partant de (1,1)
    //accesseurs en lecture aux dimensions de la matrice
    
    int getDim_l() const{return m_dim_l;}
    int getDim_c() const{return m_dim_c;}
    
    //opérations algébriques internes
    matrice& operator+=(double somme);
    matrice& operator-=(double moins);
    matrice& operator*=(double produit);
    matrice& operator/=(double divise);
    
    matrice& operator+=(const matrice& B);
    matrice& operator-=(const matrice& B);
    
    //determinant d'une matrice 3x3
    double det_3d()const;
    //inverse d'une matrice 3x3
    matrice inv_3d() const;
    //decomposition de cholesky
    matrice cholesky_3d() const;
    //transposé
    matrice transpose_3d() const;
    
    //test entre deux matrices
    bool operator==(const matrice & test);
    bool operator!=(const matrice & test);
    
    //impression
    void print() const;
    
};

//---------------------------------
//      fonctions externes à la classe
//---------------------------------


matrice operator-(const matrice & A);//-A
matrice operator+(const matrice & A);//+A
matrice operator+(const matrice & A, double x);//A+x
matrice operator-(const matrice & A, double x);//A-x
matrice operator*(const matrice & A,double x);//A*x
matrice operator/(const matrice & A, double x);//A/x
matrice operator+(double x,const matrice & A);//x+A
matrice operator-(double x,const matrice & A);//x-A
matrice operator*(double x,const matrice & A);//x*A
matrice operator+(const matrice & A,const matrice & B);//A+B
matrice operator-(const matrice & A,const matrice & B);//A-B
ostream & operator<<(ostream & os, const matrice & A); //os<<A
istream & operator<<(istream & is, const matrice & A); //is<<A

//produit matrice vecteur
vector<double> operator*(const matrice & A, const vector<double> & V);





#endif /* matrice_hpp */
