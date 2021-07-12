# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 19:23:05 2021

@author: Rafa
"""
import random
from scipy.optimize import fsolve
import numpy as np

Kw= 1E-14
class mono:
    def __init__(self,name,B):
      self.name= name
      self.B = B
      self.order= 1 
          
class di:
    def __init__(self,name,B1,B2):
        self.name = name
        self.B1= B1
        self.B2 = B2
        self.order = 2
        
class tri:
    def __init__(self, name, B1,B2,B3):
        self.name = name
        self.B1= B1
        self.B2 = B2
        self.B3=B3
        self.order= 3

class buff:
    def __init__(self, name, pka):
        self.name= name
        self.pka=pka

class metal:
    def __init__(self,name,Kf):
        self.name = name
        self.Kf= Kf
        
class electrode:
    def __init__(self, nameR,nameO, pot):
        self.nameR = nameR
        self.nameO = nameO
        self.pot = pot

class redoxtitra:
    def __init__(self,Ca,Ct,name,antes,despues):
        self.name = name 
        self.Ca= Ca
        self.Ct=Ct
        self.antes = antes
        self.despues = despues
class preci:
    def __init__(self, name, n, kps,PM):
        self.name = name
        self.n= n
        self.kps = kps 
        self.PM= PM 

def Pregunta1():
    A1,A2,A3,A4,A5,A6,A7,A8,A9=mono("Formico",3.75),mono("Acetico",4.75),di("Carbonico",10.33,16.68),di("Oxalico",3.81, 4.06),tri("Fosforico",12.32,19.53,21.69),tri("Citrico",11.29,18.05,20.31),mono("Benzoico",4.2),di("Malonico",5.69,8.54),tri("Arsenioso",11.49,18.45,20.69)
    Acids=[A1,A2,A3,A4,A5,A6,A7,A8,A9]
    Pregu=random.choice(Acids)
    C=[0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5]
    Ct=float(random.choice(C))
    O=Pregu.order
    print('###################################################')
    print('################# Pregunta 1 ######################')
    print('###################################################')
    print()
    print('Calcula el pH en equilibrio de una solución acuosa de Acido '+Pregu.name+' '+str(Ct)+'M, e ingresa tu respuesta')
    if O==1: 
        B1=10**float(Pregu.B)
        M=0
        def N(x):
            return (1+(B1*x))
        def F(x):
            return x+M
        def G(x):
            return (Kw/x)+(Ct/N(x))
        def Final(x):
            x=10**-x
            return F(x)-G(x)
        z=fsolve(Final,0.20)
    elif O==2:
        B1=10**(float(Pregu.B1))
        B2=10**(float(Pregu.B2))
        M=0.0
        def N(x):
            return (1+(B1*x)+(B2*x**2))
        def F(x):
            return M+x
        def G(x):
            return (2*Ct/N(x))+(x*B1*Ct/N(x))+(Kw/x)
        def Final(x):
            x=10**-x
            return F(x)-G(x)
        z=fsolve(Final,0.2)
        
    elif O==3:
        B1=10**(float(Pregu.B1))
        B2=10**(float(Pregu.B2))
        B3=10**(float(Pregu.B3))
        M=0.0
        def N(x):
            return (1+(B1*x)+(B2*x**2)+(B3*x**3))
        def F(x):
            return M+x
        def G(x):
            return (3*Ct/N(x))+(2*x*B1*Ct/N(x))+(x**2*B2*Ct/N(x))+(Kw/x)
        def Final(x):
            x=10**-x
            return F(x)-G(x)
        z=fsolve(Final,0.2)
    a=float(input())
    z=float(round(z[0], 2))
    err=round(abs(z-a),2) 
    print('Tu respuesta fue '+str(a)+' '+'vs calculado '+str(z)+', te pasaste por '+str(err))
    print()
def Pregunta2():
    B1,B2,B3,B4=buff("H3PO4/H2PO4",2.15),buff("H2PO4/HPO4",7.20),buff("HPO4/PO4",12.35),buff("His/His*",6.00)
    Buffers=[B1,B2,B3,B4]
    Pregu2=random.choice(Buffers)
    C=np.arange(0.1,1,0.1,dtype=float)
    Ca=round(random.choice(C),2)
    Cb=round(random.choice(C),2)
    print('###################################################')
    print('################# Pregunta 2 ######################')
    print('###################################################')
    print()
    print('Calcula el pH de una solucion buffer de'+' '+Pregu2.name+' '+'si las concentraciones de acido y su base conjugada son'+' '+str(Ca)+','+str(Cb)+' '+'respectivamente')
    Calc= float(Pregu2.pka) + np.log10(Cb/Ca)
    z=round(float(Calc),2)
    a=float(input())
    err=round(abs(z-a),2)
    print('Tu respuesta fue '+str(a)+' '+'vs calculado '+str(z)+', te pasaste por '+str(err))
    print()
    
def Pregunta3():
    M1,M2,M3,M4,M5,M6=metal("Mg",4.9E8),metal("Ca",5.0E10),metal("Mn",6.2E13),metal("Co",2.0E16),metal("Cu",6.3E18),metal("Zn",3.2E16)
    Metals=[M1,M2,M3,M4,M5,M6]
    Pregu3=random.choice(Metals)
    pH=np.arange(1,14.1,0.1,dtype=float)
    Val=round(random.choice(pH),2)
    print('###################################################')
    print('################# Pregunta 3 ######################')
    print('###################################################')
    print()
    print('Calcula la constante de formacion condicional para el complejo con EDTA de '+Pregu3.name+' a un pH de '+str(Val))
    def A4(x):
        K1,K2,K3,K4,K5,K6=1,10**-1.5,10**-2,10**-2.69,10**-6.13,10**-10.37
        return (K1*K2*K3*K4*K5*K6)/((x**6)+(K1*x**5)+(K1*K2*x**4)+(K1*K2*K3*x**3)+(K1*K2*K3*K4*x**2)+(K1*K2*K3*K4*K5*x)+(K1*K2*K3*K4*K5*K6))
    z=float(Pregu3.Kf)*A4(10**-Val)
    a=float(input())
    err=round(abs(z-a)*100/z,2)
    print('Tu respuesta fue '+str(a)+' '+'vs calculado'+ ' '+str(z)+', tuviste una desviación porcentual del'+' '+str(err)+'%')
    print()
    
def Pregunta4():
    M1,M2,M3,M4,M5,M6=metal("Mg",4.9E8),metal("Ca",5.0E10),metal("Mn",6.2E13),metal("Co",2.0E16),metal("Cu",6.3E18),metal("Zn",3.2E16)
    Metals=[M1,M2,M3,M4,M5,M6]
    Pregu4=random.choice(Metals)
    Cedta=0.1
    pH=np.arange(1,14.1,0.1,dtype=float)
    Vol=[25,50,75,100]
    C=[0.01,0.02,0.05,0.1,0.5]
    Val=round(random.choice(pH),2)
    V=random.choice(Vol)
    Ct=random.choice(C)
    Veq=Ct*V/Cedta
    print('###################################################')
    print('################# Pregunta 4 ######################')
    print('###################################################')
    print()
    print('Calcula el pM para el punto de equivalencia en una titulación con EDTA (0.1M) del metal '+Pregu4.name+' a un pH de '+str(Val)+' si el vol inicial (mL) y [M] son '+str(V)+' '+str(Ct)+' respectivamente')
    def A4(x):
        K1,K2,K3,K4,K5,K6=1,10**-1.5,10**-2,10**-2.69,10**-6.13,10**-10.37
        return (K1*K2*K3*K4*K5*K6)/((x**6)+(K1*x**5)+(K1*K2*x**4)+(K1*K2*K3*x**3)+(K1*K2*K3*K4*x**2)+(K1*K2*K3*K4*K5*x)+(K1*K2*K3*K4*K5*K6))
    z=float(Pregu4.Kf)*A4(10**-10)
    b=-(Ct*V/(V+Veq))
    coeff=[z,1,b]
    sol=np.roots(coeff)
    sol=float(sol[sol>0][0])
    pM=round(-np.log10(sol),2)
    a=float(input())
    err=round(abs(pM-a),2)
    print('Tu respuesta fue '+str(a)+' '+'vs calculado '+str(pM)+', te pasaste por '+str(err))
    print()
def Pregunta5(): 
    E1= electrode("Ag(+)|Ag","Ag|Ag(+)",0.799)
    E2= electrode("Fe(3+)|Fe(2+)","Fe(2+)|Fe(3+)",0.771)
    E3= electrode("Cd(2+)|Cd","Cd|Cd(2+)",-0.403)
    E4= electrode("Al(3+)|Al","Al|Al(3+)",-1.662)
    E5= electrode("Cu(2+)|Cu","Cu|Cu(2+)",0.337)
    E6= electrode("Co(2+)|Co","Co|Co(2+)",-0.277)
    E7= electrode("Ni(2+)|Ni","Ni|Ni(2+)",-0.250)
    E8= electrode("Zn(2+)|Zn","Zn|Zn(2+)",-0.763)
    Electrodes=[E1,E2,E3,E4,E5,E6,E7,E8]
    Right=random.choice(Electrodes)
    Left=random.choice(Electrodes)
    while Right==Left:
        Left=random.choice(Electrodes)
    else:
        Left=Left 
    print('###################################################')
    print('################# Pregunta 5 ######################')
    print('###################################################')
    print()
    print("Calcula el potencial estandar de celda para la configuración "+(Left.nameO)+'||'+Right.nameR)
    a=float(input())
    z=round(float(Right.pot-Left.pot),3)
    err=round(abs(z-a),3)
    print('Tu respuesta fue '+str(a)+' vs calculado '+str(z),', te pasaste por '+str(err))
    print()
def Pregunta6():
    S1=redoxtitra(0.02, 0.01, "Ascorbico y Fe",0.767,0.390)
    S2=redoxtitra(0.01,0.05,"Cobre y Cerio",1.7,0.161)
    S=[S1,S2]
    System=random.choice(S)
    print('###################################################')
    print('################# Pregunta 6 ######################')
    print('###################################################')
    print()
    if System== S[0]:
        Vol=[10,25]
        V=random.choice(Vol)
        Veq=System.Ca*V/(System.Ct*2)
        Avance=np.arange(0,2*Veq+1,5)
        Vc=(random.choice(Avance))
        print('Para la titulación de Fe(3+) ('+str(V)+'mL) con Acido Ascorbico ('+str(System.Ca)+', '+str(System.Ct)+' M respectivamente), calcula el potencial de celda cuando se ha añadido '+str(Vc)+' mL de titulante; considera un pH de 1 y que se realizó utilizando un electrodo de calomel')
        if Vc < Veq:
            E= System.antes -0.197 - 0.05916*np.log10(Vc/(Veq-Vc))
        elif Vc == Veq:
            E = (System.antes + 2*System.despues - 0.05916*np.log10(1/(10**-0.6)))/3 -0.197
        elif Vc > Veq: 
            E = System.despues - 0.05916*np.log10((Vc-Veq)/(Vc*10**-0.6))
        
    else:
        Vol=[75,100]
        V=random.choice(Vol)
        Veq=System.Ca*V/System.Ct
        Avance=np.arange(0,2*Veq+1,5)
        Vc=(random.choice(Avance))
        print('Para la titulación de Ce(4+) ('+str(V)+'mL) con Cu(+) ('+str(System.Ca)+', '+str(System.Ct)+' M respectivamente), calcula el potencial de celda cuando se ha añadido ' + str(Vc)+' mL de titulante; el sistema se realizó utilizando un electrodo de calomel')
        if Vc < Veq:
            E= System.antes -0.197 -0.05916*np.log10(Vc/(Veq-Vc))
        elif Vc == Veq:
            E= (System.antes + System.despues)/2 -0.197
        elif Vc > Veq:
            E= System.despues - 0.05916*np.log10((Vc-Veq)/Veq)
    z=(round(E,3))
    a=float(input())
    err=round(abs(z-a),3)
    print('Tu respuesta fue '+str(a)+' vs calculado '+str(z),', te pasaste por '+str(err))
    print()

def Pregunta7():
    P1,P2,P3,P4,P5,P6,P7,P8=preci('Al(OH)3', 3, 3E-34,78.00),preci('Ba(IO3)2',2,1.57E-9,487.13),preci('Ca(OH)2',2,6.5E-6,74.09),preci('FeCO3',1,2.1E-11,115.85),preci("NiS",1,4E-20,90.76),preci("AgCl",1,1.82E-10,143.32),preci("AgI",1, 8.3E-17,234.77),preci("ZnCO3",1,1E-10,125.39)
    Precipitados=[P1,P2,P3,P4,P5,P6,P7,P8]
    Pregu6=random.choice(Precipitados)
    Agua=np.arange(50,1001,50)
    V=random.choice(Agua)
    print('###################################################')
    print('################# Pregunta 7 ######################')
    print('###################################################')
    print()
    print('Calcula la masa (g) de '+Pregu6.name+' que precipita en '+str(V)+' mL de agua')
    z = (((Pregu6.kps/Pregu6.n**Pregu6.n)**(1/(Pregu6.n +1)))*(V/1000)*Pregu6.PM)
    a=float(input())
    err=round(abs(z-a)*100/z,2)
    print('Tu respuesta fue '+str(a)+' '+'vs calculado'+ ' '+str((z))+', tuviste una desviación porcentual del'+' '+str(err)+'%')
    print()    

def Juego():
    print('Bienvenido al juego " ", este código es únicamente una herramienta para guíar/acompañar las reglas ')
    print('expuestas; para empezar porfavor introduce el número de personas que van a participar en el juego,')
    print('dependiendo el número introducido las respuestas se repetirán la misma cantidad, siendo esto se')
    print('agradece, si se usa el codigo como juego, remitirse a las reglas compartidas; muchos exitos:')
    print()
    N=int(input('Ingrese el numero de participantes:'))
    [Pregunta1() for _ in range(N)]
    [Pregunta2() for _ in range(N)]
    [Pregunta3() for _ in range(N)]
    [Pregunta4() for _ in range(N)]
    [Pregunta5() for _ in range(N)]
    [Pregunta6() for _ in range(N)]
    [Pregunta7() for _ in range(N)]
    print('Necesita desempate? (Y/N)')
    a=str(input())
    desem=[Pregunta1,Pregunta2,Pregunta3,Pregunta4,Pregunta5,Pregunta6,Pregunta7]
    if a in {'Y'}:
        random.choice(desem)()
    elif a in {'N'}:
        print()
        print('Juego Finalizado')

Juego()         

        
        

        
        
        
        
        
        
        
        
        
    
    
