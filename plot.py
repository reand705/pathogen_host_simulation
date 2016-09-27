from __future__ import unicode_literals
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import matplotlib
plt.rcParams['text.usetex']=True
plt.rc('font', family='serif')
plt.rcParams['text.latex.unicode']=True

arquivo = open("dados_do_artigo2000.txt", 'r')
alleles = []
for linha in arquivo:
    linha = linha.split()
    for i in range(len(linha)):
        alleles.append(int(linha[i]))
plt.plot(alleles, label = 'Coevolution')
arquivo2 = open("drift2000.txt", 'r')
alleles2 = []
for linha in arquivo2:
    linha = linha.split()
    for i in range(len(linha)):
        alleles2.append(int(linha[i]))
plt.plot(alleles2, label = 'Drift')
plt.axis([0,len(alleles2),0,20])
plt.legend(loc = 'best')
plt.xlabel(r'Generations')
plt.ylabel(r'Number of alleles')
plt.title(r"$N_{host} = 1000$, $N_S = 50$, $\mu_{path} = 10^{-1}$, $\mu_{host} = 10^{-5}$")
plt.show()
