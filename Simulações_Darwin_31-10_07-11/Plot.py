from __future__ import unicode_literals
# -*- coding: utf-8 -*-
import matplotlib
import matplotlib.pyplot as plt
plt.rcParams['text.usetex']=True
plt.rc('font', family='serif')
plt.rcParams['text.latex.unicode']=True
import numpy as np

# <codecell>

arquivo = open("dados_do_artigoDARWIN[9,2].txt", 'r')
alleles = []
for linha in arquivo:
    linha = linha.split()
    for i in range(len(linha)):
        alleles.append(int(linha[i]))
alleles = np.array(alleles)

# <codecell>

plt.plot(np.arange(10000),alleles, label = "Coevolution")
plt.axis([1,10000,1,20])
plt.legend(loc = 'best')
plt.xlabel("Generations")
plt.ylabel("Number of alleles")
plt.title(r"$N_{host} = 1000$, $N_S = 50$, $\mu_{path} = 10^{-1}$, $\mu_{host} = 10^{-5}$")
plt.show()
