## Libraries:

# a) Bond graphs
import BondGraphTools as bgt

# b) Data extraction
from pyomexmeta import RDF, eUriType
import os
import sys

# c) General
import copy
import difflib
import numpy as np
import pandas as pd
import math
import operator as op
import ast
import re

# d) Plot
import matplotlib.pyplot as plt
from matplotlib import markers
import matplotlib.font_manager as font_manager
import matplotlib.colors
from pylab import rcParams

# sbml
from libsbml import*
import simplesbml

# Integration and fitting
from scipy.optimize import curve_fit, minimize, least_squares, newton_krylov, basinhopping
from scipy  import integrate
import scipy.stats as stats
from symfit import parameters, variables, log, Fit, GreaterThan, Model
from symfit.core.minimizers import BasinHopping
from gekko import GEKKO

# fluxes = pd.read_csv('fluxes.txt', delimiter=' ', index_col=False)
# species = pd.read_csv('species.txt', delimiter=' ', index_col=False)
# solutions = pd.merge(fluxes,species)
solutions = pd.read_csv('solutions.csv')

sboRef = pd.read_csv('SBO.csv')
sboID = sboRef['http://www.geneontology.org/formats/oboInOwl#id']
sboMeaning = sboRef['http://www.w3.org/2000/01/rdf-schema#label']
sboSynonyms = sboRef['Synonyms']

reader = SBMLReader()
document = reader.readSBML("MODEL1004070000_url.xml")
document.getNumErrors()

model_lib = document.getModel()
model_simple = simplesbml.loadSBMLFile("MODEL1004070000_url.xml")


print ('Num compartmetns = ', model_simple.getNumCompartments())
print ('Num parameters =', model_simple.getNumParameters())
print ('Num species =', model_simple.getNumSpecies())
print ('Num reactions = ', model_simple.getNumReactions())
print (model_simple.getListOfCompartmentIds())
print (model_simple.getListOfAllSpecies())
print ('List of reactions = ', model_simple.getListOfReactionIds())
print ('List of rules = ', model_simple.getListOfRuleIds())
print ('List of parameters = ', model_simple.getListOfParameterIds())

# Categorizing the reactions
reaction_reactants = {};
reaction_products = {};
reaction_modifiers = {};
Synthesis = [];
Dtype = []

for reaction in model_simple.getListOfReactionIds():
    reaction_reactants[reaction] = []
    reaction_products[reaction] = []
    reaction_modifiers[reaction] = []

    reactNum = model_simple.getNumReactants(reaction)
    prodNum = model_simple.getNumProducts(reaction)
    modifierNum = model_simple.getNumModifiers(reaction)

    if reactNum is 0:
        Synthesis.append(reaction)
    if prodNum is 0:
        Dtype.append(reaction)

    if modifierNum is not 0:
        for iii in range(modifierNum):
            reaction_modifiers[reaction].append((1, model_simple.getListOfModifiers(reaction)[0]))

    for i in range(reactNum):
        reaction_reactants[reaction].append(
            (model_simple.getReactantStoichiometry(reaction, i), model_simple.getReactant(reaction, i)))

    for ii in range(prodNum):
        reaction_products[reaction].append(
            (model_simple.getProductStoichiometry(reaction, ii), model_simple.getProduct(reaction, ii)))

reversibles = [];
irreversibles = [];

for (reaction, r) in zip(model_simple.getListOfReactionIds(), range(len(model_lib.getListOfReactions()))):
    #     if model_simple.getNumModifiers(reaction) is not 0:
    #         irreversibles.append(reaction)

    #     else:

    if model_lib.getListOfReactions()[r].getReversible() is True:
        reversibles.append(reaction)
    else:
        if reaction not in [synth for synth in Synthesis] and reaction not in [d for d in Dtype]:
            irreversibles.append(reaction)

speciesId = []
for species in model_simple.getListOfAllSpecies():
    for i, reacs in enumerate(reaction_reactants.values()):
        if species in [reac[1] for reac in reacs] and species not in [d for d in speciesId]:
            speciesId.append(species)

    for i, prods in enumerate(reaction_products.values()):
        if species in [prod[1] for prod in prods] and species not in [d for d in speciesId]:
            speciesId.append(species)

# Reversible reactions ==> species constants approximation

Kp_to_Kr = []

for reaction in reversibles:

    phi_Xr_1 = 1
    phi_Xp_1 = 1
    phi_Xr_e = 1
    phi_Xp_e = 1

    sigma = solutions[reaction][len(solutions[reaction]) - 1] / solutions[reaction][1]

    for stoichiometry, reactant in reaction_reactants[reaction]:
        phi_Xr_1 = phi_Xr_1 * solutions[reactant][1]
        phi_Xr_e = phi_Xr_e * solutions[reactant][len(solutions[reactant]) - 1]

    #         if reactant not in rowM:
    #             rowM.append(reactant)

    for stoichiometry, product in reaction_products[reaction]:
        phi_Xp_1 = phi_Xp_1 * solutions[product][1]
        phi_Xp_e = phi_Xp_e * solutions[product][len(solutions[product]) - 1]

    #         if product not in rowM:
    #             rowM.append(product)

    Kp_to_Kr.append((sigma * phi_Xr_1 - phi_Xr_e) / (sigma * phi_Xp_1 - phi_Xp_e))

trueValues = Kp_to_Kr

rowM = []
for s in speciesId:
    rowM.append(s)

# M initial (regarding the steady-state ratios obtained in the previous stage)

M = np.zeros([len(trueValues), len(rowM)])

for reaction, i in zip(reversibles, range(len(reversibles))):
    for stoichiometry, reactant in reaction_reactants[reaction]:
        for r in range(len(rowM)):
            if len(rowM[r]) == len(reactant) and set(rowM[r]) == set(reactant):
                M[i][r] = -stoichiometry

    for stoichiometry, product in reaction_products[reaction]:
        for r in range(len(rowM)):
            if len(rowM[r]) == len(product) and set(rowM[r]) == set(product):
                M[i][r] = stoichiometry

logTrueValues = [];

for k in trueValues:
    logTrueValues.append(math.log(k))

logTrueValues = np.array(logTrueValues)

Mpinv = np.linalg.pinv(np.array(M))
logU = Mpinv.dot(logTrueValues)

thermodynamics = [math.exp(i) for i in list(logU)]
for i in range(len(thermodynamics)):
    print(rowM[i], ' ==> ', thermodynamics[i])

diff = logTrueValues - M.dot(logU)
np.exp(diff)

reactionRates = {}
for reaction in reversibles:

    reactantsSpecies = []
    productsSpecies = []
    modifierSpecies = []

    K_reactants = 1
    K_products = 1
    K_modifier = 1

    solve = tuple()

    for stoichiometry, reactant in reaction_reactants[reaction]:
        reactantsSpecies.append(reactant)
        solve = solve + (pow(solutions[reactant], stoichiometry),)

        for i in range(len(rowM)):
            if len(rowM[i]) == len(reactant) and set(rowM[i]) == set(reactant):
                K_reactants = K_reactants * pow(thermodynamics[i], stoichiometry)

    for stoichiometry, product in reaction_products[reaction]:
        productsSpecies.append(product)
        solve = solve + (pow(solutions[product], stoichiometry),)

        for i in range(len(rowM)):
            if len(rowM[i]) == len(product) and set(rowM[i]) == set(product):
                K_products = K_products * pow(thermodynamics[i], stoichiometry)


    def func(X, r):

        reactants = np.ones((len(solutions['Time'])))
        products = np.ones((len(solutions['Time'])))
        #         modifiers = np.ones((len(solutions['time'])))

        for x in range(0, len(reactantsSpecies)):
            reactants = reactants * X[x]

        for x in range(len(reactantsSpecies), len(reactantsSpecies) + len(productsSpecies)):
            products = products * X[x]

        v = r * (K_reactants * reactants - K_products * products)

        return v


    bounds = [[0], [np.inf]]
    popt, pcov = curve_fit(func, solve, solutions[reaction], maxfev=300000, bounds=bounds)
    if popt[0] == 1:
        popt, pcov = curve_fit(func, solve, solutions[reaction], maxfev=300000, method='lm')

    reactionRates[reaction] = popt[0]

    reactants = np.ones((len(solutions['Time'])))
    products = np.ones((len(solutions['Time'])))
    for s, reac in reaction_reactants[reaction]:
        reactants = reactants * pow(solutions[reac], s)
    for s, pro in reaction_products[reaction]:
        products = products * pow(solutions[pro], s)

    plt.figure()
    plt.plot(solutions['Time'], popt[0] * (K_reactants * reactants - K_products * products), 'k--')
    plt.plot(solutions['Time'], solutions[reaction])
    plt.title(reaction)


speciesXconstant = {}
for reaction in irreversibles:
    RSS = 1; PSS=1;
    for stoichiometry,reactant in reaction_reactants[reaction]:

        diff = solutions[reactant].diff()
        diff.fillna(method='bfill', inplace=True)
        if diff.all() == 0 and solutions[reactant].all() == 0:
            RSS = RSS * 1
        else:

            length = solutions[reactant].size
            RSS = RSS * pow(solutions[reactant][length-1],stoichiometry)

    for stoichiometry,product in reaction_products[reaction]:

        diff = solutions[product].diff()
        diff.fillna(method='bfill', inplace=True)
        if diff.all() == 0 and solutions[reactant].all() == 0:
            PSS = PSS * 1
        else:
            length = solutions[product].size
            PSS = PSS * pow(solutions[product][length-1],stoichiometry)

    speciesXconstant[reaction]=1e-3*RSS/PSS

for reaction in irreversibles:

    reactantsSpecies = []
    productsSpecies = []
    modifierSpecies = []

    K_reactants = 1
    K_products = 1
    K_modifier = 1

    solve = tuple()

    for stoichiometry, reactant in reaction_reactants[reaction]:
        reactantsSpecies.append(reactant)
        solve = solve + (pow(solutions[reactant], stoichiometry),)

        for i in range(len(rowM)):
            if len(rowM[i]) == len(reactant) and set(rowM[i]) == set(reactant):
                K_reactants = K_reactants * pow(thermodynamics[i], stoichiometry)

    for stoichiometry, product in reaction_products[reaction]:
        productsSpecies.append(product)
        solve = solve + (pow(solutions[product], stoichiometry),)

        for i in range(len(rowM)):
            if len(rowM[i]) == len(product) and set(rowM[i]) == set(product):
                K_products = K_products * pow(thermodynamics[i], stoichiometry)

    K_x = speciesXconstant[reaction]


    def func(X, r):

        reactants = np.ones((len(solutions['Time'])))
        products = np.ones((len(solutions['Time'])))

        for x in range(0, len(reactantsSpecies)):
            reactants = reactants * X[x]

        for x in range(len(reactantsSpecies), len(reactantsSpecies) + len(productsSpecies)):
            products = products * X[x]

        v = r * (K_reactants * reactants - K_products * products * K_x)

        return v


    bounds = [[0], [np.inf]]
    popt, pcov = curve_fit(func, solve, solutions[reaction], maxfev=300000, bounds=bounds)
    if popt[0] == 1:
        popt, pcov = curve_fit(func, solve, solutions[reaction], maxfev=300000, method='lm')

    reactionRates[reaction] = popt[0]

    reactants = np.ones((len(solutions['Time'])))
    products = np.ones((len(solutions['Time'])))
    for s, reac in reaction_reactants[reaction]:
        reactants = reactants * solutions[reac]
    for s, pro in reaction_products[reaction]:
        products = products * solutions[pro]

    plt.figure()
    plt.plot(solutions['Time'], popt[0] * (K_reactants * reactants - K_products * products * K_x), 'k--')
    plt.plot(solutions['Time'], solutions[reaction])
    plt.title(reaction)

# BG model
BGmodel = bgt.new(name='BGmodel')
idealGassConstant = 8.314
T = 310

BGspecies = [];
BGzeroJunc = [];
BGreactions = [];
oneJunc_reac = [];
oneJunc_prod = [];
TF = [];

for reaction in reversibles + irreversibles:

    BGreactions.append(bgt.new("Re", name=reaction, library="BioChem", value={'R': idealGassConstant, 'T': T}))

    oneJunc_reac.append(bgt.new("1", name='oneReac ' + reaction))
    oneJunc_prod.append(bgt.new("1", name='oneProd ' + reaction))

    for stoichiometry, reactant in reaction_reactants[reaction]:
        if stoichiometry != 1:
            TF.append(bgt.new('TF', name='oneReac ' + reaction + ' ' + reactant, value=stoichiometry))

    for stoichiometry, product in reaction_products[reaction]:
        if stoichiometry != 1:
            TF.append(bgt.new('TF', name='oneProd ' + reaction + ' ' + product, value=1 / stoichiometry))

    if reaction in irreversibles:
        BGspecies.append(
            bgt.new("Ce_S", name=reaction + ' X', library="BioChem", value={'R': idealGassConstant, 'T': T}))
        BGzeroJunc.append(bgt.new("0", name=reaction + ' X'))

    #     for forwardStoichiometry,backwardStoichiometry in refinedStoichiometries[reaction]:
#         TF.append(bgt.new('TF', name='oneReac '+reaction, value=forwardStoichiometry))
#         TF.append(bgt.new('TF', name='oneProd '+reaction, value=1/backwardStoichiometry))


for species in speciesId:
    diff = solutions[species].diff()
    diff.fillna(method='bfill', inplace=True)
    if diff.all() != 0:
        BGspecies.append(bgt.new("Ce", name=species, library="BioChem", value={'R': idealGassConstant, 'T': T}))
        BGzeroJunc.append(bgt.new("0", name=species))
    elif diff.all() == 0:
        BGspecies.append(bgt.new("Ce_S", name=species, library="BioChem", value={'R': idealGassConstant, 'T': T}))
        BGzeroJunc.append(bgt.new("0", name=species))

for comp in TF:
    bgt.add(BGmodel, comp)

for comp1, comp2 in zip(BGspecies, BGzeroJunc):
    bgt.add(BGmodel, comp1, comp2)

for comp1, comp2, comp3 in zip(BGreactions, oneJunc_reac, oneJunc_prod):
    bgt.add(BGmodel, comp1, comp2, comp3)

# connections
connections = []

for z in BGzeroJunc:
    for s in BGspecies:
        if (z.name) == (s.name):
            connections.append((z, s))

for r in BGreactions:
    for oneR, oneP in zip(oneJunc_reac, oneJunc_prod):
        if len('oneReac ' + r.name) == len(oneR.name) and ('oneReac ' + r.name) == (oneR.name):
            connections.append((oneR, r))

        if len('oneProd ' + r.name) == len(oneP.name) and ('oneProd ' + r.name) == (oneP.name):
            connections.append((r, oneP))

for reaction in reversibles + irreversibles:
    for oneR in oneJunc_reac:
        for xR in oneR.name.split(sep=" "):
            if (reaction) == (xR):

                for z in BGzeroJunc:
                    for stoichiometry, reactant in reaction_reactants[reaction]:
                        if z.name == reactant:
                            connections.append((z, oneR))

    for oneP in oneJunc_prod:
        for xP in oneP.name.split(sep=" "):
            if (reaction) == (xP):

                for z in BGzeroJunc:
                    for stoichiometry, product in reaction_products[reaction]:
                        if z.name == product:
                            connections.append((oneP, z))

                    if reaction in irreversibles:
                        if z.name == reaction + ' X':
                            connections.append((oneP, z))

for tail, head in connections:
    bgt.connect(tail, head)

# Setting the values in BGs

for BGreaction in BGreactions:
    if BGreaction.name in [key for key in reactionRates]:
        bgt.set_param(BGreaction, 'r', reactionRates[BGreaction.name])

for bg in BGspecies:
    for i in range(len(rowM)):
        if bg.name == rowM[i]:
            bgt.set_param(bg, 'k', thermodynamics[i])

for reaction in irreversibles:
    for bg in BGspecies:
        if bg.name == reaction + ' X':
            bgt.set_param(bg, 'k', speciesXconstant[reaction])

# check if all the components have gotton the correct values
for i in range(len(BGmodel.components)):
    if BGmodel.components[i].metamodel == 'C':
        print(BGmodel.components[i], BGmodel.components[i]._params['k'])
for i in range(len(BGmodel.components)):
    if BGmodel.components[i].metamodel == 'R':
        print(BGmodel.components[i], BGmodel.components[i]._params['r'])
for i in range(len(BGmodel.components)):
    if BGmodel.components[i].metamodel == 'TF':
        print(BGmodel.components[i], BGmodel.components[i]._params['r'])

x_0 = {}
for stateVar in BGmodel.state_vars:
    for reaction in irreversibles:
        if BGmodel.state_vars[stateVar][0].name == reaction + ' X':
            x_0[stateVar] = 1

for stateVar in BGmodel.state_vars:
    for bg in BGspecies:
        if bg.name in solutions:
            if len(BGmodel.state_vars[stateVar][0].name) == len(bg.name) and set(
                    BGmodel.state_vars[stateVar][0].name) == set(bg.name):
                x_0[stateVar] = solutions[bg.name][0]

t_span = [0,solutions['Time'][len(solutions['Time'])-1]]
t, x = bgt.simulate(BGmodel, timespan=t_span, control_vars={}, x0=x_0, dt=solutions['Time'][1]-solutions['Time'][0])

# Saving the simulations as a csv file
xDict = {}
for i in range(np.size(x,1)):
    xDict[i] = x[:,i]

bgData = pd.DataFrame(xDict)
bgData.to_csv("bgData.csv",  sep=',', index=False)

for stateVar,i in zip(BGmodel.state_vars,range(len(BGmodel.state_vars))):
    if BGmodel.state_vars[stateVar][0].name in solutions:
        plt.figure()
        plt.title(BGmodel.state_vars[stateVar][0].name)
        plt.plot(solutions['Time'],solutions[BGmodel.state_vars[stateVar][0].name],label='Original')
        plt.plot(t,x[:,i],'--k',label='approximation')
        plt.legend(loc='lower right')
        plt.ylabel('Concentration(mM)')
        plt.xlabel('Time(s)')
plt.show()## Libraries:

# a) Bond graphs
import BondGraphTools as bgt

# b) Data extraction
from pyomexmeta import RDF, eUriType
import os
import sys

# c) General
import copy
import difflib
import numpy as np
import pandas as pd
import math
import operator as op
import ast
import re

# d) Plot
import matplotlib.pyplot as plt
from matplotlib import markers
import matplotlib.font_manager as font_manager
import matplotlib.colors
from pylab import rcParams

# sbml
from libsbml import*
import simplesbml

# Integration and fitting
from scipy.optimize import curve_fit, minimize, least_squares, newton_krylov, basinhopping
from scipy  import integrate
import scipy.stats as stats
from symfit import parameters, variables, log, Fit, GreaterThan, Model
from symfit.core.minimizers import BasinHopping
from gekko import GEKKO

# fluxes = pd.read_csv('fluxes.txt', delimiter=' ', index_col=False)
# species = pd.read_csv('species.txt', delimiter=' ', index_col=False)
# solutions = pd.merge(fluxes,species)
solutions = pd.read_csv('solutions.csv')

sboRef = pd.read_csv('SBO.csv')
sboID = sboRef['http://www.geneontology.org/formats/oboInOwl#id']
sboMeaning = sboRef['http://www.w3.org/2000/01/rdf-schema#label']
sboSynonyms = sboRef['Synonyms']

reader = SBMLReader()
document = reader.readSBML("MODEL1004070000_url.xml")
document.getNumErrors()

model_lib = document.getModel()
model_simple = simplesbml.loadSBMLFile("MODEL1004070000_url.xml")


print ('Num compartmetns = ', model_simple.getNumCompartments())
print ('Num parameters =', model_simple.getNumParameters())
print ('Num species =', model_simple.getNumSpecies())
print ('Num reactions = ', model_simple.getNumReactions())
print (model_simple.getListOfCompartmentIds())
print (model_simple.getListOfAllSpecies())
print ('List of reactions = ', model_simple.getListOfReactionIds())
print ('List of rules = ', model_simple.getListOfRuleIds())
print ('List of parameters = ', model_simple.getListOfParameterIds())

# Categorizing the reactions
reaction_reactants = {};
reaction_products = {};
reaction_modifiers = {};
Synthesis = [];
Dtype = []

for reaction in model_simple.getListOfReactionIds():
    reaction_reactants[reaction] = []
    reaction_products[reaction] = []
    reaction_modifiers[reaction] = []

    reactNum = model_simple.getNumReactants(reaction)
    prodNum = model_simple.getNumProducts(reaction)
    modifierNum = model_simple.getNumModifiers(reaction)

    if reactNum is 0:
        Synthesis.append(reaction)
    if prodNum is 0:
        Dtype.append(reaction)

    if modifierNum is not 0:
        for iii in range(modifierNum):
            reaction_modifiers[reaction].append((1, model_simple.getListOfModifiers(reaction)[0]))

    for i in range(reactNum):
        reaction_reactants[reaction].append(
            (model_simple.getReactantStoichiometry(reaction, i), model_simple.getReactant(reaction, i)))

    for ii in range(prodNum):
        reaction_products[reaction].append(
            (model_simple.getProductStoichiometry(reaction, ii), model_simple.getProduct(reaction, ii)))

reversibles = [];
irreversibles = [];

for (reaction, r) in zip(model_simple.getListOfReactionIds(), range(len(model_lib.getListOfReactions()))):
    #     if model_simple.getNumModifiers(reaction) is not 0:
    #         irreversibles.append(reaction)

    #     else:

    if model_lib.getListOfReactions()[r].getReversible() is True:
        reversibles.append(reaction)
    else:
        if reaction not in [synth for synth in Synthesis] and reaction not in [d for d in Dtype]:
            irreversibles.append(reaction)

speciesId = []
for species in model_simple.getListOfAllSpecies():
    for i, reacs in enumerate(reaction_reactants.values()):
        if species in [reac[1] for reac in reacs] and species not in [d for d in speciesId]:
            speciesId.append(species)

    for i, prods in enumerate(reaction_products.values()):
        if species in [prod[1] for prod in prods] and species not in [d for d in speciesId]:
            speciesId.append(species)

# Reversible reactions ==> species constants approximation

Kp_to_Kr = []

for reaction in reversibles:

    phi_Xr_1 = 1
    phi_Xp_1 = 1
    phi_Xr_e = 1
    phi_Xp_e = 1

    sigma = solutions[reaction][len(solutions[reaction]) - 1] / solutions[reaction][1]

    for stoichiometry, reactant in reaction_reactants[reaction]:
        phi_Xr_1 = phi_Xr_1 * solutions[reactant][1]
        phi_Xr_e = phi_Xr_e * solutions[reactant][len(solutions[reactant]) - 1]

    #         if reactant not in rowM:
    #             rowM.append(reactant)

    for stoichiometry, product in reaction_products[reaction]:
        phi_Xp_1 = phi_Xp_1 * solutions[product][1]
        phi_Xp_e = phi_Xp_e * solutions[product][len(solutions[product]) - 1]

    #         if product not in rowM:
    #             rowM.append(product)

    Kp_to_Kr.append((sigma * phi_Xr_1 - phi_Xr_e) / (sigma * phi_Xp_1 - phi_Xp_e))

trueValues = Kp_to_Kr

rowM = []
for s in speciesId:
    rowM.append(s)

# M initial (regarding the steady-state ratios obtained in the previous stage)

M = np.zeros([len(trueValues), len(rowM)])

for reaction, i in zip(reversibles, range(len(reversibles))):
    for stoichiometry, reactant in reaction_reactants[reaction]:
        for r in range(len(rowM)):
            if len(rowM[r]) == len(reactant) and set(rowM[r]) == set(reactant):
                M[i][r] = -stoichiometry

    for stoichiometry, product in reaction_products[reaction]:
        for r in range(len(rowM)):
            if len(rowM[r]) == len(product) and set(rowM[r]) == set(product):
                M[i][r] = stoichiometry

logTrueValues = [];

for k in trueValues:
    logTrueValues.append(math.log(k))

logTrueValues = np.array(logTrueValues)

Mpinv = np.linalg.pinv(np.array(M))
logU = Mpinv.dot(logTrueValues)

thermodynamics = [math.exp(i) for i in list(logU)]
for i in range(len(thermodynamics)):
    print(rowM[i], ' ==> ', thermodynamics[i])

diff = logTrueValues - M.dot(logU)
np.exp(diff)

reactionRates = {}
for reaction in reversibles:

    reactantsSpecies = []
    productsSpecies = []
    modifierSpecies = []

    K_reactants = 1
    K_products = 1
    K_modifier = 1

    solve = tuple()

    for stoichiometry, reactant in reaction_reactants[reaction]:
        reactantsSpecies.append(reactant)
        solve = solve + (pow(solutions[reactant], stoichiometry),)

        for i in range(len(rowM)):
            if len(rowM[i]) == len(reactant) and set(rowM[i]) == set(reactant):
                K_reactants = K_reactants * pow(thermodynamics[i], stoichiometry)

    for stoichiometry, product in reaction_products[reaction]:
        productsSpecies.append(product)
        solve = solve + (pow(solutions[product], stoichiometry),)

        for i in range(len(rowM)):
            if len(rowM[i]) == len(product) and set(rowM[i]) == set(product):
                K_products = K_products * pow(thermodynamics[i], stoichiometry)


    def func(X, r):

        reactants = np.ones((len(solutions['Time'])))
        products = np.ones((len(solutions['Time'])))
        #         modifiers = np.ones((len(solutions['time'])))

        for x in range(0, len(reactantsSpecies)):
            reactants = reactants * X[x]

        for x in range(len(reactantsSpecies), len(reactantsSpecies) + len(productsSpecies)):
            products = products * X[x]

        v = r * (K_reactants * reactants - K_products * products)

        return v


    bounds = [[0], [np.inf]]
    popt, pcov = curve_fit(func, solve, solutions[reaction], maxfev=300000, bounds=bounds)
    if popt[0] == 1:
        popt, pcov = curve_fit(func, solve, solutions[reaction], maxfev=300000, method='lm')

    reactionRates[reaction] = popt[0]

    reactants = np.ones((len(solutions['Time'])))
    products = np.ones((len(solutions['Time'])))
    for s, reac in reaction_reactants[reaction]:
        reactants = reactants * pow(solutions[reac], s)
    for s, pro in reaction_products[reaction]:
        products = products * pow(solutions[pro], s)

    plt.figure()
    plt.plot(solutions['Time'], popt[0] * (K_reactants * reactants - K_products * products), 'k--')
    plt.plot(solutions['Time'], solutions[reaction])
    plt.title(reaction)


speciesXconstant = {}
for reaction in irreversibles:
    RSS = 1; PSS=1;
    for stoichiometry,reactant in reaction_reactants[reaction]:

        diff = solutions[reactant].diff()
        diff.fillna(method='bfill', inplace=True)
        if diff.all() == 0 and solutions[reactant].all() == 0:
            RSS = RSS * 1
        else:

            length = solutions[reactant].size
            RSS = RSS * pow(solutions[reactant][length-1],stoichiometry)

    for stoichiometry,product in reaction_products[reaction]:

        diff = solutions[product].diff()
        diff.fillna(method='bfill', inplace=True)
        if diff.all() == 0 and solutions[reactant].all() == 0:
            PSS = PSS * 1
        else:
            length = solutions[product].size
            PSS = PSS * pow(solutions[product][length-1],stoichiometry)

    speciesXconstant[reaction]=1e-3*RSS/PSS

for reaction in irreversibles:

    reactantsSpecies = []
    productsSpecies = []
    modifierSpecies = []

    K_reactants = 1
    K_products = 1
    K_modifier = 1

    solve = tuple()

    for stoichiometry, reactant in reaction_reactants[reaction]:
        reactantsSpecies.append(reactant)
        solve = solve + (pow(solutions[reactant], stoichiometry),)

        for i in range(len(rowM)):
            if len(rowM[i]) == len(reactant) and set(rowM[i]) == set(reactant):
                K_reactants = K_reactants * pow(thermodynamics[i], stoichiometry)

    for stoichiometry, product in reaction_products[reaction]:
        productsSpecies.append(product)
        solve = solve + (pow(solutions[product], stoichiometry),)

        for i in range(len(rowM)):
            if len(rowM[i]) == len(product) and set(rowM[i]) == set(product):
                K_products = K_products * pow(thermodynamics[i], stoichiometry)

    K_x = speciesXconstant[reaction]


    def func(X, r):

        reactants = np.ones((len(solutions['Time'])))
        products = np.ones((len(solutions['Time'])))

        for x in range(0, len(reactantsSpecies)):
            reactants = reactants * X[x]

        for x in range(len(reactantsSpecies), len(reactantsSpecies) + len(productsSpecies)):
            products = products * X[x]

        v = r * (K_reactants * reactants - K_products * products * K_x)

        return v


    bounds = [[0], [np.inf]]
    popt, pcov = curve_fit(func, solve, solutions[reaction], maxfev=300000, bounds=bounds)
    if popt[0] == 1:
        popt, pcov = curve_fit(func, solve, solutions[reaction], maxfev=300000, method='lm')

    reactionRates[reaction] = popt[0]

    reactants = np.ones((len(solutions['Time'])))
    products = np.ones((len(solutions['Time'])))
    for s, reac in reaction_reactants[reaction]:
        reactants = reactants * solutions[reac]
    for s, pro in reaction_products[reaction]:
        products = products * solutions[pro]

    plt.figure()
    plt.plot(solutions['Time'], popt[0] * (K_reactants * reactants - K_products * products * K_x), 'k--')
    plt.plot(solutions['Time'], solutions[reaction])
    plt.title(reaction)

# BG model
BGmodel = bgt.new(name='BGmodel')
idealGassConstant = 8.314
T = 310

BGspecies = [];
BGzeroJunc = [];
BGreactions = [];
oneJunc_reac = [];
oneJunc_prod = [];
TF = [];

for reaction in reversibles + irreversibles:

    BGreactions.append(bgt.new("Re", name=reaction, library="BioChem", value={'R': idealGassConstant, 'T': T}))

    oneJunc_reac.append(bgt.new("1", name='oneReac ' + reaction))
    oneJunc_prod.append(bgt.new("1", name='oneProd ' + reaction))

    for stoichiometry, reactant in reaction_reactants[reaction]:
        if stoichiometry != 1:
            TF.append(bgt.new('TF', name='oneReac ' + reaction + ' ' + reactant, value=stoichiometry))

    for stoichiometry, product in reaction_products[reaction]:
        if stoichiometry != 1:
            TF.append(bgt.new('TF', name='oneProd ' + reaction + ' ' + product, value=1 / stoichiometry))

    if reaction in irreversibles:
        BGspecies.append(
            bgt.new("Ce_S", name=reaction + ' X', library="BioChem", value={'R': idealGassConstant, 'T': T}))
        BGzeroJunc.append(bgt.new("0", name=reaction + ' X'))

    #     for forwardStoichiometry,backwardStoichiometry in refinedStoichiometries[reaction]:
#         TF.append(bgt.new('TF', name='oneReac '+reaction, value=forwardStoichiometry))
#         TF.append(bgt.new('TF', name='oneProd '+reaction, value=1/backwardStoichiometry))


for species in speciesId:
    diff = solutions[species].diff()
    diff.fillna(method='bfill', inplace=True)
    if diff.all() != 0:
        BGspecies.append(bgt.new("Ce", name=species, library="BioChem", value={'R': idealGassConstant, 'T': T}))
        BGzeroJunc.append(bgt.new("0", name=species))
    elif diff.all() == 0:
        BGspecies.append(bgt.new("Ce_S", name=species, library="BioChem", value={'R': idealGassConstant, 'T': T}))
        BGzeroJunc.append(bgt.new("0", name=species))

for comp in TF:
    bgt.add(BGmodel, comp)

for comp1, comp2 in zip(BGspecies, BGzeroJunc):
    bgt.add(BGmodel, comp1, comp2)

for comp1, comp2, comp3 in zip(BGreactions, oneJunc_reac, oneJunc_prod):
    bgt.add(BGmodel, comp1, comp2, comp3)

# connections
connections = []

for z in BGzeroJunc:
    for s in BGspecies:
        if (z.name) == (s.name):
            connections.append((z, s))

for r in BGreactions:
    for oneR, oneP in zip(oneJunc_reac, oneJunc_prod):
        if len('oneReac ' + r.name) == len(oneR.name) and ('oneReac ' + r.name) == (oneR.name):
            connections.append((oneR, r))

        if len('oneProd ' + r.name) == len(oneP.name) and ('oneProd ' + r.name) == (oneP.name):
            connections.append((r, oneP))

for reaction in reversibles + irreversibles:
    for oneR in oneJunc_reac:
        for xR in oneR.name.split(sep=" "):
            if (reaction) == (xR):

                for z in BGzeroJunc:
                    for stoichiometry, reactant in reaction_reactants[reaction]:
                        if z.name == reactant:
                            connections.append((z, oneR))

    for oneP in oneJunc_prod:
        for xP in oneP.name.split(sep=" "):
            if (reaction) == (xP):

                for z in BGzeroJunc:
                    for stoichiometry, product in reaction_products[reaction]:
                        if z.name == product:
                            connections.append((oneP, z))

                    if reaction in irreversibles:
                        if z.name == reaction + ' X':
                            connections.append((oneP, z))

for tail, head in connections:
    bgt.connect(tail, head)

# Setting the values in BGs

for BGreaction in BGreactions:
    if BGreaction.name in [key for key in reactionRates]:
        bgt.set_param(BGreaction, 'r', reactionRates[BGreaction.name])

for bg in BGspecies:
    for i in range(len(rowM)):
        if bg.name == rowM[i]:
            bgt.set_param(bg, 'k', thermodynamics[i])

for reaction in irreversibles:
    for bg in BGspecies:
        if bg.name == reaction + ' X':
            bgt.set_param(bg, 'k', speciesXconstant[reaction])

# check if all the components have gotton the correct values
for i in range(len(BGmodel.components)):
    if BGmodel.components[i].metamodel == 'C':
        print(BGmodel.components[i], BGmodel.components[i]._params['k'])
for i in range(len(BGmodel.components)):
    if BGmodel.components[i].metamodel == 'R':
        print(BGmodel.components[i], BGmodel.components[i]._params['r'])
for i in range(len(BGmodel.components)):
    if BGmodel.components[i].metamodel == 'TF':
        print(BGmodel.components[i], BGmodel.components[i]._params['r'])

x_0 = {}
for stateVar in BGmodel.state_vars:
    for reaction in irreversibles:
        if BGmodel.state_vars[stateVar][0].name == reaction + ' X':
            x_0[stateVar] = 1

for stateVar in BGmodel.state_vars:
    for bg in BGspecies:
        if bg.name in solutions:
            if len(BGmodel.state_vars[stateVar][0].name) == len(bg.name) and set(
                    BGmodel.state_vars[stateVar][0].name) == set(bg.name):
                x_0[stateVar] = solutions[bg.name][0]

t_span = [0,solutions['Time'][len(solutions['Time'])-1]]
t, x = bgt.simulate(BGmodel, timespan=t_span, control_vars={}, x0=x_0, dt=solutions['Time'][1]-solutions['Time'][0])

# Saving the simulations as a csv file
xDict = {}
for i in range(np.size(x,1)):
    xDict[i] = x[:,i]

bgData = pd.DataFrame(xDict)
bgData.to_csv("bgData.csv",  sep=',', index=False)

for stateVar,i in zip(BGmodel.state_vars,range(len(BGmodel.state_vars))):
    if BGmodel.state_vars[stateVar][0].name in solutions:
        plt.figure()
        plt.title(BGmodel.state_vars[stateVar][0].name)
        plt.plot(solutions['Time'],solutions[BGmodel.state_vars[stateVar][0].name],label='Original')
        plt.plot(t,x[:,i],'--k',label='approximation')
        plt.legend(loc='lower right')
        plt.ylabel('Concentration(mM)')
        plt.xlabel('Time(s)')
plt.show()