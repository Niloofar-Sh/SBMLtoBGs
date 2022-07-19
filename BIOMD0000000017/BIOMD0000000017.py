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
# from symfit import parameters, variables, log, Fit, GreaterThan, Model
# from symfit.core.minimizers import BasinHopping
# from gekko import GEKKO

# Figures settings

import matplotlib as mpl
import matplotlib.ticker as ticker
mpl.rc('font',family='Times New Roman')
rcParams['figure.figsize'] = 6, 4

sboRef = pd.read_csv('SBO.csv')
sboID = sboRef['http://www.geneontology.org/formats/oboInOwl#id']
sboMeaning = sboRef['http://www.w3.org/2000/01/rdf-schema#label']
sboSynonyms = sboRef['Synonyms']

reader = SBMLReader()
document = reader.readSBML("BIOMD0000000017.xml")
document.getNumErrors()

model_lib = document.getModel()
model_simple = simplesbml.loadSBMLFile("BIOMD0000000017.xml")

solutions = pd.read_csv('solutions.csv')

print ('Num compartmetns = ', model_simple.getNumCompartments())
print ('Num parameters =', model_simple.getNumParameters())
print ('Num species =', model_simple.getNumSpecies())
print ('Num reactions = ', model_simple.getNumReactions())
print (model_simple.getListOfCompartmentIds())
print (model_simple.getListOfAllSpecies())
print ('List of reactions = ', model_simple.getListOfReactionIds())
print ('List of rules = ', model_simple.getListOfRuleIds())
print ('List of parameters = ', model_simple.getListOfParameterIds())

reaction_reactants = {};
reaction_products = {};
reaction_modifiers = {};
Synthesis = [];
Dtype = []

for reaction in model_simple.getListOfReactionIds():
    reaction_reactants[reaction] = []
    reaction_products[reaction] = []

    reactNum = model_simple.getNumReactants(reaction)
    prodNum = model_simple.getNumProducts(reaction)
    modifierNum = model_simple.getNumModifiers(reaction)

    if reactNum is 0:
        Synthesis.append(reaction)
    if prodNum is 0:
        Dtype.append(reaction)

    for i in range(reactNum):
        reaction_reactants[reaction].append(
            (model_simple.getReactantStoichiometry(reaction, i), model_simple.getReactant(reaction, i)))

    for ii in range(prodNum):
        reaction_products[reaction].append(
            (model_simple.getProductStoichiometry(reaction, ii), model_simple.getProduct(reaction, ii)))

    if modifierNum is not 0:

        reaction_modifiers[reaction] = []
        for iii in range(modifierNum):
            reaction_modifiers[reaction].append(model_simple.getListOfModifiers(reaction)[0])

reversibles = [];
irreversibles = [];
for (reaction, r) in zip(model_simple.getListOfReactionIds(), range(len(model_lib.getListOfReactions()))):

    if model_lib.getListOfReactions()[r].getReversible() is True:
        reversibles.append(reaction)
    elif reaction not in [synth for synth in Synthesis] and reaction not in [d for d in Dtype]:
        irreversibles.append(reaction)

    # Check if it's a simple irreversible reaction (v = product of the reactants)
ops = {
    '+': op.add,
    '-': op.sub,
    '*': op.mul,
    '/': op.truediv,
    '^': op.xor,
}
simpleIrreversible = [];
userDefinedIrreversible = [];
for reaction in irreversibles:
    rateLaw = model_simple.getRateLaw(reaction)
    operands = []
    for x in rateLaw.split(sep=" "):
        if x in ops:
            operands.append(x)
    if len(set(operands)) == 1 and set(operands) == {'*'}:
        simpleIrreversible.append(reaction)
    else:
        userDefinedIrreversible.append(reaction)

speciesId = []
for species in model_simple.getListOfAllSpecies():
    for i, reacs in enumerate(reaction_reactants.values()):
        if species in [reac[1] for reac in reacs] and species not in [d for d in speciesId]:
            speciesId.append(species)

    for i, prods in enumerate(reaction_products.values()):
        if species in [prod[1] for prod in prods] and species not in [d for d in speciesId]:
            speciesId.append(species)

# Using the SS concentrations & setting the K for products ==> K_p/K_R : 0.001*(PSS/RSS),
# where PSS is the product of the steaqdy-state concentration of all the products and
# RSS is the product of the steaqdy-state concentration of all the reactants in each reaction

kineticRatio = []
speciesXconstant = {}
for reaction in irreversibles:
    RSS = 1;
    PSS = 1;
    for stoichiometry, reactant in reaction_reactants[reaction]:

        diff = solutions[reactant].diff()
        diff.fillna(method='bfill', inplace=True)
        if diff.all() == 0 and solutions[reactant].all() == 0:
            RSS = RSS * 1
        else:

            length = solutions[reactant].size
            RSS = RSS * pow(solutions[reactant][length - 1], stoichiometry)

    for stoichiometry, product in reaction_products[reaction]:

        diff = solutions[product].diff()
        diff.fillna(method='bfill', inplace=True)
        if diff.all() == 0 and solutions[product].all() == 0:
            PSS = PSS * 1
        else:
            length = solutions[product].size
            PSS = PSS * pow(solutions[product][length - 1], stoichiometry)

    kineticRatio.append(1e-3 * RSS / PSS)

trueValues = copy.deepcopy(kineticRatio)

rowM = []
for s in speciesId:
    rowM.append(s)
for reaction in irreversibles:
    rowM.append(reaction + ' X')

# M initial (regarding the steady-state ratios obtained in the previous stage)

M = np.zeros([len(trueValues), len(rowM)])

for i, reaction in zip(range(len(irreversibles)), irreversibles):

    for stoichiometry, reactant in reaction_reactants[irreversibles[i]]:
        for r in range(len(rowM)):
            if rowM[r] == reactant:
                M[i][r] = -stoichiometry

    for stoichiometry, product in reaction_products[irreversibles[i]]:
        for r in range(len(rowM)):
            if rowM[r] == product:
                M[i][r] = stoichiometry

    for r in range(len(rowM)):
        if rowM[r] == reaction + ' X':
            M[i][r] = 1

logTrueValues = [];

for k in trueValues:
    logTrueValues.append(math.log(k))

logTrueValues = np.array(logTrueValues)

Mpinv = np.linalg.pinv(np.array(M))
logU = Mpinv.dot(logTrueValues)

thermodynamics = [math.exp(i) for i in list(logU)]

for i in range(len(thermodynamics)):
    print(rowM[i], ' ==> ', thermodynamics[i])

for reaction in irreversibles:
    for i in range(len(rowM)):
        if rowM[i] == reaction + ' X':
            speciesXconstant[rowM[i]] = thermodynamics[i]

diff = logTrueValues - M.dot(logU)
np.exp(diff)

reactionRates = {}
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

    K_x = speciesXconstant[reaction + ' X']


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
        reactants = reactants * pow(solutions[reac], s)
    for s, pro in reaction_products[reaction]:
        products = products * pow(solutions[pro], s)

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
                            if stoichiometry != 1:
                                for TFcomp in TF:
                                    if len(TFcomp.name) == len('oneReac ' + reaction + ' ' + z.name) and set(
                                            TFcomp.name) == set('oneReac ' + reaction + ' ' + z.name):
                                        connections.append((z, (TFcomp, 0)))
                                        connections.append(((TFcomp, 1), oneR))
                            else:
                                connections.append((z, oneR))

    for oneP in oneJunc_prod:
        for xP in oneP.name.split(sep=" "):
            if (reaction) == (xP):

                for z in BGzeroJunc:
                    for stoichiometry, product in reaction_products[reaction]:
                        if z.name == product:
                            if stoichiometry != 1:
                                for TFcomp in TF:
                                    if len(TFcomp.name) == len('oneProd ' + reaction + ' ' + z.name) and set(
                                            TFcomp.name) == set('oneProd ' + reaction + ' ' + z.name):
                                        connections.append((oneP, (TFcomp, 0)))
                                        connections.append(((TFcomp, 1), z))
                            else:
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
            bgt.set_param(bg, 'k', speciesXconstant[reaction + ' X'])

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