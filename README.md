# SBML to bond graphs: from conversion to composition
## Repository for the files related to the above-mentioned paper

In this repository I've saved the files and functions required for an automatic conversion of SBML models into Bond Graphs (BGs). This is demonstrated by converting two SBML models into bond graphs and later automatically coupling them with our previously developed framework. We have worked on 4 categories of reactions to be converted into bond graphs: reversible/irreversible mass action and reversible/irreversible Michaelis Menten. Any other type will be tried to be apprximated using mass action kinetics (might or might not represent a proper equivalent bond graph model). It's woth noting that our framework recognizes the reaction types by checking their SBO terms. We manually added SBO terms to the reactions where they were not annotated. The selected SBML models must meet certain criteria to fit into our developed framework:

a) The reactions must be described in terms of the Marcelin-de Donder kinetics (using thermodynamic parameters), reversible/irreversible mass action kinetics (or any other rate law that generates similar behaviours), or reversible/irreversible Michaelis-Menten kinetics;

b) Since the role of modifiers in SBML models is not fixed and an automatic interpretation of their function is not possible in the current SBML version, any selected SBML model must not include modifiers;

c) Considering the fact that the bond graph components in our chosen Python package (BondGraphTools) do not support non-constant and time-variant parameters, our method only works for models with constant parameters;

d) Some SBML models incorporate ``events'', which are discontinuous changes in models triggered by certain conditions or at certain times. Since these events are not particularly part of physical or biological systems, we have not covered them in our conversion framework. Some of such events can be later applied to the generated bond model, for example, removing/adding species to the system or inserting a constant flow of a concentration into the system.


Based on the manuscript, the exemplar SBML models along with their bond graph equivalents are as follows:

**1. BIOMD0000000017:** Pyruvate distribution
All the reactions have arbitrary formulations but they are all irreversible. We tried to approximate them with irreversible mass action formulation. Find the files here: [https://github.com/Niloofar-Sh/SBMLtoBGs/tree/main/BIOMD0000000017].

**2. MODEL1004070000:** The pentose phosphate pathway
The reactions were described in reversible and irreversible mass action kinetics. Find the files here: [https://github.com/Niloofar-Sh/SBMLtoBGs/tree/main/MODEL1004070000].

The composition of the two bond graph models is given in:

**3. GeneralApproachGit:** Find the files here: [https://github.com/Niloofar-Sh/SBMLtoBGs/tree/main/GeneralApproachGit].

Since we couldn't find a model that would both meet our criteria and include all the four types of supported reactions (reversible/irreversible mass action and reversible/irreversible Michaelis Menten), we created our desired SBML model: 

**4. reversibleMA_reversibleMM_irreversibleMM:** includes all four types of supported reactions. Find the files here: [https://github.com/Niloofar-Sh/SBMLtoBGs/blob/main/reversibleMA_reversibleMM_irreversibleMM/reversMAMM_irreversMMMA.ipynb].

