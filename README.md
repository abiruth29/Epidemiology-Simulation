# Epidemiology-Simulation
## Description
Simulate the spread of a disease using graph based models like the SIR model and analyze the impact of different factors.​
## Objectives
1. To simulate the spreading of dengue virus using the
SIRSI compartmental model.
2. Integrating the epidemiology model with the
population genetics of mosquitoes, considering two
alleles with the dominant trait determined by the
mortality rate of mosquitoes induced by the
insecticides.
3. To investigate whether the dominance of an
insecticide resistance gene influences the epidemic
dynamics.
## An Epidemiology Model
- A epidemiology model or an epidemic model is a mathematical model that is created to replicate the spreading of a disease over a locality. 
- It is designed to explain about the rate at which disease is spreading.​
- Once the model is created , we can use it to analyze how various factors affect the spreading of a disease. We can use it to predict the spreading rate.
- It also comes in handy when we need to discover antibiotics or vaccine for the disease.
## SIR Model
- SIR – stands for susceptible, infected and recovered.​
- It is a compartmental model – it means that the model compartmentalises people into different levels – susceptible , infected and the recovered groups​
- Infected group involves all those people who are infected by the disease at a certain point of time.​
- Recovered group involves those people who were in the infected group earlier and then got recovered.​
- Susceptible group involves all those people not belonging to these two groups – that is those people who are not infected yet. They are susceptible to get infected if they come into contact with the infected people.
- So what we do is , we calculate the rate at which susceptible people become infected people – depending on the factors infecting in the susceptible person.​
- We also calculate the rate at which an infected person recovers from the disease.​
- Based on the rates, we develop mathematical equations describing them. And this constitutes the model. And then we can simulate the model using software like MATLAB.
## What is the study about
We do three simulations in this. We divide the vector 
population into three genotypes – Aa, AA and aa. 
- The first one considers the Hardy-Weinberg equilibrium condition –
we consider all the genotypes to have same mortality rates
due to insecticides. 
- In the second simulation we consider the
genotype AA and Aa to have more mortality rate than aa.
- In the third simulation we consider aa to have a higher
mortality rate than AA or Aa.
- Since A is the dominant allele over a, the genotype Aa and AA will have the same mortality
rate.
- A model that can simulate the spreading of the dengue virus ​
- The model we are designing is going to be a compartmental model like the SIR model – but it is not exactly an SIR model. It is going to be an extension of that model​

## Differential Equation System

The following system of differential equations needs to be solved:

![image](https://github.com/user-attachments/assets/9998f93b-14cc-439d-8672-181084d04fb4)

### Abbreviations:
- Sh -susceptible human population,
- Ih - infected human population
- Rh - infected human population
- Nh -human population
- p - probability of dominant allele
- q - probability of recessive allele
- p2 - probability/proportion of homozygous dominant genotype AA
- q2 -probability/proportion of homozygous recessive genotype aa
- 2pq - probability/proportion of heterozygous genotype Aa
- γ – recovery rate of human beings
- c - birth rate of human population
- μh - death rate of human population
- βh -probability of transmission of virus from a vector to human
- Iv - infected vector population.
- SAA - susceptible vector population of genotype AA
- SAa - susceptible vector population of genotype Aa
- Saa - susceptible vector population of genotype aa
- IAA - infected vector population of genotype AA
- IAa - infected vector population of genotype Aa
- Iaa - infected vector population of genotype aa
- βv -probability of transmission of virus from a human to vector
- NT -total vector population
- 𝜽 - vector oviposition rate
- k -carrying capacity of the vectors.
- μaa, - total mortality rates for recessive homozygotes
- μAa - total mortality rates for recessive heterozygote
- μAA - total mortality rates for recessive dominant homozygote
