# Autophagic-Body-Number-Simulation
R simulation and Excel numberical methods for estimating the original size and number of autophagic bodies in the vacuole from observed cross-sections in TEM (Transmission Electron Microscopy) sections. 

The simulation creates a random (normally distributed) number of spherical bodies with random radii (lognormally distrubted) within a vacuole of random (normall distributed) size.  The bodies are (optionally) clustered together, a slice is taken at a random position (uniform distribution), and the size and number of body cross-sections captured in the slice is recorded.  The distribution of body cross-sections per vacuole cross-section for each simulation is compared to the experimentally determined distribution of body-cross sections. Each of the random parameters can be varied as appropriate to find the best fit to the experimental data - in particular, this version of the simulation is designed to find the best fit for body number.     

Citations:
Original simulation, just for size estimation: Xie Z, Nair U, Geng J, Szefler MB, Rothman ED, Klionsky DJ. Indirect estimation of the area density of Atg8 on the phagophore. Autophagy 2009; 5:217–20. Available from: https://doi.org/10.4161/auto.5.2.7201

Excel numerical methods and second published version of the simulation, for both size and number: Backues SK, Chen D, Ruan J, Xie Z, Klionsky DJ. Estimating the size and number of autophagic bodies by electron microscopy. Autophagy 2014; 10:155–64. https://doi.org/10.4161/auto.26856

Third published version, improved for number: Cawthon H, Chakraborty R, Roberts JR, Backues SK. Control of autophagosome size and number by Atg7. Biochem Biophys Res Commun. 2018 Sep 5;503(2):651-656. https://doi.org/10.1016/j.bbrc.2018.06.056

The current version includes extra code that simplifies the workflow written by Pat Wall, a former undergraduate in my lab and now a PhD student at IU Bloomington.  

How to use:
1) Measure autophagic body cross-sections - detailed instructions in the Backues 2014 protocol paper (https://doi.org/10.4161/auto.26856)
2) Use the two Excel sheets to get an initial estimate of body size and number from these measurements - detailed instructions also the Backues 2014 protocol paper (https://doi.org/10.4161/auto.26856)
3) Download the files in the "simulation" folder to your working R directory
4) Replace the data in the "experimental data" excel file with your own experimental data (number of observed body sections in each vacuole section).  This is used for comparing the results of the simulation
5) Load the simulation into R, and run it using the command resim()
6) Choose whether to use the function defaults, values from resim_input.txt, or input function values manually.  You should derive most of the function values from the initial estimates you made using the numerical methods found in the Excel sheets.  The number of simulated bodies (mean and standard deviation) are the variables to alter to try to find the best fit.  

The output of the simulation includes graphs of the mean number of simulated body cross-sections per vacuole section, the KS test (Kolmogorov–Smirnov test) results for each simulation vs your experimental data, and split-violin plots of the distribution of cross-sections per vacuole section for the simulations with the lowest KS scores.  

Other usage notes:
- "Repeats" allows you to test multiple different numbers of simulated bodies in one run.  You can either iteratively change the mean number of bodies (bm), or the standard deviation in the number of bodies (bsd)
- "Clustering" allows you to choose from one of three clustering options.  "No Clustering" is by far the fastest, and gives good results for body size and not bad results for body number.  "GLPK clustering" is the type of clustering implemented in the orginal simulation (Xie et al. 2009), and may possibly make the body number estimations more accurate.  "Type 2" clustering is an unpublished alternative developed by Pat Wall and Steven Backues at Eastern Michigan University, and gives clusters that look more realistic; however, the results may not be significantly different from those obtained with GLPK clustering.  Also, the Type 2 clustering often fails for body numbers greater than 15-20.  
- "Cores" allows you to run the simulation on multiple cores.  Note that this is somewhat buggy - it works for Linux-based systems, but often generates errors on windows-based systems, so it is best to just use 1 core when running the simulation in Windows.    
