##################################################################
#
# Program to run coupled size spectrum model using input parameters of phytoplankton and zooplankton size spectrum input, temperature and detritus
#
#
##################################################################


# There are several steps involved in running the model here:
#1. Source model code, stored in file  - 'sizebased_model_functions.r'
#3. Read in required input for model and set up matrices for storing model output.
#4. First run model from initial values in JAE 2009 paper to an equilibrium (100yrs MAY NEED TO BE LONGER) but based on new system input. 
# The model will move away from these values it is just to create new initial values of the size spectra + detritus for subsequently running the model in a new area/conditions. 
# Step 4 also requires first getting phytoplankton/zooplankton size spectrum. The input required are the slopes and intercepts
#5. Output new equilibrium values

#-------------------------------------------------------------------------------------
#
#
#
#     RUN SIZE BASED MODEL TO EQUILIBRIUM 
#
#-------------------------------------------------------------------------------------


setwd("/Users/julia/Dropbox/global f/Global model/simple example for Mike/")

source('sizebased_model_functions.r')    # Read all functions of size based model


# fixed params (these are set within the sizeparam function inside the file: 'sizebased_model_functions.r')
# read in parameters of size-based model
params<-sizeparam(equilibrium=T, dx=0.1,xmin=-12,xmax=6,xmin.consumer.u=-3,xmin.consumer.v=-3.5,tmax=300, delta_t=1/365,fmort.u =0.1,fminx.u=1, fmort.v = 0.1,fminx.v=1,depth=50,er=0.5,pp=-1,slope=-1)

params$Fmort = 0.1				     # fishing mortality rate per yr, later we will want to estimate this, but here it is just a constant

params$A.u = 640

params$A.v = 64

params$repro.on =1

params$handling = 0



# there are lots of other "default value presently set in the sizeparam() function, that can be modified 

#-----------------------------------------------------

# INPUT PARAMETERS 

# ----------------------------------------------------

# FOR EXAMPLE ONLY - THESE ONES ARE JUST MADE UP ! We will also have a few more with the new model, but similar idea.....

# log 10 intercept of plankton spectrum (units: log10 ( numbers per m^-3 at 1 g))
pp = -1

# note: in the new model we will also use slope, min and max size of plankton spectrum (4 plankton parameters including the above one)

# sea surface temperature degrees celsius ( not used here, but will affect growth if included)
sst = 20

# sea floor temperature degrees celsius (ditto but for benthic, whcih here is probably the same)
sft = 20


# detritus input, based on dafult valus form North Sea here (2009 J Anim Ecol paper), in new model we will link with depth and export functions, still only one parameter though

bio.det = W.init = 0.2

 
params$dx = 0.1
############### TEST RUN PROGRAM FOR SIMPLE INPUT

#------------------------------- NOTES: ------------------------------------------------------------------------------------

# this will run the model to "equilibrium" (actually just an abritrary number of years, here 100), can also run in "dynamic" mode, where inputs change through time (which we won't need here) after an equilibrium is reached. eash time we try to estimate parameters the model should be run to equilibrium though if possible, and we should use a better/faster way to do this 

# only gives output for final time step (can change this in sizemodel() to output all results, see "sizebased_model_functions.r")


# numerical integration: at the moment dx=0.5 and at daily time step, changing these affects numerical integration, and also runtime, should test, change them in sizeparams(), see "sizebased_model_functions.r"

#---------------------------------------------------------------------------------------------------------------------------- 

# how long does it take to run the main function for one set of input, running on my laptop?

# check runtime (on my MacBookPro: 2.7 GhZ processor, 8Gb memory) 

ptm=proc.time()
options(warn=-1)

res<-sizemodel(params=params,pp=pp,bio.det=bio.det ,sst=sst,sft=sft,U_mat=params$U.init,V_mat=params$V.init,temp.effect=T) 

         
    
print("Finished program")
print((proc.time()-ptm)/60.0)  # time is in seconds,  display minutes

                                        
           
# dx = 0.1 
# > print("Finished program")
# [1] "Finished program"
# > print((proc.time()-ptm)/60.0)  # time is in seconds,  display minutes
      # user     system    elapsed 
# 2.15286667 0.02276667 2.16591667            
  
  
# dx = 0.5            
# > print("Finished program")
# [1] "Finished program"
# > print((proc.time()-ptm)/60.0)  # time is in seconds,  display minutes
      # user     system    elapsed 
# 0.60620000 0.01193333 0.78878333           
           
           
                                       
#------------------------------------------- 

# examine outputs


# check names of output to maniplate results, plot etc:
names(res)                        

#Total biomass
TotalUbiomass<-colSums(res$U[params$ref:params$Nx,]*params$dx*10^params$x[params$ref:params$Nx]) 

TotalVbiomass<-colSums(res$V*params$dx*10^params$x) 

par(mfrow=c(3,1))
plot(TotalUbiomass,log="y")
plot(TotalVbiomass,log="y")
plot(res$W,log="y")



#Growth rates
plot(params$x[params$ref:params$Nx],res$GG.u[params$ref:params$Nx,params$Neq],log="y", type = "l", col = "blue", main = "Pelagic growth rate", xlim=c(params$x1.det,params$xmax),
     xlab = "Size", ylim = c(0.0001,100))
     
# growth rate at x=2 ( 100 g, should fall between bounds of data 0.05 and 7)
 res$GG.u[which(params$x==1)]

points(params$x[params$ref.det:params$Nx],res$GG.v[params$ref.det:params$Nx,params$Neq],log="y", type = "l", col = "red", main = "Benthic growth rate", xlim=c(params$x1.det,params$xmax),
     xlab = "Size")
res$GG.v[which(params$x==1)]

#size spectrum
plot(params$x[params$ref:params$Nx],res$U[params$ref:params$Nx,params$Neq],log="y", type = "l", col = "blue", main = "Pelagic density", xlim=c(params$x1.det,params$xmax),
     xlab = "Size")

points(params$x[params$ref.det:params$Nx],res$V[params$ref.det:params$Nx,params$Neq],log="y", type = "l", col = "red", main = "Benthic density", xlim=c(params$x1.det,params$xmax),
     xlab = "Size")



# Get the production:


# based on mortality fluxes, although here predation moratality only so not quite complete.
 totProd<-colSums(10^x[ref:end]*PM.u[ref:end,]*U[ref:end,]*dx)


#  based on growth fluxes [roughly same as the equation from Law et al. 2012]
 totProd<-colSums(10^x[ref:end]*GG.u[ref:end,]*U[ref:end,]*dx)



